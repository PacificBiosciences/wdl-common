version 1.0

import "../structs.wdl"

task split_string {
  meta {
    description: "Split a concatenated string into an array of strings"
  }

  parameter_meta {
    concatenated_string: {
      name: "Concatenated String"
    }
    delimiter: {
      name: "Delimiter"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    array: {
      name: "Array of strings"
    }
  }

  input {
    String concatenated_string
    String delimiter = ","

    RuntimeAttributes runtime_attributes
  }

  Int threads = 2
  Int mem_gb  = 2
  Int disk    = 4

  command <<<
    echo '~{sub(concatenated_string, delimiter, "\n")}'
  >>>

  output {
    Array[String] array = read_lines(stdout())
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:963c8ef4cd011cd044d5efcff2759aa37b86d0ca5739ce576f60a0ae0967292c"
    cpu: threads
    memory: "~{mem_gb} GB"
    disk: "~{disk} GB"
    disks: "local-disk ~{disk} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    queueArn: runtime_attributes.queue_arn
    zones: runtime_attributes.zones
  }
}