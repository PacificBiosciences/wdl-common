version 1.0

import "../structs.wdl"

task parabricks_pbmm2 {
  meta {
    description: "Align HiFi reads to a reference genome using NVIDIA Parabricks pbmm2."
    outputs: {
      aligned_bam: {
        description: "Aligned BAM file"
      },
      aligned_bam_index: {
        description: "Aligned BAM index file"
      }
    }
  }

  parameter_meta {
    sample_id: {
      description: "Sample ID"
    }
    hifi_reads: {
      description: "HiFi reads (BAM)"
    }
    max_reads_per_chunk: {
      description: "Maximum reads per alignment chunk"
    }
    ref_fasta: {
      description: "Reference FASTA"
    }
    ref_index: {
      description: "Reference FASTA index"
    }
    ref_name: {
      description: "Reference name"
    }
    parabricks_version: {
      name: "Parabricks container version"
    }
    runtime_attributes: {
      description: "Default runtime attribute structure"
    }
  }

  input {
    String sample_id
    File hifi_reads

    File ref_fasta
    File ref_index
    String ref_name

    String parabricks_version = "4.7.0-1.beta2"

    RuntimeAttributes runtime_attributes
  }

  String docker_image = if (runtime_attributes.backend == "AWS-HealthOmics") then runtime_attributes.container_registry else "nvcr.io/ea-nvidia-clara-parabricks" + "/clara-parabricks:~{parabricks_version}"

  Int threads                          = 128
  Int mem_gb                           = 500
  Int gpuCount                         = 2
  Int num_alignment_device_mem_buffers = 8

  Int alignment_large_pair_size = 5000
  Int max_queue_reads           = 50000
  Int max_queue_chunks          = 10000

  Int disk_size = ceil(size(hifi_reads, "GB") * 2 + size(ref_fasta, "GB") + 70)

  String movie = basename(hifi_reads, ".bam")

  command <<<
    set -euo pipefail

    # shellcheck disable=SC2034
    export TCMALLOC_MAX_TOTAL_THREAD_CACHE_BYTES=268435456

    /usr/local/parabricks/pbrun minimap2 \
      --pbmm2 --pbmm2-unmapped --eqx \
      --num-threads ~{threads} \
      --gpusort --gpuwrite \
      --num-gpus ~{gpuCount} \
      --alignment-large-pair-size ~{alignment_large_pair_size} \
      --max-queue-reads ~{max_queue_reads} \
      --max-queue-chunks ~{max_queue_chunks} \
      --num-alignment-device-mem-buffers ~{num_alignment_device_mem_buffers} \
      --ref "~{ref_fasta}" \
      --in-bam "~{hifi_reads}" \
      --out-bam "~{sample_id}.~{movie}.~{ref_name}.aligned.bam" \
  >>>

  output {
      File aligned_bam = "~{sample_id}.~{movie}.~{ref_name}.aligned.bam"
      File aligned_bam_index = "~{sample_id}.~{movie}.~{ref_name}.aligned.bam.bai"
  }

  runtime {
    docker: docker_image
    cpu: threads
    memory: mem_gb + " GiB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    gpuCount: gpuCount
    gpuType: runtime_attributes.gpuType
    acceleratorCount: gpuCount  # !UnknownRuntimeKey
    acceleratorType: runtime_attributes.gpuType  # !UnknownRuntimeKey
    nvidiaDriverVersion: "580.105.08"  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpuPlatform  # !UnknownRuntimeKey
  }
}