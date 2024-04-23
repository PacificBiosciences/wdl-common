version 1.0

import "../structs.wdl"

task pbmm2_align_wgs {
  meta {
    description: "Align HiFi reads to a reference genome with pbmm2."
  }

  parameter_meta {
    sample_id: {
      name: "Sample ID"
    }
    bam: {
      name: "HiFi reads (BAM)"
    }
    ref_fasta: {
      name: "Reference FASTA"
    }
    ref_index: {
      name: "Reference FASTA index"
    }
    ref_name: {
      name: "Reference name"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    aligned_bam: {
      name: "Aligned BAM"
    }
    aligned_bam_index: {
      name: "Aligned BAM index"
    }
    bam_stats: {
      name: "BAM stats"
    }
  }

  input {
    String sample_id
    File bam

    File ref_fasta
    File ref_index
    String ref_name

    RuntimeAttributes runtime_attributes
  }

  String movie = basename(bam, ".bam")

  Int threads   = 24
  Int mem_gb    = ceil(threads * 4)
  Int disk_size = ceil((size(bam, "GB") + size(ref_fasta, "GB")) * 4 + 20)

  command <<<
    set -euo pipefail

    pbmm2 --version

    pbmm2 align \
      --num-threads ~{threads} \
      --sort-memory 4G \
      --preset HIFI \
      --sample ~{sample_id} \
      --log-level INFO \
      --sort \
      --unmapped \
      ~{ref_fasta} \
      ~{bam} \
      ~{sample_id}.~{movie}.~{ref_name}.aligned.bam &

    # movie stats
    extract_read_length_and_qual.py \
      ~{bam} \
    > ~{sample_id}.~{movie}.read_length_and_quality.tsv

    gzip ~{sample_id}.~{movie}.read_length_and_quality.tsv

    wait
  >>>

  output {
    File aligned_bam          = "~{sample_id}.~{movie}.~{ref_name}.aligned.bam"
    File aligned_bam_index    = "~{sample_id}.~{movie}.~{ref_name}.aligned.bam.bai"
    File bam_stats            = "~{sample_id}.~{movie}.read_length_and_quality.tsv.gz"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pbmm2@sha256:ed9dcb4db98c81967fff15f50fca89c8495b1f270eee00e9bec92f46d14d7e2f"
    cpu: threads
    memory: mem_gb + " GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    queueArn: runtime_attributes.queue_arn
    zones: runtime_attributes.zones
  }
}