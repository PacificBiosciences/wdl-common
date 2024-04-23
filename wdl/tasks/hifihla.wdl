version 1.0

import "../structs.wdl"

task hifihla_call_reads {
  meta {
    description: "Call HLA alleles from HiFi reads with HiFiHLA"
  }

  parameter_meta {
    sample_id: {
      name: "Sample ID"
    }
    aligned_bam: {
      name: "Aligned BAM"
    }
    aligned_bam_index: {
      name: "Aligned BAM Index"
    }
    ref_fasta: {
      name: "Reference FASTA"
    }
    ref_index: {
      name: "Reference FASTA Index"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    summary: {
      name: "HiFiHLA summary report"
    }
    report_tsv: {
      name: "HiFiHLA report for pharmCAT"
    }
    report_json: {
      name: "HiFiHLA full report"
    }
  }

  input {
    String sample_id
    File aligned_bam
    File aligned_bam_index

    File ref_fasta  # currently requires GRCh38_no_alt_analysis_set
    File ref_index

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int mem_gb    = threads * 2 + 2
  Int disk_size = ceil(size(aligned_bam, "GB") + size(ref_fasta, "GB") + 20)

  command <<<
    set -euo pipefail
    hifihla --version

    export POLARS_MAX_THREADS=~{threads}  # Polars doesn't respect --threads # TODO: move this to the container


    mkdir -p hifihla
    hifihla call-reads \
      --threads ~{threads} \
      --preset wgs \
      --abam ~{aligned_bam} \
      --out_prefix hifihla/~{sample_id}

    mv hifihla/* .
  >>>

  output {
    File summary                 = "~{sample_id}_hifihla_summary.tsv"
    File report_tsv              = "~{sample_id}_hifihla_report.tsv"
    File report_json             = "~{sample_id}_hifihla_report.json"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/hifihla@sha256:1510a9d80686fb5c9ca24b11549c94eae92b7b1ff3577ea2f4ee9ef834f85a9b"
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
