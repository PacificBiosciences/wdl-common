version 1.0

import "../structs.wdl"

task hiphase {
  meta {
    description: "Phases VCFs and haplotags BAMs with HiPhase"
  }

  parameter_meta {
    sample_id: {
      name: "Sample ID"
    }
    vcfs: {
      name: "VCFs"
    }
    vcf_indices: {
      name: "VCF indices"
    }
    phased_vcf_names: {
      name: "Phased VCF names"
    }
    phased_vcf_index_names: {
      name: "Phased VCF index names"
    }
    bam: {
      name: "BAM"
    }
    bam_index: {
      name: "BAM index"
    }
    ref_name: {
      name: "Reference name"
    }
    ref_fasta: {
      name: "Reference FASTA"
    }
    ref_index: {
      name: "Reference index"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    phased_vcfs: {
      name: "Phased VCFs"
    }
    phased_vcf_indices: {
      name: "Phased VCF indices"
    }
    haplotagged_bam: {
      name: "Haplotagged BAM"
    }
    haplotagged_bam_index: {
      name: "Haplotagged BAM index"
    }
    phase_stats: {
      name: "Phasing statistics"
    }
    phase_blocks: {
      name: "Phase blocks"
    }
    phase_haplotags: {
      name: "Haplotag information"
    }
    stat_phased_basepairs: {
      name: "Phased basepairs"
    }
    stat_phase_block_ng50: {
      name: "Phase block N50"
    }
    stat_mapped_fraction: {
      name: "Mapped fraction"
    }
  }

  input {
    String sample_id

    Array[File] vcfs
    Array[File] vcf_indices
    Array[String] phased_vcf_names
    Array[String] phased_vcf_index_names

    File bam
    File bam_index

    String ref_name
    File ref_fasta
    File ref_index

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 16
  Int mem_gb    = threads * 5
  Int disk_size = ceil(size(vcfs, "GB") + size(ref_fasta, "GB") + size(bam, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    hiphase --version

    hiphase --threads ~{threads} \
      --sample-name ~{sample_id} \
      ~{sep=" " prefix("--vcf ", vcfs)} \
      ~{sep=" " prefix("--output-vcf ", phased_vcf_names)} \
      --bam ~{bam} \
      --output-bam ~{sample_id}.~{ref_name}.haplotagged.bam \
      --reference ~{ref_fasta} \
      --summary-file ~{sample_id}.~{ref_name}.hiphase.stats.tsv \
      --blocks-file ~{sample_id}.~{ref_name}.hiphase.blocks.tsv \
      --haplotag-file ~{sample_id}.~{ref_name}.hiphase.haplotags.tsv \
      --global-realignment-cputime 300

    gzip ~{sample_id}.~{ref_name}.hiphase.haplotags.tsv

    samtools idxstats ~{sample_id}.~{ref_name}.haplotagged.bam \
      | awk '{m+=$3;u+=4};END {print m/(m+u)}' \
      > mapped_fraction.txt
    
    awk -F '\t' -v SAMPLE="~{sample_id}" \
      '($1==SAMPLE && $2=="all") {{ print $20 }};' \
      ~{sample_id}.~{ref_name}.hiphase.stats.tsv \
      > phased_basepairs.txt

    awk -F '\t' -v SAMPLE="~{sample_id}" \
      '($1==SAMPLE && $2=="all") {{ print $21 }};' \
      ~{sample_id}.~{ref_name}.hiphase.stats.tsv \
      > phase_block_ng50.txt
  >>>

  output {
    Array[File] phased_vcfs        = phased_vcf_names
    Array[File] phased_vcf_indices = phased_vcf_index_names
    File  haplotagged_bam          = "~{sample_id}.~{ref_name}.haplotagged.bam"
    File  haplotagged_bam_index    = "~{sample_id}.~{ref_name}.haplotagged.bam.bai"
    File  phase_stats              = "~{sample_id}.~{ref_name}.hiphase.stats.tsv"
    File  phase_blocks             = "~{sample_id}.~{ref_name}.hiphase.blocks.tsv"
    File  phase_haplotags          = "~{sample_id}.~{ref_name}.hiphase.haplotags.tsv.gz"
    Int   stat_phased_basepairs    = read_int("phased_basepairs.txt")
    Int   stat_phase_block_ng50    = read_int("phase_block_ng50.txt")
    Float stat_mapped_fraction     = read_float("mapped_fraction.txt")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/hiphase@sha256:c46c8493be8b308c0433441cbafcc1b6ac999dfa6e85001d466ebd551c4a8cf0"
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
