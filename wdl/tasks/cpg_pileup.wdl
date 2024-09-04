version 1.0

import "../structs.wdl"

task cpg_pileup {
  meta {
    description: "Generate CpG methylation scores from aligned reads"
  }

  parameter_meta {
    haplotagged_bam: {
      name: "Aligned BAM"
    }
    haplotagged_bam_index: {
      name: "Aligned BAM index"
    }
    min_mapq: {
      name: "Minimum mapping quality"
    }
    min_coverage: {
      name: "Minimum coverage"
    }
    out_prefix: {
      name: "Output prefix"
    }
    ref_fasta: {
      name: "Reference FASTA"
    }
    ref_index: {
      name: "Reference FASTA index"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    combined_bed: {
      name: "5mCpG BED for all alignments"
    }
    hap1_bed: {
      name: "5mCpG BED for HP1 alignments"
    }
    hap2_bed: {
      name: "5mCpG BED for HP2 alignments"
    }
    combined_bw: {
      name: "5mCpG bigWig for all alignments"
    }
    hap1_bw: {
      name: "5mCpG bigWig for HP1 alignments"
    }
    hap2_bw: {
      name: "5mCpG bigWig for HP2 alignments"
    }
  }

  input {
    File haplotagged_bam
    File haplotagged_bam_index

    Int min_mapq     = 1
    Int min_coverage = 10

    String out_prefix

    File ref_fasta
    File ref_index

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 12
  Int mem_gb    = 12
  Int disk_size = ceil((size(haplotagged_bam, "GB") + size(ref_fasta, "GB")) * 2 + 20)

  command <<<
    set -euo pipefail

    aligned_bam_to_cpg_scores --version

    aligned_bam_to_cpg_scores \
      --threads ~{threads} \
      --bam ~{haplotagged_bam} \
      --ref ~{ref_fasta} \
      --output-prefix ~{out_prefix} \
      --min-mapq ~{min_mapq} \
      --min-coverage ~{min_coverage} \
      --model "$PILEUP_MODEL_DIR"/pileup_calling_model.v1.tflite

    # count the number of CpG sites in each bed file
    wc -l < ~{out_prefix}.hap1.bed > ~{out_prefix}.hap1.bed.count \
      || echo "0" > ~{out_prefix}.hap1.bed.count
    wc -l < ~{out_prefix}.hap2.bed > ~{out_prefix}.hap2.bed.count \
      || echo "0" > ~{out_prefix}.hap2.bed.count
    wc -l < ~{out_prefix}.combined.bed > ~{out_prefix}.combined.bed.count \
      || echo "0" > ~{out_prefix}.combined.bed.count
  >>>

  output {
    File   combined_bed            = "~{out_prefix}.combined.bed"
    File   hap1_bed                = "~{out_prefix}.hap1.bed"
    File   hap2_bed                = "~{out_prefix}.hap2.bed"
    File   combined_bw             = "~{out_prefix}.combined.bw"
    File   hap1_bw                 = "~{out_prefix}.hap1.bw"
    File   hap2_bw                 = "~{out_prefix}.hap2.bw"
    String stat_hap1_cpg_count     = read_string("~{out_prefix}.hap1.bed.count")
    String stat_hap2_cpg_count     = read_string("~{out_prefix}.hap2.bed.count")
    String stat_combined_cpg_count = read_string("~{out_prefix}.combined.bed.count")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb-cpg-tools@sha256:d6e63fe3f6855cfe60f573de1ca85fab27f4a68e24a7f5691a7a805a22af292d"
    cpu: threads
    memory: mem_gb + " GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}