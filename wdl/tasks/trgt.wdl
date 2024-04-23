version 1.0

import "../structs.wdl"

task trgt {
  meta {
    description: "Genotype tandem repeats from aligned reads using TRGT."
  }

  parameter_meta {
    sample_id: {
      name: "Sample ID"
    }
    sex: {
      name: "Sample sex",
      choices: ["MALE", "FEMALE"]
    }
    aligned_bam: {
      name: "Aligned BAM"
    }
    aligned_bam_index: {
      name: "Aligned BAM index"
    }
    ref_fasta: {
      name: "Reference FASTA"
    }
    ref_index: {
      name: "Reference FASTA index"
    }
    trgt_bed: {
      name: "TRGT tandem repeat catalog BED"
    }
    out_prefix: {
      name: "Output prefix"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    bam: {
      name: "TRGT spanning reads BAM"
    }
    bam_index: {
      name: "TRGT spanning reads BAM index"
    }
    vcf: {
      name: "TRGT repeats VCF"
    }
    vcf_index: {
      name: "TRGT repeats VCF index"
    }
  }

  input {
    String sample_id
    String? sex

    File aligned_bam
    File aligned_bam_index

    File ref_fasta
    File ref_index

    File trgt_bed

    String out_prefix

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 32
  Int mem_gb    = 16
  Int disk_size = ceil((size(aligned_bam, "GB") + size(ref_fasta, "GB")) * 2 + 20)

  String karyotype = if select_first([sex, "FEMALE"]) == "MALE" then "XY" else "XX"

  command <<<
    set -euo pipefail

    echo ~{if defined(sex) then "" else "Sex is not defined for ~{sample_id}.  Defaulting to karyotype XX for TRGT."}

    trgt --version

    trgt \
      --threads ~{threads} \
      --karyotype ~{karyotype} \
      --genome ~{ref_fasta} \
      --repeats ~{trgt_bed} \
      --reads ~{aligned_bam} \
      --output-prefix ~{out_prefix}.trgt

    bcftools --version

    bcftools sort \
      --output-type z \
      --output ~{out_prefix}.trgt.sorted.vcf.gz \
      ~{out_prefix}.trgt.vcf.gz

    bcftools index \
      ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
      --tbi \
      ~{out_prefix}.trgt.sorted.vcf.gz

    samtools --version

    samtools sort \
      ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
      -m 400M \
      -o ~{out_prefix}.trgt.spanning.sorted.bam \
      ~{out_prefix}.trgt.spanning.bam

    samtools index \
      ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
      ~{out_prefix}.trgt.spanning.sorted.bam

    bcftools view --no-header --exclude-uncalled \
      ~{out_prefix}.trgt.sorted.vcf.gz \
      | wc -l > genotyped_count.txt || echo "0" > genotyped_count.txt

    bcftools view --no-header --uncalled \
      ~{out_prefix}.trgt.sorted.vcf.gz \
      | wc -l > uncalled_count.txt || echo "0" > uncalled_count.txt
  >>>

  output {
    File bam                  = "~{out_prefix}.trgt.spanning.sorted.bam"
    File bam_index            = "~{out_prefix}.trgt.spanning.sorted.bam.bai"
    File vcf                  = "~{out_prefix}.trgt.sorted.vcf.gz"
    File vcf_index            = "~{out_prefix}.trgt.sorted.vcf.gz.tbi"
    Int  stat_genotyped_count = read_int("genotyped_count.txt")
    Int  stat_uncalled_count  = read_int("uncalled_count.txt")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/trgt@sha256:1b99363622352f31a30a7b42e3bdf6c14b6a8945d7790587577c054a6f9afb25"
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