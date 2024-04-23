version 1.0

import "../../structs.wdl"

workflow pharmcat {
  meta {
    description: "Run PharmCAT for a sample"
  }

  parameter_meta {
    haplotagged_bam: {
      name: "Haplotagged BAM"
    }
    haplotagged_bam_index: {
      name: "Haplotagged BAM index"
    }
    phased_vcf: {
      name: "Phased small variant VCF"
    }
    phased_vcf_index: {
      name: "Phased small variant VCF index"
    }
    input_tsvs: {
      name: "Pangu, StarPhase, and HiFiHLA TSVs"
    }
    ref_fasta: {
      name: "Reference FASTA"
    }
    ref_index: {
      name: "Reference FASTA index"
    }
    pharmcat_positions: {
      name: "Pharmcat positions VCF"
    }
    pharmcat_positions_index: {
      name: "Pharmcat positions VCF index"
    }
    pharmcat_min_coverage: {
      name: "Pharmcat minimum coverage"
    }
    default_runtime_attributes: {
      name: "Runtime attribute structure"
    }
    pharmcat_missing_pgx_vcf: {
      name: "Pharmcat missing PGx VCF"
    }
    pharmcat_preprocessed_filtered_vcf: {
      name: "Pharmcat preprocessed filtered VCF"
    }
    pharmcat_match_json: {
      name: "Pharmcat match JSON"
    }
    pharmcat_phenotype_json: {
      name: "Pharmcat phenotype JSON"
    }
    pharmcat_report_html: {
      name: "Pharmcat report HTML"
    }
    pharmcat_report_json: {
      name: "Pharmcat report JSON"
    }
  }

  input {
    File haplotagged_bam
    File haplotagged_bam_index
    File phased_vcf
    File phased_vcf_index
    Array[File] input_tsvs

    File ref_fasta
    File ref_index

    File pharmcat_positions
    File pharmcat_positions_index
    Int pharmcat_min_coverage

    RuntimeAttributes default_runtime_attributes
  }

  call pharmcat_preprocess {
    input:
      vcf                      = phased_vcf,
      vcf_index                = phased_vcf_index,
      ref_fasta                = ref_fasta,
      ref_index                = ref_index,
      pharmcat_positions       = pharmcat_positions,
      pharmcat_positions_index = pharmcat_positions_index,
      runtime_attributes       = default_runtime_attributes
  }

  call filter_preprocessed_vcf {
    input:
      preprocessed_vcf      = pharmcat_preprocess.preprocessed_vcf,
      haplotagged_bam       = haplotagged_bam,
      haplotagged_bam_index = haplotagged_bam_index,
      ref_index             = ref_index,
      min_coverage          = pharmcat_min_coverage,
      runtime_attributes    = default_runtime_attributes
  }

  call run_pharmcat {
    input:
      preprocessed_filtered_vcf = filter_preprocessed_vcf.filtered_vcf,
      input_tsvs                = input_tsvs,
      runtime_attributes        = default_runtime_attributes
  }

  output {
    File? pharmcat_missing_pgx_vcf           = pharmcat_preprocess.missing_pgx_vcf
    File  pharmcat_preprocessed_filtered_vcf = filter_preprocessed_vcf.filtered_vcf

    File pharmcat_match_json     = run_pharmcat.pharmcat_match_json
    File pharmcat_phenotype_json = run_pharmcat.pharmcat_phenotype_json
    File pharmcat_report_html    = run_pharmcat.pharmcat_report_html
    File pharmcat_report_json    = run_pharmcat.pharmcat_report_json
  }
}

task pharmcat_preprocess {
  meta {
    description: "Preprocess phased VCF for PharmCAT"
  }

  parameter_meta {
    vcf: {
      name: "Phased small variant VCF"
    }
    vcf_index: {
      name: "Phased small variant VCF index"
    }
    ref_fasta: {
      name: "Reference FASTA"
    }
    ref_index: {
      name: "Reference FASTA index"
    }
    pharmcat_positions: {
      name: "Pharmcat positions VCF"
    }
    pharmcat_positions_index: {
      name: "Pharmcat positions VCF index"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    preprocessed_vcf: {
      name: "Preprocessed VCF"
    }
    missing_pgx_vcf: {
      name: "Missing PGx VCF"
    }
  }

  input {
    File vcf
    File vcf_index

    File ref_fasta
    File ref_index

    File pharmcat_positions
    File pharmcat_positions_index

    RuntimeAttributes runtime_attributes
  }

  String out_prefix = basename(vcf, ".vcf.gz")
  Int    disk_size  = ceil((size(vcf, "GB") + size(ref_fasta, "GB") + size(pharmcat_positions, "GB")) * 2 + 20)

  # TODO: host the image ourselves and get the sha256
  String docker_image = if (runtime_attributes.backend == "AWS-HealthOmics") then runtime_attributes.container_registry else "pgkb" + "/pharmcat:2.3.0"

  command <<<
    set -euo pipefail

    bcftools --version

    bcftools view \
      --apply-filters PASS \
      --output-type z \
      --output ~{out_prefix}.pass_only.vcf.gz \
      ~{vcf}

    /pharmcat/pharmcat_vcf_preprocessor.py \
      --missing-to-ref \
      -vcf ~{out_prefix}.pass_only.vcf.gz \
      -refFna ~{ref_fasta} \
      -refVcf ~{pharmcat_positions} \
      -o .
  >>>

  output {
    File preprocessed_vcf = "~{out_prefix}.pass_only.preprocessed.vcf.bgz"
    File? missing_pgx_vcf = "~{out_prefix}.missing_pgx_var.vcf"
  }

  runtime {
    docker: docker_image
    cpu: 2
    memory: "4 GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    queueArn: runtime_attributes.queue_arn
    zones: runtime_attributes.zones
  }
}

# Remove ref calls with low mean coverage for sample
task filter_preprocessed_vcf {
  meta {
    description: "Filter preprocessed VCF for sample"
  }

  parameter_meta {
    preprocessed_vcf: {
      name: "Preprocessed VCF"
    }
    haplotagged_bam: {
      name: "Haplotagged BAM"
    }
    haplotagged_bam_index: {
      name: "Haplotagged BAM index"
    }
    ref_index: {
      name: "Reference FASTA index"
    }
    min_coverage: {
      name: "Minimum coverage"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    filtered_vcf: {
      name: "Filtered VCF"
    }
  }

  input {
    File preprocessed_vcf

    File haplotagged_bam
    File haplotagged_bam_index

    File ref_index

    Int min_coverage

    RuntimeAttributes runtime_attributes
  }

  String out_prefix = basename(preprocessed_vcf, ".vcf.bgz")
  Int    disk_size  = ceil((size(preprocessed_vcf, "GB") + size(haplotagged_bam, "GB")) * 2 + 20)

  command <<<
    set -euo pipefail

    bedtools coverage \
      -sorted \
      -g ~{ref_index} \
      -f 1 \
      -header \
      -mean \
      -a ~{preprocessed_vcf} \
      -b ~{haplotagged_bam} \
    | ( sed  -u '/^#CHROM/q' ; awk '$11 >= ~{min_coverage}' | cut -f1-10 ) \
    > ~{out_prefix}.filtered.vcf
  >>>

  output {
    File filtered_vcf = "~{out_prefix}.filtered.vcf"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/bedtools@sha256:f9f6ba24ebd61dbe02898097de44486691e0a337c6fd6e26f440fed5d798e321"
    cpu: 2
    memory: "4 GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    queueArn: runtime_attributes.queue_arn
    zones: runtime_attributes.zones
  }
}

task run_pharmcat {
  meta {
    description: "Run PharmCAT for sample"
  }

  parameter_meta {
    preprocessed_filtered_vcf: {
      name: "Preprocessed filtered VCF"
    }
    input_tsvs: {
      name: "Pangu, StarPhase, and/or HiFiHLA TSVs"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    pharmcat_match_json: {
      name: "Pharmcat match JSON"
    }
    pharmcat_phenotype_json: {
      name: "Pharmcat phenotype JSON"
    }
    pharmcat_report_html: {
      name: "Pharmcat report HTML"
    }
    pharmcat_report_json: {
      name: "Pharmcat report JSON"
    }
  }

  input {
    File preprocessed_filtered_vcf
    Array[File] input_tsvs

    RuntimeAttributes runtime_attributes
  }

  String out_prefix = basename(preprocessed_filtered_vcf, ".vcf")
  Int    disk_size  = ceil(size(preprocessed_filtered_vcf, "GB") * 2 + 20)

  # TODO: host the image ourselves and get the sha256
  String docker_image = if (runtime_attributes.backend == "AWS-HealthOmics") then runtime_attributes.container_registry else "pgkb" + "/pharmcat:2.3.0"

  command <<<
    set -euo pipefail

    sort -k1,1 < ~{sep=" " input_tsvs} > merged.tsv || touch merged.tsv

    # Run pharmcat
    /pharmcat/pharmcat \
      -vcf ~{preprocessed_filtered_vcf} \
      -reporterJson \
      "$([ -s merged.tsv ] && echo '-po merged.tsv')" \
      -o .
  >>>

  output {
    File pharmcat_match_json     = "~{out_prefix}.match.json"
    File pharmcat_phenotype_json = "~{out_prefix}.phenotype.json"
    File pharmcat_report_html    = "~{out_prefix}.report.html"
    File pharmcat_report_json    = "~{out_prefix}.report.json"
  }

  runtime {
    docker: docker_image
    cpu: 2
    memory: "4 GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    queueArn: runtime_attributes.queue_arn
    zones: runtime_attributes.zones
  }
}
