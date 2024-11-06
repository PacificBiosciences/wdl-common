version 1.0

import "../structs.wdl"

task sawfish_discover {
  meta {
    description: "Discover structural variant signatures with sawfish."
  }

  parameter_meta {
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
    expected_male_bed: {
      name: "Expected CN BED for sample with XY karyotype"
    }
    expected_female_bed: {
      name: "Expected CN BED for sample with XX karyotype"
    }
    out_prefix: {
      name: "Output prefix"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    discover_tar: {
      name: "Tarballed output of sawfish discover"
    }
  }

  input {
    String? sex

    File aligned_bam
    File aligned_bam_index

    File ref_fasta
    File ref_index

    File expected_male_bed
    File expected_female_bed

    String out_prefix

    RuntimeAttributes runtime_attributes
  }

  File expected_bed = if select_first([sex, "FEMALE"]) == "MALE" then expected_male_bed else expected_female_bed

  Int threads   = 16
  Int mem_gb    = threads * 8
  Int disk_size = ceil(size(aligned_bam, "GB") * 2 + size(ref_fasta, "GB") + 20)

  command <<<
    set -euo pipefail

    # sawfish stores relative filepaths of input files in the output directory
    # symlink input files into the working directory
    for i in ~{aligned_bam} ~{aligned_bam_index} ~{ref_fasta} ~{ref_index} ~{expected_bed}; do
      ln -s $i .
    done

    sawfish --version

    echo ~{if defined(sex) then "" else "Sex is not defined for ~{out_prefix}.  Defaulting to karyotype XX for sawfish."}

    sawfish discover \
      --threads ~{threads} \
      --disable-path-canonicalization \
      --ref ~{basename(ref_fasta)} \
      --expected-cn ~{basename(expected_bed)} \
      --bam ~{basename(aligned_bam)} \
      --output-dir ~{out_prefix}

    tar -cvf ~{out_prefix}.tar ~{out_prefix}

  >>>

  output {
    File discover_tar = "~{out_prefix}.tar"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/sawfish@sha256:fcd5d091908322ddeb2c86b7217b7cfdef9a103944adb3e87c76d495eb3fea5b"
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

task sawfish_call {
  meta {
    description: "Call structural variants from signatures with sawfish call."
  }

  parameter_meta {
    discover_tars: {
      name: "Tarballed output of sawfish discover"
    }
    aligned_bams: {
      name: "Aligned BAMs"
    }
    aligned_bam_indices: {
      name: "Aligned BAM indices"
    }
    ref_fasta: {
      name: "Reference FASTA"
    }
    ref_index: {
      name: "Reference FASTA index"
    }
    out_prefix: {
      name: "Output prefix"
    }
    report_supporting_reads: {
      name: "Report supporting reads"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    vcf: {
      name: "Structural variant VCF"
    }
    vcf_index: {
      name: "Structural variant VCF index"
    }
    supporting_reads: {
      name: "Supporting reads JSON"
    }
  }

  input {
    Array[File] discover_tars

    Array[File] aligned_bams
    Array[File] aligned_bam_indices

    File ref_fasta
    File ref_index

    String out_prefix

    Boolean report_supporting_reads = true

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 16
  Int mem_gb    = threads * 2
  Int disk_size = ceil(size(aligned_bams, "GB") + size(ref_fasta, "GB") + (size(discover_tars, "GB")) * 2 + 20)

  command <<<
    set -xeuo pipefail

    # sawfish stores relative filepaths of input files in the output directory
    # symlink input files into the working directory
    for i in ~{sep=" " aligned_bams} ~{sep=" " aligned_bam_indices} ~{ref_fasta} ~{ref_index}; do
      ln -s $i .
    done

    # unpack sawfish discover output
    # add the dir names to the SAMPLES list
    PREFIX="--sample "
    SAMPLES=()

    while read -r discover_tar || [[ -n "${discover_tar}" ]]; do
      sampledir=$(basename -s .tar "${discover_tar}")
      SAMPLES+=("$sampledir")
      tar xvf "${discover_tar}"
    done < ~{write_lines(discover_tars)}

    sawfish --version

    # shellcheck disable=SC2068
    sawfish joint-call \
      --threads ~{threads} \
      ~{if report_supporting_reads then "--report-supporting-reads" else ""} \
      ${SAMPLES[@]/#/$PREFIX} \
      --output-dir ~{out_prefix}

    # rename the output files to be more informative
    mv ~{out_prefix}/genotyped.sv.vcf.gz ~{out_prefix}.vcf.gz
    mv ~{out_prefix}/genotyped.sv.vcf.gz.tbi ~{out_prefix}.vcf.gz.tbi
    mv ~{out_prefix}/supporting_reads.json.gz ~{out_prefix}.supporting_reads.json.gz
  >>>

  output {
    File vcf               = "~{out_prefix}.vcf.gz"
    File vcf_index         = "~{out_prefix}.vcf.gz.tbi"
    File? supporting_reads = "~{out_prefix}.supporting_reads.json.gz"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/sawfish@sha256:fcd5d091908322ddeb2c86b7217b7cfdef9a103944adb3e87c76d495eb3fea5b"
    cpu: threads
    memory: "~{mem_gb} GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}