version 1.0

import "../structs.wdl"

task bcftools_stats_roh_small_variants {
  meta {
    description: "Run bcftools stats and bcftools roh on a small variant VCF."
  }

  parameter_meta {
    sample_id: {
      name: "Sample ID"
    }
    vcf: {
      name: "Small variant VCF"
    }
    ref_fasta: {
      name: "Reference FASTA"
    }
    ref_name: {
      name: "Reference name"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    stats: {
      name: "Small variant VCF stats"
    }
    roh_out: {
      name: "Runs of homozygosity output"
    }
    roh_bed: {
      name: "Runs of homozygosity BED"
    }
    stat_SNV_count: {
      name: "SNV count"
    }
    stat_INDEL_count: {
      name: "INDEL count"
    }
    stat_TSTV_ratio: {
      name: "Ts/Tv ratio"
    }
    stat_HETHOM_ratio: {
      name: "SNV Het/Hom ratio"
    }
  }

  input {
    String sample_id

    File vcf

    File ref_fasta
    String ref_name

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int disk_size = ceil(size(vcf, "GB") + size(ref_fasta, "GB") + 20)

  command <<<
    set -euo pipefail

    bcftools --version

    bcftools stats \
      ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
      --apply-filters PASS --samples ~{sample_id} \
      ~{"--fasta-ref " + ref_fasta} \
      ~{vcf} \
    > ~{sample_id}.~{ref_name}.vcf.stats.txt

    bcftools roh \
      ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
      --AF-dflt 0.4 \
      ~{vcf} \
    > ~{sample_id}.~{ref_name}.bcftools_roh.out

    # Convert the ROH output to BED format
    # The output format is: chrom, start, end, qual
    echo -e "#chr\\tstart\\tend\\tqual" > ~{sample_id}.~{ref_name}.roh.bed
    awk -v OFS='\t' '$1=="RG" {{ print $3, $4, $5, $8 }}' ~{sample_id}.~{ref_name}.bcftools_roh.out \
      >> ~{sample_id}.~{ref_name}.roh.bed

    # Extract the counts of SNVs and indels, Ts/Tv ratio, and SNV het/hom ratio from bcftools stats
    awk -F '\t' '($1=="SN" && $3=="number of SNPs:") {{ print $4 }}' ~{sample_id}.~{ref_name}.vcf.stats.txt \
      > snv_count.txt || echo 0 > snv_count.txt
    awk -F '\t' '($1=="SN" && $3=="number of indels:") {{ print $4 }}' ~{sample_id}.~{ref_name}.vcf.stats.txt \
      > indel_count.txt || echo 0 > indel_count.txt
    awk -F '\t' '($1=="TSTV") {{ print $5 }}' ~{sample_id}.~{ref_name}.vcf.stats.txt \
      > tstv_ratio.txt || echo 0 > tstv_ratio.txt
    awk -F '\t' '($1=="PSC") {{ print $6/$5 }}' ~{sample_id}.~{ref_name}.vcf.stats.txt \
      > hethom_ratio.txt || echo 0 > hethom_ratio.txt
  >>>

  output {
    File  stats              = "~{sample_id}.~{ref_name}.vcf.stats.txt"
    File  roh_out            = "~{sample_id}.~{ref_name}.bcftools_roh.out"
    File  roh_bed            = "~{sample_id}.~{ref_name}.roh.bed"
    Int   stat_SNV_count     = read_int("snv_count.txt")
    Int   stat_INDEL_count   = read_int("indel_count.txt")
    Float stat_TSTV_ratio    = read_float("tstv_ratio.txt")
    Float stat_HETHOM_ratio  = read_float("hethom_ratio.txt")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/bcftools@sha256:46720a7ab5feba5be06d5269454a6282deec13060e296f0bc441749f6f26fdec"
    cpu: threads
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

task concat_pbsv_vcf {
  meta {
    description: "Concatenate multiple PBSV VCFs into a single VCF."
  }

  parameter_meta {
    vcfs: {
      name: "VCFs"
    }
    vcf_indices: {
      name: "VCF indices"
    }
    out_prefix: {
      name: "Output VCF prefix"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    concatenated_vcf: {
      name: "Concatenated VCF"
    }
    concatenated_vcf_index: {
      name: "Concatenated VCF index"
    }
  }

  input {
    Array[File] vcfs
    Array[File] vcf_indices

    String out_prefix

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int disk_size = ceil(size(vcfs, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    mkdir vcfs
    while read -r input || [[ -n "${input}" ]]; do
      ln -s "${input}" vcfs
    done < ~{write_lines(flatten([vcfs,vcf_indices]))}

    find vcfs -name "*.vcf.gz" > vcf.list

    bcftools --version

    bcftools concat \
      --allow-overlaps \
      ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
      --output-type z \
      --output ~{out_prefix}.vcf.gz \
      --file-list vcf.list

    bcftools index --tbi ~{out_prefix}.vcf.gz
  >>>

  output {
    File concatenated_vcf       = "~{out_prefix}.vcf.gz"
    File concatenated_vcf_index = "~{out_prefix}.vcf.gz.tbi"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/bcftools@sha256:36d91d5710397b6d836ff87dd2a924cd02fdf2ea73607f303a8544fbac2e691f"
    cpu: threads
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

task split_vcf_by_sample {
  meta {
    description: "Split a multi-sample VCF by sample."
  }

  parameter_meta {
    sample_ids: {
      name: "Sample IDs"
    }
    vcf: {
      name: "VCF"
    }
    vcf_index: {
      name: "VCF index"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    split_vcfs: {
      name: "Split VCFs"
    }
    split_vcf_indices: {
      name: "Split VCF indices"
    }
  }

  input {
    Array[String] sample_ids
    File vcf
    File vcf_index

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int disk_size = ceil(size(vcf, "GB") * 2 + 20)

  String vcf_basename = basename(vcf, ".vcf.gz")

  command <<<
    set -euo pipefail

    bcftools --version

    for sample_id in ~{sep=" " sample_ids}; do
      bcftools view \
        ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
        --samples ${sample_id} \
        --output-type z \
        --output ${sample_id}.~{vcf_basename}.vcf.gz \
        ~{vcf}
      bcftools index --tbi ${sample_id}.~{vcf_basename}.vcf.gz
      echo ${sample_id}.~{vcf_basename}.vcf.gz >> vcf.list
      echo ${sample_id}.~{vcf_basename}.vcf.gz.tbi >> index.list
    done
  >>>

  output {
    Array[File] split_vcfs        = read_lines("vcf.list")
    Array[File] split_vcf_indices = read_lines("index.list")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/bcftools@sha256:36d91d5710397b6d836ff87dd2a924cd02fdf2ea73607f303a8544fbac2e691f"
    cpu: threads
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

task bcftools_merge {
  meta {
    description: "Merge multiple sample VCFs into a single joint VCF."
  }

  parameter_meta {
    vcfs: {
      name: "VCFs"
    }
    vcf_indices: {
      name: "VCF indices"
    }
    out_prefix: {
      name: "Output VCF name"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    merged_vcf: {
      name: "Merged VCF"
    }
    merged_vcf_index: {
      name: "Merged VCF index"
    }
  }

  input {
    Array[File] vcfs
    Array[File] vcf_indices

    String out_prefix

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int disk_size = ceil(size(vcfs, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    bcftools --version

    bcftools merge \
      ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
      --output-type z \
      --output ~{out_prefix}.vcf.gz \
      ~{sep=" " vcfs}
    bcftools index --tbi ~{out_prefix}.vcf.gz
  >>>

  output {
    File merged_vcf       = "~{out_prefix}.vcf.gz"
    File merged_vcf_index = "~{out_prefix}.vcf.gz.tbi"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/bcftools@sha256:36d91d5710397b6d836ff87dd2a924cd02fdf2ea73607f303a8544fbac2e691f"
    cpu: threads
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

task sv_stats {
  meta {
    description: "Collect statistics on structural variants in a VCF."
  }

  parameter_meta {
    vcf: {
      name: "VCF"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    stat_sv_DUP_count: {
      name: "Number of passing duplications greater than 49 bp"
    }
    stat_sv_DEL_count: {
      name: "Number of passing deletions greater than 49 bp"
    }
    stat_sv_INS_count: {
      name: "Number of passing insertions greater than 49 bp"
    }
    stat_sv_INV_count: {
      name: "Number of passing inversions greater than 49 bp"
    }
    stat_sv_BND_count: {
      name: "Number of breakends"
    }
  }

  input {
    File vcf

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int disk_size = ceil(size(vcf, "GB") + 20)

  command <<<
    # Count the number of variants of each type
    bcftools view -H -i 'FILTER="PASS" & ABS(SVLEN)>49 & SVTYPE="DUP"' \
      "~{vcf}" | wc -l \
      > stat_DUP.txt || echo 0 > stat_DUP.txt
    bcftools view -H -i 'FILTER="PASS" & ABS(SVLEN)>49 & SVTYPE="DEL"' \
      "~{vcf}" | wc -l \
      > stat_DEL.txt || echo 0 > stat_DEL.txt
    bcftools view -H -i 'FILTER="PASS" & ABS(SVLEN)>49 & SVTYPE="INS"' \
      "~{vcf}" | wc -l \
      > stat_INS.txt || echo 0 > stat_INS.txt
    bcftools view -H -i 'FILTER="PASS" & ABS(SVLEN)>49 & SVTYPE="INV"' \
      "~{vcf}" | wc -l \
      > stat_INV.txt || echo 0 > stat_INV.txt
    bcftools view -H -i 'FILTER="PASS" & SVTYPE="BND"' \
      "~{vcf}" | wc -l \
      > stat_BND.txt || echo 0 > stat_BND.txt
  >>>

  output {
    Int stat_sv_DUP_count = read_int("stat_DUP.txt")
    Int stat_sv_DEL_count = read_int("stat_DEL.txt")
    Int stat_sv_INS_count = read_int("stat_INS.txt")
    Int stat_sv_INV_count = read_int("stat_INV.txt")
    Int stat_sv_BND_count = read_int("stat_BND.txt")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/bcftools@sha256:36d91d5710397b6d836ff87dd2a924cd02fdf2ea73607f303a8544fbac2e691f"
    cpu: threads
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