version 1.0

import "../structs.wdl"

task mosdepth {
  meta {
    description: "Calculate coverage statistics using mosdepth"
  }

  parameter_meta {
    sample_id: {
      name: "Sample ID"
    }
    ref_name: {
      name: "Reference name"
    }
    aligned_bam: {
      name: "Aligned BAM"
    }
    aligned_bam_index: {
      name: "Aligned BAM index"
    }
    infer_sex: {
      name: "Infer the sex of human samples based on chrY depth"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    summary: {
      name: "Summary depth statistics"
    }
    region_bed: {
      name: "Depth BED"
    }
    inferred_sex: {
      name: "Sex inferred from chrY depth"
    }
    stat_mean_depth: {
      name: "Mean depth"
    }
  }

  input {
    String sample_id
    String ref_name

    File aligned_bam
    File aligned_bam_index

    Boolean infer_sex = false

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 4
  Int disk_size = ceil(size(aligned_bam, "GB") * 2 + 20)

  Float max_norm_female_chrY_depth = 0.1

  String out_prefix = basename(aligned_bam, ".bam")

  command <<<
    set -euo pipefail

    mosdepth --version

    mosdepth \
      --threads ~{threads - 1} \
      --by 500 \
      --no-per-base \
      --use-median \
      ~{out_prefix} \
      ~{aligned_bam}

    mv ~{out_prefix}.mosdepth.summary.txt ~{sample_id}.~{ref_name}.mosdepth.summary.txt
    mv ~{out_prefix}.regions.bed.gz ~{sample_id}.~{ref_name}.mosdepth.regions.bed.gz

    awk '($1=="total_region") {{ print $4 }}' \
      ~{sample_id}.~{ref_name}.mosdepth.summary.txt \
      > mean_depth.txt || echo 0 > mean_depth.txt

    awk -v threshold=~{max_norm_female_chrY_depth} \
      '$1 ~ /^(chr)?[[:digit:]]{{1,2}}$/ {{ acount+=1; asum+=$4 }}; $1 ~ /^(chr)?Y$/ {{ y=$4 }}; \
      END {{ y/(asum/acount) > threshold ? sex="MALE" : sex="FEMALE"; print sex }}' \
      ~{sample_id}.~{ref_name}.mosdepth.summary.txt \
      > inferred_sex.txt || echo "" > inferred_sex.txt
  >>>

  output {
    File   summary         = "~{sample_id}.~{ref_name}.mosdepth.summary.txt"
    File   region_bed      = "~{sample_id}.~{ref_name}.mosdepth.regions.bed.gz"
    String inferred_sex    = if (infer_sex) then read_string("inferred_sex.txt") else ""
    Float  stat_mean_depth = read_float("mean_depth.txt")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/mosdepth@sha256:35d5e02facf4f38742e5cae9e5fdd3807c2b431dd8d881fd246b55e6d5f7f600"
    cpu: threads
    memory: "4 GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " LOCAL"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    queueArn: runtime_attributes.queue_arn
    zones: runtime_attributes.zones
  }
}
