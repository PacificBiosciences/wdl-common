version 1.0

import "../structs.wdl"

task hificnv {
  meta {
    description: "Call copy number reads from HiFi reads based on depth with HiFiCNV"
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
    vcf: {
      name: "Small variant VCF"
    }
    vcf_index: {
      name: "Small variant VCF index"
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
    exclude_bed: {
      name: "Regions to exclude from CNV calls"
    }
    exclude_bed_index: {
      name: "Regions to exclude from CNV calls (index)"
    }
    expected_male_bed: {
      name: "Expected CN BED for sample with XY karyotype"
    }
    expected_female_bed: {
      name: "Expected CN BED for sample with XX karyotype"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    cnv_vcf: {
      name: "HiFiCNV VCF"
    }
    cnv_vcf_index: {
      name: "HiFiCNV VCF index"
    }
    copynum_bedgraph: {
      name: "HiFiCNV CN bedgraph"
    }
    depth_bw: {
      name: "Depth bigWig"
    }
    maf_bw: {
      name: "MAF bigWig"
    }
    stat_DUP_count: {
      name: "Number of passing duplications"
    }
    stat_DUP_sum: {
      name: "Total length of passing duplications (Mbp)"
    }
    stat_DEL_count: {
      name: "Number of passing deletions"
    }
    stat_DEL_sum: {
      name: "Total length of passing deletions (Mbp)"
    }
  }

  input {
    String sample_id
    String? sex

    File aligned_bam
    File aligned_bam_index

    File vcf
    File vcf_index

    File ref_fasta
    File ref_index
    String ref_name

    File exclude_bed
    File exclude_bed_index

    File expected_male_bed
    File expected_female_bed

    RuntimeAttributes runtime_attributes
  }

  File expected_bed = if select_first([sex, "FEMALE"]) == "MALE" then expected_male_bed else expected_female_bed

  Int threads   = 8
  Int mem_gb    = threads * 2
  Int disk_size = ceil((size(aligned_bam, "GB") + size(ref_fasta, "GB"))+ 20)

  command <<<
    set -euo pipefail

    echo ~{if defined(sex) then "" else "Sex is not defined for ~{sample_id}.  Defaulting to karyotype XX for HiFiCNV."}

    hificnv --version

    hificnv \
      --threads ~{threads} \
      --bam ~{aligned_bam} \
      --ref ~{ref_fasta} \
      --maf ~{vcf} \
      --exclude ~{exclude_bed} \
      --expected-cn ~{expected_bed} \
      --output-prefix hificnv

    mv hificnv.~{sample_id}.copynum.bedgraph ~{sample_id}.~{ref_name}.hificnv.copynum.bedgraph
    mv hificnv.~{sample_id}.depth.bw ~{sample_id}.~{ref_name}.hificnv.depth.bw
    mv hificnv.~{sample_id}.maf.bw ~{sample_id}.~{ref_name}.hificnv.maf.bw
    mv hificnv.~{sample_id}.vcf.gz ~{sample_id}.~{ref_name}.hificnv.vcf.gz
    bcftools index --tbi ~{sample_id}.~{ref_name}.hificnv.vcf.gz

    bcftools query \
      -i 'FILTER="PASS" & SVTYPE="DUP"' \
      -f '%INFO/SVLEN\n' \
      ~{sample_id}.~{ref_name}.hificnv.vcf.gz \
      | wc -l > stat_DUP_count.txt || echo 0 > stat_DUP_count.txt
    bcftools query \
      -i 'FILTER="PASS" & SVTYPE="DUP"' \
      -f '%INFO/SVLEN\n' ~{sample_id}.~{ref_name}.hificnv.vcf.gz \
      | awk '{{ s+=$1 }};END {{ print s/(1000000) }}' > stat_DUP_sum.txt || echo 0 > stat_DUP_sum.txt  # TODO: what units do we want on this?
    bcftools query \
      -i 'FILTER="PASS" & SVTYPE="DEL"' \
      -f '%INFO/SVLEN\n' \
      ~{sample_id}.~{ref_name}.hificnv.vcf.gz \
      | wc -l > stat_DEL_count.txt || echo 0 > stat_DEL_count.txt
    bcftools query \
      -i 'FILTER="PASS" & SVTYPE="DEL"' \
      -f '%INFO/SVLEN\n' ~{sample_id}.~{ref_name}.hificnv.vcf.gz \
      | awk '{{ s+=$1 }};END {{ print s/(1000000) }}' > stat_DEL_sum.txt || echo 0 > stat_DEL_sum.txt  # TODO: what units do we want on this?
  >>>

  output {
    File cnv_vcf          = "~{sample_id}.~{ref_name}.hificnv.vcf.gz"
    File cnv_vcf_index    = "~{sample_id}.~{ref_name}.hificnv.vcf.gz.tbi"
    File copynum_bedgraph = "~{sample_id}.~{ref_name}.hificnv.copynum.bedgraph"
    File depth_bw         = "~{sample_id}.~{ref_name}.hificnv.depth.bw"
    File maf_bw           = "~{sample_id}.~{ref_name}.hificnv.maf.bw"
    Int stat_DUP_count    = read_int("stat_DUP_count.txt")
    Int stat_DUP_sum      = read_int("stat_DUP_sum.txt")
    Int stat_DEL_count    = read_int("stat_DEL_count.txt")
    Int stat_DEL_sum      = read_int("stat_DEL_sum.txt")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/hificnv@sha256:19fdde99ad2454598ff7d82f27209e96184d9a6bb92dc0485cc7dbe87739b3c2"
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
