version 1.0

import "../../structs.wdl"

workflow parabricks_deepvariant {
  meta {
    description: "Call variants from aligned HiFi reads using Parabricks DeepVariant"
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
    regions_bed: {
      name: "Regions BED"
    }
    ref_fasta: {
      name: "Reference FASTA"
    }
    ref_index: {
      name: "Reference FAI"
    }
    ref_name: {
      name: "Reference Name"
    }
    parabricks_version: {
      name: "Parabricks Version"
    }
    default_runtime_attributes: {
      name: "Default Runtime Attributes"
    }
    vcf: {
      name: "VCF"
    }
    vcf_index: {
      name: "VCF Index"
    }
    gvcf: {
      name: "GVCF"
    }
    gvcf_index: {
      name: "GVCF Index"
    }
  }

  input {
    String sample_id
    File aligned_bam
    File aligned_bam_index

    File? regions_bed

    File ref_fasta
    File ref_index
    String ref_name

    String parabricks_version

    RuntimeAttributes default_runtime_attributes
  }

  call run_parabricks_deepvariant {
    input:
      sample_id          = sample_id,
      aligned_bam        = aligned_bam,
      aligned_bam_index  = aligned_bam_index,
      regions_bed        = regions_bed,
      ref_fasta          = ref_fasta,
      ref_index          = ref_index,
      ref_name           = ref_name,
      parabricks_version = parabricks_version,
      runtime_attributes = default_runtime_attributes
  }

  call postprocess_vcf {
    input:
      vcf                = run_parabricks_deepvariant.vcf,
      gvcf               = run_parabricks_deepvariant.gvcf,
      runtime_attributes = default_runtime_attributes
  }

  output {
    File vcf        = postprocess_vcf.vcf_gz
    File vcf_index  = postprocess_vcf.vcf_index
    File gvcf       = postprocess_vcf.gvcf_gz
    File gvcf_index = postprocess_vcf.gvcf_index
  }
}

task run_parabricks_deepvariant {
  meta {
    description: "Call variants from aligned HiFi reads using Parabricks DeepVariant"
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
    regions_bed: {
      name: "Regions BED"
    }
    ref_fasta: {
      name: "Reference FASTA"
    }
    ref_index: {
      name: "Reference FAI"
    }
    ref_name: {
      name: "Reference Name"
    }
    parabricks_version: {
      name: "Parabricks Version"
    }
    runtime_attributes: {
      name: "Runtime Attributes"
    }
    vcf: {
      name: "VCF"
    }
    gvcf: {
      name: "GVCF"
    }
  }

  input {
    String sample_id
    File aligned_bam
    File aligned_bam_index
    File? regions_bed

    File ref_fasta
    File ref_index
    String ref_name

    String parabricks_version

    RuntimeAttributes runtime_attributes
  }

  Int threads          = 48
  Int numGPUs          = 4
  Int numStreamsPerGPU = 4
  Int mem_gb           = threads * 4
  Int disk_size        = ceil((size(aligned_bam, "GB") + size(ref_fasta, "GB")) * 2 + 20)

  String docker_image = if (runtime_attributes.backend == "AWS-OMICS") then runtime_attributes.container_registry else "nvcr.io/nvidia/clara" + "/clara-parabricks:~{parabricks_version}"

  command <<<
    set -euo pipefail

    export TCMALLOC_MAX_TOTAL_THREAD_CACHE_BYTES=268435456

    /usr/local/parabricks/pbrun deepvariant \
      --num-gpus ~{numGPUs} \
      --num-streams-per-gpu ~{numStreamsPerGPU} \
      --run-partition \
      --mode pacbio \
      --gvcf \
      --ref ~{ref_fasta} \
      ~{if defined(regions_bed) then "--interval-file " + regions_bed else ""} \
      --in-bam ~{aligned_bam} \
      --out-variants ~{sample_id}.~{ref_name}.small_variants.g.vcf
  >>>

  output {
    File vcf  = "~{sample_id}.~{ref_name}.small_variants.vcf"
    File gvcf = "~{sample_id}.~{ref_name}.small_variants.g.vcf"
  }

  runtime {
    docker: docker_image
    cpu: threads
    memory: mem_gb + " GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries  # !UnknownRuntimeKey
    gpu: true
    gpuCount: numGPUs
    gpuType: runtime_attributes.gpuType
    acceleratorCount: numGPUs  # !UnknownRuntimeKey
    acceleratorType: runtime_attributes.gpuType  # !UnknownRuntimeKey
    nvidiaDriverVersion: "525.60.13" # !UnknownRuntimeKey
    zones: runtime_attributes.zones
  }
}

task postprocess_vcf {
  meta {
    description: "Filter, compress, and index VCFs"
  }

  parameter_meta {
    vcf: {
      name: "VCF"
    }
    gvcf: {
      name: "GVCF"
    }
    runtime_attributes: {
      name: "Runtime Attributes"
    }
    vcf_gz: {
      name: "Compressed VCF"
    }
    vcf_index: {
      name: "VCF Index"
    }
    gvcf_gz: {
      name: "Compressed GVCF"
    }
    gvcf_index: {
      name: "GVCF Index"
    }
  }

  input {
    File vcf
    File gvcf

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int mem_gb    = 4
  Int disk_size = ceil(size(vcf, "GB") * 2 + 20)

  String out_prefix = basename(vcf, ".vcf")

  command <<<
    set -euo pipefail

    bgzip --help
    tabix --help
    bcftools --version

    # Compress and index gVCF
    bgzip --stdout --threads ~{threads} ~{gvcf} > ~{out_prefix}.g.vcf.gz
    tabix --preset vcf ~{out_prefix}.g.vcf.gz

    # Filter VCF to remove uncalled sites
    bcftools view \
    ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
    --exclude-uncalled \
    --output-type z \
    --output-file ~{out_prefix}.vcf.gz \
    ~{vcf}

    bcftools index --tbi --force \
      ~{if threads > 1 then "--threads " + (threads - 1) else ""} \
      ~{out_prefix}.vcf.gz
  >>>

  output {
    File vcf_gz     = "~{out_prefix}.vcf.gz"
    File vcf_index  = "~{out_prefix}.vcf.gz.tbi"
    File gvcf_gz    = "~{out_prefix}.g.vcf.gz"
    File gvcf_index = "~{out_prefix}.g.vcf.gz.tbi"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87"
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
