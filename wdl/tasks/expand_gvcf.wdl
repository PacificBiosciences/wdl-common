version 1.0

import "../structs.wdl"

task expand_gvcf {
	input {
		File gvcf
		File gvcf_index

		File reference

		File expansion_targets

		RuntimeAttributes runtime_attributes
	}

	String gvcf_basename = basename(gvcf, ".g.vcf.gz")
	Int disk_size = ceil((size(gvcf, "GB") + size(reference, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		bcftools convert \
			--gvcf2vcf \
			-f ~{reference} \
			-R ~{expansion_targets} \
			-Oz \
			-o ~{gvcf_basename}.expanded.g.vcf.gz \
			~{gvcf}

		tabix ~{gvcf_basename}.expanded.g.vcf.gz
	>>>

	output {
		File expanded_gvcf = "~{gvcf_basename}.expanded.g.vcf.gz"
		File expanded_gvcf_index = "~{gvcf_basename}.expanded.g.vcf.gz.tbi"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/bcftools:1.14"
		cpu: 2
		memory: "4 GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
	}
}
