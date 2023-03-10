version 1.0

import "../structs.wdl"

task whatshap_stats {
	input {
		File phased_vcf
		File phased_vcf_index

		File reference_chromosome_lengths

		RuntimeAttributes runtime_attributes
	}

	String output_basename = basename(phased_vcf, ".vcf.gz")
	Int disk_size = ceil(size(phased_vcf, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		whatshap stats \
			--gtf ~{output_basename}.gtf \
			--tsv ~{output_basename}.tsv \
			--block-list ~{output_basename}.blocklist \
			--chr-lengths ~{reference_chromosome_lengths} \
			~{phased_vcf}
	>>>

	output {
		File gtf = "~{output_basename}.gtf"
		File tsv = "~{output_basename}.tsv"
		File blocklist = "~{output_basename}.blocklist"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/whatshap:1.4"
		cpu: 1
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
