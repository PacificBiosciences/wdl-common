version 1.0

# Concatenate and sort VCFs

import "../structs.wdl"

task concat_vcf {
	input {
		Array[File] vcfs

		String output_vcf_name

		RuntimeAttributes runtime_attributes
	}

	Int threads = 4
	Int disk_size = ceil(size(vcfs, "GB") * 2 + 20)
	
	command <<<
		set -euo pipefail

		bcftools --version

		bcftools concat \
			--threads ~{threads - 1} \
			--output-type b \
			~{sep=' ' vcfs} \
		| bcftools sort \
			--max-mem 4G \
			--output-type z \
			--output ~{output_vcf_name} \
			/dev/stdin
	>>>

	output {
		File concatenated_vcf = "~{output_vcf_name}"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/bcftools@sha256:36d91d5710397b6d836ff87dd2a924cd02fdf2ea73607f303a8544fbac2e691f"
		cpu: threads
		memory: "8 GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
	}
}
