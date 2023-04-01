version 1.0

# Call SVs using pbsv

import "../structs.wdl"

task pbsv_call {
	input {
		String sample_id
		Array[File] svsigs

		File reference
		File reference_index
		String reference_name

		RuntimeAttributes runtime_attributes
	}

	Int threads = 8
	Int disk_size = ceil((size(svsigs[0], "GB") * length(svsigs) + size(reference, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		pbsv call \
			--hifi \
			--min-sv-length 20 \
			--log-level INFO \
			--num-threads ~{threads} \
			~{reference} \
			~{sep=' ' svsigs} \
			~{sample_id}.~{reference_name}.pbsv.vcf
	>>>

	output {
		File pbsv_vcf = "~{sample_id}.~{reference_name}.pbsv.vcf"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/pbsv:2.9.0"
		cpu: threads
		memory: "64 GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		slurm_partition: runtime_attributes.slurm_partition_default
		zones: runtime_attributes.zones
	}
}
