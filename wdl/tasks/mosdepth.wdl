version 1.0

# Calculate summary stats using mosdepth

import "../structs.wdl"

task mosdepth {
	input {
		File aligned_bam
		File aligned_bam_index

		RuntimeAttributes runtime_attributes
	}

	String prefix = basename(aligned_bam, ".bam")
	Int threads = 4
	Int disk_size = ceil(size(aligned_bam, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		mosdepth \
			--threads ~{threads - 1} \
			--by 500 \
			--no-per-base \
			--use-median \
			~{prefix} \
			~{aligned_bam}
	>>>

	output {
		File summary = "~{prefix}.mosdepth.summary.txt"
		File region_bed = "~{prefix}.regions.bed.gz"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/mosdepth@sha256:3ebb896013e205db072c55a9fdcd322773f4a4bcdc7bedecc80b220b88e9b750"
		cpu: threads
		memory: "4 GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " LOCAL"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
	}
}
