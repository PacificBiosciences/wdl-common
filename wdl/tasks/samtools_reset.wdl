version 1.0

# Reset an aBAM back to uBAM _AND_ strip any kinetics tags

import "../structs.wdl"

task samtools_reset {
	input {
		File bam
		String remove_tags = "HP,PS,PC,SA,mg,rm,fi,fp,ri,rp"

		RuntimeAttributes runtime_attributes
	}

	String bam_basename = basename(bam, ".bam")
	Int threads = 4
	Int disk_size = ceil(size(bam, "GB") * 3.5 + 20)

	command <<<
		set -euo pipefail

		samtools --version

		samtools reset \
			--threads ~{threads - 1} \
			--remove-tag ~{remove_tags} \
			--reject-PG pbmm2 \
			-o {bam_basename}.reset.bam \
			~{bam}
	>>>

	output {
		File reads_fasta = "~{bam_basename}.fasta"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/samtools@sha256:cbe496e16773d4ad6f2eec4bd1b76ff142795d160f9dd418318f7162dcdaa685"
		cpu: threads
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
