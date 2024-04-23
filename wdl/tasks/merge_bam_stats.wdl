version 1.0

import "../structs.wdl"

task merge_bam_stats {
  meta {
    description: "Merge BAM stats files and create read length and quality histograms"
  }

  parameter_meta {
    sample_id: {
      name: "Sample ID"
    }
    bam_stats: {
      name: "BAM Stats"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    read_length_histogram: {
      name: "Read Length Summary"
    }
    read_quality_histogram: {
      name: "Read Quality Summary"
    }
    read_length_and_quality: {
      name: "Read Length and Quality"
    }
    stat_num_reads: {
      name: "Number of Reads"
    }
    stat_read_length_mean: {
      name: "Mean Read Length"
    }
    stat_read_length_median: {
      name: "Median Read Length"
    }
    stat_read_quality_mean: {
      name: "Mean Quality"
    }
    stat_read_quality_median: {
      name: "Median Quality"
    }
  }

  input {
    String sample_id
    Array[File] bam_stats

    RuntimeAttributes runtime_attributes
  }

  Int threads = 2
  Int disk_size = 20

  command <<<
    zcat ~{sep = " " bam_stats} > ~{sample_id}.read_length_and_quality.tsv
    # TODO: add column for movie source, adjust indices
    awk '{{ b=int($2/1000); b=(b>39?39:b); print 1000*b "\t" $2; }}' \
      ~{sample_id}.read_length_and_quality.tsv \
      | sort -k1,1g \
      | datamash -g 1 count 1 sum 2 \
      | awk 'BEGIN {{ for(i=0;i<=39;i++) {{ print 1000*i"\t0\t0"; }} }} {{ print; }}' \
      | sort -k1,1g \
      | datamash -g 1 sum 2 sum 3 \
    > ~{sample_id}.read_length_histogram.tsv

    awk '{{ print ($3>50?50:$3) "\t" $2; }}' \
      ~{sample_id}.read_length_and_quality.tsv \
      | sort -k1,1g \
      | datamash -g 1 count 1 sum 2 \
      | awk 'BEGIN {{ for(i=0;i<=60;i++) {{ print i"\t0\t0"; }} }} {{ print; }}' \
      | sort -k1,1g \
      | datamash -g 1 sum 2 sum 3 \
    > ~{sample_id}.read_quality_histogram.tsv

    datamash count 1 mean 2 median 2 mean 3 median 3 \
      < ~{sample_id}.read_length_and_quality.tsv \
      > stats.txt

    cut -f1 stats.txt > num_reads.txt || echo "0" > num_reads.txt
    cut -f2 stats.txt > read_length_mean.txt || echo "0" > read_length_mean.txt
    cut -f3 stats.txt > read_length_median.txt || echo "0" > read_length_median.txt
    cut -f4 stats.txt > read_quality_mean.txt || echo "0" > read_quality_mean.txt
    cut -f5 stats.txt > read_quality_median.txt || echo "0" > read_quality_median.txt

    gzip ~{sample_id}.read_length_and_quality.tsv
  >>>

  output {
    File  read_length_and_quality  = "~{sample_id}.read_length_and_quality.tsv.gz"
    File  read_length_histogram    = "~{sample_id}.read_length_histogram.tsv"
    File  read_quality_histogram   = "~{sample_id}.read_quality_histogram.tsv"
    Int   stat_num_reads           = read_int("num_reads.txt")
    Float stat_read_length_mean    = read_float("read_length_mean.txt")
    Float stat_read_length_median  = read_float("read_length_median.txt")
    Float stat_read_quality_mean   = read_float("read_quality_mean.txt")
    Float stat_read_quality_median = read_float("read_quality_median.txt")
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/pbmm2@sha256:ed9dcb4db98c81967fff15f50fca89c8495b1f270eee00e9bec92f46d14d7e2f"
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