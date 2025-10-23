version 1.0

struct RuntimeAttributes {
  String backend

  # The number of times to retry a task that fails due to preemption
  Int preemptible_tries
  # The number of times to retry a task that fails due a to nonzero return code
  Int max_retries

  String zones

  String gpuType

  String container_registry

  # AWS ECR registries have the format REGISTRY/NAMESPACE/CONTAINER
  # and if the namespace is not specified, HealthOmics will have permissions issues
  String? container_namespace  # Namespace within AWS ECR for HealthOmics

  # Memory override parameters
  Int? hiphase_override_mem_gb
  Int? merge_bam_stats_override_mem_gb
  Int? pbmm2_align_wgs_override_mem_gb
  Int? pbstarphase_diplotype_override_mem_gb
}
