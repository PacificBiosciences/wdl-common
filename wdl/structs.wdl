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
}