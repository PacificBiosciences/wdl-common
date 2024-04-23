version 1.0

# Set runtime attributes across environments depending on the backend in use

import "../../structs.wdl"

workflow backend_configuration {
  meta {
    description: "Set runtime attributes across environments depending on the backend in use"
  }

  parameter_meta {
    backend: {
      help: "Backend where the workflow will be executed",
      choices: ["GCP", "Azure", "AWS-AGC", "AWS-HealthOmics", "HPC"]
    }
    zones: {
      help: "Zones where compute will take place; required if backend is set to 'AWS' or 'GCP'"
    }
    aws_spot_queue_arn: {
      help: "Queue ARN for the spot batch queue; required if backend is set to 'AWS'"
    }
    aws_on_demand_queue_arn: {
      help: "Queue ARN for the on demand batch queue; required if backend is set to 'AWS'"
    }
    gpuType: {
      help: "Type of GPU/Accelerator to use"
    }
    container_registry: {
      help: "Container registry to use"
    }
  }

  input {
    String backend
    String? zones
    String? aws_spot_queue_arn
    String? aws_on_demand_queue_arn
    String? gpuType
    String container_registry = "quay.io/pacbio"
  }

  if (backend == "GCP") {
    # zones must be defined

    # preemptible_tries applies to failures due to preemption only
    # max_retries applies to failures due to a nonzero rc
    # queue_arn is not used in GCP
    # gpuCount and gpuType are optional
    # gpuType: ["nvidia-tesla-k80", "nvidia-tesla-p100", "nvidia-tesla-v100", "nvidia-tesla-p4", "nvidia-tesla-t4",
    #           "nvidia-tesla-a100", "nvidia-a100-80gb", "nvidia-l4", "nvidia-h100-80gb"]  # TODO: figure out which are compatible
    RuntimeAttributes gcp_spot_runtime_attributes = {
      "backend": "GCP",
      "preemptible_tries": 3,
      "max_retries": 0,
      "zones": select_first([zones]),
      "queue_arn": "",
      "gpuType": select_first([gpuType, ""]),
      "container_registry": container_registry
    }

    RuntimeAttributes gcp_on_demand_runtime_attributes = {
      "backend": "GCP",
      "preemptible_tries": 0,
      "max_retries": 0,
      "zones": select_first([zones]),
      "queue_arn": "",
      "gpuType": select_first([gpuType, ""]),
      "container_registry": container_registry
    }
  }

  if (backend == "Azure") {
    # Requires Cromwell on Azure v3.2+
    # preemptible_tries >= 1 will be converted to `true`; 0 will be converted to `false`
    # max_retries applies to failures due to preemption or to a nonzero rc
    # zones, queue_arn not used in Azure
    # GPUs are not available in Azure
    RuntimeAttributes azure_spot_runtime_attributes = {
      "backend": "Azure",
      "preemptible_tries": 3,
      "max_retries": 3,
      "zones": "",
      "queue_arn": "",
      "gpuType": "",
      "container_registry": container_registry
    }

    RuntimeAttributes azure_on_demand_runtime_attributes = {
      "backend": "Azure",
      "preemptible_tries": 0,
      "max_retries": 0,
      "zones": "",
      "queue_arn": "",
      "gpuType": "",
      "container_registry": container_registry
    }
  }

  if (backend == "AWS-AGC") {
    # zones must be defined
    # aws_spot_queue_arn must be defined if preemptible is set to true and engine is not miniwdl
    # aws_on_demand_queue_arn must be defined if preemptible is set to false and engine is not miniwdl
    # Using miniwdl engine, the queue ARN of the context the workflow has been submitted to will be used;
    #   the queue_arn runtime attribute will be ignored

    # max_retries applies to failures due to preemption or to a nonzero rc
    # preemptible is not used in AWS
    # gpuCount and gpuType are optional and map to acceleratorCount and acceleratorType 
    # acceleratorType: ["nvidia-tesla-a10g", "nvidia-tesla-t4", "nvidia-tesla-t4-a10g"]
    RuntimeAttributes aws_agc_spot_runtime_attributes = {
      "backend": "AWS-AGC",
      "preemptible_tries": 3,
      "max_retries": 3,
      "zones": select_first([zones]),
      "queue_arn": select_first([aws_spot_queue_arn, ""]),
      "gpuType": select_first([gpuType, ""]),
      "container_registry": container_registry
    }

    RuntimeAttributes aws_agc_on_demand_runtime_attributes = {
      "backend": "AWS-AGC",
      "preemptible_tries": 0,
      "max_retries": 0,
      "zones": select_first([zones]),
      "queue_arn": select_first([aws_on_demand_queue_arn, ""]),
      "gpuType": select_first([gpuType, ""]),
      "container_registry": container_registry
    }
  }

  if (backend == "AWS-HealthOmics") {
    # No distinction between preemptible and on-demand in AWS-HealthOmics configuration
    # zones must be defined
    # aws_spot_queue_arn must be defined if preemptible is set to true and engine is not miniwdl
    # aws_on_demand_queue_arn must be defined if preemptible is set to false and engine is not miniwdl
    # Using miniwdl engine, the queue ARN of the context the workflow has been submitted to will be used;
    #   the queue_arn runtime attribute will be ignored

    # max_retries applies to failures due to preemption or to a nonzero rc
    # preemptible is not used in AWS

    # gpuCount and gpuType are optional and map to acceleratorCount and acceleratorType 
    # acceleratorType: ["nvidia-tesla-a10g", "nvidia-tesla-t4", "nvidia-tesla-t4-a10g"]

    # AWS HealthOmics must use containers hosted on ECR and cannot use our Quay registry.
    # TODO: are zones/arns needed for AWS HealthOmics?
    RuntimeAttributes aws_healthomics_on_demand_runtime_attributes = {
      "backend": "AWS-HealthOmics",
      "preemptible_tries": 0,
      "max_retries": 0,
      "zones": select_first([zones]),
      "queue_arn": select_first([aws_on_demand_queue_arn, ""]),
      "gpuType": select_first([gpuType, ""]),
      "container_registry": container_registry
    }
  }

  if (backend == "HPC") {
    # No distinction between preemptible and on-demand in HPC configuration
    # default_hpc_partition is provided to specify the partition to use in HPC
    RuntimeAttributes hpc_runtime_attributes = {
      "backend": "HPC",
      "preemptible_tries": 0,
      "max_retries": 3,
      "zones": "",
      "queue_arn": "",
      "gpuType": select_first([gpuType, ""]),
      "container_registry": container_registry
    }
  }

  output {
    RuntimeAttributes spot_runtime_attributes = select_first([
      gcp_spot_runtime_attributes,
      azure_spot_runtime_attributes,
      aws_agc_spot_runtime_attributes,
      aws_healthomics_on_demand_runtime_attributes,
      hpc_runtime_attributes
    ])
    RuntimeAttributes on_demand_runtime_attributes = select_first([
      gcp_on_demand_runtime_attributes,
      azure_on_demand_runtime_attributes,
      aws_agc_on_demand_runtime_attributes,
      aws_healthomics_on_demand_runtime_attributes,
      hpc_runtime_attributes
    ])
  }
}
