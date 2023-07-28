# Configuring the Amazon Genomics CLI

The Amazon Genomics CLI (`agc`) allows users to orchestrate workflow execution using AWS Batch. See the [Workbench documentation](https://docs.dnastack.com/docs/cromwell-on-aws-amazon-genomics-cli) for information on installing and using the `agc` to configure and run workflows. The following section provides additional information on deploying a project using the `agc`.

## Deploying a context with `agc`

Once you have installed and authenticated with the `agc`, you can deploy a context using an agc project YAML file. This file must be named `agc-project.yaml`.

An [example agc-project.yaml file](agc-project.template.yaml) that has the workflow, reference data source, and both on-demand and spot contexts configured using Cromwell as the engine is provided here. This will create an agc project named `humanwgsAGC`, with either (or both) a `spotContext` or an `onDemandContext`. The `spotContext` will allow you to run worklfows using [AWS spot instances](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/using-spot-instances.html), which can result in substantial cost savings relative to using [on-demand instances](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-on-demand-instances.html).

Note that deploying a context **will incur costs** even if you are not actively running workflows; ensure that [contexts that are not in use are destroyed](https://aws.github.io/amazon-genomics-cli/docs/reference/agc_context_destroy/) to avoid incurring ongoing costs.

To deploy the agc project using the template file, first copy the template file to a file named `agc-project.yaml` (`cp agc-project.template.yaml agc-project.yaml`).

In the `data` section of the `agc-project.yaml` file, add any additional s3 buckets that the workflow will require access to, for example the bucket containing sample input data. Make sure that you do not remove the section granting access to the s3://dnastack-resources bucket; this is where [reference datasets are hosted](#reference-data-hosted-in-aws).

```
data:
  - location: s3://dnastack-resources
    readOnly: true
  - location: s3://<sample_data_bucket_name>
    readOnly: true
```

Then from the directory containing the `agc-project.yaml` file, run:

```bash
agc context deploy --context ${context}
```

Where `${context}` is either `spotContext` or `onDemandContext`.

If you want both spot and on-demand contexts, all contexts can be deployed at once by running:

```
agc context deploy --all
```

Note that the `miniwdl` engine run via AWS is currently not supported for this workflow.

# Checking and requesting quota in AWS

See [resources requirements](../../README.md#resource-requirements) for information on the minimum requirements for running the workflow. Typically in a new AWS environment, additional vCPU quota will be required.

## Checking current quota

1. Navigate to [the AWS console](https://console.aws.amazon.com/).
2. In the top right corner, select the region where your `agc` deployment is located.
3. Navigate to EC2.
4. In the menu on the left, select 'Limits'.
5. Filter the limits by searching for "Standard". The current limit field indicates the number of vCPUs that you currently have access to.
- Spot instance limit: `All Standard (A, C, D, H, I, M, R, T, Z) Spot Instance Requests`
- On-demand instance limit: `Running On-Demand All Standard (A, C, D, H, I, M, R, T, Z) instances`

If the number of vCPUs in the context you plan to run the workflow in is less than the limites specified in [the resources requirements](../../README.md#resource-requirements) section, you will need to request additional quota before you can run the workflow.

## Requesting additional quota

5. Continuing from the steps outlined in [checking the current quota](#checking-current-quota), select the service you want to request an increase for.
6. In the top right corner, select 'Request limit increase'.
7. Fill out the appropriate fields in the request form, ensuring that the region you select is the region where you have deployed your `agc` and where your data is located. 256 vCPUs are recommended for running trio data.

Low quota increase requests are typically fulfilled within a 1-2 hours.

# Configuring and running the workflow

## Filling out workflow inputs

Fill out any information missing in [the inputs file](inputs.aws.json). Ensure that all data files used by the workflow are at locations that have been configured in the agc-project.yaml file; see the [granting access to other data files](#granting-access-to-other-data-files) for more information.

See [the inputs section of the main README](../../README.md#workflow-inputs) for more information on the structure of the inputs.json file.

Note that you only need to fill out the queueArn corresponding to the context you are submitting the workflow to (spot or on-demand).

### Determining available zones

To determine available zones in AWS, look for the `ZoneName` attribute output by the following command:

```bash
aws ec2 describe-availability-zones --region <region>
```

For example, the zones in region us-east-2 are `"us-east-2a us-east-2b us-east-2c"`.

### Determining the AWS batch queue ARN

**Note that if you are using a `miniwdl` engine, you can skip these steps; workflows run via miniwdl will run exclusively in the job queue to which they are submitted.**

1. Visit [the AWS console](https://console.aws.amazon.com/).
2. Navigate to the Batch service.
3. In the lefthand sidebar, select "Compute environments". Note the name of the compute environment with the provisioning model SPOT (if you have deployed a context using spot instances) and the name of the compute environment with provisioning model "EC2" (if you have deployed a context that does not use spot instances).
4. In the lefthand sidebar, select "Job queues".
5. Clicking into an individual queue will show information about the compute environment ("Compute environment order"). Identify the job queue with the Compute environment name that matches the name you identified for the SPOT compute environment; copy the Amazon Resource Name (ARN) for this job queue. This is the value that should be used for the `aws_spot_queue_arn`. Repeat this process to find the ARN for the `aws_on_demand_queue_arn`.

- If `preemptible = true`, only the `aws_spot_queue_arn` is required.
- If `preemptible = false`, only the `aws_on_demand_queue_arn` is required.

## Running the workflow

### Running via `agc`

From the directory where your `agc-project.yaml` is located, run:

`agc workflow run humanwgs --context <context> --inputsFile <input_file_path.json>`

The running workflow can be monitored via [`agc workflow` commands](https://aws.github.io/amazon-genomics-cli/docs/reference/agc_workflow/), or via the AWS console.

# Reference data hosted in AWS

AWS reference data is hosted in the `us-west-2` region in the bucket `s3://dnastack-resources`.

To use AWS reference data, add the following line to the data section of your [`agc-project.yaml`](https://aws.github.io/amazon-genomics-cli/docs/concepts/projects/):

```yaml
data:
  - location: s3://dnastack-resources
    readOnly: true
```
The [AWS input file template](inputs.aws.json) has paths to the reference files in s3 prefilled. The template [agc-project.template.yaml file](agc-project.template.yaml) has this section filled out already.

### Granting access to other data files

S3 buckets outside of the reference files can be accessed by adding additional data blocks to the agc-project.yaml file. See the [agc documentation](https://aws.github.io/amazon-genomics-cli/docs/concepts/data/) for more details on adding additional data sources. All inputs referenced in the inputs.json file will need to be at locations that have been configured in the agc-project.yaml.
