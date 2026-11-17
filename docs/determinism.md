# How deterministic are variant calling results?

Most tools in the workflow produce deterministic results across most conditions. [DeepVariant](deepvariant.md) is the notable exception: it is not expected to produce identical variant calls when run on different compute hardware. This is most often observed when [comparing variant calls run on CPUs versus GPUs](https://docs.nvidia.com/clara/parabricks/tool-reference/tools/deepvariant#source-of-mismatches), but it also occurs between different CPUs. For example, if one CPU has avx512 available and the other only has avx2, they will produce different results.

## Why does this occur?

DeepVariant's optimizations and the nature of deep learning models make it sensitive to floating-point handling. Not all compute architectures handle decimal calculation in exactly the same way, and these tiny differences produce minor differences in reported variant calls. As such, the SNVs and indels produced by DeepVariant are only expected to be deterministic *when run on the same compute hardware*.

DeepVariant is one of the early steps in the workflow, so this variation can propagate to downstream steps such as phasing and star-allele calling. Independent runs of the workflow on different hardware may therefore show very minor differences in their outputs. On a compute cluster with heterogeneous nodes, tasks such as DeepVariant can land on nodes with different architectures by chance unless the node type is specified.

## How impactful is this variation?

In general, differences between runs on different hardware are very minor, and the variant calls that change have quality values near DeepVariant's calling thresholds. The workflow does not apply quality filters on top of DeepVariant's defaults; applying more stringent filters to the callset will likely remove the variants that differ.
