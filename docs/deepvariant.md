# DeepVariant subworkflow

```mermaid
flowchart TD
  aBAM[/"HiFi aBAM"/] --> make_examples["DeepVariant make_examples"]
  make_examples --> gpu{"gpu?"}
  gpu -- yes --> call_variants_gpu["DeepVariant call_variants_gpu"]
  gpu -- no --> call_variants_cpu["DeepVariant call_variants_cpu"]
  call_variants_gpu --> postprocess_variants["DeepVariant postprocess_variants"]
  call_variants_cpu --> postprocess_variants
  postprocess_variants --> vcf[/"small variant VCF"/]
  postprocess_variants --> gvcf[/"small variant gVCF"/]
```

This subworkflow runs the three steps of DeepVariant individually in order to make best use of resources.  If a GPU is available and `gpu==true`, the `call_variants` step will run on 1 GPU and 8 cpu threads, otherwise it will run on 64 CPU threads.  The `make_examples` and `postprocess_variants` steps will always run on the CPU.
