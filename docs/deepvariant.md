# DeepVariant subworkflow

```mermaid

flowchart TD
  a1[/"HiFi aBAM"/] --> a2["DeepVariant make_examples"]
  a2 --> a7{"gpu?"}
  a7 -- yes --> a8["DeepVariant call_variants_gpu"]
  a7 -- no --> a3["DeepVariant call_variants_cpu"]
  a3 --> a4["DeepVariant postprocess_variants"]
  a8 --> a4
  a4 --> a5[/"small variant VCF"/]
  a4 --> a6[/"small variant gVCF"/]
```

This subworkflow runs the three steps of DeepVariant individually in order to make best use of resources.  If a GPU is available and `gpu==true`, the `call_variants` step will run on 1 GPU and 8 cpu threads, otherwise it will run on 64 CPU threads.  The `make_examples` and `postprocess_variants` steps will always run on the CPU.
