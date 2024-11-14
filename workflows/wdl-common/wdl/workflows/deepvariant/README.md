# Custom models

Custom DeepVariant models must match the version of DeepVariant used for calling.

To prepare a `custom_deepvariant_model_tar` input:

```shell
# cd into custom model dir
tar zxvf my_custom_deepvariant_model_tar.tgz example_info.json fingerprint.pb saved_model.pb variables
```

Provide the path to this tarball as an input.
