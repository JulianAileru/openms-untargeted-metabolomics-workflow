from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from snakemake.iocontainers import Snakemake
    snakemake: Snakemake

import pyopenms as oms
from openms_utils import apply_params, flatten_params, load_feature_maps


def fit_rt_model(trafo, model_settings):
    """Fit the affine RT correction (rt' = slope * rt + intercept) to the feature pairs pose clustering found."""
    model_type = model_settings.get("type", "linear")
    if model_type != "linear":
        raise ValueError(f"Only the 'linear' (affine) RT model is supported with PoseClustering, got {model_type!r}")
    params = oms.Param()
    for name, value in flatten_params(model_settings.get("linear") or {}).items():
        params.setValue(name, value)
    trafo.fitModel("linear", params)


def alignment(feature_maps, ref_idx, aln_settings, trafo_settings):
    """Align every map's RTs onto feature_maps[ref_idx], in place."""
    aligner = oms.MapAlignmentAlgorithmPoseClustering()
    apply_params(aligner, aln_settings, skip=("reference_sample",))
    aligner.setReference(feature_maps[ref_idx])

    transformer = oms.MapAlignmentTransformer()
    store_original_rt = bool(trafo_settings.get("store_original_rt", True))
    for idx, feature_map in enumerate(feature_maps):
        if idx == ref_idx:
            continue  # the reference defines the RT scale; leave it as is
        trafo = oms.TransformationDescription()
        aligner.align(feature_map, trafo)  # fills trafo with paired RTs + a linear fit; returns None
        fit_rt_model(trafo, trafo_settings.get("model") or {})
        transformer.transformRetentionTimes(feature_map, trafo, store_original_rt)
    return feature_maps


samples = list(snakemake.params.samples)
settings = dict(snakemake.params.aln)
trafo_settings = dict(snakemake.params.trafo or {})
reference = settings["reference_sample"]
if reference not in samples:
    raise ValueError(f"reference_sample {reference!r} is not in the samplesheet: {samples}")

feature_maps = load_feature_maps(snakemake.input.featureXML)  # same order as samples
alignment(feature_maps, samples.index(reference), settings, trafo_settings)

for feature_map, out_path in zip(feature_maps, snakemake.output.aligned):
    oms.FeatureXMLFile().store(str(out_path), feature_map)
