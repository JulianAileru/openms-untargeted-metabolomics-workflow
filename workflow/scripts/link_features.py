from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from snakemake.iocontainers import Snakemake
    snakemake: Snakemake

import pyopenms as oms
from openms_utils import apply_params, load_feature_maps

# Map config names to classes; only the chosen one gets constructed.
ALGORITHMS = {
    "FeatureGroupingAlgorithmKD": oms.FeatureGroupingAlgorithmKD,
    "FeatureGroupingAlgorithmQT": oms.FeatureGroupingAlgorithmQT,
    "FeatureGroupingAlgorithmUnlabeled": oms.FeatureGroupingAlgorithmUnlabeled,
    "FeatureGroupingAlgorithmLabeled": oms.FeatureGroupingAlgorithmLabeled,
}


def make_grouper(settings):
    name = settings["algorithm"]
    if name not in ALGORITHMS:
        raise ValueError(f"Unknown feature_linking algorithm {name!r}; choose one of {sorted(ALGORITHMS)}")
    grouper = ALGORITHMS[name]()
    apply_params(grouper, settings.get(name))
    return grouper


def set_column_headers(consensus_map, feature_maps, samples):
    """Record which sample each ConsensusMap column (map index) came from."""
    headers = consensus_map.getColumnHeaders()
    for idx, (feature_map, sample) in enumerate(zip(feature_maps, samples)):
        header = headers.get(idx, oms.ColumnHeader())
        run_paths = []
        feature_map.getPrimaryMSRunPath(run_paths)
        header.filename = run_paths[0].decode() if run_paths else sample
        header.label = sample
        header.size = feature_map.size()
        header.unique_id = feature_map.getUniqueId()
        headers[idx] = header
    consensus_map.setColumnHeaders(headers)


def feature_linking(feature_maps, samples, settings):
    consensus_map = oms.ConsensusMap()
    set_column_headers(consensus_map, feature_maps, samples)
    make_grouper(settings).group(feature_maps, consensus_map)
    consensus_map.setUniqueIds()
    return consensus_map


samples = list(snakemake.params.samples)
settings = dict(snakemake.params.linking_settings)
feature_maps = load_feature_maps(snakemake.input.aligned_features)  # same order as samples
consensus_map = feature_linking(feature_maps, samples, settings)
oms.ConsensusXMLFile().store(str(snakemake.output.consensus), consensus_map)
