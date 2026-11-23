from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from snakemake.iocontainers import Snakemake
    snakemake: Snakemake

import pyopenms as oms
from openms_utils import apply_params, detect_polarity


def detect_adducts(feature_map, settings):
    """Label one sample's features with adduct, charge and neutral mass.

    Returns:
        fm_out  - the same features, labeled (this is what goes on to linking)
        groups  - features grouped by compound, within this sample
        edges   - pairs of features linked by an adduct mass difference
    """
    mfd = oms.MetaboliteFeatureDeconvolution()
    apply_params(mfd, settings)

    fm_out = oms.FeatureMap()
    groups = oms.ConsensusMap()
    edges = oms.ConsensusMap()
    mfd.compute(feature_map, fm_out, groups, edges)
    return fm_out, groups, edges


def settings_for_polarity(settings, polarity):
    """Pick the adduct list for this polarity and set negative_mode and the charge range to match.

    Config charges are written as positive numbers; in negative mode OpenMS needs them negative
    (charge_min/max 1..3 becomes -3..-1), or it fails with "implicit number of default adduct is negative".
    """
    settings = dict(settings)
    adducts = settings.pop("potential_adducts")
    if polarity not in adducts:
        raise ValueError(f"No potential_adducts listed for {polarity!r} mode in config")
    settings["potential_adducts"] = list(adducts[polarity])
    settings["negative_mode"] = "true" if polarity == "negative" else "false"
    if polarity == "negative":
        low, high = abs(settings["charge_min"]), abs(settings["charge_max"])
        settings["charge_min"], settings["charge_max"] = -high, -low
    return settings


polarity = detect_polarity(snakemake.input.mzml)
print(f"{snakemake.wildcards.sample}: {polarity} ion mode")

feature_map = oms.FeatureMap()
oms.FeatureXMLFile().load(str(snakemake.input.features), feature_map)

fm_out, groups, edges = detect_adducts(feature_map, settings_for_polarity(snakemake.params.mfd, polarity))

oms.FeatureXMLFile().store(str(snakemake.output.features), fm_out)
oms.ConsensusXMLFile().store(str(snakemake.output.groups), groups)  # for checking the assignments
