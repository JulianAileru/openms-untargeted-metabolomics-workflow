"""Shared helpers for the OpenMS scripts run by Snakemake `script:` rules."""
import pyopenms as oms


def load_experiment(path):
    exp = oms.MSExperiment()
    oms.MzMLFile().load(str(path), exp)
    exp.sortSpectra(True)
    return exp


def load_feature_maps(paths):
    """Load several featureXML files; returns FeatureMaps in the same order as `paths`."""
    feature_maps = []
    for path in paths:
        feature_map = oms.FeatureMap()  # a fresh map per file, or every file piles into one
        oms.FeatureXMLFile().load(str(path), feature_map)
        feature_maps.append(feature_map)
    return feature_maps

def apply_params(algo, settings, skip=()):
    """Set OpenMS algorithm parameters from a config dict, ignoring keys in `skip`."""
    params = algo.getDefaults()
    for name, value in flatten_params(settings or {}).items():
        if name not in skip:
            params.setValue(name, value)
    algo.setParameters(params)


def flatten_params(settings, prefix=""):
    """Turn nested YAML sections into OpenMS's colon names: {"superimposer": {"max_shift": 1}} -> {"superimposer:max_shift": 1}."""
    flat = {}
    for key, value in settings.items():
        name = f"{prefix}{key}"
        if isinstance(value, dict):
            flat.update(flatten_params(value, prefix=f"{name}:"))
        else:
            flat[name] = value
    return flat


def detect_polarity(mzml_path):
    """Return "positive" or "negative" from the MS1 spectra's polarity; fail on mixed or missing.

    Reads only spectrum metadata (no peak data), so it is fast even for large files.
    """
    od = oms.OnDiscMSExperiment()
    if not od.openFile(str(mzml_path)):
        raise OSError(f"Could not open {mzml_path}")
    polarity = oms.IonSource.Polarity
    names = {polarity.POSITIVE: "positive", polarity.NEGATIVE: "negative"}
    meta = od.getMetaData()
    found = set()
    for i in range(meta.size()):
        spec = meta.getSpectrum(i)
        if spec.getMSLevel() == 1:
            found.add(names.get(spec.getInstrumentSettings().getPolarity(), "unknown"))
    if len(found) == 1 and found <= {"positive", "negative"}:
        return found.pop()
    raise ValueError(
        f"{mzml_path}: could not determine a single polarity from MS1 spectra (found {sorted(found)}); "
        "polarity-switching runs must be split by polarity first"
    )
