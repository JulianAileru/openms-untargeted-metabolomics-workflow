"""Pair each consensus feature with one MS2 spectrum (most intense precursor across samples) and write an MSP."""
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from snakemake.iocontainers import Snakemake
    snakemake: Snakemake

import numpy as np
import pyopenms as oms
from openms_utils import detect_polarity


def load_ms2_index(mzml_path):
    """Precursor m/z, RT and precursor intensity for every MS2 scan, sorted by precursor m/z.

    Returns the open OnDiscMSExperiment (to fetch peaks later) and an array with columns
    [spectrum_index, precursor_mz, rt, precursor_intensity].
    """
    od = oms.OnDiscMSExperiment()
    if not od.openFile(str(mzml_path)):
        raise OSError(f"Could not open {mzml_path}")
    meta = od.getMetaData()
    rows = []
    for i in range(meta.size()):
        spec = meta.getSpectrum(i)
        if spec.getMSLevel() == 2 and spec.getPrecursors():
            prec = spec.getPrecursors()[0]
            rows.append((i, prec.getMZ(), spec.getRT(), prec.getIntensity()))
    ms2 = np.array(rows, dtype=float).reshape(-1, 4)
    return od, ms2[np.argsort(ms2[:, 1])]


def load_original_rts(featurexml_path):
    """Feature unique id -> RT in the raw data (before alignment)."""
    fm = oms.FeatureMap()
    oms.FeatureXMLFile().load(str(featurexml_path), fm)
    # The reference sample is never shifted, so it has no original_RT; its RT is already original.
    return {
        f.getUniqueId(): f.getMetaValue("original_RT") if f.metaValueExists("original_RT") else f.getRT()
        for f in fm
    }


def best_ms2_for_feature(ms2, mz, rt, ppm, rt_tol):
    """Row of the MS2 scan with the most intense precursor within ppm and rt_tol, or None."""
    tol = mz * ppm * 1e-6
    lo, hi = np.searchsorted(ms2[:, 1], [mz - tol, mz + tol])
    if lo == hi:
        return None
    window = ms2[lo:hi]
    window = window[np.abs(window[:, 2] - rt) <= rt_tol]
    if len(window) == 0:
        return None
    return window[np.argmax(window[:, 3])]


def select_ms2(consensus_map, ms2_by_map, rts_by_map, ppm, rt_tol):
    """For each consensus feature, the (map_index, MS2 row) with the highest precursor intensity across samples."""
    best = {}
    for cf in consensus_map:
        for handle in cf.getFeatureList():
            m = handle.getMapIndex()
            rt = rts_by_map[m].get(handle.getUniqueId())
            if rt is None:
                continue
            row = best_ms2_for_feature(ms2_by_map[m], handle.getMZ(), rt, ppm, rt_tol)
            if row is not None and (cf.getUniqueId() not in best or row[3] > best[cf.getUniqueId()][1][3]):
                best[cf.getUniqueId()] = (m, row)
    return best


def export_msp(consensus_map, best, experiments, samples, polarity, output_path):
    by_id = {cf.getUniqueId(): cf for cf in consensus_map}
    library = oms.MSExperiment()
    for cf_id, (m, row) in best.items():
        cf = by_id[cf_id]
        spec = experiments[m].getSpectrum(int(row[0]))
        native_id = spec.getNativeID()
        spec.clearMetaInfo()  # drop instrument fields; write only standard MSP fields
        spec.setName(str(cf_id))  # matches the `id` column of the feature matrix
        spec.setRT(cf.getRT())    # MSPGenericFile writes this as "Retention Time" (seconds, aligned)
        spec.setMetaValue("PRECURSORMZ", float(cf.getMZ()))
        spec.setMetaValue("RETENTIONTIME", round(cf.getRT() / 60.0, 4))  # minutes, aligned RT
        spec.setMetaValue("IONMODE", "Positive" if polarity == "positive" else "Negative")
        spec.setMetaValue("Comment", f"sample={samples[m]}; scan={native_id}; precursor_intensity={row[3]:.4g}")
        library.addSpectrum(spec)
    oms.MSPGenericFile().store(str(output_path), library)
    return library.size()


samples = list(snakemake.params.samples)
settings = dict(snakemake.params.ms2)

consensus_map = oms.ConsensusMap()
oms.ConsensusXMLFile().load(str(snakemake.input.consensus), consensus_map)

polarity = detect_polarity(snakemake.input.mzml[0])
experiments, ms2_by_map = zip(*(load_ms2_index(p) for p in snakemake.input.mzml))  # same order as samples
rts_by_map = [load_original_rts(p) for p in snakemake.input.features]
best = select_ms2(consensus_map, ms2_by_map, rts_by_map,
                  float(settings["mz_tolerance_ppm"]), float(settings["rt_tolerance"]))
n = export_msp(consensus_map, best, experiments, samples, polarity, snakemake.output.msp)
print(f"MSP: {n} of {consensus_map.size()} consensus features have an MS2 spectrum")
