import pyopenms as oms
import pandas as pd
import numpy as np 
import scipy 
import os 
import glob
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from snakemake.iocontainers import Snakemake
    snakemake: Snakemake

from openms_utils import load_experiment, apply_params


def mass_trace_detection(exp, settings):
    mtd = oms.MassTraceDetection()
    apply_params(mtd, settings, skip=("max_traces",))
    mass_traces = []
    mtd.run(exp, mass_traces, int(settings.get("max_traces", 0)))
    return mass_traces

def elution_peak_detection(mass_traces, settings):
    epd = oms.ElutionPeakDetection()
    apply_params(epd, settings)
    split_traces = []
    epd.detectPeaks(mass_traces, split_traces)
    return split_traces

def feature_finding_metabo(exp, split_traces, settings):
    ffm = oms.FeatureFindingMetabo()
    apply_params(ffm, settings)
    fm = oms.FeatureMap()
    chromatograms = []
    ffm.run(split_traces, fm, chromatograms)
    fm.setUniqueIds()
    fm.setPrimaryMSRunPath([str(snakemake.input.mzml).encode()])
    return fm

exp = load_experiment(snakemake.input.mzml)
traces = mass_trace_detection(exp, dict(snakemake.params.mtd))
split = elution_peak_detection(traces, dict(snakemake.params.epd or {}))
fm = feature_finding_metabo(exp, split, dict(snakemake.params.ffm or {}))
oms.FeatureXMLFile().store(str(snakemake.output.features), fm)

