from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from snakemake.iocontainers import Snakemake
    snakemake: Snakemake

import pyopenms as oms

fm = oms.FeatureMap()
oms.FeatureXMLFile().load(str(snakemake.input.features), fm)

# Only single-value meta values: masstrace_intensity, masstrace_centroid_rt/mz and
# isotope_distances hold one entry per isotope trace, which get_df can't put in a column.
SCALAR_META = ["label", "FWHM", "max_height", "num_of_masstraces", "legal_isotope_pattern"]

df = fm.get_df(meta_values=[m.encode() for m in SCALAR_META], export_peptide_identifications=False)
df.insert(0, "sample", snakemake.wildcards.sample)
df.to_csv(snakemake.output.table, sep="\t")
