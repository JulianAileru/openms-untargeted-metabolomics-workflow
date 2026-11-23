from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from snakemake.iocontainers import Snakemake
    snakemake: Snakemake

import pyopenms as oms


def export_feature_table(consensus_map, output_path):
    """One row per consensus feature (index `id`), one intensity column per sample."""
    df = consensus_map.get_df()
    labels = {h.filename: h.label for h in consensus_map.getColumnHeaders().values()}
    df = df.rename(columns=labels).drop(columns=["sequence"], errors="ignore")
    df.to_csv(output_path, sep="\t")
    return df


consensus_map = oms.ConsensusMap()
oms.ConsensusXMLFile().load(str(snakemake.input.consensus), consensus_map)
export_feature_table(consensus_map, snakemake.output.matrix)
