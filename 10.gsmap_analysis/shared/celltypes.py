"""Cell type grouping helpers used across gsMap analysis scripts."""

from __future__ import annotations

import pandas as pd


def assign_major_celltype(celltype: str) -> str:
    """Map fine fetal brain cell type labels to stable major groups."""
    ct = str(celltype)
    if ct.startswith("Glu_"):
        return "Glu"
    if ct.startswith("Gaba_"):
        return "Gaba"
    if ct.startswith("NPC_"):
        if "_IPC" in ct:
            return "NPC_IPC"
        if "_oRG" in ct or "_vRG" in ct or "_rg" in ct.lower():
            return "NPC_RG"
        if "_APC" in ct:
            return "NPC_APC"
        if "_precursor" in ct:
            return "NPC_precursor"
        if "_GCP" in ct:
            return "NPC_GCP"
        if "Cajal" in ct:
            return "NPC_Cajal"
        if "Crest" in ct:
            return "NPC_Crest"
        return "NPC"
    if ct == "Dopaminergic neurons":
        return "Dopaminergic"
    if ct == "choroid plexus epithelial":
        return "Choroid"
    if ct == "Immune cell":
        return "Immune"
    if ct == "OPC":
        return "OPC"
    if ct == "VLMC":
        return "VLMC"
    if ct == "unknown":
        return "unknown"
    return "other"


def add_major_celltype_columns(
    obs: pd.DataFrame,
    celltype_col: str = "celltype",
    region_col: str = "regions",
) -> pd.DataFrame:
    """Return obs with major cell type and region-celltype pair annotations."""
    if celltype_col not in obs.columns:
        raise KeyError(f"Missing required obs column: {celltype_col}")
    if region_col not in obs.columns:
        raise KeyError(f"Missing required obs column: {region_col}")

    obs = obs.copy()
    obs["major_celltype"] = obs[celltype_col].astype(str).map(assign_major_celltype)

    pair = obs[region_col].astype(str) + "|" + obs["major_celltype"].astype(str)
    obs["subregion_major_celltype"] = pair
    obs["subregion_major_celltype_pair"] = pair
    return obs

