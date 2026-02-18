# Spatial Omics Analysis with Scimap

![Python](https://img.shields.io/badge/Python-3.8%2B-blue)
![Scimap](https://img.shields.io/badge/Scimap-2.3.5-green)
![Data](https://img.shields.io/badge/Data-Demo_Set-orange)

A Python-based workflow for analyzing highly multiplexed imaging data, focusing on cell phenotyping and spatial analysis using the [Scimap](https://github.com/labsyspharm/scimap) framework.

## 📌 Overview

This repository contains scripts for **spatial proteomics analysis**. The current implementation demonstrates the workflow using **Scimap's standard demonstration dataset**, converting MCMICRO output data into `AnnData` objects and performing hierarchical cell phenotyping based on marker intensity.

**Current Capabilities:**
- Automated download and ingestion of Scimap demo data (MCMICRO format)
- Hierarchical gating strategies for cell classification
- Spatial visualization of phenotypes

**Future Modules:**
- [ ] Spatial neighborhood analysis
- [ ] Cellular interaction profiling (ligand-receptor/proximity)

---

## 🧬 Phenotyping Workflow Snippet

Below is a summary of the analysis pipeline. The workflow involves rescaling the data, defining a hierarchical gating strategy, and visualizing the resulting phenotypes [file:1].

```python
import scimap as sm
import pandas as pd

# 1. Preprocessing: Rescale data to standardize marker intensities
# This allows for consistent gating (e.g., setting a universal gate at 0.5)
adata = sm.pp.rescale(adata, gate=0.5)

# 2. Define the Hierarchical Gating Strategy
# Structure: Parent -> Child -> Marker -> Polarity (pos/neg)
phenotyping_gates = pd.DataFrame({
    'Unnamed: 0': ['all', 'all', 'all', 'Immune', 'Immune', 'Immune', 'Other myeloid cells'],
    'Unnamed: 1': ['ECAD+', 'Immune', 'SMA+', 'NK cells', 'Other myeloid cells', 'Treg', 'Dendritic cells'],
    'ECAD':       ['pos',  '',       '',     '',         '',                    '',     ''],
    'CD45':       ['',     'pos',    '',     '',         '',                    '',     ''],
    'SMA':        ['',     '',       'pos',  '',         '',                    '',     ''],
    'CD57':       ['',     '',       '',     'allpos',   '',                    '',     'allneg'],
    'CD16':       ['',     '',       '',     'allpos',   '',                    '',     'allneg'],
    'CD11B':      ['',     '',       '',     'neg',      'pos',                 '',     ''],
    'FOXP3':      ['',     '',       '',     '',         '',                    'pos',  '']
})

# 3. Run Phenotyping
adata = sm.tl.phenotype_cells(adata, phenotype=phenotyping_gates, label="phenotype")

# 4. Visualization
# Map the identified phenotypes back to spatial coordinates
sm.pl.spatial_scatter(adata, color='phenotype', s=2)

adata = sm.tl.phenotype_cells(
    adata, 
    phenotype=phenotyping_gates, 
    gate=0.5, 
    label="phenotype"  # Column name for results in adata.obs
)
