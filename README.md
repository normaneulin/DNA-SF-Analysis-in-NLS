DNA Shape Feature Analysis in Nucleosomal &amp; Linker Sequences
# DNA-SF-Analysis-in-NLS

## Overview
This project analyzes DNA shape features of nucleosome positioning datasets using Monte Carlo simulation via the DNAshapeR package in R. The analysis is based on datasets from Guo et al., organized and downloaded from Zhuo et al. The extracted shape features provide insights into nucleosomal and linker DNA regions.

## Datasets
The datasets are grouped into three categories:

- **DatasetNup1**
- **DatasetNup2**
- **DatasetNup3**

Each dataset contains DNA sequences relevant to nucleosome positioning, processed and organized for downstream analysis.

## Feature Extraction
For each file, thirteen (13) DNA shape features plus electrostatic positivity (EP) were extracted using Monte Carlo simulation via the DNAshapeR package, resulting in fourteen (14) output files per dataset. The extracted features include:

### Inter-base pair parameters
- HelT (Helical Twist)
- Rise
- Roll
- Shift
- Slide
- Tilt

### Intra-base pair parameters
- Buckle
- Opening
- ProT (Propeller Twist)
- Shear
- Stagger
- Stretch

### Other features
- EP (Electrostatic Positivity)
- MGW (Minor Groove Width)

## Visualization
The shape features are visualized using schematic representations inspired by Figure 1 from Li et al., but generated using the provided Python scripts in this repository. These visualizations help interpret the structural properties of nucleosomal and linker DNA.

## Code and Usage
The repository includes Python scripts for:
- Distribution analysis
- Heatmap generation
- Position and variance profiling
- Statistical testing
- Shape profile visualization

Refer to the individual scripts for usage instructions and details on input/output formats.

## References
- Guo et al. (source of nucleosome positioning datasets)
- Zhuo et al. (dataset organization and download)
- Li et al. (schematic representation of DNA shape features)
- DNAshapeR package (Monte Carlo simulation and feature extraction)
