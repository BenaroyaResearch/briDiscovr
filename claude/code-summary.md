# briDiscovr Code Summary

## Overview

**briDiscovr** is an R package that implements the DISCOV-R (Distribution analysis across clusters of a parent population overlaid with a rare subpopulation) analysis method for flow cytometry and mass cytometry (CyTOF) data. The method was originally published in Wiedeman, Muir, et al. 2020 (PMID 31815738) for analyzing autoreactive CD8+ T cell exhaustion in type 1 diabetes.

**Version:** 0.4.2
**License:** Custom (see LICENSE file)
**Authors:** Mario Rosasco, Virginia Muir, Matt Dufort
**Repository:** https://github.com/BenaroyaResearch/briDiscovr

## Purpose

The package enables researchers to:
1. Define phenotypic clusters from a parent cell population (e.g., all CD8+ T cells)
2. Map rare cell populations (e.g., antigen-specific memory cells) onto the parent population's phenotypic landscape
3. Robustly characterize differences in rare cell populations through metaclustering and visualization

## Package Structure

### Core Object Class

**discovrExperiment (S3 Class)**
- Central data structure that holds all experiment information
- Tracks status through analysis phases: initialized → clustered → normalized → metaclustered
- Contains marker information, FCS file data, clustering results, and visualization parameters

### Main Analysis Workflow

The package follows a sequential workflow:

1. **Setup Phase** (`setupDiscovrExperiment`)
   - Loads FCS files and marker information
   - Validates data integrity
   - Transforms data using arcsinh transformation
   - Creates discovrExperiment object

2. **Clustering Phase** (`clusterDiscovrExperiment`)
   - Performs PhenoGraph clustering on each sample
   - Uses k-nearest neighbors and Louvain community detection
   - Identifies phenotypic clusters within parent populations

3. **Normalization Phase** (`normalizeDiscovrExperiment`)
   - Normalizes marker expression across samples
   - Supports multiple normalization methods: z-score, warpSet, or none
   - Essential for meaningful metaclustering

4. **Metaclustering Phase** (`metaclusterDiscovrExperiment`)
   - Groups similar clusters across samples
   - Applies hierarchical clustering with configurable linkage/distance metrics
   - Filters low-abundance clusters
   - Maps rare populations onto metaclusters

5. **Visualization Phase** (`makeMetaclusterHeatmaps`)
   - Generates comprehensive heatmaps
   - Shows normalized counts and marker intensities
   - Visualizes metacluster phenotypes

6. **Dimensionality Reduction** (`runUmapDiscovrExperiment`)
   - Projects high-dimensional data to 2D using UMAP
   - Enables visualization of cell populations

## Directory Structure

```
briDiscovr/
├── R/                          # R source code
│   ├── briDiscovr.R           # Package imports and dynamic library loading
│   ├── discovrExperiment.R    # S3 class definition and methods
│   ├── clusteringPhase.R      # Setup and clustering functions
│   ├── metaclusteringPhase.R  # Normalization and metaclustering
│   ├── heatmapGenerators.R    # Visualization functions
│   ├── exportedUtils.R        # Utility functions for users
│   ├── internalUtils.R        # Internal helper functions
│   ├── phenograph.R           # PhenoGraph clustering implementation
│   └── RcppExports.R          # R wrappers for C++ functions
├── src/                        # C++ source code
│   ├── jaccard_coeff.cpp      # Jaccard coefficient computation
│   └── RcppExports.cpp        # Rcpp interface (auto-generated)
├── man/                        # R documentation files
├── DESCRIPTION                 # Package metadata and dependencies
├── NAMESPACE                   # Exported functions and imports
└── README.md                   # User documentation
```

## Key Components

### R/briDiscovr.R
- Minimal file that loads the compiled C++ library
- Imports Rcpp for C++ integration

### R/discovrExperiment.R
Core S3 class implementation:
- `is.discovrExperiment()` - Class validation
- `print.discovrExperiment()` - Display experiment summary
- `getSampleNames()` - Extract sample identifiers

### R/clusteringPhase.R
Main analysis setup (~350 lines):
- `setupDiscovrExperiment()` - Initialize experiment from FCS files
  - Validates marker and FCS information files
  - Checks memory requirements
  - Loads and processes FCS data
  - Applies arcsinh transformation
  - Supports optional downsampling
- `clusterDiscovrExperiment()` - Perform PhenoGraph clustering
  - Clusters each sample independently
  - Computes cluster statistics and means

### R/metaclusteringPhase.R
Metaclustering and normalization (~600 lines):
- `normalizeDiscovrExperiment()` - Normalize marker expression
  - Supports z-score normalization
  - Supports warpSet normalization (with configurable peak numbers)
  - Can apply different methods to different markers
- `metaclusterDiscovrExperiment()` - Create metaclusters
  - Hierarchical clustering of normalized cluster means
  - Filters low-abundance clusters
  - Maps events to metaclusters
  - Computes occupancy statistics
- `recutMetaclusters()` - Re-cut dendrogram to change metacluster count
- `testNormalizationMethodByMarker()` - Test different normalization approaches

### R/exportedUtils.R
User-facing utility functions (~1900 lines):
- File validation and quality control:
  - `checkFcsFiles()` - Validate FCS file structure
  - `checkForFcsByteOffsetIssue()` - Detect corrupted FCS headers
  - `getFcsNEvents()` - Count events in FCS files
- Data extraction:
  - `getSubjectCounts()` - Event counts per sample
  - `getSubjectClusters()` - Cluster counts per sample
  - `getMetaclusterOccupancy()` - Cell distribution across metaclusters
  - `getEventIntensities()` - Extract marker intensities for specific events
  - `getEventZScores()` - Extract normalized z-scores
- Downsampling:
  - `downsampleFcsList()` - Reduce dataset size for large experiments
- Visualization helpers:
  - `plotDensityNormalizedExprsDiscovrExperiment()` - Assess normalization quality
  - `runUmapDiscovrExperiment()` - UMAP dimensionality reduction

### R/heatmapGenerators.R
Visualization functions (~700 lines):
- `makeMetaclusterHeatmaps()` - Main heatmap generation
  - Creates multiple heatmap types:
    - Normalized count heatmaps showing cluster phenotypes
    - Intensity heatmaps showing marker expression
    - Separate heatmaps for clustering vs non-clustering markers
  - Highly customizable color schemes and layouts
  - Exports high-resolution PNG files

### R/internalUtils.R
Internal helper functions (~190 lines):
- `buildFcsList()` - Load FCS files into list structure
- `processFcsList()` - Clean and validate marker names
- `getColorList()` - Generate distinct colors for visualization
- `checkMarkerInfoNormalizationMethod()` - Validate normalization specifications
- `normalizeWarpSetMergedExpr()` - Apply warpSet normalization
- `normalizeWarpSetSingleMarker()` - Normalize individual markers

### R/phenograph.R
PhenoGraph clustering implementation (~120 lines):
- `Rphenograph()` - Main clustering function
  - Implements graph-based clustering for single-cell data
  - Creates k-nearest neighbor graph
  - Computes Jaccard coefficients between neighbor sets
  - Applies Louvain community detection
- `find_neighbors()` - K-nearest neighbor search using RANN library
- Based on the algorithm from Levine et al., Cell 2015
- Licensed under Artistic-2.0

### src/jaccard_coeff.cpp
C++ implementation for performance (~40 lines):
- `jaccard_coeff()` - Fast computation of Jaccard coefficients
  - Computes similarity between nearest-neighbor sets
  - Symmetrizes the undirected graph
  - Critical for PhenoGraph performance on large datasets
- Uses Rcpp for R integration
- Original code from Chen Hao, used under Artistic-2.0 license

### src/RcppExports.cpp
Auto-generated Rcpp interface:
- R-to-C++ function wrappers
- Automatically maintained by Rcpp

## Key Features

### Data Input
- Accepts FCS files from flow cytometry or CyTOF
- Requires two CSV configuration files:
  1. **Marker information**: Maps FCS channel names to common marker names, specifies clustering markers
  2. **FCS information**: Lists sample IDs, cell subsets, and file paths

### Normalization Methods
1. **Z-score normalization**: Standard normalization to mean 0, SD 1 per sample
2. **warpSet normalization**: Aligns marker distributions across samples using landmark registration
   - Can specify number of peaks for better alignment
   - Applied marker-by-marker for reproducibility
3. **No normalization**: For markers that shouldn't be normalized

### Memory Management
- Built-in memory checks before loading large datasets
- Optional downsampling to reduce memory footprint
- Garbage collection at strategic points
- Efficient data structures to minimize duplication

### Clustering Algorithm
- **PhenoGraph**: Graph-based clustering using:
  - K-nearest neighbors (default k=30)
  - Jaccard coefficient for edge weights
  - Louvain method for community detection
- Per-sample clustering preserves biological variation
- Metaclustering groups similar clusters across samples

### Quality Control
- FCS file validation (byte offsets, event counts)
- Marker name consistency checks
- Memory usage warnings
- Low-abundance cluster filtering
- Normalization quality visualization

### Visualization
- Multiple heatmap types with customizable color schemes
- Dendrogram-based metacluster ordering
- Separate displays for clustering and non-clustering markers
- UMAP projections for intuitive data exploration
- High-resolution export suitable for publication

## Dependencies

### Core Dependencies
- **Rcpp**: C++ integration for performance
- **flowCore, flowStats**: FCS file handling and processing
- **igraph**: Graph-based clustering
- **RANN**: Fast nearest-neighbor search
- **ComplexHeatmap**: Advanced heatmap visualization
- **umap**: Dimensionality reduction

### Data Manipulation
- **dplyr, tidyr, tibble**: Data wrangling
- **reshape2**: Data reshaping
- **magrittr**: Pipe operators

### Visualization Support
- **ggplot2**: Plotting framework
- **RColorBrewer, circlize, viridisLite**: Color palettes
- **grid**: Graphics system

### Other
- **stringr**: String manipulation
- **memuse**: Memory monitoring
- **rlang, methods**: R programming utilities

## Technical Considerations

### Performance Optimizations
- C++ implementation for computationally intensive Jaccard coefficient calculation
- Efficient nearest-neighbor search using kd-trees (RANN library)
- Marker-by-marker normalization for reproducibility
- Strategic memory cleanup

### Data Transformation
- **arcsinh transformation**: Applied to raw intensity values
  - Configurable parameters (a, b, c) for different cytometry platforms
  - Default: a=0, b=0.2 (1/5 for CyTOF), c=0
  - Recommended: b=1/150 for traditional flow cytometry

### Workflow Design
- Sequential processing ensures data integrity
- Each major function returns updated discovrExperiment object
- Status tracking prevents out-of-order operations
- Intermediate results can be saved/loaded using standard R data handling

## Use Cases

### Primary Applications
1. **Immunology research**: Characterizing rare antigen-specific T cell populations
2. **Biomarker discovery**: Identifying phenotypic signatures associated with disease
3. **Longitudinal studies**: Tracking immune cell phenotypes over time
4. **Comparative studies**: Comparing immune responses across conditions

### Typical Workflow Example
```r
# 1. Setup
expt <- setupDiscovrExperiment(
  markerInfoFile = "markers.csv",
  fcsInfoFile = "samples.csv",
  parentPopulation = "CD8_T_cells"
)

# 2. Cluster
expt <- clusterDiscovrExperiment(expt)

# 3. Normalize
expt <- normalizeDiscovrExperiment(expt)

# 4. Metacluster
expt <- metaclusterDiscovrExperiment(
  expt,
  nMetaclusters = 12
)

# 5. Visualize
makeMetaclusterHeatmaps(expt)

# 6. Extract results
occupancy <- getMetaclusterOccupancy(expt)
```

## Version History Context

The current version (0.4.2) represents a major evolution:
- **v0.4**: Introduced flexible normalization methods (breaking change from v0.3)
  - Previously all markers were z-score normalized automatically during metaclustering
  - Now normalization is a separate, configurable step
- Recent updates focused on:
  - More reproducible warpSet normalization
  - Better memory handling
  - Enhanced diagnostic functions
  - Pre-normalization mean calculations for metaclustering

## Summary

briDiscovr provides a complete, well-structured solution for analyzing rare cell populations in high-dimensional cytometry data. The package balances ease of use with flexibility, offering sensible defaults while allowing customization at each analysis step. Its modular design, comprehensive documentation, and robust error handling make it suitable for both routine analysis and advanced research applications.
