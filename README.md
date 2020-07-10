# briDiscovr

Functions implementing and supporting the "Distribution analysis across clusters of a parent population overlaid with a rare subpopulation" (DISCOV-R) analysis method (Wiedeman, Muir, et al. 2020).

**NB** - This code has been tested with R version 4.0.1 and 4.0.2. The change from R version 3.x to 4.x came with a number of substantial differences, and so while there's a possibility this package may work with earlier versions of R, it should be considered untested.

---

## Installation

1. Install R v4.0+ from https://cloud.r-project.org/bin/windows/base/
  * If you're using Windows, you'll need to install Rtools40 from https://cran.r-project.org/bin/windows/Rtools/ as well
2. Open R and install devtools using ```install.packages("devtools")```
3. Install this package using ```devtools::install_github("BenaroyaResearch/briDiscovr")```
  * If you encounter difficulties with any dependency packages, please see the "Dependencies" section at the bottom of this page.
  
---

## Usage

This package is based around the `discovrExperiment` object class, which is a list that holds all of the information about your experiment. You can save, load, and copy `discovrExperiment` object to facilitate different phases of the analysis.

### Preparing to start your analysis

Before running DISCOV-R, you'll need to prepare a marker information file and an. fcs information file. Both of these can be created in Excel, and should be saved in the 'comma-separated values (.csv)' format.

#### marker information file

This file should contain 3 columns:

* A column for the common name of your marker (eg: "CD45RA") - by default this column is expected to be named 'fixed'
* A column for the marker name as it's represented in the .fcs data (eg: "143Nd_CD45RA") - by default this column is expected to be named 'desc'
* A column named 'useToCluster' that indicates whether the marker will be used for clustering. This column can only contain the values "TRUE"" or "FALSE".

#### .fcs information file

This file should contain 3 columns:

* A column named 'subject' indicating the subject/sample ID (eg: "SubjXpb997364476"). 
  * Note that this should NOT start with a number.
* A column named 'cellSubset' indicating the parent and child populations to be used in the analysis (eg: "A2_CD8_T_Cells")
* A column named 'filename' indicating where to find the .fcs file for that patient and cell subset (eg: "/Users/mrosasco/Documents/projects/briUtils/briDiscovr/testFiles/Xpb997364476_A2_CD8_T_Cells.fcs"). 
  * Note that it's possible to use paths relative to your working directory, but it's generally preferred to provide an absolute path.
  
### Setting up your experiment

Once you have the information about your study design in your marker info and .fcs info files, you're ready to begin your analysis. This is done by using the `setupDiscovrExperiment()` function as shown below. More information about how to use this function can be accessed using `?setupDiscovrExperiment`.

```R
library(briDiscovr)

myExpt <- setupDiscovrExperiment(
    markerInfoFile = "discovrTestData/markerInformation.csv",
    fcsInfoFile = "discovrTestData/fcsManifest.csv",
    parentPopulation = "A2_CD8_T_Cells"
)
```

This will load in the .fcs files, check for any issues with the data, clean the marker names, and transform the event data using an arcsinh transform. For details on the appropriate transformation parameters for CyTOF and FACS data, please see `?setupDiscovrExperiment`.

Once the data are loaded, you can access summary information by using the functions `print()` and `getSubjectCounts()`. This is useful to do before clustering as a 'sanity check' to make sure that your markers look as expected, and that you have reasonable numbers of event counts for all individuals in your dataset.

```R
print(myExpt)
## Prints out:
# An object of class 'discovrExperiment'
# All markers: CD45, GRZMA, CD57, CD45RA, CD38, CD8, CD4, KLRG1, CD14, CD127, CD226, HELIOS, CD2, NKG2C, CD3, TIGIT, CD25, CD27, CD161, TBET, CD39, EOMES, CXCR3, CD95, CD19, NKG2A, CCR7, CD122, CD103, KI67, PD1, CD56, CD16
# Clustering markers: CD45, GRZMA, CD57, CD45RA, CD38, CD8, CD4, KLRG1, CD14, CD127, CD226, HELIOS, CD2, NKG2C, CD3, TIGIT, CD25, CD27, CD161, TBET, CD39, EOMES, CXCR3, CD95, CD19, NKG2A, CCR7, CD122, CD103, KI67, PD1, CD56, CD16
# Experiment status: initialized

getSubjectCounts(myExpt)
## Prints out:
#             sample nEvents
# 1 SubjXpb997364476  172634
```

The experiment can then be saved using normal R data handling utilities. For example:

```R
saveRDS(myExpt, "200709-textExpt.RDS")
```

### Clustering

Clustering is handled through the function `clusterDiscovrExperiment()`. Note that this step takes the longest, and may require several hours depending on the size of the dataset. Additional information about this function can be accessed using `?clusterDiscovrExperiment`. You can use the utility function `getSubjectClusters()` to assess the results of the clustering.

Note that this function returns a `discovrExperiment` object, and so you'll need to assign this to a variable. This can be the same variable you use to store your initialized `discovrExperiment`.

```R
myExpt <- clusterDiscovrExperiment(myExpt)
# May run for several hours

print(myExpt)
## Prints:
# An object of class 'discovrExperiment'
# All markers: CD45, GRZMA, CD57, CD45RA, CD38, CD8, CD4, KLRG1, CD14, CD127, CD226, HELIOS, CD2, NKG2C, CD3, TIGIT, CD25, CD27, CD161, TBET, CD39, EOMES, CXCR3, CD95, CD19, NKG2A, CCR7, CD122, CD103, KI67, PD1, CD56, CD16
# Clustering markers: CD45, GRZMA, CD57, CD45RA, CD38, CD8, CD4, KLRG1, CD14, CD127, CD226, HELIOS, CD2, NKG2C, CD3, TIGIT, CD25, CD27, CD161, TBET, CD39, EOMES, CXCR3, CD95, CD19, NKG2A, CCR7, CD122, CD103, KI67, PD1, CD56, CD16
# Experiment status: clustered

getSubjectClusters(myExpt)
## Prints:
#             sample kClusters
# 1 SubjXpb997364476        26
```

### Metaclustering

Metaclustering is handled throught the function `metaclusterDiscovrExperiment()`. The documentation at `?metaclusterDiscovrExperiment` contains information about how to change cluster occupancy thresholds, target numbers of metaclusters, metaclustering distance and linkage metrics, and how to exclude markers or samples from the metacluster analysis. Note that dropping markers or samples will *not* affect the clustered results, and to exclude samples altogether from the clustering a new experiment should be started from the initialization step. 

Note that as with the clustering function, `metaclusterDiscovrExperiment` will return a new `discovrExperiment` object that needs to be assigned to a variable.

```R
myExpt <- metaclusterDiscovrExperiment(myExpt)

print(myExpt)
## Prints:
# An object of class 'discovrExperiment'
# All markers: CD45, GRZMA, CD57, CD45RA, CD38, CD8, CD4, KLRG1, CD14, CD127, CD226, HELIOS, CD2, NKG2C, CD3, TIGIT, CD25, CD27, CD161, TBET, CD39, EOMES, CXCR3, CD95, CD19, NKG2A, CCR7, CD122, CD103, KI67, PD1, CD56, CD16
# Clustering markers: CD45, GRZMA, CD57, CD45RA, CD38, CD8, CD4, KLRG1, CD14, CD127, CD226, HELIOS, CD2, NKG2C, CD3, TIGIT, CD25, CD27, CD161, TBET, CD39, EOMES, CXCR3, CD95, CD19, NKG2A, CCR7, CD122, CD103, KI67, PD1, CD56, CD16
# Experiment status: metaclustered
```

### Heatmaps

Once you've performed metaclustering, you can use the function 'makeMetaclusterHeatmaps()` to produce and save heatmaps to summarize your data. By default this function assumes that the first subset listed in the `fcsInfoFile` used to initialize the experimnet is the parent population (eg: all CD8+s), and that all the remaining cell subsets are the rare/child populations. Details on explicitly assigning these populations for plotting are available in the docs at `?makeMetaclusterHeatmaps`. By default this function will save all heatmaps in the current working directory with the current date as the prefix, but this can be changed as in the example below.

```R
makeMetaclusterHeatmaps(
  myExpt,
  filenamePrefix = file.path(
    "Users", "mrosasco", "discovrOutput", format(Sys.Date(), "%y%m%d-")
  )
)
```

---

## Dependencies

This package depends on bioconductor packages, which may run into issues when installing automatically. If you encounter any difficulties during installation, you can manually install the bioconductor dependencies using the following code:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install("flowCore")
BiocManager::install("ComplexHeatmap")
```

The remaining depndencies should be automatically installed via CRAN, but if you encounter any issues you can manually install using the following code:

```R
install.packages(pkgs=c("rlang", "igraph", "RANN", "Rcpp", "methods", "tidyverse", "reshape2", "RColorBrewer", "circlize", "viridisLite")
```
