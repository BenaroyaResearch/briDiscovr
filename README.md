# briDiscovr

## Dependencies

This package depends on a number of bioconductor packages, which may not install automatically. If you encounter any difficulties during installation, you can manually install the bioconductor dependencies using the following code:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install("flowStats")
BiocManager::install("Rtsne")
BiocManager::install("umap")
```
