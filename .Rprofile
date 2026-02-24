local({
  # Bioconductor version must match the "bioconductor.version" in renv/settings.json
  bioc_version <- "3.18"
  bioc_repos <- c(
    BioCsoft = sprintf("https://bioconductor.org/packages/%s/bioc", bioc_version),
    BioCann = sprintf("https://bioconductor.org/packages/%s/data/annotation", bioc_version),
    BioCexp = sprintf("https://bioconductor.org/packages/%s/data/experiment", bioc_version),
    BioCworkflows = sprintf("https://bioconductor.org/packages/%s/workflows", bioc_version),
    BioCbooks = sprintf("https://bioconductor.org/packages/%s/books", bioc_version),
    CRAN = "https://cloud.r-project.org"
  )
  existing <- getOption("repos", default = character(0))
  options(repos = c(bioc_repos, existing[!names(existing) %in% names(bioc_repos)]))
})

source("renv/activate.R")
