## Internal helpers for clusterDiscovrExperiment(). Not exported.

# Cluster the events of one sample, returning one cluster ID per row of `data`.
# Steps: exact k-nearest-neighbor search -> PhenoGraph graph (Jaccard-weighted
# shared nearest neighbors) -> Louvain (method = "phenograph") or Leiden
# (method = "leiden") community detection -> isolated events given the majority
# cluster of their neighbors.
# The caller is responsible for setting the random seed beforehand.
clusterEvents <- function(
  data, method, k, nThreads, leidenResolution, leidenIterations, verbose
){
  nDigits <- 3 # number of digits for time display
  data <- as.matrix(data)
  if(nrow(data) < k + 2){
    stop("k must be smaller than the total number of points!")
  }
  algorithmName <- if(method == "leiden") "Leiden" else "Louvain"

  if(verbose){
    message(
      "Run clustering starts:", "\n",
      "  -Input data of ", nrow(data), " rows and ", ncol(data), " columns", "\n",
      "  -k is set to ", k, "\n",
      "  Finding nearest neighbors...")
  }
  t1 <- system.time(neighborMatrix <- findNearestNeighbors(data, k = k, nThreads = nThreads))

  if(verbose){
    message(
      "DONE ~", round(t1[3], nDigits), "s\n",
      "Compute jaccard coefficient between nearest-neighbor sets and build graph...")
  }
  # builds the graph with one vertex per row, in row order
  t2 <- system.time(g <- buildPhenographGraph(neighborMatrix))

  if(verbose){
    message("DONE ~", round(t2[3], nDigits), "s\n", "Run ", algorithmName, " clustering on the graph...")
  }
  t3 <- system.time(
    community <- if(method == "leiden"){
      runLeiden(g, resolution = leidenResolution, nIterations = leidenIterations)
    } else {
      igraph::cluster_louvain(g)
    }
  )
  if(verbose){
    message("DONE ~", round(t3[3], nDigits), "s\n")
  }

  # give events with no graph edges the majority cluster of their nearest neighbors
  isolatedResult <- assignIsolatedCells(
    as.integer(igraph::membership(community)), g, neighborMatrix
  )
  clusterLabels <- isolatedResult$membership

  if(verbose){
    message("Run clustering DONE, took a total of ", round(sum(c(t1[3], t2[3], t3[3])), nDigits), "s.")
    message("  -Modularity value: ", igraph::modularity(community), "\n")
    message("  -Number of clusters: ", length(unique(clusterLabels)))
  }
  if(isolatedResult$nIsolated > 0){
    message(
      "  -", isolatedResult$nIsolated, " of ", nrow(data), " events shared no ",
      "nearest neighbors with any of their own nearest neighbors; they were ",
      "assigned to the most common cluster among their ", k, " nearest neighbors. ",
      "Consider increasing k if this number is large."
    )
  }
  return(clusterLabels)
}

# Find the k nearest neighbors of every event, excluding the event itself, in
# order of increasing distance. Uses exact Euclidean search (the k-means for
# k-nearest neighbors algorithm, KMKNN), so the result is the same for any
# number of threads. Returns an n x k matrix of row indices.
findNearestNeighbors <- function(data, k, nThreads){
  BiocNeighbors::findKNN(
    data,
    k = k,
    BNPARAM = BiocNeighbors::KmknnParam(),
    num.threads = nThreads,
    get.distance = FALSE
  )$index
}

# Leiden community detection, optimizing modularity (as Louvain does).
# Newer versions of igraph renamed the argument `resolution_parameter` to
# `resolution`; pass whichever name the installed version uses.
runLeiden <- function(g, resolution, nIterations){
  leidenArgs <- list(
    graph = g,
    objective_function = "modularity",
    n_iterations = nIterations
  )
  if("resolution" %in% names(formals(igraph::cluster_leiden))){
    leidenArgs$resolution <- resolution
  } else {
    leidenArgs$resolution_parameter <- resolution
  }
  do.call(igraph::cluster_leiden, leidenArgs)
}
