#' Build the PhenoGraph shared-nearest-neighbor graph with one vertex per cell
#'
#' Computes Jaccard coefficients between nearest-neighbor sets and builds an
#' undirected weighted graph. Every cell is declared as a vertex, in row order,
#' so that vertex \code{i} always corresponds to row \code{i} of the input data.
#'
#' Without declaring the vertices explicitly, \code{graph_from_data_frame()}
#' orders vertices by first appearance in the edge list. At small k, some cells
#' share no neighbors with any of their own neighbors and have no edge of their
#' own; such cells are then either moved out of row order (silently misaligning
#' cluster labels for all subsequent rows) or dropped from the graph entirely
#' (causing a row-count error). Declaring vertices prevents both.
#'
#' \code{jaccard_coeff()} records a pair of cells once from each cell's
#' neighbor list, each time with half the Jaccard coefficient. These two
#' half-weight edges are merged into one edge carrying the full weight. The
#' merged graph has the same weighted connections, and therefore the same
#' modularity for any partition, but has fewer edges, so clustering is faster.
#' Self-loops (possible when cells are duplicated) are kept, as before.
#'
#' @param neighborMatrix An integer matrix with one row per cell, giving the
#' row indices of each cell's k nearest neighbors (excluding the cell itself)
#' @return An igraph graph with \code{nrow(neighborMatrix)} vertices named
#' \code{1..nrow(neighborMatrix)} in row order
#' @importFrom igraph make_graph simplify vcount V E
#' @keywords internal
buildPhenographGraph <- function(neighborMatrix){
  nCells <- nrow(neighborMatrix)

  links <- jaccard_coeff(neighborMatrix)
  # keep only rows that were filled in (cell pairs that share at least one neighbor)
  links <- links[links[, 1] > 0, , drop = FALSE]

  # build from integer vertex IDs; setting n declares every cell as a vertex,
  # including cells with no edges, so vertex i is cell i
  edgeList <- as.vector(t(links[, 1:2]))
  g <- igraph::make_graph(edgeList, n = nCells, directed = FALSE)
  igraph::E(g)$weight <- links[, 3]
  igraph::V(g)$name <- as.character(seq_len(nCells))

  # merge the two half-weight edges recorded for each pair into one full-weight edge
  g <- igraph::simplify(
    g,
    remove.multiple = TRUE,
    remove.loops = FALSE,
    edge.attr.comb = list(weight = "sum")
  )

  # confirm that vertex i is cell i; cluster labels are matched to cells by position
  if(igraph::vcount(g) != nCells ||
     !identical(as.integer(igraph::V(g)$name), seq_len(nCells))){
    stop(
      "Internal error: PhenoGraph graph vertices are not aligned with input rows. ",
      "Please contact the package developer for assistance."
    )
  }

  return(g)
}


#' Assign isolated cells to the majority cluster of their nearest neighbors
#'
#' Cells with no edges in the PhenoGraph graph (they share no neighbors with
#' any of their own nearest neighbors, which happens mostly at small k) are
#' placed by Louvain into one-cell clusters. This function instead gives each
#' such cell the most common cluster label among its k nearest neighbors,
#' ignoring neighbors that are themselves isolated. Ties are broken in favor of
#' the label of the closest neighbor among the tied clusters. Cluster IDs are
#' then renumbered to be consecutive (1, 2, 3, ...), preserving their order.
#'
#' An isolated cell whose neighbors are all isolated keeps its own cluster.
#'
#' Reassigning isolated cells does not change the modularity of the partition,
#' because a vertex with no edges contributes nothing to modularity.
#'
#' @param membershipVect An integer vector of cluster labels, one per cell, in row order
#' @param g The graph from \code{buildPhenographGraph()}
#' @param neighborMatrix The nearest-neighbor matrix used to build \code{g},
#' with neighbors in order of increasing distance
#' @return A list with elements \code{membership} (the updated, renumbered
#' cluster labels) and \code{nIsolated} (the number of isolated cells)
#' @importFrom igraph degree
#' @keywords internal
assignIsolatedCells <- function(membershipVect, g, neighborMatrix){
  # self-loops (possible with duplicated cells) do not connect a cell to anything else
  isIsolated <- igraph::degree(g, loops = FALSE) == 0
  isolatedCells <- which(isIsolated)

  newMembership <- membershipVect
  for(cell in isolatedCells){
    neighbors <- neighborMatrix[cell, ]
    connectedNeighbors <- neighbors[!isIsolated[neighbors]]
    if(length(connectedNeighbors) == 0){
      next # nothing to borrow a label from; keep the singleton cluster
    }
    # labels of connected neighbors, in order of increasing distance
    neighborLabels <- membershipVect[connectedNeighbors]
    labelCounts <- table(neighborLabels)
    tiedLabels <- names(labelCounts)[labelCounts == max(labelCounts)]
    # first (closest) neighbor whose label is among the most common
    newMembership[cell] <- neighborLabels[match(TRUE, as.character(neighborLabels) %in% tiedLabels)]
  }

  # renumber so cluster IDs are consecutive, keeping their original order
  newMembership <- match(newMembership, sort(unique(newMembership)))

  return(list(membership = newMembership, nIsolated = length(isolatedCells)))
}
