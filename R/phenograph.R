#' RphenoGraph clustering
#'
#' R implementation of the PhenoGraph algorithm
#'
#' A simple R implementation of the [PhenoGraph](http://www.cell.com/cell/abstract/S0092-8674(15)00637-6) algorithm,
#' which is a clustering method designed for high-dimensional single-cell data analysis. It works by creating a graph ("network") representing
#' phenotypic similarities between cells by calclating the Jaccard coefficient between nearest-neighbor sets, and then identifying communities
#' using the well known [Louvain method](https://sites.google.com/site/findcommunities/) in this graph.
#'
#' Author: Chen Hao, 2015
#' Copied from https://github.com/JinmiaoChenLab/Rphenograph/blob/master/R/phenograph.R
#' for use in briDiscovr on 200611 under the Artistic-2.0 license
#'
#' Modified in briDiscovr 0.4.7: the graph now has one vertex per input row, in
#' row order, so \code{membership()} always lines up with the rows of
#' \code{data}. Rows that share no nearest neighbors with any of their own k
#' nearest neighbors (and therefore have no edges in the graph) are assigned
#' the most common cluster among their k nearest neighbors instead of forming
#' one-row clusters. Previously, such rows could cause an error or silently
#' misalign cluster labels with rows, mostly at small k.
#'
#' @param data matrix; input data matrix
#' @param k integer; number of nearest neighbours (default:30)
#'
#' @return a list contains an igraph graph object for \code{graph_from_data_frame} and a communities object, the operations of this class contains:
#' \item{print}{returns the communities object itself, invisibly.}
#' \item{length}{returns an integer scalar.}
#' \item{sizes}{returns a numeric vector.}
#' \item{membership}{returns a numeric vector, one number for each vertex in the graph that was the input of the community detection.}
#' \item{modularity}{returns a numeric scalar.}
#' \item{algorithm}{returns a character scalar.}
#' \item{crossing}{returns a logical vector.}
#' \item{is_hierarchical}{returns a logical scalar.}
#' \item{merges}{returns a two-column numeric matrix.}
#' \item{cut_at}{returns a numeric vector, the membership vector of the vertices.}
#' \item{as.dendrogram}{returns a dendrogram object.}
#' \item{show_trace}{returns a character vector.}
#' \item{code_len}{returns a numeric scalar for communities found with the InfoMAP method and NULL for other methods.}
#' \item{plot}{for communities objects returns NULL, invisibly.}
#'
#' @references Jacob H. Levine and et.al. Data-Driven Phenotypic Dissection of AML Reveals Progenitor-like Cells that Correlate with Prognosis. Cell, 2015.
#' @examples
#' \dontrun{
#' iris_unique <- unique(iris) # Remove duplicates
#' data <- as.matrix(iris_unique[,1:4])
#' Rphenograph_out <- Rphenograph(data, k = 45)
#' modularity(Rphenograph_out[[2]])
#' membership(Rphenograph_out[[2]])
#' iris_unique$phenograph_cluster <- factor(membership(Rphenograph_out[[2]]))
#'  ggplot(
#'   iris_unique,
#'   aes(x=Sepal.Length, y=Sepal.Width, col=Species, shape=phenograph_cluster)
#'  ) +
#'  geom_point(size = 3)+theme_bw()
#' }
#' @importFrom igraph cluster_louvain modularity membership
#' @useDynLib briDiscovr
#'
#' @export
Rphenograph <- function(data, k=30){
  if(is.data.frame(data))
    data <- as.matrix(data)

  if(!is.matrix(data))
    stop("Wrong input data, should be a data frame of matrix!")

  if(k<1){
    stop("k must be a positive integer!")
  }else if (k > nrow(data)-2){
    stop("k must be smaller than the total number of points!")
  }

  message("Run Rphenograph starts:","\n",
    "  -Input data of ", nrow(data)," rows and ", ncol(data), " columns","\n",
    "  -k is set to ", k)

  cat("  Finding nearest neighbors...")
  t1 <- system.time(neighborMatrix <- find_neighbors(data, k=k+1)[,-1])
  cat("DONE ~",t1[3],"s\n", " Compute jaccard coefficient between nearest-neighbor sets...")
  # Builds the graph with one vertex per row, in row order (see buildPhenographGraph)
  t2 <- system.time(g <- buildPhenographGraph(neighborMatrix))

  # Other community detection algorithms:
  #    cluster_walktrap, cluster_spinglass,
  #    cluster_leading_eigen, cluster_edge_betweenness,
  #    cluster_fast_greedy, cluster_label_prop
  cat("DONE ~",t2[3],"s\n", " Run louvain clustering on the graph ...")
  t3 <- system.time(community <- cluster_louvain(g))
  cat("DONE ~",t3[3],"s\n")

  # Give cells with no graph edges the majority cluster of their nearest neighbors.
  # Only $membership is updated; modularity is unchanged by this reassignment.
  isolatedResult <- assignIsolatedCells(as.integer(membership(community)), g, neighborMatrix)
  community$membership <- isolatedResult$membership

  message("Run Rphenograph DONE, totally takes ", sum(c(t1[3],t2[3],t3[3])), "s.")
  cat("  Return a community class\n  -Modularity value:", modularity(community),"\n")
  cat("  -Number of clusters:", length(unique(membership(community))), "\n")
  if(isolatedResult$nIsolated > 0){
    message(
      "  -", isolatedResult$nIsolated, " of ", nrow(data), " rows shared no ",
      "nearest neighbors with any of their own nearest neighbors; they were ",
      "assigned to the most common cluster among their ", k, " nearest neighbors."
    )
  }

  return(list(g, community))
}


#' K Nearest Neighbour Search
#'
#' Uses a kd-tree to find the p number of near neighbours for each point in an input/output dataset.
#'
#' Use the nn2 function from the RANN package, utilizes the Approximate Near Neighbor (ANN) C++ library,
#' which can give the exact near neighbours or (as the name suggests) approximate near neighbours
#' to within a specified error bound. For more information on the ANN library please
#' visit http://www.cs.umd.edu/~mount/ANN/.
#'
#' @param data matrix; input data matrix
#' @param k integer; number of nearest neighbours
#'
#' @return a n-by-k matrix of neighbor indices
#'
#' @examples
#' iris_unique <- unique(iris) # Remove duplicates
#' data <- as.matrix(iris_unique[,1:4])
#' neighbors <- find_neighbors(data, k=10)
#'
#' @importFrom RANN nn2
#' @export
find_neighbors <- function(data, k){
  nearest <- nn2(data, data, k, searchtype = "standard")
  return(nearest[[1]])
}
