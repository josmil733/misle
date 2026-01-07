#' @export mvc

mvc <- function(adj, adj.build, lattice){
  # adj <- A
  n <- ifelse(length(unique(dim(adj)))==1, unique(dim(adj)), stop("Adjacency matrix must be a square matrix."))

  if(!lattice) return(sample(n, n/2))

  # edges <- which(adj!=0, arr.ind = TRUE)
  # for now, use a dummy value

  # 2/6/25: building output only for the 2D connected lattice

  if(lattice){
    if(nrow(adj.build)%%2){
      return(2*(1:(floor(n/2))) ) #if odd number of rows, select even numbered points
    }
    if(!nrow(adj.build)%%2){
      points <- c()
      for(j in 1:ncol(adj.build)){
        position <- (j-1)*nrow(adj.build)
        if(j%%2){
          points <- c(points, position + which(!((position+1:nrow(adj.build))%%2) ) )
        } else {
          points <- c(points, position + which(as.logical((position+1:nrow(adj.build))%%2) ) )
        }
      }
      return(points)
    }
  }

}
