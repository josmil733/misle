#' @export decompose

decompose <- function(A, A.build, y, X, lattice) {
    n <- ifelse(
        length(unique(dim(A))) == 1,
        unique(dim(A)),
        stop("Adjacency matrix must be a square matrix.")
    )
    n.sim <- ncol(y)

    mvc <- mvc(A, A.build, lattice)
    mvc.c <- setdiff(1:n, mvc)

    A.mvc <- A[mvc, mvc]
    A.m <- A[mvc, mvc.c]

    A.ordered <- matrix(0, nrow = n, ncol = n)

    A.ordered[1:length(mvc), 1:length(mvc)] <- A.mvc
    A.ordered[1:length(mvc), (length(mvc) + 1):n] <- A.m
    A.ordered[(length(mvc) + 1):n, 1:length(mvc)] <- t(A.m)
    A.ordered[(length(mvc) + 1):n, (length(mvc) + 1):n] <- A[mvc.c, mvc.c]

    rownames(A) <- 1:n
    colnames(A) <- 1:n
    rownames(A.ordered) <- c(mvc, mvc.c)
    colnames(A.ordered) <- rownames(A.ordered)

    # y <- sample(c(-1,1), n, replace=TRUE)
    # y <- rep(-1, n)
    # y <- sample(1:5, 5) #for testing partitioning method
    y.ordered <- matrix(0, nrow = n, ncol = n.sim)
    y.ordered[1:length(mvc), ] <- y[mvc, ]
    y.ordered[(length(mvc) + 1):n, ] <- y[mvc.c, ]

    rownames(y) <- 1:n
    rownames(y.ordered) <- c(mvc, mvc.c)

    X.ordered <- X[c(mvc, mvc.c), , ]
    rownames(X.ordered) <- c(mvc, mvc.c)

    # will need to add arguments to get this part to work

    return(list(
        A.ordered = A.ordered,
        y.ordered = y.ordered,
        X.ordered = X.ordered,
        mvc = mvc
    ))
}
