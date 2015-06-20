#' SCS - Splitting Conic Solver 
#'
#' @param A data matrix
#' @param b primal constraint vector
#' @param c primal objective vector
#' @param cone list of cone sizes
#' @param params solver parameters
#' @return list of solution vectors x,y,s and information about run
#' @seealso \code{\link{nchar}} which this function wraps
#' @export
#' @examples
#' scs(matrix(c(0.5,2),2,1), c(3, 1), 1, list(l=2), list(max_iters=5000))
scs <- function(A, b, c, cone, params) {
        Ax <- c(as.matrix(A)) # values
        data <- list(m = dim(A)[1], n = dim(A)[2], Ax = Ax, b = b, c = c)
        ret <- .Call("scsr", data, cone, params, PACKAGE="scs")
}
