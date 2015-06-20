\name{scs}
\alias{scs}
\title{SCS - Splitting Conic Solver}
\usage{
scs(A, b, c, cone, params, use_indirect = false)
}
\arguments{
\item{A}{data matrix}

\item{b}{primal constraint vector}

\item{c}{primal objective vector}

\item{cone}{list of cone sizes}

\item{params}{solver parameters}

\item{use_indirect}{boolean to use the indirect solver}
}
\value{
list of solution vectors x,y,s and information about run
}
\description{
Solves a convex cone program
}
\examples{
sol <- scs(matrix(c(0.5,2),2,1), c(3, 1), 1, list(l=2), list(max_iters=5000))
}
\seealso{
\code{\link{nchar}} which this function wraps
}

