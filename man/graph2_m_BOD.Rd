% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph2_m_BOD.R
\name{graph2_m_BOD}
\alias{graph2_m_BOD}
\title{Graph to select m}
\usage{
graph2_m_BOD(output, mseries, B, RTS, ORIENTATION, print = TRUE)
}
\arguments{
\item{output}{matrix (or vector) of indicators along which the units are evaluated.}

\item{mseries}{vector containing the different values of f that needed to be tested.}

\item{B}{number of bootstrap replicates}

\item{RTS}{For more details see the dea function in the package Benchmarking. Text string or a number defining the underlying DEA technology / returns to scale assumption.
0    fdh    Free disposability hull, no convexity assumption
1    vrs    Variable returns to scale, convexity and free disposability
2    drs    Decreasing returns to scale, convexity, down-scaling and free disposability
3    crs    Constant returns to scale, convexity and free disposability
4    irs    Increasing returns to scale, (up-scaling, but not down-scaling), convexity and free disposability
5    irs2    Increasing returns to scale (up-scaling, but not down-scaling), additivity, and free disposability
6    add    Additivity (scaling up and down, but only with integers), and free disposability; also known af replicability and free disposability, the free disposability and replicability hull (frh) -- no convexity assumption
7    fdh+    A combination of free disposability and restricted or local constant return to scale
10    vrs+    As vrs, but with restrictions on the individual lambdas via param}

\item{ORIENTATION}{For more details see the dea function in the package Benchmarking. Input efficiency "in" (1), output efficiency "out" (2), and graph efficiency "graph" (3). For use with DIRECT, an additional option is "in-out" (0).}

\item{print}{If print = TRUE the number of the unit under evaluation is printed. In case of large sample the function could require some time, so it could be useful to control how many units have already been evaluated and which one still have to be evaluated. By default print = FALSE.}
}
\description{
This function allows to draw a graph that relates the average efficiency score and the choice of m
}
\examples{
graph2_m_BOD()

}
\keyword{BOD,}
\keyword{graph}
\keyword{m,}
