#' Graph to select m
#'
#' This function allows to draw a graph that relates the number of super efficient units and the choice of m
#'
#' @param output matrix (or vector) of indicators along which the units are evaluated.
#' @param mseries vector containing the different values of f that needed to be tested.
#' @param check vector containing the values of the thresholds to be considered to define the superefficient units
#' @param col vector containing the colors. the vector col must contain the same number of element of the vector check.
#' @param legend legend to be added to the graph.
#' @param B number of bootstrap replicates
#' @param RTS For more details see the dea function in the package Benchmarking. Text string or a number defining the underlying DEA technology / returns to scale assumption.
#' 0	fdh	Free disposability hull, no convexity assumption
#' 1	vrs	Variable returns to scale, convexity and free disposability
#' 2	drs	Decreasing returns to scale, convexity, down-scaling and free disposability
#' 3	crs	Constant returns to scale, convexity and free disposability
#' 4	irs	Increasing returns to scale, (up-scaling, but not down-scaling), convexity and free disposability
#' 5	irs2	Increasing returns to scale (up-scaling, but not down-scaling), additivity, and free disposability
#' 6	add	Additivity (scaling up and down, but only with integers), and free disposability; also known af replicability and free disposability, the free disposability and replicability hull (frh) -- no convexity assumption
#' 7	fdh+	A combination of free disposability and restricted or local constant return to scale
#' 10	vrs+	As vrs, but with restrictions on the individual lambdas via param
#' @param ORIENTATION For more details see the dea function in the package Benchmarking. Input efficiency "in" (1), output efficiency "out" (2), and graph efficiency "graph" (3). For use with DIRECT, an additional option is "in-out" (0).
#' @param print If print = TRUE the number of the unit under evaluation is printed. In case of large sample the function could require some time, so it could be useful to control how many units have already been evaluated and which one still have to be evaluated. By default print = FALSE.
#' @keywords BOD, m, graph
#' @export
#' @examples
#' graph1_m_BOD()
#'

graph1_m_BOD <- function(output, mseries, B, RTS, ORIENTATION, check, col, legend, print = TRUE) {

  meff <- matrix(0, nrow = length(mseries), ncol = length(check))

  for (m in 1:length(mseries)) {
    if (print == TRUE) {
      sprintf("The code is now computing the Robust BOD scores for m = %s", mseries[m])
    }
    BOD <-  robust_BOD(output, m = mseries[m], B, RTS = RTS, ORIENTATION = ORIENTATION)

    for(l in 1:length(check)) {
      meff[m,l] <- length(BOD[BOD>check[l]])
    }

  }

  plot(x=mseries, y=meff[,1], type = "b", lwd = 2,
       main = "Number of super-efficient units",
       xlab = c("m"),
       ylab = c("number of super-efficient units"),
       ylim=c(0, max(meff[,1])) )

  for(l in 2:length(check)) {
    lines(x=mseries, y=meff[,l], type = "b", lwd = 2, col = col[l])
  }

  legend("topright", legend, col=col, lwd=2)

}
