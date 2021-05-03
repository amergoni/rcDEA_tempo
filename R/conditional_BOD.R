#' Conditional BOD function
#'
#' This function allows to compute Robust and Conditional BOD scores.
#'
#' @param output matrix (or vector) of indicators along which the units are evaluated.
#' @param exogenous matrix (or vector) of exogenous variables involved in the conditional analysis. The similarity among the units is determined according to the exogeneous variable(s) using the function npudensbw and npudens (from the package np) with epanechnikov kernel.
#' @param similarity matrix of similarities. In alternative to provide the exogenous variables, the matrix of similarities can be directly provided. This allow to customize the estimation of the similarities.
#' @param m number of unit to be included in the reference set
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
#' @param inclusion If inclusion = TRUE the unit under analysis is included in the reference set. So, no super efficient scores are allowed. By default inclusion = FALSE.
#' @param print If print = TRUE the number of the unit under evaluation is printed. In case of large sample the function could require some time, so it could be useful to control how many units have already been evaluated and which one still have to be evaluated. By default print = FALSE.
#' @keywords BOD, Conditional
#' @export
#' @examples
#' conditional_BOD()
#'

conditional_BOD <- function(output, exogenous = FALSE, m, B, RTS, ORIENTATION, similarity = FALSE, inclusion = FALSE, print=FALSE) {

  #define preliminary variables
  n <- nrow(output)
  k <- ncol(output)
  exogenous <- as.data.frame(exogenous)

  BOD <- rep(0, n)
  BOD_B<- rep(0, B) #bootstrap BOD

  if(similarity == FALSE) {
  print(c("R is now computing the bandwidth using the function npudensbw in the package np"))
  bw = npudensbw(dat=exogenous)
  }

  if(print == TRUE) {
    print(c("R is now computing the conditional BOD"))
  }

    for (i in 1:n) {
      if(print == TRUE) {
      print(i)
      }

      #in case the code needs to compute the similarity for each obs i
      if(similarity == FALSE) {
         kerz <- npudens(bws=bw,
                          cykertype="epanechnikov",cxkertype="epanechnikov",
                          tdat=exogenous[i,],edat=exogenous)
          similarity_i <- cbind(kerz$dens)
      }

      #in case the matrix of similarities has already been computed
      if(similarity != FALSE) { #equivalent of asking if the values for the similarity are inserted
        similarity_i <- similarity[i,]
      }

      #consider only the units that perform at least as good as unit i
      y <- output[i, ]
      x <- 1
      Y_Rob <- output
      similarity_Rob <- similarity_i
      if (ORIENTATION == "out") {
        for (l in 1:k) {
          similarity_Rob <- similarity_Rob[Y_Rob[, l] >= y[l]]
          Y_Rob <- as.data.frame(Y_Rob[Y_Rob[, l] >= y[l], ])
        }
      }
      Y_Rob <- as.data.frame(Y_Rob)
      n_sample <- nrow(Y_Rob)

      #pick a sample of random unit in the reference set if there are at least 2 units in the ref
      if (n_sample < 2) {
        BOD[i] <- 1
      }
      else {
        for (j in 1:B) {
          #select m random units
          m_sample <- sample(n_sample, m, prob = similarity_Rob, replace = TRUE)
          Y_ref <- as.data.frame(Y_Rob[m_sample,])

          if (inclusion == TRUE) {#I add the same observation to be sure it is included
            Y_ref <- as.data.frame(rbind(y, Y_ref))
          }

          n_m <- nrow(Y_ref)
          ones_ref <- rep(1, n_m)

          #compute the BOD for unit i
          BOD_B[j] <- dea(X = x, Y = y, RTS = RTS, ORIENTATION = ORIENTATION, XREF = ones_ref, YREF = Y_ref)$eff
        }
        BOD[i] <- mean(BOD_B)
      }
    }
  c_BOD <- BOD
}







