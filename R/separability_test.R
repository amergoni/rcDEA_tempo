#' Separability condition test
#'
#' This function test the separability condition and return the value of t and the p-value.  The main reference is the paper by Daraio and Simar (2015)
#'
#' @param rDEA vector containing the robust DEA (or BOD or FDH) scores
#' @param cDEA vector containing the conditional DEA (or BOD or FDH) scores
#' @param p    number input used in the analysys (in case of BOD analysis set p=1)
#' @param q    number output used in the analysis
#' @param B    number of bootstrap replicates
#' @param model model implemented to obtain the efficiency score. Possible values for this parameter are "DEA" or "FDH"
#' @param alpha  significance level. By default alpha = 0.05
#' @param print if print = TRUE the function print the result of the test. By default print = TRUE.
#' @keywords test, separability condition,
#' @export
#' @examples
#' separability_test()
#'

separability_test <- function(rDEA, cDEA, p, q, B, model, alpha = 0.05, print = TRUE) {

  n = length(rDEA)

  #define the convergence rate
  if(model == "FDH") {
    k = 1/(p+q)  }
  if(model == "DEA") {
    k = 2/(p+q+1)
  }


  if (model == "DEA" & (p+q) <= 4 | model == "FDH" & (p+q) <= 3) {
  #pick two independent sample  from the robust and the conditional
  n1 <- floor(n/2)
  n2 <- n - n1

  sample1 <- sample(n, n1)
  sample2 <- rep(1:n)[-sample1]

  rDEA_test <- rDEA[sample1]
  cDEA_test <- cDEA[sample2]

  mean_rDEA_test <- mean(rDEA_test)
  mean_cDEA_test <- mean(cDEA_test)
  sd_rDEA_test <-   sd(rDEA_test)
  sd_cDEA_test <-   sd(cDEA_test)
  sd2_rDEA_test <- (sd_rDEA_test)^2
  sd2_cDEA_test <- (sd_cDEA_test)^2

  B_rDEA_test <- rep(0, B)
  B_cDEA_test <- rep(0, B)
  B_Bias_r <- rep(0, B)
  B_Bias_c <- rep(0, B)

  for (i in 1:B) {
    n11 <- floor(n1/2)
    n12 <- n1 - n11
    sample11 <- sample(n1, n11)
    sample12 <- rep(1:n1)[-sample11]
    B_rDEA_test[i] <- (mean(rDEA_test[sample11]) +
                       mean(rDEA_test[sample12])  )/2
    B_Bias_r[i] <- (1/(2^k - 1))*(B_rDEA_test[i] - mean_rDEA_test)

    n21 <- floor(n2/2)
    n22 <- n2 - n21
    sample21 <- sample(n2, n21)
    sample22 <- rep(1:n2)[-sample21]
    B_cDEA_test[i] <- (mean(cDEA_test[sample21]) +
                       mean(cDEA_test[sample22])  )/2
    B_Bias_c[i] <- (1/(2^k - 1))*(B_cDEA_test[i] - mean_cDEA_test)
  }

  Bias_r <- mean(B_Bias_r)
  B_Bias_c <- mean(B_Bias_c)

  #t distribution with n degree of freedom
  t <- ((mean_rDEA_test - mean_cDEA_test) - (Bias_r - B_Bias_c) ) /
        sqrt((sd2_rDEA_test/n1) + (sd2_cDEA_test/n2))

  p_value <- pt(abs(t),df=n)
  }


  else {

    B_rDEA_test <- rep(0, B)
    B_cDEA_test <- rep(0, B)
    B_Bias_r <- rep(0, B)
    B_Bias_c <- rep(0, B)
    mean_rDEA_test <- rep(0, B)
    mean_cDEA_test <- rep(0, B)
    sd_rDEA_test <- rep(0, B)
    sd_cDEA_test <- rep(0, B)
    sd2_rDEA_test <- rep(0, B)
    sd2_cDEA_test <- rep(0, B)

    for (i in 1:B) {

      n1 <- floor(n/2)
      n2 <- n - n1

      sample1 <- sample(n, n1)
      sample2 <- rep(1:n)[-sample1]

      rDEA_test <- rDEA[sample1]
      cDEA_test <- cDEA[sample2]

      mean_rDEA_test[i] <- mean(rDEA_test)
      mean_cDEA_test[i] <- mean(cDEA_test)
      sd_rDEA_test[i] <-   sd(rDEA_test)
      sd_cDEA_test[i] <-   sd(cDEA_test)
      sd2_rDEA_test[i] <- (sd_rDEA_test[i])^2
      sd2_cDEA_test[i] <- (sd_cDEA_test[i])^2

      n11 <- floor(n1/2)
      n12 <- n1 - n11
      sample11 <- sample(n1, n11)
      sample12 <- rep(1:n1)[-sample11]
      B_rDEA_test[i] <- (mean(rDEA_test[sample11]) +
                           mean(rDEA_test[sample12])  )/2
      B_Bias_r[i] <- (1/(2^k - 1))*(B_rDEA_test[i] - mean_rDEA_test[i])

      n21 <- floor(n2/2)
      n22 <- n2 - n21
      sample21 <- sample(n2, n21)
      sample22 <- rep(1:n2)[-sample21]
      B_cDEA_test[i] <- (mean(cDEA_test[sample21]) +
                          mean(cDEA_test[sample22])  )/2
      B_Bias_c[i] <- (1/(2^k - 1))*(B_cDEA_test[i] - mean_cDEA_test[i])
    }
  }
    mean_rDEA_test <- mean(mean_rDEA_test)
    mean_cDEA_test <- mean(mean_cDEA_test)
    Bias_r <- mean(B_Bias_r)
    B_Bias_c <- mean(B_Bias_c)
    sd2_rDEA_test <- mean(sd2_rDEA_test)
    sd2_cDEA_test <- mean(sd2_cDEA_test)

    #t distribution with n degree of freedom
    t <- ((mean_rDEA_test - mean_cDEA_test) - (Bias_r - B_Bias_c) ) /
      sqrt((sd2_rDEA_test/n1) + (sd2_cDEA_test/n2))

    p_value <- pt(abs(t),df=n)

  if(print == TRUE) {
    if (alpha > p_value) {
      print("The null hypothesis of separability is rejected")
    }
    else {
      print("The null hypothesis of separability cannot be rejected")
    }
  }

  return(data.frame(p_value, t))

}
