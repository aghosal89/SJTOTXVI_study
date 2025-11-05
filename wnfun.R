
## This is the cost function for estimating the index parameter

## Inputs: 

# 1) et            - index parameter in polar coordinates
# 2) Y             - a list object of length n whose each element is an connectivity
#                    matrix of order mxm as the response 
# 3) linear_vars   - the names of the p covariates for the linear part of model 
#                    included in datosx
# 4) formula.      - a single character denoting the formula for the variables in the 
#                    linear part to be considered for the model, e.g. interactions, 
#                    higher order terms etc.
# 5) si_vars       - the names of the q covariates for the SI part of model 
#                    included in X.
# 7) sp            - the degree of polynomial considered for spline regression.
# 8) dfs           - degrees of freedom as an alternative to specifying the knots. 
# 9) X             - a data frame with n rows whose columns include the covariates 
#                    for the model.

## Output: the mean square prediction error. 

# to run this function 'lnr.R' file has to be sourced first
# linear_vars=NA, formula=NA, si_vars=NA

wnfun <- function(et=NULL, Y=NULL, x=NULL, h=NULL, kern= NULL, metric=NULL,
                  alpha= NULL) {
  
  th <- as.matrix(frechet:::pol2car(c(1, et))) # convert the polar coordinates to cartesian
  # compute the single index
  z <- x %*% th
  
  ft <- lnr(gl=Y, x= z, optns = list(metric= metric, alpha=alpha, bw= h, kernel= kern))
  #Yhat <- ft$fit
  
  #get Frobenius distance between the response and its prediction
  RSS <- sapply(1:length(ft$fit), function(i) distance_between_gl(ft$fit[[i]], ft$gl[[i]]))
  
  return(mean(RSS))
  
}


distance_between_gl <- function(GL_1,
                                GL_2,
                                euc = TRUE,
                                sqrt = FALSE,
                                proc = FALSE) {
  if ((euc == TRUE) && (sqrt == FALSE) && (proc == FALSE)) {
    return(norm(GL_1 - GL_2, type = 'F'))
  }
  
  if ((euc == FALSE) && (sqrt == TRUE) && (proc == FALSE)) {
    return(norm(shapes::rootmat(GL_1) - shapes::rootmat(GL_2), type = 'F'))
  }
  
  if ((euc == FALSE) && (sqrt == FALSE) && (proc == TRUE)) {
    return(procdist(
      shapes::rootmat(GL_1),
      shapes::rootmat(GL_2),
      type = 'sizeandshape',
      reflect = TRUE
    ))
  }
  
  else
    print('Error, one and only one metric must be TRUE')
}



is.laplacian <- function(L, tol = 1e-08) {
  L <- as.matrix(L)
  if (!matrixcalc::is.square.matrix(L)) {
    print("argument L is not a square matrix")
    return(FALSE)
  }
  
  if (!matrixcalc::is.symmetric.matrix(L)) {
    print("argument L is not a symmetric matrix")
    return(FALSE)
  }
  
  if (!is.numeric(L)) {
    print("argument L is not a numeric matrix")
    return(FALSE)
  }
  upper_diag <- L[upper.tri(L)]
  row_sums <- rowSums(L)
  if (any(row_sums < -tol)) {
    # print("Rows do no sum to zero")
    return(FALSE)
  }
  if (any(row_sums > tol)) {
    # print("Rows do no sum to zero")
    return(FALSE)
  }
  if (any(upper_diag > tol)) {
    # print("Off diagonal elements are not all non-positive")
    return(FALSE)
  }
  return(TRUE)
}



X <- mrn1_riskitt[,c(4,5,6)]
nsp <- 3

p <- ncol(X)
# specified spacing between staring points in each coordinate
spc <- pi/nsp  
# equally spaced starting points in polar coordinates
f <- lapply(1:(p - 1), function(j) seq(-pi/2 + spc/2, pi/2 - spc/2, by = spc)) 
# create grid of starting values
etaStart <- as.matrix(expand.grid(f))   

## To provide information about optimization as output
optInf <- list()

# provide criteria for termination of the algorithm
optim_optns <- list(factr = 1e11, maxit = 100)

WnMin <- rep(NA, nrow(etaStart))
etaMin <- matrix(NA, nrow = nrow(etaStart), ncol = p - 1)

# main optimization loop over starting values
for (k in 1:nrow(etaStart)) {
  
  WnOpt <- optim(par = etaStart[k, ], fn = WnCost, method = "L-BFGS-B",
                 lower = -pi/2, upper = pi/2, control = optim_optns,
                 X = X, tt = tt, Y = Y, kern = kern, h = h)
  
  optInf[[k]] <- WnOpt
  
  WnMin[k] <- WnOpt$value
  etaMin[k, ] <- WnOpt$par

}

