
##############################################################################
# Simulation for the paper "Network Regression with Graph Laplacians" #
# Within the set-up of the Partially Linear Frechet Regression (PL-FSI) model 
##############################################################################

#### Clear the environment first ####
rm(list= ls())

## Set working directory
setwd("/Users/aghosal/Documents/MatchedT16vsCon_MRTrix")

source("gnr.R")
source("lnr.R")
source("kerFctn.R")
source("functions_needed.R")

Q <- 1000  # number of simulations Monte-Carlo
nVec <- c(50, 100, 200, 500) # number of observations per simulation 
N <- 9 
m <- c(5, 5)
theta <- c(0.5, 0.2, 0.5) 
d <- c(m[1]*(m[1] - 1)/2, m[1]*m[2], m[2]*(m[2] - 1)/2)

alpha <- 0.5 # for the square-root metric as alternative to the Frobenius metric
# initialization of OSQP solver
W <- 2^32 # bound on the weights of edges in the graph
nConsts <- sum(m)^2 # number of constraints
l <- c(rep.int(0, sum(m) * (sum(m) + 1)/ 2), rep.int(-W, sum(m) * (sum(m) - 1) / 2))
u <- rep.int(0, nConsts)
q <- rep.int(0, sum(m)^2)
P <- diag(sum(m)^2)

consts <- matrix(0, nrow = nConsts, ncol = sum(m)^2)
k <- 0
for (i in 1:(sum(m) - 1)) {
  for (j in (i + 1):sum(m)) {
    k <- k + 1
    consts[k, (j - 1) * sum(m) + i] <- 1
    consts[k, (i - 1) * sum(m) + j] <- -1
  }
}

for (i in 1:sum(m)) {
  consts[k + i, ((i - 1) * sum(m) + 1):(i * sum(m))] <- rep(1, sum(m))
}
k <- k + sum(m)
for (i in 1:(sum(m) - 1)) {
  for (j in (i + 1):sum(m)) {
    k <- k + 1
    consts[k, (j - 1) * sum(m) + i] <- 1
  }
}

m <- osqp::osqp(P, q, consts, l, u, osqp::osqpSettings(verbose = FALSE))


################################################################################
# Partially Linear Frechet Single Index model for Graph Laplacian matrices
# endowed with the Frobenius/ square root metric
################################################################################

# this function computes the euclidean norm of the vector x:
e_norm <- function(x) {sqrt(sum(x^2))}

value_threshold <- function(x=NA, thresh=NA) {
  
  if( min(x)> thresh) {
    return("Y") 
  } else {
    return("N")
  }
}

icc_fun <- function(x1=NA, x2=NA, l=NA) {
  
  x1n<- x1[l]
  x2n<- x2[l]
  
  # compute the intra-class correlation
  c<- irr::icc(cbind(x1, x2),  model="oneway")
  
  return(c$value)
}  

perc_edge <- function(node_x=NULL, node_y=NULL, cm=NULL, cntrl_array=NULL, 
                      pos=NULL) {
  
  ay <- cm[node_x, node_y, cntrl_array]
  #aycdf <- ecdf(ay)
  
  l_1<- length(which(ay <= cm[node_x, node_y, pos] )) / length(cn)
  
  l_2<- length(which(ay > cm[node_x, node_y, pos] )) / length(cn)
  
  if(l_1 < l_2) {
    return(list(pvalue= l_1, tail="lower"))
  } else {
    return(list(pvalue= l_2, tail="upper"))
  }
}


# read the file names for the connectivity matrices
mrnc <- read.csv("MRN_control.csv", header = TRUE)
mrnp <- read.csv("MRN_patient.csv", header = TRUE)

node_labels <- read.csv("connectome_labels.csv", header =TRUE)

#head(node_labels)[,-3]

#tail(node_labels)[,-3]

# the dimension of the connectivity matrix 
d <- 379

# number of control participants
nc<- 59

# number of patient participants
np<- 335

# measurement of connectivity of the control participants 
# first time point
y_c_1 <- array(NA, dim=c(d, d, nc))

# second time point
y_c_2 <- array(NA, dim=c(d, d, nc))

# measurement of connectivity of patients 
# first time point
y_p_1 <- array(NA, dim=c(d, d, np))

# second time point
y_p_2 <- array(NA, dim=c(d, d, np))

for (i in 1:nc) {
  y_c_1[,,i] <-as.matrix(read.csv(paste(c("Con_", mrnc[i,], "-1_hcpmmp3",".csv"), collapse =""), header = TRUE)[,-1])
  y_c_2[,,i] <- as.matrix(read.csv(paste(c("Con_", mrnc[i,], "-3_hcpmmp3",".csv"), collapse =""), header = TRUE)[,-1])
}

for (i in 1:np) {
  y_p_1[,,i] <- as.matrix(read.csv(paste(c("Pat_", mrnp[i,], "-1_hcpmmp3", ".csv"), collapse =""), header = TRUE)[,-1])
  y_p_2[,,i] <- as.matrix(read.csv(paste(c("Pat_", mrnp[i,], "-2_hcpmmp3", ".csv"), collapse =""), header = TRUE)[,-1])
}

y_1 <- array(NA, dim= c(d, d, nc+np))

y_1[,,1:nc] <- y_c_1
y_1[,,c(60:394)] <- y_p_1

y_2 <- array(NA, dim= c(d,d, nc+np))

y_2[,,1:nc] <- y_c_2
y_2[,,c(60:394)] <- y_p_2

# then obtain the matrix of connectivity differences
y_diff <- y_2 - y_1

# Assign names to the 3rd dimension of the connectivity matrices
dimnames(y_1)[[3]] <- c(mrnc$MRN, mrnp$MRN)
dimnames(y_2)[[3]] <- c(mrnc$MRN, mrnp$MRN)
dimnames(y_diff)[[3]] <- c(mrnc$MRN, mrnp$MRN)

# Read the covariate data for analysis of the connectivity matrices
Xdat <- read.csv("Covariate_data_forregression.csv", header=TRUE)
Xdat$Sex <- as.factor(ifelse(Xdat$Sex=="2", "M", "F"))
Xdat <- within(Xdat, Risk <- relevel(as.factor(Risk), ref = "C"))
Xdat$Risk

Xdat$Group.ITT <- NA

Xdat$Group.ITT <- sapply(1:nrow(Xdat), function(i) {
  
  #TeX("SHR$\\leq 26$", output="character"))
  if(Xdat$Risk[i] == "LR" & Xdat$IT.total[i] <= 20 ) {"LR<20"}
  else if (Xdat$Risk[i] == "LR" & Xdat$IT.total[i] > 20) {"LR>=20"}
  else if (Xdat$Risk[i] == "SHR" & Xdat$IT.total[i] <= 26) {"SHR<27"}
  else if (Xdat$Risk[i] == "SHR" & Xdat$IT.total[i] > 26) {"SHR>=27"}
  else if (Xdat$Risk[i] == "C") {"C"}
  
})

# get the indexes corresponding to the risk categories in our study
cn <- which(Xdat$Risk=="C")
lr <- which(Xdat$Risk=="LR")
shr <- which(Xdat$Risk=="SHR")
l <- length(y_1[1,1,])

# get the MRNs for the patients in the LR and SHR groups 
mrn_lr <- Xdat$MRN[lr]
mrn_shr <- Xdat$MRN[shr]

# read the neurocognitive data first~
ndat <- read.csv("Neurocog_scores.csv", header = TRUE)

cnames_ndat <- colnames(ndat)

roi_diff21 <- read.csv("Edges_Connectivity_difference.csv")

wech_edges <- c(403,150,159,210,401,400,75,203,254,600,61,286,551,525,524,287,526,
                253,118,454,553,554)

roi_diff <- subset(roi_diff21, Edges %in% wech_edges & ExtremePatients >= 80)

temp <- ndat[, c(1,8,9)]  # choose the Wech working memory score and the MRN variable

# choose the neurocog scores that are non-empty
temp1 <- which(temp[,2] != "." & temp[,3] != ".")

mrn1 <- ndat$MRN[temp1]  # get MRNs for patients whose neurocog scores are non-empty

###############################
# Prepare the covariate data

#mrn1_riskitt <- cbind.data.frame(subset(Xdat, MRN %in% mrn1), 
#                                 Scores= as.numeric(temp[temp1,2]))

Xdat2 <- cbind.data.frame(subset(Xdat, MRN %in% temp$MRN), 
                          Wechsler= as.numeric(temp$wech.totalds.scs),
                          AudWM= as.numeric(temp$wj3.audwm.ss))

Xdat3 <- na.omit(Xdat2)

############################################################
# Choice of bandwidth range for the Total 16 covariate data
############################################################

# scale the covariate columns to standardize them
X_ctr <- cbind.data.frame(Xdat3[,c(1,2,3,6)], scale(Xdat3[, c(4,5,7,8)]))

# to find the bandwidth range for analysis~
h_max = max(apply(X_ctr[,c(5,6,7,8)], 1, e_norm))

# number of observations
n <- nrow(X_ctr)

## to find the lowest possible value for bandwidth h
metric_v_temp <- matrix(NA, n ,1) 

for (j in 1:n) {
  
  mv <-matrix(NA, n ,1)   
  
  for (i in 1:n) {
    if(i!=j)
      # computing the euclidean distance between rows the standardized X matrix
      mv[i] <- e_norm(X_ctr[j, c(5,6,7,8)] - X_ctr[i, c(5,6,7,8)]) 
  }
  metric_v_temp[j] = min(mv[-j]) # taking minimum distance 
}

h_min <- (min(metric_v_temp)+1)* 1.5

# the sequence of bandwidths to optimize over ~
h <- exp(seq(log(h_min), log(h_max), length.out= 10))


################################################
# Creating the response data for our regression
################################################

# gathering the edges attributed to working memory ~
rois_nr <- unique(c(roi_diff$Row, roi_diff$Column))
ydn <- y_diff[rois_nr, rois_nr, temp1]

gnreg <- gnr(gl= ydn, x= X_ctr[,c(5,6,7,8)]) 

lnreg <- lnr(gl= ydn, x= X_ctr[,5]) 

##############################################
# Estimate the Single-Index parameter for 

## To provide information about optimization as output
optInf <- list()

# provide criteria for termination of the algorithm
optim_optns <- list(factr = 1e11, maxit = 100)

WnMin <- rep(NA, nrow(etaStart))
etaMin <- matrix(NA, nrow = nrow(etaStart), ncol = p - 1)

# main optimization loop over starting values
for (k in 1:nrow(etaStart)) {
  WnOpt <- optim(par = etaStart[k, ], fn = wnfun, method = "L-BFGS-B",
                 lower = -pi/2, upper = pi/2, control = optim_optns,
                 Y=Y, x=x, h=2, kern= "gauss", metric="frobenius",
                 alpha= 1)
  
  optInf[[k]] <- WnOpt
  
  WnMin[k] <- WnOpt$value
  etaMin[k, ] <- WnOpt$par
}

# the op etimizer, i.e. thetaHat    
thetaHat <- frechet:::pol2car(c(1,etaMin[which.min(WnMin),]) )

optvalue <- min(WnMin)  # updated to find the minimized Wn in training set



# exploratory analysis to understand the relationship among the response and the 
# covariates

plot(mrn1_riskitt[,4], ydn[1,4, ])
abline(lm(ydn[1,4, ]~ mrn1_riskitt[,4]))


ggraph_plot <- function(GL, label_nodes=0, quantile_edges=0, delete_edges=0, 
                        tol=-10^-10, colour_node='slateblue', colour_line='springgreen')#delete nodes less then degree of label nodes,

#set anything belew quantile edges to this quantile, delete edges less then delete edges quantile
{GL[GL> tol]<-0
graph<-igraph_from_lapl(GL, TRUE)
igraph::V(graph)$node_label <- unname(ifelse(strength(graph)[igraph::V(graph)] >label_nodes, names(igraph::V(graph)), "")) 
igraph::V(graph)$node_size <- unname(ifelse(strength(graph)[igraph::V(graph)] > label_nodes, strength(graph), 0))

graph<-delete_edges(graph,E(graph)[E(graph)$weight[E(graph)]< quantile(E(graph)$weight,delete_edges)]) 
E(graph)$test<-E(graph)$weight
E(graph)$test[E(graph)$test<quantile(E(graph)$test,quantile_edges)]<-quantile(E(graph)$test,quantile_edges)

ggraph(graph, layout = 'linear', circular = TRUE) + 
  geom_edge_arc(edge_width=0.125, aes(alpha=(E(graph)$test)^1)) +
  geom_node_label(aes(label=node_label, size=node_size),
                  label.size=0, fill="#ffffff66", segment.colour=colour_line,
                  color=colour_node, repel=TRUE,  fontface="bold") +
  coord_fixed() +
  scale_size_area(trans="sqrt") +
  labs(title="", subtitle="") +
  theme_graph() +
  theme(legend.position="none")}


###################################################################################
# function to compute the statistical mean of an array of graph-Laplacian matrices

Mean_GL <- function(Array,
                    euc = TRUE,
                    sqrt = FALSE,
                    proc = FALSE, project=FALSE)
{
  if (is.list(Array) == TRUE) {
    rownms <- rownames(Array[[1]])
    AArray <-
      array(0, dim = c(dim(Array[[1]])[1], dim(Array[[1]])[2], length(Array)))
    for (i in 1:length(Array))
      AArray[, , i] <- as.matrix(Array[[i]])
    Array <- AArray
    rownames(Array) <-rownms
  }
  n <- dim(Array)[3]
  k <- dim(Array)[2]
  rownms <- rownames(Array[, , 1])
  if (euc == TRUE) {
    a <- Array[,,1]
    for ( i in 2:n)
      a<- a+Array[,,i]
    a<- a/n
    h <- Array[, , 1]
    h[] <- a
    rownames(h) <- rownms
    colnames(h) <- rownms
    return(h)
  }
  if (sqrt == TRUE) {
    #square root euclidean mean
    a <- array(rep(0, k * k * n), dim = c(k, k, n))
    for (i in 1:n)
      a[, , i] <-  shapes::rootmat(Array[, , i])
    b <- a[,,1]
    for (i in 2:n)
      b<- b+a[,,i]
    b<-b/n
    b <- b %*% t(b)
    h <- Array[, , 1]
    h[] <- b
    rownames(b) <- rownms
    colnames(b) <- rownms
    if (project==FALSE){
      return(b)
    }
    if(project==TRUE){
      b <-proj_sparse(b)
      rownames(b) <- rownms
      colnames(b) <- rownms
      return(b)
    }
  }
  if (proc == TRUE) {
    #size-and-shape mean
    c <-estSS(Array)
    h <- Array[, , 1]
    h[] <- c
    rownames(c) <- rownms
    colnames(c) <- rownms
    if (project==FALSE){
      return(c)
    }
    if(project==TRUE){
      c <- proj_sparse(c)
      rownames(c) <- rownms
      colnames(c) <- rownms
      return(c)
    }
  }
}

Mean_GL(gnreg$gl)
















# gather the connectivity brain edges involved in working memory
wmem_edges <- c(403,150,159,210,401,400,75,203,254,600,61,286,551,525,524,287,526,
                253,118,454,553,554)

gl <- y_diff[,, which(dimnames(y_diff)[[3]] %in% as.character(Xdat3$MRN))]

#################################################
# R code for the global network regression model
#################################################

gnr <- function(gl = NULL, x = NULL, xOut = NULL, optns = list()) {
  if (is.null(gl) | is.null(x)) {
    stop("requires the input of both gl and x")
  }
  if (is.null(optns$metric)) {
    optns$metric <- "frobenius"
  }
  if (!(optns$metric %in% c("frobenius", "power"))) {
    stop("metric choice not supported")
  }
  if (is.null(optns$alpha)) {
    optns$alpha <- 1
  }
  if (optns$alpha < 0) {
    stop("alpha must be non-negative")
  }
  if (is.null(optns$digits)) {
    optns$digits <- NA
  }
  if (!is.matrix(x)) {
    if (is.data.frame(x) | is.vector(x)) {
      x <- as.matrix(x)
    } else {
      stop("x must be a matrix or a data frame or a vector")
    }
  }
  n <- nrow(x) # number of observations
  p <- ncol(x) # number of covariates
  if (!is.list(gl)) {
    if (is.array(gl)) {
      gl <- lapply(seq(dim(gl)[3]), function(i) gl[, , i])
    } else {
      stop("gl must be a list or an array")
    }
  }
  if (length(gl) != n) {
    stop("the number of rows in x must be the same as the number of graph Laplacians in gl")
  }
  if (!is.null(xOut)) {
    if (!is.matrix(xOut)) {
      if (is.data.frame(xOut)) {
        xOut <- as.matrix(xOut)
      } else if (is.vector(xOut)) {
        if (p == 1) {
          xOut <- as.matrix(xOut)
        } else {
          xOut <- t(xOut)
        }
      } else {
        stop("xOut must be a matrix or a data frame or a vector")
      }
    }
    if (ncol(xOut) != p) {
      stop("x and xOut must have the same number of columns")
    }
    nOut <- nrow(xOut) # number of predictions
  } else {
    nOut <- 0
  }
  nodes <- colnames(gl[[1]])
  m <- ncol(gl[[1]]) # number of nodes
  if (is.null(nodes)) nodes <- 1:m
  glVec <- matrix(unlist(gl), ncol = m^2, byrow = TRUE) # n by m^2
  if (substr(optns$metric, 1, 1) == "p") {
    glAlpha <- lapply(gl, function(gli) {
      eigenDecom <- eigen(gli)
      Lambda <- pmax(Re(eigenDecom$values), 0) # exclude 0i
      U <- eigenDecom$vectors
      U %*% diag(Lambda^optns$alpha) %*% t(U)
    })
    glAlphaVec <- matrix(unlist(glAlpha), ncol = m^2, byrow = TRUE) # n by m^2
  }
  
  # initialization of OSQP solver
  W <- 2^32 # bound on the weights of edges in the graph
  nConsts <- m^2 # number of constraints
  l <- c(rep.int(0, m * (m + 1) / 2), rep.int(-W, m * (m - 1) / 2))
  u <- rep.int(0, nConsts)
  q <- rep.int(0, m^2)
  P <- diag(m^2)
  consts <- matrix(0, nrow = nConsts, ncol = m^2)
  k <- 0
  for (i in 1:(m - 1)) {
    for (j in (i + 1):m) {
      k <- k + 1
      consts[k, (j - 1) * m + i] <- 1
      consts[k, (i - 1) * m + j] <- -1
    }
  }
  for (i in 1:m) {
    consts[k + i, ((i - 1) * m + 1):(i * m)] <- rep(1, m)
  }
  k <- k + m
  for (i in 1:(m - 1)) {
    for (j in (i + 1):m) {
      k <- k + 1
      consts[k, (j - 1) * m + i] <- 1
    }
  }
  model <- osqp::osqp(P, q, consts, l, u, osqp::osqpSettings(verbose = FALSE))
  
  xMean <- colMeans(x)
  invVa <- solve(var(x) * (n - 1) / n)
  fit <- vector(mode = "list", length = n)
  residuals <- rep.int(0, n)
  wc <- t(apply(x, 1, function(xi) t(xi - xMean) %*% invVa)) # n by p
  totVa <- sum((scale(glVec, scale = FALSE))^2)
  if (nrow(wc) != n) wc <- t(wc) # for p=1
  if (substr(optns$metric, 1, 1) == "f") {
    for (i in 1:n) {
      w <- apply(wc, 1, function(wci) 1 + t(wci) %*% (x[i, ] - xMean))
      qNew <- apply(glVec, 2, weighted.mean, w) # m^2
      model$Update(q = -qNew)
      temp <- matrix(model$Solve()$x, ncol = m, dimnames = list(nodes, nodes))
      temp <- (temp + t(temp)) / 2 # symmetrize
      if (!is.na(optns$digits)) temp <- round(temp, digits = optns$digits) # round
      temp[temp > 0] <- 0 # off diagonal should be negative
      diag(temp) <- 0
      diag(temp) <- -colSums(temp)
      fit[[i]] <- temp
      residuals[i] <- sqrt(sum((gl[[i]] - temp)^2))
    }
    resVa <- sum(residuals^2)
    RSquare <- 1 - resVa / totVa
    AdjRSquare <- RSquare - (1 - RSquare) * p / (n - p - 1)
    if (nOut > 0) {
      predict <- vector(mode = "list", length = nOut)
      for (i in 1:nOut) {
        w <- apply(wc, 1, function(wci) 1 + t(wci) %*% (xOut[i, ] - xMean))
        qNew <- apply(glVec, 2, weighted.mean, w) # m^2
        model$Update(q = -qNew)
        temp <- matrix(model$Solve()$x, ncol = m, dimnames = list(nodes, nodes))
        temp <- (temp + t(temp)) / 2 # symmetrize
        if (!is.na(optns$digits)) temp <- round(temp, digits = optns$digits) # round
        temp[temp > 0] <- 0 # off diagonal should be negative
        diag(temp) <- 0
        diag(temp) <- -colSums(temp)
        predict[[i]] <- temp
      }
      res <- list(fit = fit, predict = predict, RSquare = RSquare, AdjRSquare = AdjRSquare, residuals = residuals, gl = gl, x = x, xOut = xOut, optns = optns)
    } else {
      res <- list(fit = fit, RSquare = RSquare, AdjRSquare = AdjRSquare, residuals = residuals, gl = gl, x = x, optns = optns)
    }
  } else if (substr(optns$metric, 1, 1) == "p") {
    # bAlpha <- matrix(colMeans(glAlphaVec), ncol = m)# m by m
    # eigenDecom <- eigen(bAlpha)
    # Lambda <- pmax(Re(eigenDecom$values), 0)# projection to M_m
    # U <- eigenDecom$vectors
    # qNew <- as.vector(U%*%diag(Lambda^(1/optns$alpha))%*%t(U))# inverse power
    # model$Update(q = -qNew)
    # omegaPlus <- matrix(model$Solve()$x, ncol = m, dimnames = list(nodes, nodes))
    # eigenDecom <- eigen(omegaPlus)
    # Lambda <- pmax(Re(eigenDecom$values), 0)
    # U <- eigenDecom$vectors
    # omegaPlusAlpha <- U%*%diag(Lambda^optns$alpha)%*%t(U)
    # totVa <- sum((t(glAlphaVec)-as.vector(omegaPlusAlpha))^2)# using Euclidean power metric
    for (i in 1:n) {
      w <- apply(wc, 1, function(wci) 1 + t(wci) %*% (x[i, ] - xMean))
      bAlpha <- matrix(apply(glAlphaVec, 2, weighted.mean, w), ncol = m) # m by m
      eigenDecom <- eigen(bAlpha)
      Lambda <- pmax(Re(eigenDecom$values), 0) # projection to M_m
      U <- eigenDecom$vectors
      qNew <- as.vector(U %*% diag(Lambda^(1 / optns$alpha)) %*% t(U)) # inverse power
      model$Update(q = -qNew)
      temp <- matrix(model$Solve()$x, ncol = m, dimnames = list(nodes, nodes))
      temp <- (temp + t(temp)) / 2 # symmetrize
      if (!is.na(optns$digits)) temp <- round(temp, digits = optns$digits) # round
      temp[temp > 0] <- 0 # off diagonal should be negative
      diag(temp) <- 0
      diag(temp) <- -colSums(temp)
      fit[[i]] <- temp
      residuals[i] <- sqrt(sum((gl[[i]] - temp)^2))
      # eigenDecom <- eigen(fit[[i]])
      # Lambda <- pmax(Re(eigenDecom$values), 0)
      # U <- eigenDecom$vectors
      # fitiAlpha <- U%*%diag(Lambda^optns$alpha)%*%t(U)
      # residuals[i] <- sqrt(sum((glAlpha[[i]]-fitiAlpha)^2))# using Euclidean power metric
    }
    resVa <- sum(residuals^2)
    RSquare <- 1 - resVa / totVa
    AdjRSquare <- RSquare - (1 - RSquare) * p / (n - p - 1)
    if (nOut > 0) {
      predict <- vector(mode = "list", length = nOut)
      for (i in 1:nOut) {
        w <- apply(wc, 1, function(wci) 1 + t(wci) %*% (xOut[i, ] - xMean))
        bAlpha <- matrix(apply(glAlphaVec, 2, weighted.mean, w), ncol = m) # m^2
        eigenDecom <- eigen(bAlpha)
        Lambda <- pmax(Re(eigenDecom$values), 0) # projection to M_m
        U <- eigenDecom$vectors
        qNew <- as.vector(U %*% diag(Lambda^(1 / optns$alpha)) %*% t(U)) # inverse power
        model$Update(q = -qNew)
        temp <- matrix(model$Solve()$x, ncol = m, dimnames = list(nodes, nodes))
        temp <- (temp + t(temp)) / 2 # symmetrize
        if (!is.na(optns$digits)) temp <- round(temp, digits = optns$digits) # round
        temp[temp > 0] <- 0 # off diagonal should be negative
        diag(temp) <- 0
        diag(temp) <- -colSums(temp)
        predict[[i]] <- temp
      }
      res <- list(fit = fit, predict = predict, RSquare = RSquare, AdjRSquare = AdjRSquare, residuals = residuals, gl = gl, x = x, xOut = xOut, optns = optns)
    } else {
      res <- list(fit = fit, RSquare = RSquare, AdjRSquare = AdjRSquare, residuals = residuals, gl = gl, x = x, optns = optns)
    }
  }
  class(res) <- "nr"
  res
}


#################################################
# R code for the local network regression model
#################################################

lnr <- function(gl = NULL, x = NULL, xOut = NULL, optns = list()) {
  if (is.null(gl) | is.null(x)) {
    stop("requires the input of both gl and x")
  }
  if (is.null(optns$metric)) {
    optns$metric <- "frobenius"
  }
  if (!(optns$metric %in% c("frobenius", "power"))) {
    stop("metric choice not supported")
  }
  if (is.null(optns$alpha)) {
    optns$alpha <- 1
  }
  if (optns$alpha < 0) {
    stop("alpha must be non-negative")
  }
  if (is.null(optns$kernel)) {
    optns$kernel <- "gauss"
  }
  if (is.null(optns$bw)) {
    optns$bw <- NA
  }
  if (is.null(optns$digits)) {
    optns$digits <- NA
  }
  if (!is.matrix(x)) {
    if (is.data.frame(x) | is.vector(x)) {
      x <- as.matrix(x)
    } else {
      stop("x must be a matrix or a data frame or a vector")
    }
  }
  n <- nrow(x) # number of observations
  p <- ncol(x) # number of covariates
  if (p > 2) {
    stop("local method is designed to work in low dimensional case (p is either 1 or 2)")
  }
  if (!is.na(sum(optns$bw))) {
    if (sum(optns$bw <= 0) > 0) {
      stop("bandwidth must be positive")
    }
    if (length(optns$bw) != p) {
      stop("dimension of bandwidth does not agree with x")
    }
  }
  if (!is.list(gl)) {
    if (is.array(gl)) {
      gl <- lapply(seq(dim(gl)[3]), function(i) gl[, , i])
    } else {
      stop("gl must be a list or an array")
    }
  }
  if (length(gl) != n) {
    stop("the number of rows in x must be the same as the number of graph Laplacians in gl")
  }
  if (!is.null(xOut)) {
    if (!is.matrix(xOut)) {
      if (is.data.frame(xOut)) {
        xOut <- as.matrix(xOut)
      } else if (is.vector(xOut)) {
        if (p == 1) {
          xOut <- as.matrix(xOut)
        } else {
          xOut <- t(xOut)
        }
      } else {
        stop("xOut must be a matrix or a data frame or a vector")
      }
    }
    if (ncol(xOut) != p) {
      stop("x and xOut must have the same number of columns")
    }
    nOut <- nrow(xOut) # number of predictions
  } else {
    nOut <- 0
  }
  nodes <- colnames(gl[[1]])
  m <- ncol(gl[[1]]) # number of nodes
  if (is.null(nodes)) nodes <- 1:m
  glVec <- matrix(unlist(gl), ncol = m^2, byrow = TRUE) # n by m^2
  if (substr(optns$metric, 1, 1) == "p") {
    glAlpha <- lapply(gl, function(gli) {
      eigenDecom <- eigen(gli)
      Lambda <- pmax(Re(eigenDecom$values), 0) # exclude 0i
      U <- eigenDecom$vectors
      U %*% diag(Lambda^optns$alpha) %*% t(U)
    })
    glAlphaVec <- matrix(unlist(glAlpha), ncol = m^2, byrow = TRUE) # n by m^2
  }
  
  # initialization of OSQP solver
  W <- 2^32 # bound on the weights of edges in the graph
  nConsts <- m^2 # number of constraints
  l <- c(rep.int(0, m * (m + 1) / 2), rep.int(-W, m * (m - 1) / 2))
  u <- rep.int(0, nConsts)
  q <- rep.int(0, m^2)
  P <- diag(m^2)
  consts <- matrix(0, nrow = nConsts, ncol = m^2)
  k <- 0
  for (i in 1:(m - 1)) {
    for (j in (i + 1):m) {
      k <- k + 1
      consts[k, (j - 1) * m + i] <- 1
      consts[k, (i - 1) * m + j] <- -1
    }
  }
  for (i in 1:m) {
    consts[k + i, ((i - 1) * m + 1):(i * m)] <- rep(1, m)
  }
  k <- k + m
  for (i in 1:(m - 1)) {
    for (j in (i + 1):m) {
      k <- k + 1
      consts[k, (j - 1) * m + i] <- 1
    }
  }
  model <- osqp::osqp(P, q, consts, l, u, osqp::osqpSettings(verbose = FALSE))
  
  # select kernel
  Kern <- kerFctn(optns$kernel)
  K <- function(x, h) {
    k <- 1
    for (i in 1:p) {
      k <- k * Kern(x[i] / h[i])
    }
    return(as.numeric(k))
  }
  
  # choose bandwidth by cross-validation
  if (is.na(sum(optns$bw))) {
    hs <- matrix(0, p, 20)
    for (l in 1:p) {
      hs[l, ] <- exp(seq(
        from = log(n^(-1 / (1 + p)) * (max(x[, l]) - min(x[, l])) / 10),
        to = log(5 * n^(-1 / (1 + p)) * (max(x[, l]) - min(x[, l]))),
        length.out = 20
      ))
    }
    cv <- array(0, 20^p)
    for (k in 0:(20^p - 1)) {
      h <- array(0, p)
      for (l in 1:p) {
        kl <- floor((k %% (20^l)) / (20^(l - 1))) + 1
        h[l] <- hs[l, kl]
      }
      for (j in 1:n) {
        a <- x[j, ]
        if (p > 1) {
          mu1 <- rowMeans(apply(as.matrix(x[-j, ]), 1, function(xi) K(xi - a, h) * (xi - a)))
          mu2 <- matrix(rowMeans(apply(as.matrix(x[-j, ]), 1, function(xi) K(xi - a, h) * ((xi - a) %*% t(xi - a)))), ncol = p)
        } else {
          mu1 <- mean(apply(as.matrix(x[-j, ]), 1, function(xi) K(xi - a, h) * (xi - a)))
          mu2 <- mean(apply(as.matrix(x[-j, ]), 1, function(xi) K(xi - a, h) * ((xi - a) %*% t(xi - a))))
        }
        skip <- FALSE
        tryCatch(solve(mu2), error = function(e) skip <<- TRUE)
        if (skip) {
          cv[k + 1] <- Inf
          break
        }
        wc <- t(mu1) %*% solve(mu2) # 1 by p
        w <- apply(as.matrix(x[-j, ]), 1, function(xi) {
          K(xi - a, h) * (1 - wc %*% (xi - a))
        }) # weight
        if (substr(optns$metric, 1, 1) == "f") {
          qNew <- apply(glVec[-j, ], 2, weighted.mean, w) # m^2
          model$Update(q = -qNew)
          fitj <- matrix(model$Solve()$x, ncol = m, dimnames = list(nodes, nodes))
          fitj <- (fitj + t(fitj)) / 2 # symmetrize
          if (!is.na(optns$digits)) fitj <- round(fitj, digits = optns$digits) # round
          fitj[fitj > 0] <- 0 # off diagonal should be negative
          diag(fitj) <- 0
          diag(fitj) <- -colSums(fitj)
          cv[k + 1] <- cv[k + 1] + sum((gl[[j]] - fitj)^2) / n
        } else if (substr(optns$metric, 1, 1) == "p") {
          bAlpha <- matrix(apply(glAlphaVec[-j, ], 2, weighted.mean, w), ncol = m) # m by m
          eigenDecom <- eigen(bAlpha)
          Lambda <- pmax(Re(eigenDecom$values), 0) # projection to M_m
          U <- eigenDecom$vectors
          qNew <- as.vector(U %*% diag(Lambda^(1 / optns$alpha)) %*% t(U)) # inverse power
          model$Update(q = -qNew)
          fitj <- matrix(model$Solve()$x, ncol = m, dimnames = list(nodes, nodes))
          fitj <- (fitj + t(fitj)) / 2 # symmetrize
          if (!is.na(optns$digits)) fitj <- round(fitj, digits = optns$digits) # round
          fitj[fitj > 0] <- 0 # off diagonal should be negative
          diag(fitj) <- 0
          diag(fitj) <- -colSums(fitj)
          cv[k + 1] <- cv[k + 1] + sum((gl[[j]] - fitj)^2) / n
          # eigenDecom <- eigen(fitj)
          # Lambda <- pmax(Re(eigenDecom$values), 0)
          # U <- eigenDecom$vectors
          # fitjAlpha <- U%*%diag(Lambda^optns$alpha)%*%t(U)
          # cv[k+1] <- cv[k+1] + sum((glAlpha[[j]]-fitjAlpha)^2)/n# using Euclidean power metric
        }
      }
    }
    bwi <- which.min(cv)
    optns$bw <- array(0, p)
    for (l in 1:p) {
      kl <- floor((bwi %% (20^l)) / (20^(l - 1))) + 1
      optns$bw[l] <- hs[l, kl]
    }
  }
  
  fit <- vector(mode = "list", length = n)
  residuals <- rep.int(0, n)
  if (substr(optns$metric, 1, 1) == "f") {
    for (i in 1:n) {
      a <- x[i, ]
      if (p > 1) {
        mu1 <- rowMeans(apply(x, 1, function(xi) K(xi - a, optns$bw) * (xi - a)))
        mu2 <- matrix(rowMeans(apply(x, 1, function(xi) K(xi - a, optns$bw) * ((xi - a) %*% t(xi - a)))), ncol = p)
      } else {
        mu1 <- mean(apply(x, 1, function(xi) K(xi - a, optns$bw) * (xi - a)))
        mu2 <- mean(apply(x, 1, function(xi) K(xi - a, optns$bw) * ((xi - a) %*% t(xi - a))))
      }
      wc <- t(mu1) %*% solve(mu2) # 1 by p
      w <- apply(x, 1, function(xi) {
        K(xi - a, optns$bw) * (1 - wc %*% (xi - a))
      }) # weight
      qNew <- apply(glVec, 2, weighted.mean, w) # m^2
      model$Update(q = -qNew)
      temp <- matrix(model$Solve()$x, ncol = m, dimnames = list(nodes, nodes))
      temp <- (temp + t(temp)) / 2 # symmetrize
      if (!is.na(optns$digits)) temp <- round(temp, digits = optns$digits) # round
      temp[temp > 0] <- 0 # off diagonal should be negative
      diag(temp) <- 0
      diag(temp) <- -colSums(temp)
      fit[[i]] <- temp
      residuals[i] <- sqrt(sum((gl[[i]] - temp)^2))
    }
    if (nOut > 0) {
      predict <- vector(mode = "list", length = nOut)
      for (i in 1:nOut) {
        a <- xOut[i, ]
        if (p > 1) {
          mu1 <- rowMeans(apply(x, 1, function(xi) K(xi - a, optns$bw) * (xi - a)))
          mu2 <- matrix(rowMeans(apply(x, 1, function(xi) K(xi - a, optns$bw) * ((xi - a) %*% t(xi - a)))), ncol = p)
        } else {
          mu1 <- mean(apply(x, 1, function(xi) K(xi - a, optns$bw) * (xi - a)))
          mu2 <- mean(apply(x, 1, function(xi) K(xi - a, optns$bw) * ((xi - a) %*% t(xi - a))))
        }
        wc <- t(mu1) %*% solve(mu2) # 1 by p
        w <- apply(x, 1, function(xi) {
          K(xi - a, optns$bw) * (1 - wc %*% (xi - a))
        }) # weight
        qNew <- apply(glVec, 2, weighted.mean, w) # m^2
        model$Update(q = -qNew)
        temp <- matrix(model$Solve()$x, ncol = m, dimnames = list(nodes, nodes))
        temp <- (temp + t(temp)) / 2 # symmetrize
        if (!is.na(optns$digits)) temp <- round(temp, digits = optns$digits) # round
        temp[temp > 0] <- 0 # off diagonal should be negative
        diag(temp) <- 0
        diag(temp) <- -colSums(temp)
        predict[[i]] <- temp
      }
      res <- list(fit = fit, predict = predict, residuals = residuals, gl = gl, x = x, xOut = xOut, optns = optns)
    } else {
      res <- list(fit = fit, residuals = residuals, gl = gl, x = x, optns = optns)
    }
  } else if (substr(optns$metric, 1, 1) == "p") {
    for (i in 1:n) {
      a <- x[i, ]
      if (p > 1) {
        mu1 <- rowMeans(apply(x, 1, function(xi) K(xi - a, optns$bw) * (xi - a)))
        mu2 <- matrix(rowMeans(apply(x, 1, function(xi) K(xi - a, optns$bw) * ((xi - a) %*% t(xi - a)))), ncol = p)
      } else {
        mu1 <- mean(apply(x, 1, function(xi) K(xi - a, optns$bw) * (xi - a)))
        mu2 <- mean(apply(x, 1, function(xi) K(xi - a, optns$bw) * ((xi - a) %*% t(xi - a))))
      }
      wc <- t(mu1) %*% solve(mu2) # 1 by p
      w <- apply(x, 1, function(xi) {
        K(xi - a, optns$bw) * (1 - wc %*% (xi - a))
      }) # weight
      bAlpha <- matrix(apply(glAlphaVec, 2, weighted.mean, w), ncol = m) # m by m
      eigenDecom <- eigen(bAlpha)
      Lambda <- pmax(Re(eigenDecom$values), 0) # projection to M_m
      U <- eigenDecom$vectors
      qNew <- as.vector(U %*% diag(Lambda^(1 / optns$alpha)) %*% t(U)) # inverse power
      model$Update(q = -qNew)
      temp <- matrix(model$Solve()$x, ncol = m, dimnames = list(nodes, nodes))
      temp <- (temp + t(temp)) / 2 # symmetrize
      if (!is.na(optns$digits)) temp <- round(temp, digits = optns$digits) # round
      temp[temp > 0] <- 0 # off diagonal should be negative
      diag(temp) <- 0
      diag(temp) <- -colSums(temp)
      fit[[i]] <- temp
      residuals[i] <- sqrt(sum((gl[[i]] - temp)^2))
      # eigenDecom <- eigen(fit[[i]])
      # Lambda <- pmax(Re(eigenDecom$values), 0)
      # U <- eigenDecom$vectors
      # fitiAlpha <- U%*%diag(Lambda^optns$alpha)%*%t(U)
      # residuals[i] <- sqrt(sum((glAlpha[[i]]-fitiAlpha)^2))# using Euclidean power metric
    }
    if (nOut > 0) {
      predict <- vector(mode = "list", length = nOut)
      for (i in 1:nOut) {
        a <- xOut[i, ]
        if (p > 1) {
          mu1 <- rowMeans(apply(x, 1, function(xi) K(xi - a, optns$bw) * (xi - a)))
          mu2 <- matrix(rowMeans(apply(x, 1, function(xi) K(xi - a, optns$bw) * ((xi - a) %*% t(xi - a)))), ncol = p)
        } else {
          mu1 <- mean(apply(x, 1, function(xi) K(xi - a, optns$bw) * (xi - a)))
          mu2 <- mean(apply(x, 1, function(xi) K(xi - a, optns$bw) * ((xi - a) %*% t(xi - a))))
        }
        wc <- t(mu1) %*% solve(mu2) # 1 by p
        w <- apply(x, 1, function(xi) {
          K(xi - a, optns$bw) * (1 - wc %*% (xi - a))
        }) # weight
        bAlpha <- matrix(apply(glAlphaVec, 2, weighted.mean, w), ncol = m) # m^2
        eigenDecom <- eigen(bAlpha)
        Lambda <- pmax(Re(eigenDecom$values), 0) # projection to M_m
        U <- eigenDecom$vectors
        qNew <- as.vector(U %*% diag(Lambda^(1 / optns$alpha)) %*% t(U)) # inverse power
        model$Update(q = -qNew)
        temp <- matrix(model$Solve()$x, ncol = m, dimnames = list(nodes, nodes))
        temp <- (temp + t(temp)) / 2 # symmetrize
        if (!is.na(optns$digits)) temp <- round(temp, digits = optns$digits) # round
        temp[temp > 0] <- 0 # off diagonal should be negative
        diag(temp) <- 0
        diag(temp) <- -colSums(temp)
        predict[[i]] <- temp
      }
      res <- list(fit = fit, predict = predict, residuals = residuals, gl = gl, x = x, xOut = xOut, optns = optns)
    } else {
      res <- list(fit = fit, residuals = residuals, gl = gl, x = x, optns = optns)
    }
  }
  class(res) <- "nr"
  res
}


