
# Get the items in the global R environment
#ls(envir = .GlobalEnv)

#### Clear the environment first ####
rm(list= ls())

# read the libraries of functions
library('gdata')
library('sampling')
library('matrixcalc')
library("ggplot2")
library('reshape')
library('viridis')
library('MASS')
library('ggpubr')
library('transport')
library('latex2exp')
library('patchwork')
library('irr')
library('statip')
library('cubature')

# Study of Brain Connectivity (TOTXVI, St. Jude Children's Research Hospital)
# author: Aritra Ghosal, St Jude Children's Research Hospital
# date: 03/19/2024

# functions to read before running the codes
#boxcox_fun <- function(x, lambda=NA) {
  
#  if(lambda==0) {
#    y<- lox(x+1)
#  } else {
#    y<- ((x^(lambda))-1)/(lambda)
#  }
  
#  return(y)
#} 

# The following is the function to run the regression for every brain edge

# Input items:

# x      : 
# xin    :
# covs   :
# formula :
# optns   :

#Nreg_fun <- function(x=NULL, xin=NULL, covs=NULL, formula=NULL, optn= NULL) {
  
#  reg_temp <- data.frame(Y= x, xin[,covs])
#  if(!is.na(formula)) {
#    ftemp.1 <- (lm(formula, data = reg_temp))
#  } else {
#    ftemp.1 <- (lm(Y~., data = reg_temp))
#  }
#  
#  ftemp <- summary(ftemp.1)
#  
#  if(optn == "Coeffs") {
#    return(ftemp.1$coefficients)
#  } else if(optn == "p-values") {
#    return(ftemp$coefficients[,4])
#  } else if(optn == "Residuals") {
#    return(ftemp$residuals)
#  } else if(optn== "Fitted") {
#    return(ftemp.1$fitted.values)
#  } else if(optn =="Leverages") {
#    return(hatvalues(ftemp.1)) 
#  } else if(optn == "Cook's D") {
#    return(cooks.distance(ftemp.1))
#  }
#}

# function to compute the minimum connectivity value in the edge

value_threshold <- function(x=NA, thresh=NA) {
  
  if( min(x)> thresh) {
    return("Y") 
  } else {
    return("N")
  }
}


hdist<- function (x, y, lower = 0, upper = Inf) {
  fx <- densityfun(x)
  fy <- densityfun(y)
  
  g <- function(z) (fx(z)^0.5 - fy(z)^0.5)^2
  h2 <- cubature::adaptIntegrate(g, lowerLimit= lower, upperLimit= upper)$integral/2
  
  return(sqrt(h2))
}


# function to compute the Intra-Class Correlation of the connectivity edges

icc_fun <- function(x1=NA, x2=NA, l=NA) {
  
  x1n<- x1[l]
  x2n<- x2[l]
  
  # compute the intra-class correlation
  c<- irr::icc(cbind(x1, x2),  model="oneway")
  
  return(c$value)
}  

###############################################################
# The following is the description of the function 'perc_edge' 
###############################################################

## Inputs

# node_x      = the x-coordinate of the brain ROI. 
# node_y      = the y-coordinate of the brain ROI. 
# cm          = the 379 x 379 connectivity matrix from either the time 1 or 2.
# cntrl_array = the array position vector for the control group.
# pos         = the position of the patient in the array of risk groups.

## Output Values

# pvalue = the probability of observing a more extreme observation on either the 
#          upper or the lower tail of the distribution of control weights.

# tail   = A bivariate indicator with values 'upper' or 'lower', denoting if the 
#          weight whose p-value is computed is situated at either the upper or the
#          lower tail of the control weight distribution.

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


############################################################
# reading the connectivity data for Object data regression:
############################################################

# set working directory
setwd("/Users/aghosal/Documents/MatchedT16vsCon_MRTrix")

# read the file names for the connectivity matrices
mrnc <- read.csv("MRN_control.csv", header = TRUE)
mrnp <- read.csv("MRN_patient.csv", header = TRUE)

# read the dataset containing the labels, areas of the ROIs in our study
node_labels <- read.csv("connectome_labels.csv", header =TRUE)

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

for(i in 1:nrow(Xdat)) {
  
  if(Xdat$Risk[i] == "LR" & Xdat$IT.total[i] <= 19){
    Xdat$Group.ITT[i] <- "LR(LD)"
  }
  else if (Xdat$Risk[i] == "LR" & Xdat$IT.total[i] > 19) {
    Xdat$Group.ITT[i] <- "LR(HD)"
  }
  else if (Xdat$Risk[i] == "SHR" & Xdat$IT.total[i] <= 26) {
    Xdat$Group.ITT[i] <- "SHR(LD)"
  }
  else if (Xdat$Risk[i] == "SHR" & Xdat$IT.total[i] > 26) {
    Xdat$Group.ITT[i] <- "SHR(HD)"
  }
  else {
    Xdat$Group.ITT[i] <- "C"
  }
}

Xdat$Group.ITT <- factor(Xdat$Group.ITT, levels=c("C", "LR(LD)", "LR(HD)", 
                                                  "SHR(LD)", "SHR(HD)"))

# this file contains the neuro-cognitive scores of the patients
Xdat2 <- read.csv("Xdat2.csv", header = T)
X3 <- merge(Xdat, Xdat2, by="MRN")
X3 <- X3[,-c(2:12)]

tsq <- 0.1

# get the indexes corresponding to the risk categories in our study
cn <- which(Xdat$Risk=="C")
lr <- which(Xdat$Risk=="LR")
shr <- which(Xdat$Risk=="SHR")

lrld21 <- which(Xdat$Risk=="LR" & Xdat$IT.total <= 19)
lrhd21 <- which(Xdat$Risk=="LR" & Xdat$IT.total > 19)
shrld27 <- which(Xdat$Risk=="SHR" & Xdat$IT.total <= 26)
shrhd27 <- which(Xdat$Risk=="SHR" & Xdat$IT.total > 26)

lp <- length(y_1[1,1,])

# get the MRNs for the patients in the LR and SHR groups 
mrn_lr <- Xdat$MRN[lr]
mrn_shr <- Xdat$MRN[shr]

## 1st, 3rd quartiles and median of every edge
#y_2_1stq <- apply(y_2[,,cn], c(1,2), FUN = quantile, prob= 0.25)
#y_2_3rdq <- apply(y_2[,,cn], c(1,2), FUN = quantile, prob= 0.75)
#y_2_med <- apply(y_2[,,cn], c(1,2), FUN = median)
#y_2_min <- apply(y_2[,,], c(1,2), FUN = min)

# we shortlist the edges with minimum value at time 2, since the patients have recovered 
# more and are more stable at that time.

value_count <- apply(y_2[,,cn], c(1,2), FUN =value_threshold, thresh =tsq)

value_count1 <- data.frame(value= c(value_count), X1 = rep(seq(1, 379), 379), 
                           X2= rep(seq(1,379), each= 379) ) 

vc <- subset(value_count1, value=="Y" & X1> X2)[,c("X1","X2")]

vc_symm <- subset(value_count1, value=="Y")[,c("X1","X2")]

vc1_tab <- subset(value_count1, X1> X2)
a <- table(vc1_tab$value)
a

Risk_id <- rep(NA, len= lp)
Risk_id[cn]<- "C"
Risk_id[lr]<- "LR"
Risk_id[shr]<- "SHR"

# As a starting example first we consider the weight distribution between the regions: 
# Left Accumbens (ROI $8$) & Left Ventral Diencephalon (ROI $9$), we consider the weight 
# distribution both at the start and at the 120 weeks.

#We also study the boxplot of these weights at $(8, 9)$ in the same order as follows:
  
par(mfrow=c(1,2))

dat_t1 <- data.frame(Risk=Risk_id, Weight=y_1[8,9,])

p <- ggplot(dat_t1, aes(x= Risk, y=Weight, fill = Risk))
p1<- p + geom_boxplot() +
  labs(y = "Connectivity Strength", x = "", title = "Brain connectivity at time 1", 
       fill= "Correlation") +
  theme(axis.text.x=element_text(size=12),
        axis.title.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.y=element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        title = element_text(size=12),
        strip.text = element_text(size=12))+
  
geom_jitter() +
  guides(guide_colorbar="none")

#ggsave(file="Rplot_connectivity_weights_8_9_t1.png", height =4.5, width =4.5,
#       units ="in", dpi= 300)

dat_t2 <- data.frame(Risk=Risk_id, Weight=y_2[8,9,])

p <- ggplot(dat_t2, aes(x= Risk, y=Weight, fill = Risk))
p2<- p + geom_boxplot() +
  labs(y = "Connectivity Strength",x = "",title = "Brain connectivity at time 2") +
  theme(axis.text.x=element_text(size=12),
        axis.title.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.y=element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        title = element_text(size=12),
        strip.text = element_text(size=12)) +
  
  geom_jitter()

p1+p2

#ggsave(file="Rplot_connectivity_weights_8_9_t2.png", height =4.5, width =4.5,
#       units ="in", dpi= 300)


############################################################################
# Plot of the connectivity distributions:

edens_plotdat1 <- cbind.data.frame(Risk= rep(Risk_id, 3), 
                                  Connectivity= c(y_1[208, 210,], y_2[208, 210,],
                                                  y_diff[208,210,]),
                                  Time= rep(c("RE-1", "W-120", "Difference"), 
                                            each= length(Risk_id)) )

edens_plotdat2 <- cbind.data.frame(Risk= rep(Risk_id, 3), 
                                   Connectivity= c(y_1[3, 5,], y_2[3, 5,],
                                                   y_diff[3, 5,]),
                                   Time= rep(c("RE-1", "W-120", "Difference"), 
                                             each= length(Risk_id)) )

edens_plotdat3 <- cbind.data.frame(Risk= rep(Risk_id, 3), 
                                   Connectivity= c(y_1[228, 359,], y_2[228, 359,],
                                                   y_diff[228, 359,]),
                                   Time= rep(c("RE-1", "W-120", "Difference"), 
                                             each= length(Risk_id)) )

edens_plotdat <- data.frame(rbind.data.frame(edens_plotdat1, edens_plotdat2, edens_plotdat3),
                            Edge= rep(c("RFrontalEF vs RArea55b","LCaudate vs LPallidum","RArea7m vs RArea31pd"),
                                      each= nrow(edens_plotdat1)))

edens_plotdat$Time <- factor(edens_plotdat$Time, levels=c("RE-1", "W-120", 
                                                    "Difference"))

ggplot(edens_plotdat, aes(x= Connectivity, fill= Risk)) +
  geom_density(alpha = 0.6)+
  facet_grid( Edge ~ Time, scales = "free_x")+
  scale_fill_hue()+
  #scale_fill_manual(values = cols)+
  labs(title ="Connectivity distributions of the Risk Groups",
       fill= "Risk groups")+
  #geom_rug(aes(color= Risk, alpha = 0.1))+
  theme_bw()

ggsave(file="Rplot_connectivity_densities_t1t2diff.png", height =8, width =9.5,
       units ="in", dpi= 300)

df <- data.frame(Connectivity= c(y_diff[y_cor_control2$Row[1], y_cor_control2$Column[1],],
                                y_diff[y_cor_control2$Row[2], y_cor_control2$Column[2],],
                                y_diff[y_cor_control2$Row[3], y_cor_control2$Column[3],],
                                y_diff[y_cor_control2$Row[4], y_cor_control2$Column[4],],
                                y_diff[y_cor_control2$Row[5], y_cor_control2$Column[5],],
                                y_diff[y_cor_control2$Row[6], y_cor_control2$Column[6],],
                                y_diff[y_cor_control2$Row[7], y_cor_control2$Column[7],],
                                y_diff[y_cor_control2$Row[8], y_cor_control2$Column[8],],
                                y_diff[y_cor_control2$Row[9], y_cor_control2$Column[9],]), 
                Edge= c("N1","N2","N3","N4","N5","N6","N7","N8","N9"),
                Groups= rep(Xdat$Risk, 9) )

df$membership <- NA

for(i in 1:nrow(df)) {
  
  if(df$Groups[i] == "C"){
    df$membership[i] <- 1
  }
  else if (df$Groups[i] == "LR") {
    df$membership[i] <- 2
  }
  else {
    df$membership[i] <- 3
  }
}


qnts<- quantile(df$Connectivity[1:59], probs = c(0.05, 0.95) )

# Helper variables
limits <- range(df$Connectivity)
step   <- diff(limits) * 0.05
size   <- 0.45 * step

ggplot(df, aes(x = Connectivity, fill = Groups)) + 
  geom_density(aes(alpha = 0.3))+
  
  geom_segment(
    aes(colour = Groups, xend = Connectivity,
      y    = 0 - as.numeric(membership) * step + size,
      yend = 0 - as.numeric(membership) * step - size
    )
  ) +
  facet_wrap(~ Edge)+
  theme_bw()+
  geom_vline(xintercept = qnts[1], col="black", linetype= "dashed")+
  geom_vline(xintercept = qnts[2], col="black", linetype= "dashed")+
  scale_y_continuous(limits = c(-0.3, 3))
  


##########################################################################
# Next, we consider the weight distribution between the regions: Right Area 6m Anterior & Right Area 6 Anterior, here as well we consider the weights distribution both at the start and at the 120 weeks.

par(mfrow=c(1,2))

dat_t1 <- data.frame(Risk=Risk_id, Weight=y_1[294,242,])

p <- ggplot(dat_t1, aes(x= Risk, y=Weight))
p + geom_boxplot() +
  labs(y = "Connectivity Strength", x = "", title = "Brain connectivity at time 1", 
       fill = "Correlation") +
  theme(axis.text.x=element_text(size=12),
        axis.title.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.y=element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        title = element_text(size=12),
        strip.text = element_text(size=12))

#ggsave(file="Rplot_connectivity_weights_294_242_t1.png", height =4.5, width =4.5,
#       units ="in", dpi= 300)

dat_t2<- data.frame(Risk=Risk_id, Weight=y_2[294,242,])

p <- ggplot(dat_t2, aes(x= Risk, y=Weight))
p + geom_boxplot() +
  labs(y = "Connectivity Strength", x = "", title = "Brain connectivity at time 2", 
       fill= "Correlation") +
  theme(axis.text.x=element_text(size=12),
        axis.title.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.y=element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        title = element_text(size=12),
        strip.text = element_text(size=12))

#ggsave(file="Rplot_connectivity_weights_294_242_t2.png",, height =4.5, width =4.5,
#       units ="in", dpi= 300)

age_min <- min(Xdat$Age)
age_max <- max(Xdat$Age)

dn_age_c<- density(x= Xdat$Age[cn], from =age_min, to= age_max, n=512)
dn_age_lr<- density(Xdat$Age[lr], from =age_min, to= age_max, n=512)
dn_age_shr<- density(x=Xdat$Age[shr], from =age_min, to= age_max, n=512)

p1<- ggplot(Xdat, aes(y= Age, x=Group.ITT, fill= Group.ITT)) +
  geom_boxplot(alpha = 0.7) +
  labs(title ="Age distribution of the subjects")+
  scale_fill_brewer(palette = "Set2")+
  guides(x = guide_axis(title = NULL, position = "NULL"))+
  labs(fill= "Groups")
  #geom_jitter(size=0.5)

p2<- ggplot(Xdat[-c(1:59),], aes(x= IT.total, fill= Risk)) +
  geom_density(alpha = 0.4)+
  labs(title ="ITT distribution of the patients", x= "ITT dosage")+
  geom_vline(xintercept = 19, col="black", linetype= "dashed")+
  geom_vline(xintercept = 26, col="black", linetype= "dashed")+
  theme_bw()

p1 + p2

ggsave(file="Rplot_Age_ITT_distributions_Risk.png", height =4, width =9,
       units ="in", dpi= 300)


## What motivates us to gather the extreme edges and people in this study:
par(mfrow= c(1,1))

# connectivity between Right Caudate and Right Thalamus

# Time = 1
png(filename = "Rplot_extrmedges_294_242_t1.png", width = 420, height = 450,
    units = "px", pointsize = 12)
qntl <- quantile(y_1[294,242, cn], probs = c(0.05, 0.95))

plot(density(y_1[294, 242, cn]), ylim=c(-0.1, 1.1), xlim= c(-1, 5),
     main = "R Area 6A vs R Area 6mA", xlab ="Connectivity, Time = 1",
     cex.lab=1.4, cex.axis=1.3, cex.main=1.4)
abline(v=qntl[2], col= "red")
abline(v=qntl[1], col= "red")
points(y_1[294, 242, shr], rep(0, length(shr)), col="blue", pch= 16, cex=0.6)
points(y_1[294, 242, lr], rep(-0.055, length(lr)), col="green3", pch= 16,
       cex=0.6)
dev.off()

# Time = 2
png(filename = "Rplot_extrmedges_294_242_t2.png", width = 420, height = 450,
    units = "px", pointsize = 12)
qntl <- quantile(y_2[294,242, cn], probs = c(0.05, 0.95))

plot(density(y_2[294, 242, cn]), ylim=c(-.09,1.1), xlim= c(-1, 5),
     main = "R Area 6A vs R Area 6mA", xlab ="Connectivity, Time = 2",
     cex.lab=1.4, cex.axis=1.3, cex.main=1.4)
abline(v=qntl[2], col= "red")
abline(v=qntl[1], col= "red")
points(y_2[294, 242, shr], rep(-0.045, length(shr)), col="blue", pch= 16, cex=0.6)
points(y_2[294, 242, lr], rep(-0.055, length(lr)), col="green3", pch= 16,
       cex=0.6)
dev.off()

# Time difference
png(filename = "Rplot_extrmedges_294_242_t21.png", width = 420, height = 450,
    units = "px", pointsize = 12)
y_diff1 <- y_2[294, 242,] - y_1[294, 242,]
qntl <- quantile(y_diff1[cn], probs = c(0.05, 0.95))

plot(density(y_diff1[cn]), ylim=c(-.085,1.1), xlim =c(-2.5, 2.5),
     main = "R Area 6A vs R Area 6mA", xlab ="Connectivity difference",
     cex.lab=1.4, cex.axis=1.3, cex.main=1.4)
abline(v=qntl[2], col= "red")
abline(v=qntl[1], col= "red")
points(y_diff1[shr], rep(-0.04, length(shr)), col="blue", pch= 16, cex=0.6)
points(y_diff1[lr], rep(-0.055, length(lr)), col="green3", pch= 16,
       cex=0.6)
dev.off()


# connectivity between Left Caudate and Left Thalamus

# Time = 1
png(filename = "Rplot_extrmedges_8_9_t1.png", width = 420, height = 450,
    units = "px", pointsize = 12)
qntl <- quantile(y_1[8,9, cn], probs = c(0.05, 0.95))

plot(density(y_1[8,9, 1:59]), xlim=c(-1, 4.4), ylim=c(-.1, 1.1),
     main = "L Accumbens vs L V Diencephalon", xlab ="Connectivity, time=1",
     cex.lab=1.4, cex.axis=1.3, cex.main=1.4)
abline(v=qntl[2], col= "red")
abline(v=qntl[1], col= "red")
points(y_1[8,9, shr], rep(0, length(shr)), col="blue", pch= 16, cex=0.6)
points(y_1[8,9, lr], rep(-0.055, length(lr)), col="green3", pch= 16,
       cex=0.6)
dev.off()

# Time = 2
png(filename ="Rplot_extrmedges_8_9_t2.png", width = 420, height = 450,
    units = "px", pointsize = 12)
qntl <- quantile(y_2[8, 9, cn], probs = c(0.05, 0.95))

plot(density(y_2[8, 9, 1:59]), ylim=c(-.1, 1.1), xlim= c(-1, 4.4),
     main = "L Accumbens vs L V Diencephalon", xlab ="Connectivity, time=2",
     cex.lab=1.4, cex.axis=1.3, cex.main=1.4)
abline(v=qntl[2], col= "red")
abline(v=qntl[1], col= "red")
points(y_2[8, 9, shr], rep(0, length(shr)), col="blue", pch= 16, cex=0.6)
points(y_2[8, 9, lr], rep(-0.055, length(lr)), col="green3", pch= 16,
       cex=0.6)
dev.off()

# Time difference
png(filename = "Rplot_extrmedges_8_9_t21.png", width = 420, height = 450,
    units = "px", pointsize = 12)
y_diff1 <- y_2[8, 9,] - y_1[8, 9,]
qntl <- quantile(y_diff1[cn], probs = c(0.05, 0.95))

plot(density(y_diff1[cn]), ylim=c(-.1, 1), xlim= c(-2.5, 2.5), 
     main = "L Accumbens vs L V Diencephalon", xlab ="Connectivity difference",
     cex.lab=1.4, cex.axis=1.3, cex.main=1.4)
abline(v=qntl[2], col= "red")
abline(v=qntl[1], col= "red")
points(y_diff1[shr], rep(0, length(shr)), col="blue", pch= 16, cex=0.6)
points(y_diff1[lr], rep(-0.055, length(lr)), col="green3", pch= 16,
       cex=0.6)
dev.off()

# for those edges with minimum value for the control participants is > 0.1, we compute
# the Intra-Class Correlations (ICC) as follows

y_cor_control <- data.frame( ICC= sapply(1:nrow(vc_symm), 
                  function(i) icc_fun(y_1[vc_symm$X1[i],
                  vc_symm$X2[i],] , y_2[vc_symm$X1[i], vc_symm$X2[i],], l=cn) ), 
                  Row= vc_symm$X1, Column= vc_symm$X2)


ggplot(y_cor_control, aes(y= Row, x= Column, color= ICC)) +
  geom_point(size=0.3) +
  scale_color_gradientn(colors = terrain.colors(10)) +
  labs(y = "Row", x = "Column", title = "ICC of the 1093 connectivity edges")   +
  theme(axis.text.x=element_text(size=12),
        axis.title.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.y=element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        title = element_text(size=12),
        strip.text = element_text(size=12))

ggsave(file="Rplot_ICC_control_grp_b01.png", height =4, width =5.5,
       units ="in", dpi= 300)

#Next, the following table provides the count of edges that contain the ICC at least $50 \%$ and weights $\ge 0.1$ out of the 875 edges chosen earlier. The plot provides the 576 edges satisfying the criteria mentioned earlier. Notice again that the plot will contain $576 \times 2= 1152$ edges due to symmetry.



##################################################################################
# Connectivity plot of a representative participant, means over various risk groups

set.seed(1323)
rp <- sample(seq(60, 394), 1)
rp<- 19

# for this representative sample rp is 323, a female, SHR with IT dosage 27, age 12.7
l <- c(c(y_1[,, rp]), c(y_2[,, rp]))
Connectivity_rp<- sapply(1:length(l), function(i) {
  if(l[i]==0) {return("0")}
  else if (l[i]> 0 & l[i]<=0.1) {return("(0, 0.1]")}
  else if(l[i]> 0.1 & l[i]<=0.5) {return("(0.1, 0.5]")}
  else if(l[i]> 0.5 & l[i]<=1) {return("(0.5, 1]")}
  else if(l[i]> 1 & l[i]<=5) {return("(1, 5]")}
  else if(l[i]> 5) {return("> 5")}
} )

y_cor_control_rep <- data.frame( Connectivity= Connectivity_rp,  
                               Row= rep(c(1:379),2), 
                               Column= rep(rep(c(1:379), each= 379), 2),
                               Time= rep(c("RE-1", "W-120"), 
                                         each= length(c(y_1[,,rp]))) )


# get the mean connectivity at each pair of ROIs for the Low Risk patients

y_1_meanP <- apply(y_1[,,-cn], c(1,2), FUN = mean)
y_2_meanP <- apply(y_2[,,-cn], c(1,2), FUN = mean)

l_mean <- c(c(y_1_meanP), c(y_2_meanP))

Connectivity_mean <- sapply(1:length(l_mean), function(i) {
  if(l_mean[i]==0) {return("0")}
  else if (l_mean[i]> 0 & l_mean[i]<=0.1) {return("(0, 0.1]")}
  else if(l_mean[i]> 0.1 & l_mean[i]<=0.5) {return("(0.1, 0.5]")}
  else if(l_mean[i]> 0.5 & l_mean[i]<=1) {return("(0.5, 1]")}
  else if(l_mean[i]> 1 & l_mean[i]<=5) {return("(1, 5]")}
  else if(l_mean[i]> 5) {return("> 5")}
})

y_cor_control_meanP <- data.frame(Connectivity = Connectivity_mean,  
                                 Row= rep(c(1:379),2), 
                                 Column= rep(rep(c(1:379), each= 379), 2),
                                 Time= rep(c("RE-1", "W-120"), 
                                           each= length(c(y_1_meanP))) )

# get the mean connectivity at each pair of ROIs for the control cohort

y_1_meanC <- apply(y_1[,,cn], c(1,2), FUN = mean)
y_2_meanC <- apply(y_2[,,cn], c(1,2), FUN = mean)

l_mean <- c(c(y_1_meanC), c(y_2_meanC))
Connectivity_mean <- sapply(1:length(l_mean), function(i) {
  if(l_mean[i]==0) {return("0")}
  else if (l_mean[i]> 0 & l_mean[i]<=0.1) {return("(0, 0.1]")}
  else if(l_mean[i]> 0.1 & l_mean[i]<=0.5) {return("(0.1, 0.5]")}
  else if(l_mean[i]> 0.5 & l_mean[i]<=1) {return("(0.5, 1]")}
  else if(l_mean[i]> 1 & l_mean[i]<=5) {return("(1, 5]")}
  else if(l_mean[i]> 5) {return("> 5")}
} )

y_cor_control_meanC <- data.frame(Connectivity = Connectivity_mean,  
                                 Row= rep(c(1:379),2), 
                                 Column= rep(rep(c(1:379), each= 379), 2),
                                 Time= rep(c("RE-1", "W-120"), 
                                           each= length(c(y_1_meanC))) )


dat_combine<- data.frame( rbind.data.frame(y_cor_control_rep, y_cor_control_meanC, 
                                y_cor_control_meanP),
                          pid= rep(c("Sample", "Mean: Control", "Mean: Patients"),
                                   each= nrow(y_cor_control_rep)) )

dat_combine$pid <- factor(dat_combine$pid, levels=c("Sample", "Mean: Control", 
                                             "Mean: Patients"))

ggplot(dat_combine, aes(x = Row, y = Column)) +
  geom_tile(aes(fill=Connectivity), size=0.4) +
  scale_fill_manual(breaks= c("0", "(0, 0.1]","(0.1, 0.5]", "(0.5, 1]", "(1, 5]",
                              "> 5"),
                    values = c("white","green3","yellow","orange","red","black")) +
  
  facet_grid(pid ~ Time, labeller = label_parsed) +
  labs(y = "Row", x = "Column", title= "Connectivity matrices: means and sample")   +
  theme(axis.text.x=element_text(size=12),
        axis.title.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.y=element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        title = element_text(size=12),
        strip.text = element_text(size=12))

ggsave(file="Rplot_Connectivity_means.png", height =10, width =9,
       units ="in", dpi= 300)


# we further select those edges with ICC at least 0.5 for the control cohort
y_cor_control1 <- subset(y_cor_control[,c(1:3)], ICC >= 0.5)

ggplot(y_cor_control1, aes(y = Row, x = Column)) +
  geom_point(size=0.4) +
  labs(y = "Row", x = "Column", 
       title = "Edges with connectivity constraints") +
  theme(axis.text.x=element_text(size=12),
        axis.title.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.y=element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        title = element_text(size=12),
        strip.text = element_text(size=12))

ggsave(file="Rplot_ICC_control_grp_a01.png", height =4, width =5,
       units ="in", dpi= 300)

# skim the duplicate edges counted due to symmetry of the connectivity matrix
y_cor_control2 <- subset(y_cor_control1, Row > Column) # remove the repeated points

y_cor_control4 <- y_cor_control2
y_cor_control4$Row <- as.factor(y_cor_control4$Row)
y_cor_control4$Column <- as.factor(y_cor_control4$Column)


#################################################################
# Massive multiple regression analysis for the 671 selected edges
# For sex, we get the coefficient for the male
# For risk we set the control group as the baseline. 
########

cn1 <- c(cn, lr)

# subset the covariate data for control and low risk groups
xdat_cr1 <- Xdat[cn1, ]
xdat_cr1$Age <- scale(xdat_cr1$Age)
colnames(xdat_cr1$Age) <- "Age"

# for Time = 1

# p-value collection
reg1_CLR_pvalues <- data.frame(sapply(1:nrow(y_cor_control2), 
                                      function(i) {
                                        a1 <- y_cor_control2$Row[i]
                                        b1 <- y_cor_control2$Column[i]
                                        Y <- y_1[a1, b1, cn1]
                                        fit<- lm(Y ~ xdat_cr1$Sex + xdat_cr1$Risk + xdat_cr1$Age)
                                        pv<- matrix(summary(fit)$coefficients[,4], nrow =1)
                                        colnames(pv) <- c("Intercept", "Sex", "Risk","Age")
                                        return(pv)
                                      } ) 
)

# coefficient estimate collection
reg1_CLR_coeffs <- data.frame(sapply(1:nrow(y_cor_control2), 
                                     function(i) {
                                       a1 <- y_cor_control2$Row[i]
                                       b1 <- y_cor_control2$Column[i]
                                       Y <- y_1[a1, b1, cn1]
                                       fit<- lm(Y ~ xdat_cr1$Sex + xdat_cr1$Risk + xdat_cr1$Age)
                                       coeffs<- matrix(summary(fit)$coefficients[,1], nrow =1)
                                       colnames(coeffs) <- c("Intercept", "Sex", "Risk","Age")
                                       return(coeffs)
                                     } ) 
)

reg1_CLR_pvalues_adjcoeff <- apply(reg1_CLR_pvalues, 1, p.adjust, method= "fdr")

reg1_CLR_pvalues_adj <- which(reg1_CLR_pvalues_adjcoeff[,2] < 0.05)
reg1_CLR_sex <- reg1_CLR_coeffs[2, reg1_CLR_pvalues_adj]

reg1_CLR_pvalues_adj <- which(reg1_CLR_pvalues_adjcoeff[,3] < 0.05)
reg1_CLR_risk <- reg1_CLR_coeffs[3, reg1_CLR_pvalues_adj]

reg1_CLR_pvalues_adj <- which(reg1_CLR_pvalues_adjcoeff[,4] < 0.05)
reg1_CLR_age <- reg1_CLR_coeffs[4, reg1_CLR_pvalues_adj]


# for Time = 2

# p-values collection
reg2_CLR_pvalues <- data.frame(sapply(1:nrow(y_cor_control2), 
                                      function(i) {
                                        a1 <- y_cor_control2$Row[i]
                                        b1 <- y_cor_control2$Column[i]
                                        Y <- y_2[a1, b1, cn1]
                                        fit<- lm(Y ~ xdat_cr1$Sex + xdat_cr1$Risk + xdat_cr1$Age)
                                        pv<- matrix(summary(fit)$coefficients[,4], nrow =1)
                                        colnames(pv) <- c("Intercept", "Sex", "Risk","Age")
                                        return(pv)
                                      } ) 
)

# coefficient estimates collection
reg2_CLR_coeffs <- data.frame(sapply(1:nrow(y_cor_control2), 
                                     function(i) {
                                       a1 <- y_cor_control2$Row[i]
                                       b1 <- y_cor_control2$Column[i]
                                       Y <- y_2[a1, b1, cn1]
                                       fit<- lm(Y ~ xdat_cr1$Sex + xdat_cr1$Risk + xdat_cr1$Age)
                                       coeffs<- matrix(summary(fit)$coefficients[,1], nrow =1)
                                       colnames(coeffs) <- c("Intercept", "Sex", "Risk","Age")
                                       return(coeffs)
                                     } ) 
)

reg2_CLR_pvalues_adjcoeff <- apply(reg2_CLR_pvalues, 1, p.adjust, method= "fdr")

reg2_CLR_pvalues_adjcoeff <- apply(reg2_CLR_pvalues, 1, p.adjust, method= "fdr")
reg2_CLR_pvalues_adj <- which(reg2_CLR_pvalues_adjcoeff[,2] < 0.05)
reg2_CLR_sex <- reg2_CLR_coeffs[2, reg2_CLR_pvalues_adj]

reg2_CLR_pvalues_adj <- which(reg2_CLR_pvalues_adjcoeff[,3] < 0.05)
reg2_CLR_risk <- reg2_CLR_coeffs[3, reg2_CLR_pvalues_adj]

reg2_CLR_pvalues_adj <- which(reg2_CLR_pvalues_adjcoeff[,4] < 0.05)
reg2_CLR_age <- reg2_CLR_coeffs[4, reg2_CLR_pvalues_adj]


# for Time = 2 subtract Time = 1
reg21_CLR_pvalues <- data.frame(sapply(1:nrow(y_cor_control2), 
                                       function(i) {
                                         a1 <- y_cor_control2$Row[i]
                                         b1 <- y_cor_control2$Column[i]
                                         Y <- y_2[a1, b1, cn1] - y_1[a1, b1, cn1]
                                         fit<- lm(Y ~ xdat_cr1$Sex + xdat_cr1$Risk + xdat_cr1$Age)
                                         coeffs<- matrix(summary(fit)$coefficients[,4], nrow =1)
                                         colnames(coeffs) <- c("Intercept", "Sex", "Risk","Age")
                                         return(coeffs)
                                       } ) 
)

reg21_CLR_coeffs <- data.frame(sapply(1:nrow(y_cor_control2), 
                                      function(i) {
                                        a1 <- y_cor_control2$Row[i]
                                        b1 <- y_cor_control2$Column[i]
                                        Y <- y_2[a1, b1, cn1] - y_1[a1, b1, cn1]
                                        fit<- lm(Y ~ xdat_cr1$Sex + xdat_cr1$Risk + xdat_cr1$Age)
                                        coeffs<- matrix(summary(fit)$coefficients[,1], nrow =1)
                                        colnames(coeffs) <- c("Intercept", "Sex", "Risk","Age")
                                        return(coeffs)
                                      } ) 
)

reg21_CLR_pvalues_adjcoeff <- apply(reg21_CLR_pvalues, 1, p.adjust, method= "fdr")

reg21_CLR_pvalues_adj <- which(reg21_CLR_pvalues_adjcoeff[,2] < 0.05)
reg21_CLR_sex <- reg21_CLR_coeffs[2, reg21_CLR_pvalues_adj]

reg21_CLR_pvalues_adj <- which(reg21_CLR_pvalues[,3] < 0.05)
reg21_CLR_risk <- reg21_CLR_coeffs[3, reg21_CLR_pvalues_adj]

reg21_CLR_pvalues_adj <- which(reg21_CLR_pvalues[,4] < 0.05)
reg21_CLR_age <- reg21_CLR_coeffs[4, reg21_CLR_pvalues_adj]


######## For Control vs Standard-High Risk
cn2 <- c(cn, shr)

# subset the covariate data for control and low risk groups
xdat_cr2 <- Xdat[cn2, ]
#xdat_cr2$Risk <- relevel(xdat_cr2$Risk, ref = "")
xdat_cr2$Age <- scale(xdat_cr2$Age)
colnames(xdat_cr2$Age) <- "Age"

# for Time = 1

# p-values collection
reg1_CSHR_pvalues <- data.frame(sapply(1:nrow(y_cor_control2), 
                                       function(i) {
                                         a1 <- y_cor_control2$Row[i]
                                         b1 <- y_cor_control2$Column[i]
                                         Y <- y_1[a1, b1, cn2]
                                         fit<- lm(Y ~ xdat_cr2$Sex + xdat_cr2$Risk + xdat_cr2$Age)
                                         pvalues<- matrix(summary(fit)$coefficients[,4], nrow =1)
                                         colnames(pvalues) <- c("Intercept", "Sex", "Risk","Age")
                                         return(pvalues)
                                       } ) 
)

#reg1_CSHR_pvalues_adj <- apply(reg1_CLR_pvalues, 1, p.adjust, method= "fdr")
# coefficient estimate collection
reg1_CSHR_coeffs <- data.frame(sapply(1:nrow(y_cor_control2), 
                                      function(i) {
                                        a1 <- y_cor_control2$Row[i]
                                        b1 <- y_cor_control2$Column[i]
                                        Y <- y_1[a1, b1, cn2]
                                        fit <- lm(Y ~ xdat_cr2$Sex + xdat_cr2$Risk + xdat_cr2$Age)
                                        coeffs<- matrix(summary(fit)$coefficients[,1], nrow =1)
                                        colnames(coeffs) <- c("Intercept", "Sex", "Risk","Age")
                                        return(coeffs)
                                      }) 
)

reg1_CSHR_pvalues_adjcoeff <- apply(reg1_CSHR_pvalues, 1, p.adjust, method= "fdr")

reg1_CSHR_pvalues_adj <- which(reg1_CSHR_pvalues_adjcoeff[,2] < 0.05)
reg1_CSHR_sex <- reg1_CSHR_coeffs[2, reg1_CSHR_pvalues_adj]

reg1_CSHR_pvalues_adj <- which(reg1_CSHR_pvalues_adjcoeff[,3] < 0.05)
reg1_CSHR_risk <- reg1_CSHR_coeffs[3, reg1_CSHR_pvalues_adj]

reg1_CSHR_pvalues_adj <- which(reg1_CSHR_pvalues_adjcoeff[,4] < 0.05)
reg1_CSHR_age <- reg1_CSHR_coeffs[4, reg1_CSHR_pvalues_adj]


# for Time = 2
# p-values collection
reg2_CSHR_pvalues <- data.frame(sapply(1:nrow(y_cor_control2), 
                                       function(i) {
                                         a1 <- y_cor_control2$Row[i]
                                         b1 <- y_cor_control2$Column[i]
                                         Y <- y_2[a1, b1, cn2]
                                         fit<- lm(Y ~ xdat_cr2$Sex + xdat_cr2$Risk + xdat_cr2$Age)
                                         pvalues<- matrix(summary(fit)$coefficients[,4], nrow =1)
                                         colnames(pvalues) <- c("Intercept", "Sex", "Risk","Age")
                                         return(pvalues)
                                       } ) 
)

# coefficient estimate collection
reg2_CSHR_coeffs <- data.frame(sapply(1:nrow(y_cor_control2), 
                                      function(i) {
                                        a1 <- y_cor_control2$Row[i]
                                        b1 <- y_cor_control2$Column[i]
                                        Y <- y_2[a1, b1, cn2]
                                        fit<- lm(Y ~ xdat_cr2$Sex + xdat_cr2$Risk + xdat_cr2$Age)
                                        coeffs<- matrix(summary(fit)$coefficients[,1], nrow =1)
                                        colnames(coeffs) <- c("Intercept", "Sex", "Risk","Age")
                                        return(coeffs)
                                      } ) 
)

reg2_CSHR_pvalues_adjcoeff <- apply(reg2_CSHR_pvalues, 1, p.adjust, method= "fdr")

reg2_CSHR_pvalues_adj <- which(reg2_CSHR_pvalues_adjcoeff[,2] < 0.05)
reg2_CSHR_sex<- reg2_CSHR_coeffs[2,reg2_CSHR_pvalues_adj]

reg2_CSHR_pvalues_adj <- which(reg2_CSHR_pvalues_adjcoeff[,3] < 0.05)
reg2_CSHR_risk <- reg2_CSHR_coeffs[3, reg2_CSHR_pvalues_adj]

reg2_CSHR_pvalues_adj <- which(reg2_CSHR_pvalues_adjcoeff[,4] < 0.05)
reg2_CSHR_age<- reg2_CSHR_coeffs[4,reg2_CSHR_pvalues_adj]

# for connectivity differences: Time = 2 - Time = 1
# p-values collection
reg21_CSHR_pvalues <- data.frame(sapply(1:nrow(y_cor_control2), 
                                        function(i) {
                                          a1 <- y_cor_control2$Row[i]
                                          b1 <- y_cor_control2$Column[i]
                                          Y <- y_2[a1, b1, cn2] - y_1[a1, b1, cn2]
                                          fit<- lm(Y ~ xdat_cr2$Sex + xdat_cr2$Risk + xdat_cr2$Age)
                                          pv<- matrix(summary(fit)$coefficients[,4], nrow =1)
                                          colnames(pv) <- c("Intercept", "Sex", "Risk","Age")
                                          return(pv)
                                        } ) 
)

# coeffcient estimate collection
reg21_CSHR_coeffs <- data.frame(sapply(1:nrow(y_cor_control2), 
                                       function(i) {
                                         a1 <- y_cor_control2$Row[i]
                                         b1 <- y_cor_control2$Column[i]
                                         Y <- y_2[a1, b1, cn2] - y_1[a1, b1, cn2]
                                         fit <- lm(Y ~ xdat_cr2$Sex + xdat_cr2$Risk + xdat_cr2$Age)
                                         coeffs<- matrix(summary(fit)$coefficients[,1], nrow =1)
                                         colnames(coeffs) <- c("Intercept", "Sex", "Risk","Age")
                                         return(coeffs)
                                       } ) 
)


#reg21_CSHR_pvalues_adj <- which(p.adjust(reg21_CSHR_pvalues[2,], method= "fdr") < 0.05)
#reg21_CSHR_sex<- reg21_CSHR_coeffs[2,reg21_CSHR_pvalues_adj]

reg21_CSHR_pvalues_adjcoeff <- apply(reg21_CSHR_pvalues, 1, p.adjust, method= "fdr")

reg21_CSHR_pvalues_adj <- which(reg21_CSHR_pvalues_adjcoeff[,3] < 0.05)
reg21_CSHR_risk<- reg21_CSHR_coeffs[3, reg21_CSHR_pvalues_adj]

reg21_CSHR_pvalues_adj<- which(reg21_CSHR_pvalues_adjcoeff[,4] < 0.05)
reg21_CSHR_age<- reg21_CSHR_coeffs[4,reg21_CSHR_pvalues_adj]

# histograms of age coefficients for comparison of Control vs Standard-High Risk
age_CSHR <- data.frame(Coefficient= c(as.numeric(reg1_CSHR_age), 
                                      as.numeric(reg2_CSHR_age), 
                                      as.numeric(reg21_CSHR_age)),
                       Design = c(rep("Time 1, C vs. SHR", length(reg1_CSHR_age)),
                                  rep("Time 2, C vs. SHR", length(reg2_CSHR_age)),
                                  rep("Difference, C vs. SHR", 
                                      length(reg21_CSHR_age))))

age_CLR <- data.frame(Coefficient= c(as.numeric(reg1_CLR_age), 
                                     as.numeric(reg2_CLR_age), 
                                     as.numeric(reg21_CLR_age)),
                      Design = c(rep("Time 1, C vs. LR", length(reg1_CLR_age)),
                                 rep("Time 2, C vs. LR", length(reg2_CLR_age)),
                                 rep("Difference, C vs. LR", 
                                     length(reg21_CLR_age))))

age_plotdat<- rbind.data.frame(age_CLR, age_CSHR)

ggplot(age_plotdat, aes(Coefficient)) + geom_histogram() + 
  facet_wrap(~ factor(Design, c("Time 1, C vs. LR", "Time 2, C vs. LR", "Difference, C vs. LR",
                                "Time 1, C vs. SHR", "Time 2, C vs. SHR", "Difference, C vs. SHR")) )+
  xlab("Coefficient of Age")+
  ylab("Count of significant coefficients after FDR correction")

ggsave(file="Rplot_coeff_massiveregression_age.png", height = 6, width =8.5,
       units ="in", dpi= 300)

# histograms of risk coefficients for comparison of Control vs Standard-High Risk
risk_CSHR <- data.frame(Coefficient= c(as.numeric(reg1_CSHR_risk), 
                                       as.numeric(reg2_CSHR_risk), 
                                       as.numeric(reg21_CSHR_risk)),
                        Design = c(rep("Time 1, C vs. SHR", length(reg1_CSHR_risk)),
                                   rep("Time 2, C vs. SHR", length(reg2_CSHR_risk)),
                                   rep("Difference, C vs. SHR", 
                                       length(reg21_CSHR_risk))))

risk_CLR <- data.frame(Coefficient= c(as.numeric(reg1_CLR_risk), 
                                      as.numeric(reg2_CLR_risk), 
                                      as.numeric(reg21_CLR_risk)),
                       Design = c(rep("Time 1, C vs. LR", length(reg1_CLR_risk)),
                                  rep("Time 2, C vs. LR", length(reg2_CLR_risk)),
                                  rep("Difference, C vs. LR", 
                                      length(reg21_CLR_risk))))

risk_plotdat <- rbind.data.frame(risk_CLR, risk_CSHR)

ggplot(risk_plotdat, aes(Coefficient)) + geom_histogram() + 
  facet_wrap(~ factor(Design, c("Time 1, C vs. LR", "Time 2, C vs. LR", 
                                "Difference, C vs. LR",
                                "Time 1, C vs. SHR", "Time 2, C vs. SHR", 
                                "Difference, C vs. SHR")) )+
  xlab("Coefficient of Risk")+
  ylab("Count of significant coefficients after FDR correction") 

ggsave(file="Rplot_coeff_massiveregression_risk.png", height = 6, width =8.5,
       units ="in", dpi= 300)

# histograms of sex coefficients for comparison of Control vs Standard-High Risk
sex_CSHR <- data.frame(Coefficient= c(as.numeric(reg1_CSHR_sex), 
                                      as.numeric(reg2_CSHR_sex)),
                       Design = c(rep("Time 1, C vs. SHR", length(reg1_CSHR_sex)),
                                  rep("Time 2, C vs. SHR", length(reg2_CSHR_sex))))

sex_CLR <- data.frame(Coefficient= c(as.numeric(reg1_CLR_sex), 
                                     as.numeric(reg2_CLR_sex)),
                      Design = c(rep("Time 1, C vs. LR", length(reg1_CLR_sex)),
                                 rep("Time 2, C vs. LR", length(reg2_CLR_sex))))

sex_plotdat <- rbind.data.frame(sex_CLR, sex_CSHR)

ggplot(sex_plotdat, aes(Coefficient)) + geom_histogram() + facet_wrap(~Design) + 
  xlab("Coefficient of Sex") +
  ylab("Count of significant coefficients after FDR correction")

ggsave(file="Rplot_coeff_massiveregression_sex.png", height = 5, width = 6.5,
       units ="in", dpi= 300)

# finally a large dataset containing all the 671 edges with the coefficients, the p-values and the 
# adjusted p-values with 'fdr' correction:

dat_final<- data.frame(Sex_CLR_t1_coeff= as.numeric(reg1_CLR_coeffs[2,]), 
                       Sex_CLR_t1_pvalue= as.numeric(reg1_CLR_pvalues[2,]),
                       Sex_CLR_t1_pvadj= as.numeric(reg1_CLR_pvalues_adjcoeff[,2]),
                       Sex_CLR_t2_coeff= as.numeric(reg2_CLR_coeffs[2,]), 
                       Sex_CLR_t2_pvalue= as.numeric(reg2_CLR_pvalues[2,]),
                       Sex_CLR_t2_pvadj= as.numeric(reg2_CLR_pvalues_adjcoeff[,2]),
                       Sex_CLR_t21_coeff= as.numeric(reg21_CLR_coeffs[2,]), 
                       Sex_CLR_t21_pvalue= as.numeric(reg21_CLR_pvalues[2,]),
                       Sex_CLR_t21_pvadj= as.numeric(reg21_CLR_pvalues_adjcoeff[,2]),
                       
                       Sex_CSHR_t1_coeff= as.numeric(reg1_CSHR_coeffs[2,]), 
                       Sex_CSHR_t1_pvalue= as.numeric(reg1_CSHR_pvalues[2,]),
                       Sex_CSHR_t1_pvadj= as.numeric(reg1_CSHR_pvalues_adjcoeff[,2]),
                       Sex_CSHR_t2_coeff= as.numeric(reg2_CSHR_coeffs[2,]), 
                       Sex_CSHR_t2_pvalue= as.numeric(reg2_CSHR_pvalues[2,]),
                       Sex_CSHR_t2_pvadj= as.numeric(reg2_CSHR_pvalues_adjcoeff[,2]),
                       Sex_CSHR_t21_coeff= as.numeric(reg21_CSHR_coeffs[2,]), 
                       Sex_CSHR_t21_pvalue= as.numeric(reg21_CSHR_pvalues[2,]),
                       Sex_CSHR_t21_pvadj= as.numeric(reg21_CSHR_pvalues_adjcoeff[,2]),
                       
                       Risk_CLR_t1_coeff= as.numeric(reg1_CLR_coeffs[3,]), 
                       Risk_CLR_t1_pvalue= as.numeric(reg1_CLR_pvalues[3,]),
                       Risk_CLR_t1_pvadj= as.numeric(reg1_CLR_pvalues_adjcoeff[,3]),
                       Risk_CLR_t2_coeff= as.numeric(reg2_CLR_coeffs[3,]), 
                       Risk_CLR_t2_pvalue= as.numeric(reg2_CLR_pvalues[3,]),
                       Risk_CLR_t2_pvadj= as.numeric(reg2_CLR_pvalues_adjcoeff[,3]),
                       Risk_CLR_t21_coeff= as.numeric(reg21_CLR_coeffs[3,]), 
                       Risk_CLR_t21_pvalue= as.numeric(reg21_CLR_pvalues[3,]),
                       Risk_CLR_t21_pvadj= as.numeric(reg21_CLR_pvalues_adjcoeff[,3]),
                       
                       Risk_CSHR_t1_coeff= as.numeric(reg1_CSHR_coeffs[3,]), 
                       Risk_CSHR_t1_pvalue= as.numeric(reg1_CSHR_pvalues[3,]),
                       Risk_CSHR_t1_pvadj= as.numeric(reg1_CSHR_pvalues_adjcoeff[,3]),
                       Risk_CSHR_t2_coeff= as.numeric(reg2_CSHR_coeffs[3,]), 
                       Risk_CSHR_t2_pvalue= as.numeric(reg2_CSHR_pvalues[3,]),
                       Risk_CSHR_t2_pvadj= as.numeric(reg2_CSHR_pvalues_adjcoeff[,3]),
                       Risk_CSHR_t21_coeff= as.numeric(reg21_CSHR_coeffs[3,]), 
                       Risk_CSHR_t21_pvalue= as.numeric(reg21_CSHR_pvalues[3,]),
                       Risk_CSHR_t21_pvadj= as.numeric(reg21_CSHR_pvalues_adjcoeff[,3]),
                       
                       Age_CLR_t1_coeff= as.numeric(reg1_CLR_coeffs[3,]), 
                       Age_CLR_t1_pvalue= as.numeric(reg1_CLR_pvalues[3,]),
                       Age_CLR_t1_pvadj= as.numeric(reg1_CLR_pvalues_adjcoeff[,3]),
                       Age_CLR_t2_coeff= as.numeric(reg2_CLR_coeffs[2,]), 
                       Age_CLR_t2_pvalue= as.numeric(reg2_CLR_pvalues[2,]),
                       Age_CLR_t2_pvadj= as.numeric(reg2_CLR_pvalues_adjcoeff[,2]),
                       Age_CLR_t21_coeff= as.numeric(reg21_CLR_coeffs[2,]), 
                       Age_CLR_t21_pvalue= as.numeric(reg21_CLR_pvalues[2,]),
                       Age_CLR_t21_pvadj= as.numeric(reg21_CLR_pvalues_adjcoeff[,2]),
                       
                       Age_CSHR_t1_coeff= as.numeric(reg1_CSHR_coeffs[2,]), 
                       Age_CSHR_t1_pvalue= as.numeric(reg1_CSHR_pvalues[2,]),
                       Age_CSHR_t1_pvadj= as.numeric(reg1_CSHR_pvalues_adjcoeff[,2]),
                       Age_CSHR_t2_coeff= as.numeric(reg2_CSHR_coeffs[2,]), 
                       Age_CSHR_t2_pvalue= as.numeric(reg2_CSHR_pvalues[2,]),
                       Age_CSHR_t2_pvadj= as.numeric(reg2_CSHR_pvalues_adjcoeff[,2]),
                       Age_CSHR_t21_coeff= as.numeric(reg21_CSHR_coeffs[2,]), 
                       Age_CSHR_t21_pvalue= as.numeric(reg21_CSHR_pvalues[2,]),
                       Age_CSHR_t21_pvadj= as.numeric(reg21_CSHR_pvalues_adjcoeff[,2]))

write.csv(dat_final, "Massive_regression_output.csv")

################################################################################
# We gather the extreme brain edges for every patient under study and extreme 
# patients for every brain edge for every brain connectivity edge under study

# rows represent the persons in study, columns represent the edges in the study. 
#m_t1 <- matrix(NA, nrow= l, ncol= nrow(y_cor_control2))  # at the start
#m_t2 <- matrix(NA, nrow= l, ncol= nrow(y_cor_control2))  # at the 120 weeks study
#m_t1_tail <- matrix(NA, nrow= l, ncol= nrow(y_cor_control2))  # at the start
#m_t2_tail <- matrix(NA, nrow= l, ncol= nrow(y_cor_control2))  # at the 120 weeks study

#for (i1 in 1:nrow(m_t1)) {
#  for (i2 in 1:ncol(m_t1)) {
#    temp <- perc_edge(node_x = y_cor_control2$Row[i2], 
#                              node_y = y_cor_control2$Column[i2], cm = y_1, 
#                              cntrl_array = cn, pos = i1)
#    m_t1[i1, i2] <- temp$pvalue
#    m_t1_tail[i1, i2] <- temp$tail
#  }
#}

#for (i1 in 1:nrow(m_t2)) {
#  for (i2 in 1:ncol(m_t2)) {
#    temp <- perc_edge(node_x= y_cor_control2$Row[i2], 
#                              node_y= y_cor_control2$Column[i2],
#                              cm = y_2, cntrl_array = cn, pos=i1)
#    
#    m_t2[i1, i2] <- temp$pvalue
    
#    m_t2_tail[i1, i2] <- temp$tail
    
#  }
#}

# add risk designation for each patient, at time 1 and 2
#m_t1<- cbind(Risk_id, m_t1) 
#m_t2<- cbind(Risk_id, m_t2) 

# get the tail where the patient/participant belongs, at time 2: 120 weeks
#m_t1_risk <- data.frame(subset(m_t1, Risk_id != "C"))

# get the tail where the patient/participant belongs, at time 2: 120 weeks
#m_t2_risk <- data.frame(subset(m_t2, Risk_id != "C"))

# get the tail where the patient/participant belongs, at time 1
#m_t1_tail <- cbind(Risk_id, m_t1_tail)  

# get the tail where the patient/participant belongs, at time 2: 120 weeks
#m_t2_tail <- cbind(Risk_id, m_t2_tail) 

# get the tail where the patient belongs, at time 1
#m_t1_risktail <- data.frame(subset(m_t1_tail, Risk_id != "C"))

# get the tail where the patient belongs, at time 2: 120 weeks
#m_t2_risktail <- data.frame(subset(m_t2_tail, Risk_id != "C"))

# Each column summarizes the number of people, hence corresponding to each column summary, 
# we have the # of people that are significant in each column.

#############
# Now consider the extreme patients and edges with respect to connectivity differences 2 - 1
####

# rows represent the persons in study, columns represent the edges in the study. 

# get the p-value for every connectivity differences
m21 <- matrix(NA, nrow= lp, ncol= nrow(y_cor_control2))   
# get the tail position w.r.t. the control distribution for every connectivity
m21_tail <- matrix(NA, nrow= lp, ncol= nrow(y_cor_control2))

# difference of connectivity matrix at difference
for (i1 in 1:nrow(m21)) {
  for (i2 in 1:ncol(m21)) {
    temp <- perc_edge(node_x = y_cor_control2$Row[i2], 
                              node_y = y_cor_control2$Column[i2], cm = y_diff, 
                              cntrl_array = cn, pos = i1)
    
    m21[i1, i2] <- temp$pvalue
    
    m21_tail[i1, i2] <- temp$tail
  }
}

m21 <- cbind.data.frame(MRN= Xdat$MRN,Risk_id, m21)

m21_risk <- data.frame(subset(m21, Risk_id != "C"))

m21_tail <- cbind.data.frame(MRN= Xdat$MRN, Risk_id, m21_tail)

m21_rtail <- data.frame(subset(m21_tail, Risk_id != "C"))

# x1 is the number of significant patients in the time 1.
x21 <- sapply(3:ncol(m21_risk), function(i) {
  
  # find out for that edge which patients are extreme
  a1 <- which(m21_risk[,i] < 0.05)
  
  # find out the LR among those extreme patients
  a2 <- length(which(m21_risk[a1, 2] == "LR"))
  
  # find out the SHR among those extreme patients
  a3 <- length(which(m21_risk[a1, 2] == "SHR")) 
  
  # find out among the LR extreme patients, those are in the upper tail of the control dist
  a2_utail <- length(which(m21_risk[a1, 2] == "LR" & m21_rtail[a1, i] == "upper")) 
  
  # find out among the LR extreme patients, those are in the lower tail of the control dist
  a2_ltail <- length(which(m21_risk[a1, 2] == "LR" & m21_rtail[a1, i] == "lower")) 
  
  # find out among the LR extreme patients, those are in the upper tail of the control dist
  a3_utail <- length(which(m21_risk[a1, 2] == "SHR" & m21_rtail[a1, i] == "upper")) 
  
  # find out among the LR extreme patients, those are in the lower tail of the control dist
  a3_ltail <- length(which(m21_risk[a1, 2] == "SHR" & m21_rtail[a1, i] == "lower")) 
  
  a <- matrix(c(length(a1), a2, a2_utail, a2_ltail, a3, a3_utail, a3_ltail), nrow = 1)
  
  return(a)
})

# For connectivity difference we perform this operation:

# Now we shortlist those edges that have extreme patients count 70 or more 
roi_diff21 <- cbind.data.frame(Edges= seq(1: nrow(y_cor_control2)), y_cor_control2[, -1]) 

# collect the number of extreme patients in the edge
#which(x21[1,] >=70)
roi_diff21$ExtremePatients <- x21[1, ] 
roi_diff21$LR <- x21[2, ] 
roi_diff21$LR_UT <- x21[3, ]
roi_diff21$LR_LT <- x21[4, ]
roi_diff21$SHR <- x21[5, ] 
roi_diff21$SHR_UT <- x21[6, ]
roi_diff21$SHR_LT <- x21[7, ]

roi_diff21$Lobe_ROI2 <- roi_diff21$Lobe_ROI1 <- roi_diff21$ROI2 <- roi_diff21$ROI1 <- NA 

roi_diff21$Region_ROI2 <- roi_diff21$Region_ROI1 <- NA

# collect the ROI pair in the edge
for (i in 1:nrow(roi_diff21)) {
  roi_diff21$ROI1[i] <- paste(node_labels[roi_diff21[i,"Row"], "Side"], 
                            node_labels[roi_diff21[i,"Row"], "Name"], sep =" ")
  roi_diff21$ROI2[i] <- paste(node_labels[roi_diff21[i,"Column"], "Side"], 
                            node_labels[roi_diff21[i,"Column"], "Name"], sep =" ")
  
  roi_diff21$Lobe_ROI1[i] <- node_labels[roi_diff21[i,"Row"], "Lobe"]
  roi_diff21$Lobe_ROI2[i] <- node_labels[roi_diff21[i,"Column"], "Lobe"]
  
  roi_diff21$Region_ROI1[i] <- node_labels[roi_diff21[i,"Row"], "Region.Name"]
  roi_diff21$Region_ROI2[i] <- node_labels[roi_diff21[i,"Column"], "Region.Name"]
  
  #print(i)
  
}

# next we order the edges from highest to the lowest count of extreme patients
roi_diff21 <- roi_diff21[order(-roi_diff21$ExtremePatients),]
rownames(roi_diff21) <- NULL

write.csv(roi_diff21, "Edges_Connectivity_difference.csv")

################################################################################
# Regression where Neurocognitive scores are the response variables 
# Brain connectivity are the covariates. In another regression, we consider the 
# Risk categories classified by the ITT scores as a categorical covariate also
################################################################################

# subset the data to include only those edges with most extreme patients at least 66
#roi_diff <- subset(roi_diff21, ExtremePatients >= 66)

# read the neurocognitive data first~
ndat <- read.csv("Neurocog_scores.csv", header = TRUE)

########################################################
# first we try to check if the IQ scores were the same

temp20 <- ndat[, c(1,7)]  # choose the neurocog score and the MRN variable

# choose the neurocog scores that are non-empty
temp21 <- temp20[which(temp20[,2] != "."), ]

mrn21 <- temp21$MRN  # get MRNs for patients whose neurocog scores are non-empty

temp22 <- X3[, c(1,3)] # get the IQ scores at re-induction 1 phase

temp23 <- merge(temp21, temp22, by= "MRN")
temp24 <- merge(temp23, Xdat, by= "MRN")
temp25 <- na.omit(temp24)

t.test(temp25$re1.sb5.estiq.ss, as.numeric(temp25$sb5.estIQ.ss), paired = T, 
       alternative = "greater")

plotdat1<- cbind.data.frame(IQ = c(temp25$re1.sb5.estiq.ss, as.numeric(temp25$sb5.estIQ.ss),
                                    as.numeric(temp25$sb5.estIQ.ss) - temp25$re1.sb5.estiq.ss),
                             Time= rep(c("RE-1","W-120","Difference"), 
                                       each= length(temp25$re1.sb5.estiq.ss)),
                             Groups= rep(temp25$Group.ITT, 3) )

plotdat12<- subset(plotdat1, Time== "Difference")
plotdat13<- subset(plotdat1, Time!="Difference")

p1<- ggplot(plotdat12, aes(y=IQ, fill= Groups))+
  geom_boxplot(alpha=0.8)+
  facet_wrap(~Time)+ 
  scale_fill_brewer(palette = "Set2") +
  guides(x = guide_axis(position = "NULL"), fill="none")+
  geom_hline(yintercept = 0, col="red")+
  labs(y= "IQ differences")+theme_bw()

p2<- ggplot(plotdat13, aes(y=IQ, fill= Groups))+
  geom_boxplot(alpha=0.8)+
  facet_wrap(~Time) + 
  #scale_fill_brewer(labels = unname(TeX(c("$LR,ITT < 20$", "$LR,ITT \\geq 20$",
  #                                         "$SHR,ITT<27$","$LR,ITT \\geq 27$"))),
  #                                     palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  labs(fill= "Groups")+
guides(x = guide_axis(position = "NULL"))+theme_bw()

p1 + p2 + plot_layout(widths = c(3, 6)) 

ggsave(filename ="Rplot_IQdifferences_byRGITT.png", height =3, width = 10, 
       units = "in", dpi =300)
dev.off()


t.test(`IQ difference` ~ as.character(Risk), data = plotdat12, var.equal=T)

t.test(`IQ difference` ~ as.character(Sex), data = plotdat2, var.equal=T)

t1<- subset(plotdat2, `Risk & ITT` == "LR<20")[,1]
t2<- subset(plotdat2, `Risk & ITT` == "LR>=20")[,1]
t3<- subset(plotdat2, `Risk & ITT` == "SHR<27")[,1]
t4<- subset(plotdat2, `Risk & ITT` == "SHR>=27")[,1]

t.test(t1, t2, var.equal = T)
t.test(t1, t4, var.equal = T)
t.test(t3, t4, var.equal = T, alternative ="less")
t.test(t2, t4, var.equal = T)
t.test(t2, t3, var.equal = T, alternative= "greater")
t.test(t1, t3, var.equal = T, alternative ="greater")

roi_diff_reg <- NA

for (j in 1:nrow(roi_diff21)) {

y1d <- y_diff[roi_diff21$Row[j], roi_diff21$Column[j], ]
yd <- subset(y1d, names(y1d) %in% temp25$MRN)

dat10 <- cbind.data.frame(temp25, "Connectivity difference"= yd, 
                "IQ difference"= as.numeric(temp25$sb5.estIQ.ss)-
                  as.numeric(temp25$re1.sb5.estiq.ss) )
  
ft <- (lm(`IQ difference` ~ yd*Group.ITT, data = dat10))
ft1 <- (lm(`IQ difference` ~ yd+ Group.ITT, data = dat10))

aoft <- anova(ft, ft1)
ft_summ <- summary(lm(`IQ difference` ~ yd*Group.ITT, data = dat10))
ft1_summ <- summary(lm(`IQ difference` ~ yd+ Group.ITT, data = dat10))

Fpv <- pf(ft_summ$fstatistic[1], ft_summ$fstatistic[2], ft_summ$fstatistic[3], lower.tail = F)

# 
if(Fpv< 0.05 & aoft$`Pr(>F)`[2] <0.05) {
print(j)
roi_diff_reg <- rbind.data.frame(roi_diff_reg, data.frame("Edge"= roi_diff21$Edges[j],
                                                          "Row"= roi_diff21$Row[j],
                                                          "Column"= roi_diff21$Column[j],
                                                          "Intercept"= ft_summ$coefficients[1,1],
                                                          "Connectivity"= ft_summ$coefficients[2,1],
                                                          "LR,ITTge20"= ft_summ$coefficients[3,1],
                                                          "SHR,ITTle26"= ft_summ$coefficients[4,1],
                                                          "SHR,ITTge27"= ft_summ$coefficients[5,1],
                                                          "Yd:LR,ITTge20"= ft_summ$coefficients[6,1],
                                                          "Yd:SHR,ITTle26"= ft_summ$coefficients[7,1],
                                                          "Yd:SHR,ITTge27"=ft_summ$coefficients[8,1],
                                                          "Adj.R2"= ft_summ$adj.r.squared,
                                                          "Fpvalue"=Fpv))

}
}

roi_diff_reg <- roi_diff_reg[-1,]

write.csv(roi_diff_reg, "ROI_regression_IQ.csv")

# we order the connectivity edges by the decreasing order the adjusted R-squared of the regression
roi_diff_reg2 <- roi_diff_reg[order(-roi_diff_reg$Adj.R2), ]

yd_vec <- NA

for (i in 1:6) {
  y1d <- y_diff[roi_diff_reg2$Row[i], roi_diff_reg2$Column[i], ]
  yd <- subset(y1d, names(y1d) %in% temp25$MRN)
  
  yd_vec <- c(yd_vec, yd)
  
}

yd_vec <- yd_vec[-1]

dat10 <- cbind.data.frame("Connectivity difference"= yd_vec, 
          "IQ difference"= rep((as.numeric(temp25$sb5.estIQ.ss)- as.numeric(temp25$re1.sb5.estiq.ss)), 6),
          "Edge"= rep(c("R PFcm vs R OP2-3/VS", "R Dorsal 24d vs R 5m Ventral",
                        "L PHT vs L TE1 posterior", "L TPOccipital Junction 1 vs 2",
                        "L 7m vs L 7p", "L 31p vs L 7m"), each= nrow(temp25)),
          "Groups" = rep(temp25$Group.ITT, each= 6 ) ) 

#ft <- (lm(`IQ difference` ~ `Connectivity difference`*Group.ITT, data = dat10))

#ft_summ <- summary(ft)

#interactions::interact_plot(model =ft,  pred = `Connectivity difference`, 
#                            modx = Group.ITT, interval =TRUE)

p1<- ggplot(dat10, aes(x=`Connectivity difference`, y=`IQ difference`, color = Groups)) +
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~ Edge, scale= "free_x")+
  scale_color_discrete(labels = unname(TeX(c("$LR,ITT < 20$", "$LR,ITT \\geq 20$",
                                          "$SHR,ITT<27$","$LR,ITT \\geq 27$"))))+
  theme_bw()

  
  #labs( subtitle = TeX( "R Area PFcm vs. R Area OP2-3/VS, $\\, \\bar{R}^2 \\approx 9\\%$" )) +
  #guides(color = "none")
p1

ggsave(filename ="Rplot_IQvsConn_differences_byRGITT.png", height =7, width = 8, 
       units = "in", dpi =300)
dev.off()


##########################################################################################
# We try to check if the digit span/ auditory WM scores were the same between time points

# the following are the connectivity edges involved in the working memory of brain
wmem_edges <- c(403,150,159,210,401,400,75,203,254,600,61,286,551,525,524,287,526,
                253,118,454,553,554)

i <- 8 # obtain the wechsler scores: higher is better
#i <- 9 # obtain the WJ3 auditory memory score: higher is better

temp30 <- ndat[, c(1,i)]  # choose the neurocog score and the MRN variable

# choose the neurocog scores that are non-empty
temp31 <- temp30[which(temp30[,2] != "."), ]

#mrn31 <- temp31$MRN  # get MRNs for patients whose neurocog scores are non-empty

temp32 <- X3[, c(1,4,5)]

temp33 <- merge(temp31, temp32, by= "MRN")[,c(1,2)]
temp34 <- merge(temp33, Xdat, by= "MRN")
#temp35 <- na.omit(temp34)

t.test(as.numeric(wech.totalds.scs) ~ as.character(Risk), data=temp34)
#t.test(as.numeric(wj3.audwm.ss) ~ as.character(Risk), data=temp34)

p1<- ggplot(temp34, aes(y=as.numeric(wech.totalds.scs), fill = Risk)) +
  geom_boxplot(alpha=0.8)+
  scale_fill_hue() +
  guides(x = guide_axis(position = "NULL"))+
  labs(y= "Wechsler's score",
       title = "By Risk Groups")

p1

p2<- ggplot(temp34, aes(y=as.numeric(wech.totalds.scs), fill = Group.ITT)) +
  geom_boxplot(alpha=0.8) +
  scale_fill_brewer(labels = unname(TeX(c("$LR,ITT < 20$", "$LR,ITT \\geq 20$",
                                          "$SHR,ITT<27$","$LR,ITT \\geq 27$"))),
                    palette = "Set2")+
  guides(x = guide_axis(position = "NULL"))+
  labs(y= "Wechsler's score",
       title = "By Risk Groups and ITT",
       fill= "Groups & ITT")

p2
  
p1 + p2 + plot_layout(widths = c(3, 4))

ggsave(filename ="Rplot_wechWM_byRGITT.png", height =4, width = 9, 
       units = "in", dpi =300)
dev.off()


t1 <- as.numeric(subset(temp34, Group.ITT == "SHR, ITT<=26")[,2])
t2 <- as.numeric(subset(temp34, Group.ITT == "SHR, ITT>=27")[,2])
t3 <- as.numeric(subset(temp34, Group.ITT == "LR, ITT<=26")[,2])

t.test(t1, t2, var.equal = T, alternative ="greater")
t.test(t2, t3, var.equal = T)
t.test(t1, t3, var.equal = T, alternative ="greater")

roi_diff_reg <- NA

for (j in 1:nrow(roi_diff21)) {
  
  y1d <- y_diff[roi_diff21$Row[j], roi_diff21$Column[j], ]
  yd <- subset(y1d, names(y1d) %in% temp34$MRN)
  
  dat10 <- cbind.data.frame(temp34, "Connectivity difference"= yd, 
                            "WJ3AudWM"= temp34$wj3.audwm.ss)
  #dat10 <- cbind.data.frame(temp34, "Connectivity difference"= yd, 
  #                          "Wechsler"= temp34$wech.totalds.scs)
  
  ft <- summary(lm(WJ3AudWM ~ yd*Group.ITT, data = dat10))
  
  Fpv <- pf(ft$fstatistic[1], ft$fstatistic[2], ft$fstatistic[3], lower.tail = F)
  
  wmem_edgeID <- ifelse(roi_diff21$Edges[j] %in% wmem_edges, "TRUE", "FALSE")
  
  if(Fpv<= 0.05) {
    roi_diff_reg <- rbind.data.frame(roi_diff_reg, data.frame("Edge"= roi_diff21$Edges[j],
                                                              "Edge_WM"=wmem_edgeID,
                                                              "Row"= roi_diff21$Row[j],
                                                              "Column"= roi_diff21$Column[j],
                                                              "Connectivity"= ft$coefficients[2,1],
                                                              "pvalue"= ft$coefficients[2,4],
                                                              "SHR,ITTle26"=ft$coefficients[3,1],
                                                              "pvalue"= ft$coefficients[3,4],
                                                              "SHR,ITTge27"= ft$coefficients[4,1],
                                                              "pvalue"= ft$coefficients[4,4],
                                                              "Yd:SHR,ITTle26"=ft$coefficients[5,1],
                                                              "pvalue"= ft$coefficients[5,4],
                                                              "Yd:SHR,ITTge27"=ft$coefficients[6,1],
                                                              "pvalue"= ft$coefficients[6,4],
                                                              "Adj.R2"= ft$adj.r.squared,
                                                              "Fpvalue"=Fpv))
    
  }
}

roi_diff_reg <- roi_diff_reg[-1, ]

#write.csv( roi_diff_reg, "ROI_regression_WechslerWM.csv")
write.csv(roi_diff_reg, "ROI_regression_WJ3AWM.csv")



################################################################################
# We check for the neurocognitive scores on ommissions based on the attention edges

# the following are the connectivity edges involved in the working memory of brain
att_edges <- c(172,125,122,67,137,456,495,463,74,169,345,458,504,410,498,126,130,496)

#roi_diff <- subset(roi_diff21, Edges %in% att_edges)

i <- 13 # obtain the CPT omissions scores: higher is worse
#i<- 14 # obtain the CPT hit reaction time score: higher is worse
#i<- 15 # obtain the CPT variability score: higher is worse
#i<- 14 # obtain the CPT detectability score: higher is worse
#i<- 14 # obtain the CPT response style score: higher is worse

temp40 <- ndat[, c(1,i)]  # choose the neurocog score and the MRN variable

# choose the neurocog scores that are non-empty
temp41 <- temp40[which(temp40[,2] != "."), ]

#mrn31 <- temp31$MRN  # get MRNs for patients whose neurocog scores are non-empty

temp42 <- X3[, c(1,8)]

temp43 <- merge(temp41, temp42, by= "MRN")
temp44 <- merge(temp43, Xdat, by= "MRN")
temp45 <- na.omit(temp44)

#t.test(as.numeric(wech.totalds.scs) ~ as.character(Risk), data=temp34)
t.test(as.numeric(cpt.omm.t) ~ as.character(Risk), data=temp45)

boxplot(as.numeric(cpt.omm.t) ~ Group.ITT, data = temp45, pch=16, cex=0.6,
        ylab ="CPT Ommissions score")

t1 <- as.numeric(subset(temp45, Group.ITT == "SHR, ITT<=26")[,2])
t2 <- as.numeric(subset(temp45, Group.ITT == "SHR, ITT>=27")[,2])
t3 <- as.numeric(subset(temp45, Group.ITT == "LR, ITT<=26")[,2])

t.test(t1, t2, var.equal = T, alternative ="greater")
t.test(t2, t3, var.equal = T)
t.test(t1, t3, var.equal = T, alternative ="greater")


roi_diff_reg <- NA

for (j in 1:nrow(roi_diff21)) {
  
  y1d <- y_diff[roi_diff21$Row[j], roi_diff21$Column[j], ]
  yd <- subset(y1d, names(y1d) %in% temp45$MRN)
  
  dat10 <- cbind.data.frame(temp45, "Connectivity difference"= yd, 
                            "CPTOmmissions"= temp45$cpt.omm.t)
  #dat10 <- cbind.data.frame(temp34, "Connectivity difference"= yd, 
  #                          "Wechsler"= temp34$wech.totalds.scs)
  
  ft <- summary(lm(CPTOmmissions ~ yd*Group.ITT, data = dat10))
  
  Fpv <- pf(ft$fstatistic[1], ft$fstatistic[2], ft$fstatistic[3], lower.tail = F)
  
  att_edgeID <- ifelse(roi_diff21$Edges[j] %in% att_edges, "TRUE", "FALSE")
  
  if(Fpv<= 0.05) {
    roi_diff_reg <- rbind.data.frame(roi_diff_reg, data.frame("Edge"= roi_diff21$Edges[j],
                                                              "Edge_Attention"=att_edgeID,
                                                              "Row"= roi_diff21$Row[j],
                                                              "Column"= roi_diff21$Column[j],
                                                              "Connectivity"= ft$coefficients[2,1],
                                                              "pvalue"= ft$coefficients[2,4],
                                                              "SHR,ITTle26"=ft$coefficients[3,1],
                                                              "pvalue"= ft$coefficients[3,4],
                                                              "SHR,ITTge27"= ft$coefficients[4,1],
                                                              "pvalue"= ft$coefficients[4,4],
                                                              "Yd:SHR,ITTle26"=ft$coefficients[5,1],
                                                              "pvalue"= ft$coefficients[5,4],
                                                              "Yd:SHR,ITTge27"=ft$coefficients[6,1],
                                                              "pvalue"= ft$coefficients[6,4],
                                                              "Adj.R2"= ft$adj.r.squared,
                                                              "Fpvalue"=Fpv))
    
  }
}

roi_diff_reg <- roi_diff_reg[-1, ]

#write.csv( roi_diff_reg, "ROI_regression_WechslerWM.csv")
write.csv(roi_diff_reg, "ROI_regression_WJ3AWM.csv")



###############################################################################
# We try to understand the neuro-cognitive score on retrieval fluency

temp52 <- X3[, c(1,6)]

temp54 <- merge(temp52, Xdat, by= "MRN")
temp55 <- na.omit(temp54)

plotdat2 <- cbind.data.frame("Retrieval fluency" = temp55$wk120.retfl.ss,
                             Risk= temp55$Risk, "Risk & ITT"= temp55$Group.ITT,
                             Age= temp55$Age, Sex= temp55$Sex)

boxplot(`Retrieval fluency` ~ `Risk & ITT`, data = plotdat2, pch=16, cex=0.6)
#abline(h=0, col="red")

boxplot(`Retrieval fluency` ~ as.character(Risk), data = plotdat2, pch=16, cex=0.6)
#abline(h=0, col="red")

boxplot(`IQ difference` ~ as.character(Sex), data = plotdat2, pch=16, cex=0.6)
abline(h=0, col="red")

t.test(`Retrieval fluency` ~ as.character(Risk), data = plotdat2, var.equal=T)

t1<- subset(plotdat2, `Risk & ITT` == "SHR, ITT<=26")[,1]
t2<- subset(plotdat2, `Risk & ITT` == "SHR, ITT>=27")[,1]
t3<- subset(plotdat2, `Risk & ITT` == "LR, ITT<=26")[,1]

t.test(t1, t2, var.equal = T)
t.test(t2, t3, var.equal = T)
t.test(t1, t3, var.equal = T)

# add columns for the regression outputs: 
roi_diff <- cbind.data.frame(roi_diff, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)

nc1 <- ncol(roi_diff)

colnames(roi_diff)[c(nc1 -5, nc1 -4, nc1 -3, nc1-2, nc1-1, nc1)] <- c(paste(cnames_ndat[i], "coeff_connectivity", sep ="_"),
                                           paste(cnames_ndat[i], "coeff_shr_ITTle26", sep = "_"),
                                           paste(cnames_ndat[i], "coeff_shr_ITTge27", sep ="_"),
                                           paste(cnames_ndat[i], "pvalue_connectivity", sep = "_"),
                                           paste(cnames_ndat[i], "pvalue_shr_ITTle26", sep = "_"),
                                           paste(cnames_ndat[i], "pvalue_shr_ITTge27", sep = "_"))

mrn1_riskitt <- cbind.data.frame(subset(Xdat, MRN %in% mrn1), 
                                  Scores= as.numeric(temp1[,2]))

png(filename = "Rplot_WMWechslerscore_RiskGrp.png", height = 5, width = 4.5, res = 300, 
    units = "in")
boxplot(Scores ~ as.character(Risk), data = mrn1_riskitt, ylab ="Wechsler scores", pch=16, cex=0.6, 
        main= "Working memory (digit span) scores", xlab = "Risk Groups")
dev.off()

png(filename = "Rplot_WMWechslerscore_RiskGrpITT.png", height = 5, width = 5.5, res = 300, 
    units = "in")
boxplot(Scores ~ Group.ITT, data = mrn1_riskitt, ylab ="Wechsler scores", pch=16, cex=0.6, 
        main= "Working memory (digit span) scores", xlab = "Risk Groups & ITT")
dev.off()


# for 5,7,19
j <- 1

y1d <- y_diff[roi_diff$Row[j], roi_diff$Column[j], ]
yd <- subset(y1d, names(y1d) %in% temp25$MRN)

# Boxplot for the Wech: working memory scores
par(mfrow=c(1,1))
png(filename = "Rplot_wech_scores_boxplot.png", height = 6, width = 6, 
     units ="in", res = 300)
boxplot(Scores ~ Group.ITT, data = mrn1_riskitt, ylab ="Wech scores", pch=16, cex=0.6, 
        main= "Working memory scores", xlab = "Risk Groups & ITT")
dev.off()

# Boxplot for the WJ3: auditory working memory scores
par(mfrow=c(1,1))
png(filename = "Rplot_wj3.audwm_scores_boxplot.png", height = 6, width = 6, 
    units ="in", res = 300)
boxplot(Scores ~ Group.ITT, data = mrn1_riskitt, ylab ="WJ3: Auditory memory scores", pch=16, cex=0.6, 
        main= "Working memory scores", xlab="Risk Groups & ITT")
dev.off()

# 
ggplot(dat1, aes(x=Connectivity, y=Scores, colour = Group.ITT)) +
  geom_point()+
  labs(y = "Wech scores", 
       title = "Edge: Left MedialArea 7p vs. Left Parieto-Occipital SulcusArea2")

ggsave(filename ="Rplot_wech_connectivity_LMedialA7p_LParietoOccipitalSA2.png", 
       height =5, width = 7, units = "in", dpi =300)

for (i2 in 1:nrow(roi_diff)) {
  
  #Y <- y_diff[roi_diff$Row[i2], roi_diff$Column[i2], ]
  #Y_shr_1 <- subset(Y, names(Y) %in% mrn1_shr_ittle26)
  #Y_shr_2 <- subset(Y, names(Y) %in% mrn1_shr_ittge27)
  
  y1d <- y_diff[roi_diff$Row[i2], roi_diff$Column[i2], ]
  yd <- subset(y1d, names(y1d) %in% mrn1_riskitt$MRN)
  
  dat1 <- cbind.data.frame(mrn1_riskitt, Connectivity= yd)
  
  #Y_new <- scale(Y[mdt])
  #ft1 <- summary(lm(score1_shr_ittle26 ~ Y_shr_1))
  #ft2 <- summary(lm(score1_shr_ittge27 ~ Y_shr_2))
  
  ft1 <- summary(lm(Scores ~ Connectivity+Group.ITT, data = dat1))
  #ft2 <- summary(lm(NScores ~ Group, data = dat1))
  
  roi_diff[i2, c(nc1 -5, nc1 -4, nc1 -3)] <- ft1$coefficients[c(2,3,4), 1]
  roi_diff[i2, c(nc1 -2, nc1 -1, nc1)] <- ft1$coefficients[c(2,3,4), 4]
  
}

write.csv(roi_diff, "Neurocog_regression3.csv")

################################################################################
# Another plotting exercise for the Neurocognitive scores and the connectivity

m21_risk_nc <- subset(m21_risk, MRN %in% mrn1_riskitt$MRN)
m21_risk_nc2 <- cbind.data.frame(m21_risk_nc[,c(1,2)], 
                                 m21_risk_nc[, paste("X", wmem_edges, sep = "")])

ncp_extr <- sapply(1:nrow(m21_risk_nc2), function(i) {
  length(which(m21_risk_nc2[i,-c(1,2)] < 0.05))
})

ncp_extr2 <- cbind.data.frame(mrn1_riskitt, ExtremeCount=ncp_extr)

ggplot(ncp_extr2, aes(x=ExtremeCount, y=Scores, colour = Risk)) +
  geom_point()+
  labs(x = "Extreme Count", y= "Wech Score",
       title= "Wech score: higher is better")

ggsave(filename ="Rplot_WM_vs_ExtremeCt.png", 
       height =5, width = 7, units = "in", dpi =300)

ggplot(ncp_extr2, aes(x=ExtremeCount, y=Scores, colour = Group.ITT)) +
  geom_point()+
  labs(x = "Extreme Count", y= "Wech Score",
       title= "Wech: higher is better")

summary(lm(Scores ~ Group.ITT, data = ncp_extr2))

for (i2 in 1:nrow(roi_diff)) {
  
  Y <- y_diff[roi_diff$Row[i2], roi_diff$Column[i2], ]
  Y_shr_1 <- subset(ncp_extr2, Group.ITT %in% mrn1_shr_ittle26)
  Y_shr_2 <- subset(Y, names(Y) %in% mrn1_shr_ittge27)
  
  #Y_new <- scale(Y[mdt])
  ft1 <- summary(lm(score1_shr_ittle26 ~ Y_shr_1))
  ft2 <- summary(lm(score1_shr_ittge27 ~ Y_shr_2))
  
  roi_diff[i2, c(nc1 -3, nc1 -2)] <- ft1$coefficients[2, c(1,4)]
  roi_diff[i2, c(nc1 -1, nc1)] <- ft2$coefficients[2, c(1,4)]
  
}


###########################
# We run the regression of the Connectivity on the Dosage scores of the SHR group of patients

xdat_shr <- subset(Xdat, MRN %in% mrn_shr)

Groups <- ifelse(xdat_shr$IT.total >= 27, ">= 27", "<= 26")
roi_diff21_cshr <- cbind.data.frame(roi_diff21, NA, NA)

for (i2 in 1:nrow(roi_diff21_cshr)) {
  
  Y <- y_diff[roi_diff21_cshr$Row[i2], roi_diff21_cshr$Column[i2], ]
  Y_new <- subset(Y, names(Y) %in% mrn_shr)
  ft1 <- summary(lm(Y_new ~ Groups))
  
  roi_diff21_cshr[i2, c(ncol(roi_diff21_cshr)-1, ncol(roi_diff21_cshr))] <- ft1$coefficients[2, c(1,4)]
  
}

roi_diff21_cshr <- subset(roi_diff21_cshr, Edges %in% roi_diff$Edges)

colnames(roi_diff21_cshr)
write.csv(roi_diff, "neurocog_regression2.csv")




################################################################################
# We try to understand the connection of the neurocognitive scores to measure 
# basal ganglia and the edges functioning on attention
################################################################################

# the following are the connectivity edges involved in the basal ganglia part of the brain
bg_edges <- c(7, 6, 3, 350, 349, 4)

# get the corresponding rows and columns of the connectivity matrix for the edges above
roi_diff <- subset(roi_diff21, Edges %in% bg_edges)

i <- 22 # neurocognitive score variable 

temp <- ndat[, c(1,i)]  # choose the neurocog score and the MRN variable

# choose the neurocog scores that are non-empty
temp1 <- temp[which(temp[,2] != "."), ]

mrn1 <- temp1$mrn  # get MRNs for patients whose neurocog scores are non-empty

roi_diff <- cbind.data.frame(roi_diff, NA, NA, NA, NA, NA, NA)

nc1 <- ncol(roi_diff)

colnames(roi_diff)[c(nc1 -5, nc1 -4, nc1 -3, nc1-2, nc1-1, nc1)] <- c(paste(cnames_ndat[i], 
                                           "coeff_connectivity", sep ="_"),
                          paste(cnames_ndat[i], "coeff_shr_ITTle26",sep = "_"),
                          paste(cnames_ndat[i], "coeff_shr_ITTge27", sep ="_"),
                          paste(cnames_ndat[i], "pvalue_connectivity", sep = "_"),
                          paste(cnames_ndat[i], "pvalue_shr_ITTle26", sep = "_"),
                          paste(cnames_ndat[i], "pvalue_shr_ITTge27", sep = "_"))

mrn1_riskitt <- cbind.data.frame(subset(Xdat, MRN %in% mrn1), 
                                 Scores= as.numeric(temp1[,2]))

mrn1_riskitt$Group.ITT <- NA

mrn1_riskitt$Group.ITT <- sapply(1:nrow(mrn1_riskitt), function(i) {
  
  if(mrn1_riskitt$Risk[i] == "LR" & mrn1_riskitt$IT.total[i] <= 26 ) {"LR, ITT<=26" }
  else if (mrn1_riskitt$Risk[i] == "SHR" & mrn1_riskitt$IT.total[i] <= 26) {"SHR, ITT<=26"}
  else if (mrn1_riskitt$Risk[i] == "SHR" & mrn1_riskitt$IT.total[i] >= 27) {"SHR, ITT>=27"}
})

# boxplots of the NC scores
boxplot(Scores ~ as.character(Risk), data = mrn1_riskitt, ylab ="Wech scores", 
        pch=16, cex=0.6, main= "Working memory scores", xlab = "Risk Groups")

boxplot(Scores ~ Group.ITT, data = mrn1_riskitt, ylab ="Wech scores", pch=16, cex=0.6, 
        main= "Working memory scores", xlab = "Risk Groups & ITT")

t1<- subset(mrn1_riskitt, Group.ITT == "SHR, ITT<=26")$Scores
t2<- subset(mrn1_riskitt, Group.ITT == "LR, ITT<=26")$Scores
t3<- subset(mrn1_riskitt, Group.ITT == "SHR, ITT>=27")$Scores

t1.1 <- subset(mrn1_riskitt, Risk == "SHR")$Scores
t2.1 <- subset(mrn1_riskitt, Risk == "LR")$Scores

# t-tests to know if there are significant differences
t.test(t1, t2, var.equal = T)
t.test(t1, t3, var.equal = T)
t.test(t2, t3, var.equal = T)

t.test(t2.1, t1.1)

m21_risk_nc <- subset(m21_risk, MRN %in% mrn1_riskitt$MRN)
m21_risk_nc2 <- cbind.data.frame(m21_risk_nc[,c(1,2)], 
                                 m21_risk_nc[, paste("X", wmem_edges, sep = "")])

ncp_extr <- sapply(1:nrow(m21_risk_nc2), function(i) {
  length(which(m21_risk_nc2[i,-c(1,2)] < 0.05))
})

ncp_extr2 <- cbind.data.frame(mrn1_riskitt, ExtremeCount=ncp_extr)

ggplot(ncp_extr2, aes(x=ExtremeCount, y=Scores, colour = Risk)) +
  geom_point()+
  labs(x = "Extreme Count", y= "Wech Score",
       title= "Wech score: higher is better")

ggplot(ncp_extr2, aes(x=ExtremeCount, y=Scores, colour = Group.ITT)) +
  geom_point()+
  labs(x = "Extreme Count", y= "Wech Score",
       title= "Wech: higher is better")


for (i2 in 1:nrow(roi_diff)) {
  
  #Y <- y_diff[roi_diff$Row[i2], roi_diff$Column[i2], ]
  #Y_shr_1 <- subset(Y, names(Y) %in% mrn1_shr_ittle26)
  #Y_shr_2 <- subset(Y, names(Y) %in% mrn1_shr_ittge27)
  
  y1d <- y_diff[roi_diff$Row[i2], roi_diff$Column[i2], ]
  yd <- subset(y1d, names(y1d) %in% mrn1_riskitt$MRN)
  
  dat1 <- cbind.data.frame(mrn1_riskitt, Connectivity= yd)
  
  #Y_new <- scale(Y[mdt])
  #ft1 <- summary(lm(score1_shr_ittle26 ~ Y_shr_1))
  #ft2 <- summary(lm(score1_shr_ittge27 ~ Y_shr_2))
  
  ft1 <- summary(lm(Scores ~ Connectivity+Group.ITT, data = dat1))
  #ft2 <- summary(lm(NScores ~ Group, data = dat1))
  
  roi_diff[i2, c(nc1 -5, nc1 -4, nc1 -3)] <- ft1$coefficients[c(2,3,4), 1]
  roi_diff[i2, c(nc1 -2, nc1 -1, nc1)] <- ft1$coefficients[c(2,3,4), 4]
  
}

j<- 1
j <- j+1

y1d <- y_diff[roi_diff$Row[j], roi_diff$Column[j], ]
yd <- subset(y1d, names(y1d) %in% mrn1_riskitt$MRN)

dat1 <- cbind.data.frame(mrn1_riskitt, Connectivity= yd)

ggplot(dat1, aes(x=Connectivity, y=Scores, colour = Risk)) +
  geom_point()+
  labs(y = "Neurocog scores", x="Connectivity difference", 
       title = "Edge: Right Area 8C vs. Right Premotor Eye Field")

ggplot(dat1, aes(x=Connectivity, y=Scores, colour = Group.ITT)) +
  geom_point()+
  labs(y = "Neurocog scores", x="Connectivity difference", 
       title = "Edge: Right Area 8C vs. Right Premotor Eye Field")


# x1 is the number of significant patients in the time 1.
x1 <- sapply(2:ncol(m_t1_risk), function(i) {

  # find out for that edge which patients are extreme
  a1 <- which(m_t1_risk[,i] < 0.05)
  
  # find out the LR among those extreme patients
  a2 <- length(which(m_t1_risk[a1, 1] == "LR"))
  
  # find out among the LR extreme patients, those are in the upper tail of the control dist
  a2_utail <- length(which(m_t1_risk[a1, 1] == "LR" & m_t1_risktail[a1, i] == "upper")) 
  
  # find out among the LR extreme patients, those are in the lower tail of the control dist
  a2_ltail <- length(which(m_t1_risk[a1, 1] == "LR" & m_t1_risktail[a1, i] == "lower")) 
  
  # find out the SHR among those extreme patients
  a3 <- length(which(m_t1_risk[a1, 1] == "SHR")) 
  
  # find out among the LR extreme patients, those are in the upper tail of the control dist
  a3_utail <- length(which(m_t1_risk[a1, 1] == "SHR" & m_t1_risktail[a1, i] == "upper")) 
  
  # find out among the LR extreme patients, those are in the lower tail of the control dist
  a3_ltail <- length(which(m_t1_risk[a1, 1] == "SHR" & m_t1_risktail[a1, i] == "lower")) 
  
  a <- matrix(c(length(a1), a2, a2_utail, a2_ltail, a3, a3_utail, a3_ltail), nrow = 1)
  
  return(a)
})

# For time = 1 we perform this operation:

# Now we shortlist those edges that have extreme patients count 70 or more 
roi_t1 <- cbind.data.frame(Edges= seq(1: nrow(y_cor_control2)), y_cor_control2[, -1])

# collect the number of extreme patients in the edge
roi_t1$ExtremePatients <- x1[1, ] 
roi_t1$LR <- x1[2, ] 
roi_t1$LR_UT <- x1[3, ]
roi_t1$LR_LT <- x1[4, ]
roi_t1$SHR <- x1[5, ] 
roi_t1$SHR_UT <- x1[6, ]
roi_t1$SHR_LT <- x1[7, ]

roi_t1$ROI2 <- roi_t1$ROI1 <- NA
roi_t1$Lobe_ROI2 <- roi_t1$Lobe_ROI1 <- NA 

roi_t1$Region_ROI2 <- roi_t1$Region_ROI1 <- NA

# collect the ROI pair in the edge
for (i in 1:nrow(roi_t1)) {
  roi_t1$ROI1[i] <- paste(node_labels[roi_t1[i,"Row"], "Side"], 
                           node_labels[roi_t1[i,"Row"], "Name"], sep =" ")
  roi_t1$ROI2[i] <- paste(node_labels[roi_t1[i,"Column"], "Side"], 
                           node_labels[roi_t1[i,"Column"], "Name"], sep =" ")
  
  roi_t1$Lobe_ROI1[i] <- node_labels[roi_t1[i,"Row"], "Lobe"]
  roi_t1$Lobe_ROI2[i] <- node_labels[roi_t1[i,"Column"], "Lobe"]
  
  roi_t1$Region_ROI1[i] <- node_labels[roi_t1[i,"Row"], "Region.Name"]
  roi_t1$Region_ROI2[i] <- node_labels[roi_t1[i,"Column"], "Region.Name"]
  
  print(i)
  #print(roi70_t1$ROI2[i])
  
}

# Now we fetch the information of what region of the brain they are in
#rownames(roi70_t1) <- NULL
roi_t1 <- roi_t1[order(-roi_t1$ExtremePatients), ]
rownames(roi_t1)<- NULL

write.csv(roi_t1, "Edges_time1.csv")

# Same as for time 1 we perform this operation for time = 2:

# x2 is the number of significant patients at the time 2.
x2 <- sapply(2:ncol(m_t2_risk), function(i) {
  
  # find out for that edge which patients are extreme
  a1 <- which(m_t2_risk[,i] < 0.05)
  
  # find out the LR among those extreme patients
  a2 <- length(which(m_t2_risk[a1, 1] == "LR"))
  
  # find out among the LR extreme patients, those are in the upper tail of the control dist
  a2_utail <- length(which(m_t2_risk[a1, 1] == "LR" & m_t2_risktail[a1, i] == "upper")) 
  
  # find out among the LR extreme patients, those are in the lower tail of the control dist
  a2_ltail <- length(which(m_t2_risk[a1, 1] == "LR" & m_t2_risktail[a1, i] == "lower")) 
  
  # find out the SHR among those extreme patients
  a3 <- length(which(m_t2_risk[a1, 1] == "SHR")) 
  
  # find out among the LR extreme patients, those are in the upper tail of the control dist
  a3_utail <- length(which(m_t2_risk[a1, 1] == "SHR" & m_t2_risktail[a1, i] == "upper")) 
  
  # find out among the LR extreme patients, those are in the lower tail of the control dist
  a3_ltail <- length(which(m_t2_risk[a1, 1] == "SHR" & m_t2_risktail[a1, i] == "lower")) 
  
  a <- matrix(c(length(a1), a2, a2_utail, a2_ltail, a3, a3_utail, a3_ltail), nrow = 1)
  
  return(a)
})

# For time = 2 we perform this operation:

roi_t2 <- cbind.data.frame(Edges= seq(1: nrow(y_cor_control2)), y_cor_control2[, -1])

# collect the number of extreme patients in the edge
roi_t2$ExtremePatients <- x2[1, ] 
roi_t2$LR <- x2[2, ] 
roi_t2$LR_UT <- x2[3, ]
roi_t2$LR_LT <- x2[4, ]
roi_t2$SHR <- x2[5, ] 
roi_t2$SHR_UT <- x2[6, ]
roi_t2$SHR_LT <- x2[7, ]

roi_t2$ROI2 <- roi_t2$ROI1 <- NA
roi_t2$Lobe_ROI2 <- roi_t2$Lobe_ROI1 <- NA

roi_t2$Region_ROI2 <- roi_t2$Region_ROI1 <- NA

# collect the ROI pair in the edge
for (i in 1:nrow(roi_t2)) {
  roi_t2$ROI1[i] <- paste(node_labels[roi_t2[i,"Row"], "Side"], 
                          node_labels[roi_t2[i,"Row"], "Name"], sep =" ")
  roi_t2$ROI2[i] <- paste(node_labels[roi_t2[i,"Column"], "Side"], 
                          node_labels[roi_t2[i,"Column"], "Name"], sep =" ")
  
  roi_t2$Lobe_ROI1[i] <- node_labels[roi_t2[i,"Row"], "Lobe"]
  roi_t2$Lobe_ROI2[i] <- node_labels[roi_t2[i,"Column"], "Lobe"]
  
  roi_t2$Region_ROI1[i] <- node_labels[roi_t2[i,"Row"], "Region.Name"]
  roi_t2$Region_ROI2[i] <- node_labels[roi_t2[i,"Column"], "Region.Name"]
  
  print(i)
  
}

# Now we fetch the information of what region of the brain they are in
#rownames(roi70_t1) <- NULL
roi_t2 <- roi_t2[order(-roi_t2$ExtremePatients), ]
rownames(roi_t2)<- NULL

write.csv(roi_t2, "Edges_time2.csv")

####################################################################
# Regression where Neurocognitive scores are the response variables 
# Brain connectivity are the covariates
####################################################################

# subset the data to include only those edges with most extreme patients at least 66
#roi_diff <- subset(roi_diff21, ExtremePatients >= 66)

# read the neurocognitive data first~
ndat <- read.csv("Neurocog_scores.csv", header = TRUE)

cnames_ndat <- colnames(ndat)

# wmem_rois <- read.csv("Working_Memory_ROI_GReddick.csv", header =TRUE)
# the following are the connectivity edges involved in the working memory of brain
wmem_edges <- c(284,39,228,157,117,55,244,315,179,105)












x1_r1prop <- sapply(2:ncol(m_t1_risk), function(i) {
  length(which((m_t1_risk$Risk_id[which(m_t1_risk[,i]<0.05)]=="Risk1")== TRUE))/length(lr)
}) 

x1_r2prop <- sapply(2:ncol(m_t1_risk), function(i) {
  length(which((m_t1_risk$Risk_id[which(m_t1_risk[,i]<0.05)]== "Risk2")== TRUE))/length(shr)
}) 

x1_prop <- data.frame(Risk= rep(c("Risk1", "Risk2"), each = length(x1_r1prop)), 
                      Proportion= c(x1_r1prop, x1_r2prop), Edge = c(1:671))

ggplot(data=x1_prop, aes(x=Edge, y=Proportion, fill= Risk)) +
  geom_bar(stat="identity") +
  labs(xlab="Brain connectivity edge", 
       title ="Extreme patients by Risk group, Time = 1")+
  theme(axis.text.x=element_text(size=12),
        axis.title.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.y=element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        title = element_text(size=12),
        strip.text = element_text(size=12))
  #ylim(c(0,1.2))
ggsave(file="Rplot_extremepatients_t1_risk.png", height = 4, width =5.5,
       units ="in", dpi= 300)


x1_r1prop <- sapply(2:ncol(m_t1_risk), function(i) {
  length(which((m_t1_risk$Risk_id[which(m_t1_risk[,i] < 0.05)]== "Risk1") == TRUE))/length(lr)
})

x2_r1prop <- sapply(2:ncol(m_t2_risk), function(i) {
  length(which((m_t2_risk$Risk_id[which(m_t2_risk[,i]<0.05)]=="Risk1")== TRUE))/length(lr)
}) 

x2_r2prop <- sapply(2:ncol(m_t2_risk), function(i) {
  length(which((m_t2_risk$Risk_id[which(m_t2_risk[,i]<0.05)]== "Risk2")== TRUE))/length(shr)
}) 


x2_prop <- data.frame(Risk= rep(c("Risk1","Risk2"), each=length(x2_r1prop)), 
                      Proportion= c(x2_r1prop, x2_r2prop), Edge= c(1:571))

ggplot(data=x2_prop, aes(x = Edge, y = Proportion, fill = Risk)) +
  geom_bar(stat="identity") +
  labs(xlab="Brain connectivity edge", 
       title ="Extreme patients by Risk group, Time = 2")+
  theme(axis.text.x=element_text(size=12),
        axis.title.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.y=element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        title = element_text(size=12),
        strip.text = element_text(size=12))+
  ylim(c(0,1.2))
ggsave(file="Rplot_extremepatients_t2_risk.png", height =4, width =5.5,
       units ="in", dpi= 300)

# Each row summarizes the number of edges for the given person, hence corresponding 
# to each column summary, we have the # of edges that are significant in each row.

# x1_pat is the number of people significant/extreme at the time 1.
x1_pat <- sapply(1:nrow(m_t1_risk), function(i) {
  length(which(as.numeric(m_t1_risk[i,-1]) <0.05 ))
})

x1_edge <- sapply(2:ncol(m_t1_risk), function(i) {
  length(which(as.numeric(m_t1_risk[,i]) <0.05 ))
})

# x21_pat is the number of people significant/extreme at the time difference 2 - 1.
x21_pat <- sapply(1:nrow(m21_risk), function(i) {
  length(which(as.numeric(m21_risk[i,-1]) < 0.05 ))
})

# x21_edge is the number of people significant/extreme at the time difference 2 - 1.
x21_edge <- sapply(2:ncol(m21_risk), function(i) {
  length(which(as.numeric(m21_risk[,i]) < 0.05 ))
})

# x2_pat is the number of people significant/extreme at the time 2.
x2_pat <- sapply(1:nrow(m_t2_risk), function(i) {
  length(which(as.numeric(m_t2_risk[i,-1]) <0.05))
})

x2_edge <- sapply(2:ncol(m_t2_risk), function(i) {
  length(which(as.numeric(m_t2_risk[,i]) <0.05 ))
})

par(m)

par(mfrow= c(1,2))
# regress the count of extreme patients on the covariates, Time = 1
fit1 <- glm(x1_pat ~ Xdat$IT.total[-c(1:59)], family = poisson(link ="log"))
pred_glm <- predict(fit1, type = "response")

summary(fit1)

plot(x1_pat ~ Xdat$IT.total[-c(1:59)], xlab ="ITT", ylab= "extreme patients at time 1",
     pch=16, cex=0.6, ylim= c(20, 250), xlim= c(8, 30))
lines(pred_glm ~ Xdat$IT.total[-c(1:59)], col = 2)
#lines(x1_pat ~ Xdat$IT.total[-c(1:59)], col=4)


# regress the count of extreme patients on the covariates, Time = 2 
fit2 <- glm(x2_pat ~ Xdat$IT.total[-c(1:59)], family = poisson(link ="log"))
pred_glm <- predict(fit2, type = "response")

summary(fit2)

plot(x2_pat ~ Xdat$IT.total[-c(1:59)], xlab ="ITT", ylab= "extreme patients at time 2",
     pch=16, cex=0.6, ylim= c(20, 250), xlim= c(8, 30))
lines(pred_glm ~ Xdat$IT.total[-c(1:59)], col = 2)
#lines(x1_pat ~ Xdat$IT.total[-c(1:59)], col=4)
#legend("topleft", legend = c("poisson GLM"))


# fitted for the covariate Age:

fit1 <- glm(x1_pat ~ Xdat$Age[-c(1:59)], family = poisson(link ="log"))
pred_glm <- predict(fit1, type = "response")

summary(fit1)

plot(x1_pat ~ Xdat$Age[-c(1:59)], xlab ="Age", ylab= "extreme patients at time 1",
     pch=16, cex=0.6, ylim= c(27, 250), xlim= c(0, 22))
lines(pred_glm ~ Xdat$Age[-c(1:59)], col = 2)
#lines(x1_pat ~ Xdat$IT.total[-c(1:59)], col=4)


# regress the count of extreme patients on the covariates, Time = 2 
fit2 <- glm(x2_pat ~ Xdat$Age[-c(1:59)], family = poisson(link ="log"))
pred_glm <- predict(fit2, type = "response")

summary(fit2)

plot(x2_pat ~ Xdat$Age[-c(1:59)], xlab ="Age", ylab= "extreme patients at time 2",
     pch=16, cex=0.6, ylim= c(27, 250), xlim= c(0, 22))
lines(pred_glm ~ Xdat$Age[-c(1:59)], col = 2)
#lines(x1_pat ~ Xdat$IT.total[-c(1:59)], col=4)


# fit for the covariate Risk
fit1 <- glm(x1_pat ~ Xdat$Risk[-c(1:59)], family = poisson(link ="log"))
pred_glm1 <- predict(fit1, type = "response")

summary(fit1)

# regress the count of extreme patients on the covariates, Time = 2 
fit2 <- glm(x2_pat ~ Xdat$Risk[-c(1:59)], family = poisson(link ="log"))
pred_glm2 <- predict(fit2, type = "response")

summary(fit2)

plot(x1_pat ~ Xdat$Risk[-c(1:59)], xlab ="Risk", ylab= "Extreme patients at time 1",
     pch=16, cex=0.6)
lines(pred_glm1 ~ Xdat$Risk[-c(1:59)], col = 2)
#lines(x1_pat ~ Xdat$IT.total[-c(1:59)], col=4)
#legend("topleft", legend = c("poisson GLM"))

plot(x2_pat ~ Xdat$Risk[-c(1:59)], xlab ="Risk", ylab= "Extreme patients at time 2",
     pch=16, cex=0.6)
lines(pred_glm2 ~ Xdat$Risk[-c(1:59)], col = 2)
#lines(x1_pat ~ Xdat$IT.total[-c(1:59)], col=4)
#legend("topleft", legend = c("poisson GLM"))


Xdat$Sex <- as.factor(ifelse(Xdat$Sex==2, "M", "F"))

fit1 <- glm(x1_pat ~ Xdat$Sex[-c(1:59)], family = poisson(link ="log"))
pred_glm1 <- predict(fit1, type = "response")

summary(fit1)

# regress the count of extreme patients on the covariates, Time = 2 
fit2 <- glm(x2_pat ~ Xdat$Sex[-c(1:59)], family = poisson(link ="log"))
pred_glm2 <- predict(fit2, type = "response")

summary(fit2)

plot(x1_pat ~ Xdat$Sex[-c(1:59)], xlab ="Sex", ylab= "Extreme patients at time 1",
     pch=16, cex=0.6)
lines(pred_glm1 ~ Xdat$Sex[-c(1:59)], col = 2)
#lines(x1_pat ~ Xdat$IT.total[-c(1:59)], col=4)

plot(x2_pat ~ Xdat$Sex[-c(1:59)], xlab ="Sex", ylab= "Extreme patients at time 2",
     pch=16, cex=0.6)
lines(pred_glm2 ~ Xdat$Risk[-c(1:59)], col = 2)
#lines(x1_pat ~ Xdat$IT.total[-c(1:59)], col=4)


edge_names <- sapply(1:nrow(y_cor_control2), function(i) {
  paste("V_",paste(y_cor_control2$Row[i], y_cor_control2$Column[i], sep ="_"), 
        sep="")
})

colnames(m_t1_risk) <- c("Risk.ID", edge_names)
colnames(m_t2_risk) <- c("Risk.ID", edge_names)

#write.csv(m_t1_risk, "Risk_edges_matrix_T1.csv")
#write.csv(m_t2_risk, "Risk_edges_matrix_T2.csv")

# Each row summarizes the number of edges for the given person, hence corresponding 
# to each column summary, we have the # of edges that are significant in each row.

aed_t1 <- NA
aed_t2 <- NA

for (i in 1:nrow(m_t1_risk)) {
  aed_t1 <- c(aed_t1, c(y_cor_control2[which(m_t1_risk[i,-1] < 0.05), c("Row")], 
                       y_cor_control2[which(m_t1_risk[i,-1] < 0.05), c("Column")]))
}

for (i in 1:nrow(m_t2_risk)) {
  aed_t2 <- c(aed_t2, c(y_cor_control2[which(m_t2_risk[i, -1] < 0.05), c("Row")], 
                       y_cor_control2[which(m_t2_risk[i, -1] < 0.05), c("Column")]))
}

aed_t1 <- data.frame(table(aed_t1))
aed_t2 <- data.frame(table(aed_t2))

# for the patients participating in the studies
paed_t1 <- NA
paed_t2 <- NA

for (i in 2:ncol(m_t1_risk)) {
  paed_t1 <- c(paed_t1, which(m_t1_risk[,i] < 0.05))
}

for (i in 2:ncol(m_t2_risk)) {
  paed_t2 <- c(paed_t2, which(m_t2_risk[,i] < 0.05))
}

paed_t1 <- data.frame(table(paed_t1), Risk= Xdat$Risk[-c(1:59)],
                      Sex= Xdat$Sex[-c(1:59)], Age= Xdat$Age[-c(1:59)],
                      "ITTotal"= Xdat$IT.total[-c(1:59)])

paed_t2 <- data.frame(table(paed_t2), Risk= Xdat$Risk[-c(1:59)],
                      Sex= Xdat$Sex[-c(1:59)], Age= Xdat$Age[-c(1:59)],
                      "ITTotal"= Xdat$IT.total[-c(1:59)])

paed_t1$Risk<- as.factor(paed_t1$Risk)
paed_t2$Risk<- as.factor(paed_t2$Risk)

paed_t1$Sex <- as.factor(ifelse(paed_t1$Sex==1, 'F', 'M'))
paed_t2$Sex <- as.factor(ifelse(paed_t2$Sex==1, 'F', 'M'))

paed_t1$paed_t1 <- as.numeric(paed_t1$paed_t1)
paed_t2$paed_t2 <- as.numeric(paed_t2$paed_t2)

paed <- cbind.data.frame(Time= as.factor(rep(c(1,2), each= nrow(paed_t1))) , 
                  rbind.data.frame(paed_t1[,-1], paed_t2[,-1]))

#ggplot(paed_t1, aes(x=paed_t1, y=Freq)) + 
#geom_segment( aes(x=paed_t1, xend=paed_t1, y=0, yend=Freq)) + 
#geom_point(aes(color = Sex),size=1.5) +
#  labs(y = "Count", x = "Patient", 
#       title = "Number of extreme edges for each patient at Time = 1") +
#  theme(axis.text.x=element_text(size=12),
#        axis.title.x=element_text(size=12),
#        axis.text.y=element_text(size=12),
#        axis.title.y=element_text(size=12),
#        legend.text=element_text(size=12),
#        legend.title=element_text(size=12),
#        title = element_text(size=12),
#        strip.text = element_text(size=12))

#ggsave(file="Rplot_extremeedges_t1_gender.png", height =4.5, width =4.5,
#       units ="in", dpi= 300)

#ggplot(paed_t2, aes(x=paed_t2, y=Freq)) + 
#geom_segment( aes(x=paed_t2, xend=paed_t2, y=0, yend=Freq)) + 
#geom_point(aes(color = Sex),size=1.5) +
#  labs(y = "Count", x = "Patient", 
#       title = "Number of extreme edges for each patient at Time = 2") +
#  theme(axis.text.x=element_text(size=12),
#        axis.title.x=element_text(size=12),
#        axis.text.y=element_text(size=12),
#        axis.title.y=element_text(size=12),
#        legend.text=element_text(size=12),
#        legend.title=element_text(size=12),
#        title = element_text(size=12),
#        strip.text = element_text(size=12))

#ggsave(file="Rplot_extremeedges_t2_gender.png", height =4.5, width =4.5,
#       units ="in", dpi= 300)


ggplot(paed_t1, aes(x=paed_t1, y=Freq)) + 
  geom_segment( aes(x=paed_t1, xend=paed_t1, y=0, yend=Freq), alpha=0.5) + 
  geom_point(aes(color = Risk),size=1) +
  labs(y = "Count", x = "Patient", 
       title = "Extreme brain edges, Time = 1") +
  ylim(c(0,240))+
  theme(axis.text.x=element_text(size=12),
        axis.title.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.y=element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        title = element_text(size=12),
        strip.text = element_text(size=12))

ggsave(file="Rplot_extremeedges_t1_risk.png", height =4, width = 5,
       units ="in", dpi= 300)

r1t1<- subset(paed_t1, Risk==1)[,2]
r2t1<- subset(paed_t1, Risk==2)[,2]

t.test(r1t1, r2t1)

ggplot(paed_t2, aes(x=paed_t2, y=Freq)) + 
  geom_segment( aes(x=paed_t2, xend=paed_t2, y=0, yend=Freq), alpha= 0.5) + 
  geom_point(aes(color = Risk),size=1) +
  labs(y = "Count", x ="Patient", 
       title = "Extreme brain edges, Time = 2") +
  ylim(c(0,240))+
  theme(axis.text.x=element_text(size=12),
        axis.title.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.y=element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        title = element_text(size=12),
        strip.text = element_text(size=12))

ggsave(file="Rplot_extremeedges_t2_risk.png", height = 4, width = 5,
       units ="in", dpi= 300)

r1t2<- subset(paed_t2, Risk==1)[,2]
r2t2<- subset(paed_t2, Risk==2)[,2]

t.test(r1t2, r2t2)

aed_t1 <- cbind.data.frame(aed_t1, Lobe= node_labels$Lobe[aed_t1$aed_t1],
                           Region= node_labels$Region.Name[aed_t1$aed_t1])

aed_t2 <- cbind.data.frame(aed_t2, Lobe= node_labels$Lobe[aed_t2$aed_t2],
                           Region= node_labels$Region.Name[aed_t2$aed_t2])

####
# Plot of the extreme edges per patient classified by the IT total score of treatment
####

## Time = 1
ggplot(paed_t1, aes(x=paed_t1, y=Freq)) + 
  geom_segment( aes(x=paed_t1, xend=paed_t1, y=0, yend=Freq),
                alpha=0.5) + 
  geom_point(aes(color = ITTotal),size=1) +
  scale_colour_gradientn(colors=terrain.colors(10))+
  labs(y = "Count", x = "Patient", 
       title = "Extreme brain edges, Time = 1") +
  ylim(c(0,240))+
  theme(axis.text.x=element_text(size=12),
        axis.title.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.y=element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        title = element_text(size=12),
        strip.text = element_text(size=12)) 

ggsave(file="Rplot_extremeedges_t1_ITTotal.png", height =4, width =5,
       units ="in", dpi= 300)

## Time = 2
ggplot(paed_t2, aes(x=paed_t2, y=Freq)) + 
  geom_segment( aes(x=paed_t2, xend=paed_t2, y=0, yend=Freq),
                alpha=0.5) + 
  geom_point(aes(color = ITTotal), size=1) +
  scale_colour_gradientn(colors = terrain.colors(10))+
  labs(y = "Count", x = "Patient", 
       title = "Extreme brain edges, Time = 2") +
  ylim(c(0,240))+
  theme(axis.text.x=element_text(size=12),
        axis.title.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.y=element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        title = element_text(size=12),
        strip.text = element_text(size=12)) 

ggsave(file="Rplot_extremeedges_t2_ITTotal.png", height =4, width =5,
       units ="in", dpi= 300)


################################################################################
# study of relation of the number of brain connectivity edges on ITTotal score
##

dev.off()
par(mfrow = c(1,1))

# for time = 1
png(filename = "Rplot_extremeedges_byITT_t1.png", 
    width = 600, height = 500, units ="px", pointsize =16)

plot(paed_t1$ITTotal, paed_t1$Freq, pch=16, cex= 0.6,
     xlab= "IT Total", ylab ="Count of extreme edges",
     main = "Connectivity edges vs ITTotal, Time 1",
     ylim = c(10, 250))
ft <- lm(paed_t1$Freq ~ paed_t1$ITTotal)
abline(ft, col= "red")
abline(a=ft$coefficients[1], b=1)

legend("topright", legend = c("Simple Regression", "y=x"), 
       col = c("red", "black"), 
       lty = 1)

dev.off()


# for time = 1, with age in the x-axis
png(filename = "Rplot_extremeedges_byAge_t1.png", 
    width = 600, height = 500, units ="px", pointsize =16)

plot(paed_t1$Age, paed_t1$Freq, pch=16, cex= 0.6,
     xlab= "IT Total", ylab ="Count of extreme edges",
     main = "Connectivity edges vs ITTotal, Time 1",
     ylim = c(10, 250))
ft <- lm(paed_t1$Freq ~ paed_t1$Age)
abline(ft, col= "red")
abline(a=ft$coefficients[1], b=1)

legend("topright", legend = c("Simple Regression", "y=x"), 
       col = c("red", "black"), 
       lty = 1)

dev.off()


# for time = 2
png(filename = "Rplot_extremeedges_byITT_t2.png", 
    width = 600, height = 500, units ="px", pointsize =16)

plot(paed_t2$ITTotal, paed_t2$Freq, pch=16, cex= 0.6,
     xlab= "IT Total", ylab ="Count of extreme brain edges",
     main = "Connectivity edges by ITTotal, Time 2",
     ylim = c(10, 250))
ft <- lm(paed_t2$Freq ~ paed_t2$ITTotal)
abline(ft, col= "red")
abline(a=ft$coefficients[1], b=1)

legend("topright", legend = c("Simple Regression", "y=x"), 
       col = c("red", "black"), 
       lty = 1)
dev.off()

node_seq<- seq(1,379)
ndls1 <- node_seq[!(node_seq %in% aed_t1$aed_t1)]
View(node_labels[ndls1,])

ndls2 <- node_seq[!(node_seq %in% aed_t2$aed_t2)]
View(node_labels[ndls2,])

png(filename = "Rplot_ROIs_inextremeedges_t1.png", 
    width = 600, height = 500, units ="px", pointsize =16)
barplot(Freq~ aed_t1, data = aed_t1, xlab= "ROIs at time 1", ylab="Frequency",
        ylim=c(0,550))
abline(h=300, col="red")
dev.off()


png(filename = "Rplot_ROIs_inextremeedges_t2.png", 
    width = 600, height = 500, units ="px", pointsize =16)
barplot(Freq~ aed_t2, data = aed_t2, xlab= "ROIs at time 2", ylab="Frequency",
        ylim=c(0, 550))
abline(h=300, col="red")

intg_t1<- as.integer(subset(aed_t1, Freq > 250)[,1])
node_labels[intg_t1, -c(1,3)]

intg_t2 <- as.integer(subset(aed_t2, Freq > 250)[,1])
node_labels[intg_t2, -c(1,3)]


png(filename = "Rplot_summary_bypatients.png", 
    width = 500, height = 500, units ="px", pointsize =16)
plot(x1_edge, x2_edge, pch=20, cex=.4, xlab="Extreme patients at time 1", 
     ylab="Extreme patients at time 2", xlim =c(0, 200), ylim =c(0,200),
     main= "Extreme patients in all edges")
abline(lm(x2_edge ~ x1_edge))
abline(a=0, b=1, col="red")
legend("topleft", legend = c("Simple Regression", "y=x"), 
       col = c("black", "red"), 
       lty = 1)
dev.off()

# t.test(x1, x2, paired = TRUE, alternative ="greater")

png(filename = "Rplot_summary_byedges.png", 
    width = 500, height = 500, units ="px", pointsize =16)
plot(x1_pat, x2_pat, pch=20, cex=.4, 
     xlab="Extreme edges at time 1", 
     ylab="Extreme edges at time 2", 
     main="Extreme edges in all patients", xlim = c(20, 250),
     ylim =c(20, 250))
abline(lm(x2_pat ~ x1_pat))
abline(a=0, b=1, col="red")
legend("topleft", legend = c("Simple Regression", "y=x"), 
       col = c("black", "red"), 
       lty = 1)
dev.off()

t.test(x1_pat, x2_pat, paired = TRUE, alternative ="greater")

## Running the logistic regression on the weights:

m_t1_risk_bin<- apply(m_t1_risk[,-1], c(1,2), function(x) {
  ifelse(x< 0.05, 1, 0)
})

m_t2_risk_bin<- apply(m_t2_risk[,-1], c(1,2), function(x) {
  ifelse(x< 0.05, 1, 0)
})

logistic_intercepts_coeff<- apply(m_t1_risk_bin, 2, function(x) {
  summary(glm(x ~ Xdat$Sex[-c(1:59)] + scale(Xdat$Age[-c(1:59)])+ 
                scale(Xdat$IT.total[-c(1:59)]), 
              family ="binomial"))$coefficients[1,1]
  
})

logistic_intercepts_pv<- apply(m_t1_risk_bin, 2, function(x) {
  summary(glm(x ~ Xdat$Sex[-c(1:59)] + scale(Xdat$Age[-c(1:59)])+ 
                scale(Xdat$IT.total[-c(1:59)]), 
              family ="binomial"))$coefficients[1,4]
  
})

par(mfrow=c(1,2))
hist(logistic_intercepts_coeff, breaks =15,
     col="lightgreen",
     border="black",
     #prob = TRUE,
     xlab = "Intercept",
     main = "Coefficient of Intercept, Time 1" )
#abline(v=-log(mean(y_cor_control$Correlation)+2, base=10), col="red")
grid(nx = NULL, ny = NULL, col = "lightgrey", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)

hist(log(logistic_intercepts_pv, base=10), breaks =15,
     col="lightgreen",
     border="black",
     #prob = TRUE,
     xlab = "Intercept",
     main = "p-value of intercept, Time 1")
abline(v=log(0.05, base=10), col="red")
grid(nx = NULL, ny = NULL, col = "lightgrey", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)


logistic_age_coeff<- apply(m_t1_risk_bin, 2, function(x) {
  summary(glm(x ~ Xdat$Sex[-c(1:59)] + scale(Xdat$Age[-c(1:59)])+ 
                scale(Xdat$IT.total[-c(1:59)]), 
              family ="binomial"))$coefficients[3,1]
  
})

logistic_age_pv <- apply(m_t1_risk_bin, 2, function(x) {
  summary(glm(x ~ Xdat$Sex[-c(1:59)] + scale(Xdat$Age[-c(1:59)])+ 
                scale(Xdat$IT.total[-c(1:59)]), 
              family ="binomial"))$coefficients[3,4]
  
})

hist(logistic_age_coeff, breaks =15,
     col="lightgreen",
     border="black",
     #prob = TRUE,
     xlab = "Intercept",
     main = "Coefficient of Age, Time 1" )
#abline(v=-log(mean(y_cor_control$Correlation)+2, base=10), col="red")
grid(nx = NULL, ny = NULL, col = "lightgrey", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)

hist(log(logistic_age_pv, base=10), breaks =15,
     col="lightgreen",
     border="black",
     #prob = TRUE,
     xlab = "Intercept",
     main = "p-value of Age, Time 1")
abline(v=log(0.05, base=10), col="red")
grid(nx = NULL, ny = NULL, col = "lightgrey", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)


logistic_sex_coeff<- apply(m_t1_risk_bin, 2, function(x) {
  summary(glm(x ~ Xdat$Sex[-c(1:59)] + scale(Xdat$Age[-c(1:59)])+ 
                scale(Xdat$IT.total[-c(1:59)]), 
              family ="binomial"))$coefficients[2,1]
})

logistic_sex_pv <- apply(m_t1_risk_bin, 2, function(x) {
  summary(glm(x ~ Xdat$Sex[-c(1:59)] + scale(Xdat$Age[-c(1:59)])+ 
                scale(Xdat$IT.total[-c(1:59)]), 
              family ="binomial"))$coefficients[2,4]
})

hist(logistic_sex_coeff, breaks =15,
     col="lightgreen",
     border="black",
     #prob = TRUE,
     xlab = "Intercept",
     main = "Coefficient of Sex, Time 1" )
#abline(v=-log(mean(y_cor_control$Correlation)+2, base=10), col="red")
grid(nx = NULL, ny = NULL, col = "lightgrey", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)

hist(log(logistic_sex_pv, base=10), breaks =15,
     col="lightgreen",
     border="black",
     #prob = TRUE,
     xlab = "Intercept",
     main = "p-value of Sex, Time 1")
abline(v=log(0.05, base=10), col="red")
grid(nx = NULL, ny = NULL, col = "lightgrey", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)

logistic_itt_coeff<- apply(m_t1_risk_bin, 2, function(x) {
  summary(glm(x ~ Xdat$Sex[-c(1:59)] + scale(Xdat$Age[-c(1:59)])+ 
                scale(Xdat$IT.total[-c(1:59)]), 
              family ="binomial"))$coefficients[4,1]
})

logistic_itt_pv <- apply(m_t1_risk_bin, 2, function(x) {
  summary(glm(x ~ Xdat$Sex[-c(1:59)] + scale(Xdat$Age[-c(1:59)])+ 
                scale(Xdat$IT.total[-c(1:59)]), 
              family ="binomial"))$coefficients[4,4]
})

hist(logistic_itt_coeff, breaks =15,
     col="lightgreen",
     border="black",
     #prob = TRUE,
     xlab = "Intercept",
     main = "Coefficient of IT Total, Time 1" )
#abline(v=-log(mean(y_cor_control$Correlation)+2, base=10), col="red")
grid(nx = NULL, ny = NULL, col = "lightgrey", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)

hist(log(logistic_itt_pv, base=10), breaks =15,
     col="lightgreen",
     border="black",
     #prob = TRUE,
     xlab = "Intercept",
     main = "p-value of IT Total, Time 1")
abline(v=log(0.05, base=10), col="red")
grid(nx = NULL, ny = NULL, col = "lightgrey", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)

# No evaluation for the code section below:
mydat1 <- cbind.data.frame( m_t1_risk, Xdat[-cn,])  # at the first time 
mydat2 <- cbind.data.frame( m_t2_risk, Xdat[-cn,])  # at 120 weeks time

Intercept<- rep(NA, 576)
Sex2<- rep(NA, 576)
Age<- rep(NA, 576)
Risk2<- rep(NA, 576)
IT.total<- rep(NA, 576)

Intercept_pv<- rep(NA, 576)
Sex2_pv<- rep(NA, 576)
Age_pv<- rep(NA, 576)
Risk2_pv<- rep(NA, 576)
IT.total_pv<- rep(NA, 576)

for (i in 2:ncol(m_t1_risk)) {
  
  mydat<- cbind.data.frame(Y= as.numeric(m_t1_risk[,i]), Xdat[-cn, c("Sex","Risk","IT.total","Age")])
  mydat$Risk<-as.factor(mydat$Risk)
  m<- glm(Y ~ Sex+ Age+ Risk+ IT.total, data =mydat, family = binomial(link ="logit"))
  
  Intercept[i]<- m$coefficients[1]
  Sex2[i]<- m$coefficients[2]
  Age[i]<- m$coefficients[3]
  Risk2[i]<- m$coefficients[4]
  IT.total[i]<- m$coefficients[5]
  
  Intercept_pv[i]<- m$coefficients[1]
  Sex2_pv[i]<- m$coefficients[2]
  Age_pv[i]<- m$coefficients[3]
  Risk2_pv[i]<- m$coefficients[4]
  IT.total_pv[i]<- m$coefficients[5]
  
}

par(mfrow=c(1,2))
hist(Intercept[-1], breaks =15,
     col="lightgreen",
     border="black",
     #prob = TRUE,
     xlab = "Intercept",
     main = "" )
#abline(v=-log(mean(y_cor_control$Correlation)+2, base=10), col="red")
grid(nx = NULL, ny = NULL, col = "lightgrey", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)

hist(Sex2[-1], breaks =15,
     col="lightgreen",
     border="black",
     #prob = TRUE,
     xlab = "Sex",
     main = "")
#abline(v=-log(mean(y_cor_control$Correlation)+2, base=10), col="red")
grid(nx = NULL, ny = NULL, col = "lightgrey", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)

hist(Age[-1], breaks =15,
     col="lightgreen",
     border="black",
     #prob = TRUE,
     xlab = "Age",
     main = "" )
#abline(v=-log(mean(y_cor_control$Correlation)+2, base=10), col="red")
grid(nx = NULL, ny = NULL, col = "lightgrey", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)

hist(Risk2[-1], breaks =15,
     col="lightgreen",
     border="black",
     #prob = TRUE,
     xlab = "Risk Groups",
     main = "" )
#abline(v=-log(mean(y_cor_control$Correlation)+2, base=10), col="red")
grid(nx = NULL, ny = NULL, col = "lightgrey", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)

hist(IT.total[-1], breaks =15,
     col="lightgreen",
     border="black",
     #prob = TRUE,
     xlab = "IT.Total",
     main = "" )
#abline(v=-log(mean(y_cor_control$Correlation)+2, base=10), col="red")
grid(nx = NULL, ny = NULL, col = "lightgrey", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)

### Now we perform the clustering analysis for the connectivity data as follows:

# the following function computes the euclidean norm of a vector
euclidean_norm <- function(x) {
  x<- matrix(x, ncol =1)
  nr<- t(x) %*% x
  return(sqrt(nr))
}

mw_t1 <- matrix(NA, nrow= dim(y_1)[3], ncol= nrow(y_cor_control2))
mw_t2 <- mw_t1

for (i2 in 1:ncol(mw_t1)) {
  mw_t1[, i2] <- y_1[y_cor_control2$Row[i2], y_cor_control2$Column[i2],]
}

for (i2 in 1:ncol(mw_t1)) {
  mw_t2[, i2] <- y_2[y_cor_control2$Row[i2], y_cor_control2$Column[i2],]
}

mw_t1_risk <- mw_t1[-c(1:59),]
mw_t2_risk <- mw_t2[-c(1:59),]

kclust_patient_t1 <- kmeans(x= mw_t1_risk, centers = 2, nstart =25, iter.max = 30)
#totwss_patient <- kclust_patient$tot.withinss

kclust_patient_t2 <- kmeans(x= mw_t2_risk, centers = 2, nstart =25, iter.max = 30)

aov_df <- data.frame(ITTotal= Xdat$IT.total[-c(1:59)], Risk= Xdat$Risk[-c(1:59)], Cluster=kclust_patient_t1$cluster, Sex= Xdat$Sex[-c(1:59)], Age= Xdat$Age[-c(1:59)])

fit<- lm(ITTotal ~ Risk+ Cluster+ Sex+ Age, data= aov_df)
anova(fit)
summary(fit)



##################################################################################
# Compute the Hellinger distances between the connectivity distributions among the
# groups Control, LR and SHR with their further subdivisions based on the ITT score.
###################################################################################

# we create a function to compute the Hellinger distance between the distributions:

# Hellinger distances computation between C and LR groups
# at time = 2
hdist2_C_LR <- sapply(1:nrow(y_cor_control2), function(i) {
  l11 <- y_2[y_cor_control2$Row[i], y_cor_control2$Column[i], ]
  x11 <- as.numeric(l11[cn])
  y11 <- as.numeric(l11[lr])
  statip::hellinger(x11, y11, lower = 0, upper = Inf)
})

# at time = 1
hdist1_C_LR <- sapply(1:nrow(y_cor_control2), function(i) {
  l11<- y_1[y_cor_control2$Row[i], y_cor_control2$Column[i], ]
  x11<- as.numeric(l11[cn])
  y11<- as.numeric(l11[lr])
  statip::hellinger(x11, y11, lower = 0, upper = Inf)
})

# create data frame for the plots
hdat_C_LR <- data.frame(HDist = c(hdist1_C_LR, 
                                       hdist2_C_LR),
                          Time = c(rep("1", nrow(y_cor_control2)), 
                                   rep("2", nrow(y_cor_control2))) )

ggplot(hdat_C_LR, aes(y= HDist, x= Time)) +
  geom_boxplot(aes(fill = Time), color= "black", 
                 position = "identity") +
  labs(y="Hellinger Distance",
       title = "Hellinger distance between C & SHR") +
  scale_fill_brewer(labels=unname(TeX(c("$T_1$", "$T_2$")) ), palette = "Pastel1")+
  guides(x = guide_axis(title = NULL, position = "NULL"))+
  theme_bw()


# Hellinger distances computation between C and SHR groups
# at time = 2
hdist2_C_SHR <- sapply(1:nrow(y_cor_control2), function(i) {
  l11 <- y_2[y_cor_control2$Row[i], y_cor_control2$Column[i], ]
  x11 <- as.numeric(l11[cn])
  y11 <- as.numeric(l11[shr])
  statip::hellinger(x11, y11, lower = 0, upper = Inf)
})

# at time = 1
hdist1_C_SHR <- sapply(1:nrow(y_cor_control2), function(i) {
  l11<- y_1[y_cor_control2$Row[i], y_cor_control2$Column[i], ]
  x11<- as.numeric(l11[cn])
  y11<- as.numeric(l11[shr])
  hdist(x11, y11, lower = 0, upper = Inf)
})

# create data frame for the plots
hdat_C_SHR <- data.frame(HDist = c(hdist1_C_SHR, 
                                  hdist2_C_SHR),
                        Time = c(rep("1", nrow(y_cor_control2)), 
                                 rep("2", nrow(y_cor_control2))) )

ggplot(hdat_C_SHR, aes(y= HDist, x= Time, alpha =0.6)) +
  geom_boxplot(aes(fill = Time), 
                 position = "identity") +
  labs(x="Hellinger Distance",
       title = "Hellinger distance between C & SHR") +
  scale_fill_brewer(labels= )



# Hellinger distances computation between LR and SHR groups
# at time = 2
hdist2_LR_SHR <- sapply(1:nrow(y_cor_control2), function(i) {
  l11<- y_2[y_cor_control2$Row[i], y_cor_control2$Column[i], ]
  x11<- l11[lr]
  y11<- l11[shr]
  statip::hellinger(x11, y11, lower = 0, upper = Inf)
})

# at time = 1
hdist1_LR_SHR <- sapply(1:nrow(y_cor_control2), function(i) {
  l11<- y_1[y_cor_control2$Row[i], y_cor_control2$Column[i], ]
  x11<- l11[lr]
  y11<- l11[shr]
  statip::hellinger(x11, y11)
})

# create data frame for the plots
hdat_LR_SHR <- data.frame(HDist = c(hdist1_LR_SHR, 
                                       hdist2_LR_SHR),
                          Time = c(rep("1", nrow(y_cor_control2)), 
                                   rep("2", nrow(y_cor_control2))) )

ggplot(hdat_LR_SHR, aes(x= HDist)) +
  geom_histogram(aes(fill = Time, color= Time, alpha =0.5), color= "black",
                 position = "identity") +
  labs(x="Hellinger Distance",
       title = "Distance between LR & SHR")


# Hellinger distance between Control LR with Low Dose <= 21 & Control groups
hdist2_C_LRLD21 <- sapply(1:nrow(y_cor_control2), function(i) {
  l11 <- y_2[y_cor_control2$Row[i], y_cor_control2$Column[i], ]
  x11 <- l11[cn]
  y11 <- l11[lrld21]
  statip::hellinger(x11, y11)
})

hdist1_C_LRLD21 <- sapply(1:nrow(y_cor_control2), function(i) {
  l11<- y_1[y_cor_control2$Row[i], y_cor_control2$Column[i], ]
  x11<- l11[cn]
  y11<- l11[lrld21]
  statip::hellinger(x11, y11)
})

hdat_C_LRLD21 <- data.frame(HDist= c(hdist1_C_LRLD21, 
                                             hdist2_C_LRLD21),
                                Time= c(rep("1", nrow(y_cor_control2)), 
                                        rep("2", nrow(y_cor_control2))) )

ggplot(hdat_C_LRLD21, aes(x= HDist)) +
  geom_histogram(aes(fill = Time, color= Time, alpha =0.5), color= "black",
               position = "identity") +
  labs(y="2-Wasserstein Distance",
       title = "Distance between C & LR with LD")

# Hellinger distance between Risk LR with Low Dose >= 21 & Control groups
hdist2_C_LRHD21 <- sapply(1:nrow(y_cor_control2), function(i) {
  l11<- y_2[y_cor_control2$Row[i], y_cor_control2$Column[i], ]
  x11<- l11[cn]
  y11<- l11[lrhd21]
  statip::hellinger(x11, y11, lower = 0, upper = Inf)
})

hdist1_C_LRHD21 <- sapply(1:nrow(y_cor_control2), function(i) {
  l11<- y_1[y_cor_control2$Row[i], y_cor_control2$Column[i], ]
  x11<- l11[cn]
  y11<- l11[lrhd21]
  hdist(x11, y11, lower = 0, upper = Inf)
})

hdat_C_LRHD21 <- data.frame(HDist= c(hdist1_C_LRHD21, 
                                                  hdist2_C_LRHD21),
                                     Time= c(rep("1", nrow(y_cor_control2)), 
                                             rep("2", nrow(y_cor_control2))) )

# Hellinger distance between Risk LR with Low Dose >= 21 & LR with high dose groups
hdist2_LRLD21_LRHD21 <- sapply(1:nrow(y_cor_control2), function(i) {
  l11<- y_2[y_cor_control2$Row[i], y_cor_control2$Column[i], ]
  x11<- l11[lrld21]
  y11<- l11[lrhd21]
  statip::hellinger(x11, y11, lower = 0, upper = Inf)
})

hdist1_LRLD21_LRHD21 <- sapply(1:nrow(y_cor_control2), function(i) {
  l11<- y_1[y_cor_control2$Row[i], y_cor_control2$Column[i], ]
  x11<- l11[lrld21]
  y11<- l11[lrhd21]
  statip::hellinger(x11, y11, lower = 0, upper = Inf)
})

hdat_LRLD21_LRHD21 <- data.frame(HDist= c(hdist1_LRLD21_LRHD21, 
                                     hdist2_LRLD21_LRHD21),
                            Time= c(rep("1", nrow(y_cor_control2)), 
                                    rep("2", nrow(y_cor_control2))) )


boxplot(HDist ~ Time, data = hdat_C_LRHD21, ylab="2-Wasserstein Distance",
        main= "Distance: C vs LR at HD", pch=20, cex=.5)


# Hellinger distance between Control and SHR with Low Dose >= 21
hdist2_C_SHRLD27 <- sapply(1:nrow(y_cor_control2), function(i) {
  l11<- y_2[y_cor_control2$Row[i], y_cor_control2$Column[i], ]
  x11<- l11[cn]
  y11<- l11[shrld27]
  statip::hellinger(x11, y11, lower = 0, upper = Inf)
})

hdist1_C_SHRLD27 <- sapply(1:nrow(y_cor_control2), function(i) {
  l11<- y_1[y_cor_control2$Row[i], y_cor_control2$Column[i], ]
  x11<- l11[cn]
  y11<- l11[shrld27]
  statip::hellinger(x11, y11, lower = 0, upper = Inf)
})

hdat_C_SHRLD27 <- data.frame(HDist= c(hdist1_C_SHRLD27, 
                                                  hdist2_C_SHRLD27),
                                     Time= c(rep("1", nrow(y_cor_control2)), 
                                             rep("2", nrow(y_cor_control2))) )

ggplot(hdat_C_SHRLD27, aes(x= HDist)) +
  geom_histogram(aes(fill = Time, color= Time, alpha =0.5), color= "black",
               position = "identity") +
  labs(x="Hellinger Distance",
       title = "Distance between C & SHR at LD")

boxplot(HDist ~ Time, data = hdat_C_SHRLD27, ylab="Hellinger Distance",
        main= "Distance: C vs LR at HD", pch=20, cex=.5)

# Hellinger distance between SHR with High Dose >= 27 & Control groups
# at time = 2
hdist2_C_SHRHD27 <- sapply(1:nrow(y_cor_control2), function(i) {
  l11<- y_2[y_cor_control2$Row[i], y_cor_control2$Column[i], ]
  x11<- l11[cn]
  y11<- l11[shrhd27]
  statip::hellinger(x11, y11, lower =0, upper = Inf)
})

# at time = 1 (hdist function was used because round-off error was reported while 
# using the statip::hellinger function)
hdist1_C_SHRHD27 <- sapply(1:nrow(y_cor_control2), function(i) {
  l11<- y_1[y_cor_control2$Row[i], y_cor_control2$Column[i], ]
  x11<- l11[cn]
  y11<- l11[shrhd27]
  hdist(x11, y11, lower =0, upper = Inf)
})

hdat_C_SHRHD27<- data.frame(HDist = c(hdist1_C_SHRHD27, 
                                                   hdist2_C_SHRHD27),
                                      Time= c(rep("1", nrow(y_cor_control2)), 
                                              rep("2", nrow(y_cor_control2))) )

ggplot(hdat_C_SHRHD27, aes(x= HDist)) +
  geom_histogram(aes(fill = Time, color= Time, alpha =0.5), color= "black",
               position = "identity") +
  labs(y="Hellinger Distance",
       title = "Distance between C & SHR at HD")

ggplot(hdat_C_SHRHD27, aes(x= Time, y= HDist, fill = Time)) +
  geom_boxplot()+
  labs(y ="Hellinger Distance", title = "Distance: C vs SHR at HD")

# Hellinger distance between SHR LD group & SHR HD group
# at time = 2
hdist2_SHRLD_SHRHD27 <- sapply(1:nrow(y_cor_control2), function(i) {
  l11 <- y_2[y_cor_control2$Row[i], y_cor_control2$Column[i], ]
  x11 <- l11[shrld27]
  y11 <- l11[shrhd27]
  statip::hellinger(x11, y11, lower =0, upper = Inf)
})

# at time = 1 (hdist function was used because round-off error was reported while 
# using the statip::hellinger function)
hdist1_SHRLD_SHRHD27 <- sapply(1:nrow(y_cor_control2), function(i) {
  l11<- y_1[y_cor_control2$Row[i], y_cor_control2$Column[i], ]
  x11<- l11[shrld27]
  y11<- l11[shrhd27]
  statip::hellinger(x11, y11, lower =0, upper = Inf)
})

hdat_SHRLD_SHRHD27 <- data.frame(HDist = c(hdist1_SHRLD_SHRHD27, 
                                      hdist2_SHRLD_SHRHD27),
                            Time= c(rep("1", nrow(y_cor_control2)), 
                                    rep("2", nrow(y_cor_control2))) )

ggplot(hdat_SHRLD_SHRHD27, aes(x= HDist)) +
  geom_histogram(aes(fill = Time, color= Time, alpha =0.5), color= "black",
                 position = "identity") +
  labs(y="Hellinger Distance",
       title = "Distance between C & SHR(HD)")

ggplot(hdat_SHRLD_SHRHD27, aes(x= Time, y= HDist, fill = Time)) +
  geom_boxplot()+
  labs(y ="Hellinger Distance", title = "Distance: C vs SHR(HD)")+
  theme_bw()


# perform t-test
t.test(hdist1_C_LR, hdist2_C_LR, paierd=TRUE, alternative ="greater")

t.test(hdist1_C_SHR, hdist2_C_SHR, paierd=TRUE, alternative = "greater")

t.test(hdist1_LR_SHR, hdist2_LR_SHR, paierd=TRUE, alternative ="greater")

t.test(hdist1_C_LRLD21, hdist2_C_LRLD21, paired =TRUE, alternative ="greater")

t.test(hdist1_LRLD21_LRHD21, hdist2_LRLD21_LRHD21, paired =TRUE, alternative ="greater")

t.test(hdist1_C_LRHD21, hdist2_C_LRHD21, paired =TRUE, alternative ="greater")

t.test(hdist1_C_SHRLD27, hdist2_C_SHRLD27, paired =TRUE, alternative ="greater")

t.test(hdist1_C_SHRHD27, hdist2_C_SHRHD27, paired =TRUE,
       alternative= "greater")

t.test(hdist1_SHRLD_SHRHD27, hdist2_SHRLD_SHRHD27, paired =TRUE,
       alternative= "greater")


hdat_all <- data.frame(rbind.data.frame(hdat_C_LR, hdat_C_SHR, hdat_LR_SHR,
                             hdat_C_SHRLD27, hdat_C_SHRHD27,
                             hdat_SHRLD_SHRHD27), 
                       Index= rep(c("C vs LR", "C vs SHR","LR vs SHR",
                                    "C vs SHR(LD)", "C vs SHR(HD)",
                                    "SHR(LD) vs SHR(HD)"), each= nrow(hdat_C_LR)),
                       pvals= rep(c("paired t-test p= 0.79","paired t-test p << 0.05", "paired t-test p << 0.05",
                                  "paired t-test p << 0.05","paired t-test p << 0.05","paired t-test p= 0.28"), 
                                  each= nrow(hdat_C_LR)) )

hdat_all$Index <- factor(hdat_all$Index, levels= c("C vs LR", "C vs SHR","LR vs SHR",
                                                   "C vs SHR(LD)", "C vs SHR(HD)",
                                                   "SHR(LD) vs SHR(HD)"))

ggplot(hdat_all, aes(x= Time, y= HDist, fill = Time)) +
  geom_boxplot(outlier.size = 0.5)+
  labs(y= "Hellinger distance")+
  facet_wrap(~ Index+pvals)+
  scale_fill_brewer(labels=unname(TeX(c("$T_1$", "$T_2$")) ), palette= "Pastel1") +
  guides(x = guide_axis(title = NULL, position = "NULL"))+
  theme_bw()
  
ggsave(file="Rplot_hellinger_distance.png", height =6, width =7, 
       units ="in", dpi= 300)


df <- data.frame(Gene, Age, Group)
df$Group <- as.factor(df$Group)

mybreaks <- seq(min(df$Age)-1, to=max(df$Age)+10, by=10)
df$groups_age <- cut(df$Age, breaks = mybreaks, by=10)

bp <- ggplot(df, aes(x=groups_age, y=Gene, group=groups_age)) + 
  geom_boxplot(aes(fill=groups_age)) + 
  facet_grid(. ~ Group)

bp

pval <- df %>%
  group_by(Group) %>%
  summarize(Kruskal_pvalue = kruskal.test(Gene ~ groups_age)$p.value)

# This is to create new labels for the facetgrid where we can show the phenotype and the KW pvalue.
labels <- c(paste('Group 1\n KW p-val:', signif(subset(pval$Kruskal_pvalue, pval$Group=="Group1"), digits = 3)),
            paste('Group 2\n  KW p-val:', signif(subset(pval$Kruskal_pvalue, pval$Group=="Group2"), digits = 3)),
            paste('Group 3\n  KW p-val:', signif(subset(pval$Kruskal_pvalue, pval$Group=="Group3"), digits = 3)))

df$KW <- factor(df$Group, levels = levels(df$Group), labels = labels)

bp <- ggplot(df, aes(x=groups_age, y=Gene, group=groups_age)) + 
  geom_boxplot(aes(fill=groups_age)) + 
  facet_grid(. ~ KW) +
  theme(legend.position="none")
bp


library(dplyr)
pval <- df %>%
  group_by(Group) %>%
  summarize(Kruskal_pvalue = kruskal.test(Gene ~ groups_age)$p.value)

# This is to create new labels for the facetgrid where we can show the phenotype and the KW pvalue.
labels <- c(paste('Group 1\n KW p-val:', signif(subset(pval$Kruskal_pvalue, pval$Group=="Group1"), digits = 3)),
            paste('Group 2\n  KW p-val:', signif(subset(pval$Kruskal_pvalue, pval$Group=="Group2"), digits = 3)),
            paste('Group 3\n  KW p-val:', signif(subset(pval$Kruskal_pvalue, pval$Group=="Group3"), digits = 3)))

df$KW <- factor(df$Group, levels = levels(df$Group), labels = labels)



# No evaluation for the code section below:
#y_cor_control2$SCT_RG1_t1<- sapply(1:length(lst_condt_RG1_t1_pv), function(i) {
#  length(which(lst_condt_RG1_t1_pv[[i]]<0.05))/length(lst_condt_RG1_t1_pv[[i]]) }
#)

#y_cor_control2$SCT_RG1_t2<- sapply(1:length(lst_condt_RG1_t2_pv), function(i) {
#  length(which(lst_condt_RG1_t2_pv[[i]]<0.05))/length(lst_condt_RG1_t2_pv[[i]]) }
#)

#y_cor_control2$SCT_RG2_t1<- sapply(1:length(lst_condt_RG2_t1_pv), function(i) {
#  length(which(lst_condt_RG2_t1_pv[[i]]<0.05))/length(lst_condt_RG2_t1_pv[[i]]) }
#)

#y_cor_control2$SCT_RG2_t2<- sapply(1:length(lst_condt_RG2_t2_pv), function(i) {
#  length(which(lst_condt_RG2_t2_pv[[i]]<0.05))/length(lst_condt_RG2_t2_pv[[i]]) }
#)

#y_cor_control2$Tail_RG1_t1<- sapply(1:length(lst_condt_RG1_t1_tail), function(i) {
#  length(which(lst_condt_RG1_t1_tail[[i]]=="upper"))/length(lst_condt_RG1_t1_tail[[i]]) }
#)

#y_cor_control2$Tail_RG1_t2<- sapply(1:length(lst_condt_RG1_t2_tail), function(i) {
#  length(which(lst_condt_RG1_t2_tail[[i]]=="upper"))/length(lst_condt_RG1_t2_tail[[i]]) }
#)

#y_cor_control2$Tail_RG2_t1<- sapply(1:length(lst_condt_RG2_t1_pv), function(i) {
#  length(which(lst_condt_RG2_t1_tail[[i]]=="upper"))/length(lst_condt_RG2_t1_tail[[i]]) }
#)

#y_cor_control2$Tail_RG2_t2<- sapply(1:length(lst_condt_RG2_t2_pv), function(i) {
#  length(which(lst_condt_RG2_t2_tail[[i]]=="upper"))/length(lst_condt_RG2_t2_tail[[i]]) }
#)

#y_cor_control2<- y_cor_control2[,-1]


#ggplot(df_plot_SCT_t1, aes(y = Row, x = Column)) +
#  geom_point(aes(color=SCT), size=0.4) +
  #scale_size_manual(values=c(1,1))+
  #scale_color_manual(breaks= c("0< to <10%", "No 0 weight"),
  #                  values = c("lightgreen", "red3")) +
#  scale_color_gradientn(colors = viridis(10)) +
#  facet_wrap(~ Risk) +
#  labs(y = "Row", x = "Column", title = "P-values of patients at Start", color="Signifcant Percentage") +
#  theme(axis.text.x=element_text(size=12),
#        axis.title.x=element_text(size=12),
#        axis.text.y=element_text(size=12),
#        axis.title.y=element_text(size=12),
#        legend.text=element_text(size=12),
#        legend.title=element_text(size=12),
#        title = element_text(size=12),
#        strip.text = element_text(size=12))

#ggplot(df_plot_SCT_t2, aes(y = Row, x = Column)) +
#  geom_point(aes(color=SCT), size=0.4) +
  #scale_size_manual(values=c(1,1))+
  #scale_color_manual(breaks= c("0< to <10%", "No 0 weight"),
  #                  values = c("lightgreen", "red3")) +
#  scale_color_gradientn(colors = viridis(10)) +
#  facet_wrap(~ Risk) +
#  labs(y = "Row", x = "Column", title = "P-values of patients at Week 120", color="Significant Percentage") +
#  theme(axis.text.x=element_text(size=12),
#        axis.title.x=element_text(size=12),
#        axis.text.y=element_text(size=12),
#        axis.title.y=element_text(size=12),
#        legend.text=element_text(size=12),
#        legend.title=element_text(size=12),
#        title = element_text(size=12),
#        strip.text = element_text(size=12))


# No evaluation for the code section below:
#pvalue_table_t1<- sapply(1:nrow(y_cor_control2), function(i) {
#  x_1<- lst_condt_RG1_t1_sig[i]
#  y_1<- lst_condt_RG2_t1_sig[i]
#  dat<- chisq.test(matrix(c(x_1, y_1, length(lr)-x_1, length(shr)-y_1), nrow=2))
#  dat$p.value
#})

#pvalue_table_t2<- sapply(1:nrow(y_cor_control2), function(i) {
#  x_1<- lst_condt_RG1_t2_sig[i]
#  y_1<- lst_condt_RG2_t2_sig[i]
#  dat<- chisq.test(matrix(c(x_1, y_1, length(lr)-x_1, length(shr)-y_1), nrow=2))
#  dat$p.value
#})

#l11<- (which(pvalue_table_t1=="NaN"))
#l22<- length(which(pvalue_table_t2=="NaN"))

#pvalue_table_t1 <- data.frame(pvalue_table_t1)
#pvalue_table_t1$Significance <- ifelse(pvalue_table_t1[,1] < 0.05, "Yes","No")

#pvalue_table_t2 <- data.frame(pvalue_table_t2)
#pvalue_table_t2$Significance <- ifelse(pvalue_table_t2[,1] < 0.05, "Yes","No")

#ggdensity(-log(pvalue_table_t1,10), 
          #x = "weight",
#          add = "mean", rug = TRUE,
          #color = "Risk", fill = "Risk", 
#          title = "Weight distributions classified w.r.t. Risk at (8,9)",
#          xlab= "-log10(p-value)"
          #palette = c("#00AFBB", "#E7B800","Red")
#)

#ggdensity(-log(pvalue_table_t2,10), 
          #x = "weight",
#          add = "mean", rug = TRUE,
          #color = "Risk", fill = "Risk", 
#          title = "Weight distributions classified w.r.t. Risk at (8,9) at 120 weeks",
#          xlab= "-log10(p-value)"
          #palette = c("#00AFBB", "#E7B800","Red")
#)


# No evaluation for the code section below:
#df_plot_tail_t1<- data.frame(Row = rep(c(y_cor_control1$Row),2) , 
#                             Column= rep(c(y_cor_control1$Column),2),       
#                             Tail= c(rep(c(y_cor_control2$Tail_RG1_t1), 2), rep(y_cor_control2$Tail_RG2_t1, 2)), 
#                             Risk= rep(c("Risk Group 1", "Risk Group 2"), each= nrow(y_cor_control1) ))

#df_plot_tail_t2<- data.frame(Row = rep(c(y_cor_control1$Row),2) , 
#                             Column= rep(c(y_cor_control1$Column),2),       
#                             Tail= c(rep(c(y_cor_control2$Tail_RG1_t2), 2), rep(y_cor_control2$Tail_RG2_t2, 2)), 
#                             Risk= rep(c("Risk Group 1", "Risk Group 2"), each= nrow(y_cor_control1) ))

#ggplot(df_plot_tail_t1, aes(y = Row, x = Column)) +
#  geom_point(aes(color=Tail), size=0.4) +
  #scale_size_manual(values=c(1,1))+
  #scale_color_manual(breaks= c("0< to <10%", "No 0 weight"),
  #                  values = c("lightgreen", "red3")) +
#  scale_color_gradientn(colors = viridis(10)) +
#  facet_wrap(~ Risk) +
#  labs(y = "Row", x = "Column", title = "Percentage of patients on the upper tail at Start", color="Percentage") +
#  theme(axis.text.x=element_text(size=12),
#        axis.title.x=element_text(size=12),
#        axis.text.y=element_text(size=12),
#        axis.title.y=element_text(size=12),
#        legend.text=element_text(size=12),
#        legend.title=element_text(size=12),
#        title = element_text(size=12),
#        strip.text = element_text(size=12))

#ggplot(df_plot_tail_t2, aes(y = Row, x = Column)) +
#  geom_point(aes(color=Tail), size=0.4) +
  #scale_size_manual(values=c(1,1))+
  #scale_color_manual(breaks= c("0< to <10%", "No 0 weight"),
  #                  values = c("lightgreen", "red3")) +
#  scale_color_gradientn(colors = viridis(10)) +
#  facet_wrap(~ Risk) +
#  labs(y = "Row", x = "Column", title = "Percentage of patients on the upper tail at Week 120", color="Percentage") +
#  theme(axis.text.x=element_text(size=12),
#        axis.title.x=element_text(size=12),
#        axis.text.y=element_text(size=12),
#        axis.title.y=element_text(size=12),
#        legend.text=element_text(size=12),
#        legend.title=element_text(size=12),
#        title = element_text(size=12),
#        strip.text = element_text(size=12))

#full_ut_r2_t1<- which(y_cor_control2$Tail_RG2_t1> 0.9)
#full_ut_r1_t1<- which(y_cor_control2$Tail_RG1_t1> 0.9)
#full_ut_r1_t2<- which(y_cor_control2$Tail_RG1_t2> 0.9)
#full_ut_r2_t2<- which(y_cor_control2$Tail_RG2_t2> 0.8)

#Risk<- rep(NA, dim(y_2)[3])
#Risk[cn]<- "Control"; Risk[lr] <- "R1"; Risk[shr] <- "R2"
#y_cor_control2<- cbind.data.frame(y_cor_control2, Risk)

# Now we consider the nodes where majority of the patients ($>80\%$) in the Risk 
# groups 1 and 2 are in the upper/lower tail of the control distribuion, i.e. the 
# patients are significantly different from the distribution. This we consider for 
# both the start and the week 120 of the experiment. 

# No evaluation for the code section below:
#a11<- which(y_cor_control2$SCT_RG1_t1>.6 & y_cor_control2$Tail_RG1_t1> .8)
#a12<- which(y_cor_control2$SCT_RG1_t2> .6 & y_cor_control2$Tail_RG1_t2> .8)
#a21<- which(y_cor_control2$SCT_RG2_t1> .6 & y_cor_control2$Tail_RG2_t1> .8)
#a22<- which(y_cor_control2$SCT_RG2_t2> .6 & y_cor_control2$Tail_RG2_t2< .4)
#a11
#a12
#a21
#a22
#aec<- intersect(c(c(a11,a12),a21), a22)

#y_cor_control3<- y_cor_control2[aec,c("Row","Column")]
#y_cor_control3<- data.frame(Row= c(y_cor_control3$Row, y_cor_control3$Column),
#                            Column= c(y_cor_control3$Column, y_cor_control3$Row))

#ggplot(y_cor_control3, aes(y = Row, x = Column)) +
#  geom_point(size=0.4) +
#  labs(y = "Row", x = "Column", 
#       title = "Nodes all weights significantly on the upper tail of the control distribution") +
#  theme(axis.text.x=element_text(size=12),
#        axis.title.x=element_text(size=12),
#        axis.text.y=element_text(size=12),
#        axis.title.y=element_text(size=12),
#        legend.text=element_text(size=12),
#        legend.title=element_text(size=12),
#        title = element_text(size=12),
#        strip.text = element_text(size=12))

# Now we compute the p-values of the chi-square tests for the 2x2 contingency tables described as follows: In rows we have the patient counts at the low (top row) and high (bottom row) risk repectively, on the columns we have count of non-significant (left column), significant (right column) respectively.

#dat11<- cbind.data.frame(y_cor_control2[,c("Row","Column")], pValue=-log(pvalue_table_t1,10) ) 
#dat11<- data.frame(Row= c(dat11$Row, dat11$Column), Column= c(dat11$Column, dat11$Row),
#                   pValue= c(dat11$pValue, dat11$pValue) )
#dat22<- cbind.data.frame(y_cor_control2[,c("Row","Column")], pValue=-log(pvalue_table_t2,10) )
#dat22<- data.frame(Row= c(dat22$Row, dat22$Column), Column= c(dat22$Column, dat22$Row),
#                   pValue= c(dat22$pValue, dat22$pValue) )

#ggplot(dat11, aes(y = Row, x = Column)) +
#  geom_point(aes(color=pValue), size=0.4) +
#  scale_color_gradientn(colors = viridis(10)) +
#  labs(y = "Row", x = "Column", title = "P-values of Chi-Square test for the patients at Start", color="P-Value") +
#  theme(axis.text.x=element_text(size=12),
#        axis.title.x=element_text(size=12),
#        axis.text.y=element_text(size=12),
#        axis.title.y=element_text(size=12),
#        legend.text=element_text(size=12),
#        legend.title=element_text(size=12),
#        title = element_text(size=12),
#        strip.text = element_text(size=12))

#ggplot(dat22, aes(y = Row, x = Column)) +
#  geom_point(aes(color=pValue), size=0.4) +
#  scale_color_gradientn(colors = viridis(10)) +
#  labs(y = "Row", x = "Column", title = "P-values of Chi-Square test for the patients at Week 120", color="P-Value") +
#  theme(axis.text.x=element_text(size=12),
#        axis.title.x=element_text(size=12),
#        axis.text.y=element_text(size=12),
#        axis.title.y=element_text(size=12),
#        legend.text=element_text(size=12),
#        legend.title=element_text(size=12),
#        title = element_text(size=12),
#        strip.text = element_text(size=12))





