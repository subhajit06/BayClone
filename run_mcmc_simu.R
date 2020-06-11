#############################################################################
  
  # BayClone :  A Bayesian Feature Allocation Model for Inference of Tumor 
  # Subclones Using Next-Generation Sequencing Data
  # Copyright (C) <2014>  <Jin Wang>
  # This file is part of BayClone.

  #  BayClone is free software: you can redistribute it and/or modify
  #  it under the terms of the GNU General Public License as published by
  #  the Free Software Foundation, either version 3 of the License, or
  #  (at your option) any later version.

  #  BayClone is distributed in the hope that it will be useful,
  #  but WITHOUT ANY WARRANTY; without even the implied warranty of
  #  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  #  GNU General Public License for more details.

  #  You should have received a copy of the GNU General Public License
  #  along with BayClone. If not, see <http://www.gnu.org/licenses/>.
  ############################################################################  
  
### Run MCMC on simulation data for a list of C
### return posterior estimates as well as heatmaps of Z and w
### Based on "BayClone" [Sengupta, 2014] paper

require("gplots")
source("mcmc_cIBP.R")
source("compute_CPO.R")

###############################################
### data loading
load("simu_data.RData")
S = dim(data$n)[1]; T = dim(data$n)[2]

###############################################
### initialize hyperparameters
hparam = NULL
hparam$alpha0 = 1; 
hparam$beta0 = c(2,2);
hparam$a00 = 1;
hparam$b00 = 100;    

hparam$sigma_p0 = 0.2;
hparam$sigma_theta = 0.1;

# mcmc variables
mcmcVar = NULL
mcmcVar$n_burnin = 3000;
mcmcVar$n_samples = 1000;
mcmcVar$n_thin = 3;

###############################################
### apply MCMC on different C 
CList = c(3,4,5,6,7);
LPsML = NULL; pList = NULL; wList = NULL; ZList = NULL; p0List = NULL; i=0;

seed = 26675 # set seed for simulation retrieve

for(C in CList){
    i = i+1;
    hparam$alpha = rep(1,C);

    # main function
    mchain = mcmc_cIBP(data, mcmcVar, hparam, seed);
    
    # compute LPsML
    cpo = compute_CPO(mchain,data);    
    LPsML[i] = sum(log(cpo))

    # posterior mean
    p_all = sapply(mchain, function(x) x$p)[,-1]
    pList[[i]] = matrix(rowMeans(p_all),S,T)    

    w_all = sapply(mchain, function(x) x$w)[,-1]
    wList[[i]] = matrix(rowMeans(w_all),T,C)

    p0List[i] = mean(sapply(mchain, function(x) x$p0)[-1])

    # closest Z matrix to posterior mean
    Z_all = sapply(mchain, function(x) x$Z)[,-1]
    tmp = apply(Z_all, 1, function(x){
        dist = c(sum(abs(x)), sum(abs(x-0.5)), sum(abs(x-1)))
        (which(dist==min(dist))[1]-1)*0.5
        })
    ZList[[i]] = matrix(tmp,S,C)    
}    

print(LPsML)


###############################################
### heatmap of Z
C_len = length(CList)
mycol = rev(terrain.colors(3))

par(mfrow=c(1,C_len+1),mar=c(1,1,1,1),oma=c(0,0,2,0))
image(t(data$Z)[-1,S:1],axes=FALSE,col=mycol)
title("True Z")
for(i in 1:C_len){
    image(t(ZList[[i]])[-1,S:1],axes=FALSE,col=mycol)
    title(paste("C =",CList[i]-1))
}
mtext("Estimated Z", outer = TRUE, cex = 1)

# reorder Z
ZList.new = ZList
ZList.new[[1]] = ZList[[1]][,c(1,3,2)]
ZList.new[[2]] = ZList[[2]][,c(1,2,3,4)]
ZList.new[[3]] = ZList[[3]][,c(1,2,3,5,4)]
ZList.new[[4]] = ZList[[4]][,c(1,2,3,4,5,6)]
ZList.new[[5]] = ZList[[5]][,c(1,2,4,3,6,7,5)]

dev.new()
par(mfrow=c(1,C_len+1),mar=c(1,1,1,1),oma=c(0,0,2,0))
image(t(data$Z)[-1,S:1],axes=FALSE,col=mycol)
title("True Z")
for(i in 1:C_len){
    image(t(ZList.new[[i]])[-1,S:1],axes=FALSE,col=mycol)
    title(paste("C =",CList[i]-1))
}
mtext("Reordered Z", outer = TRUE, cex = 1)

###############################################
### heatmap of w
mycol = colorRampPalette(c("white","steelblue"))(256)
w.new = wList[[3]]

dev.new()
par(cex.main=1)
heatmap.2(w.new[,c(2,3,5,4)],Rowv=NA,Colv=NA,dendrogram="none", 
    key=TRUE,keysize=0.5,density.info="none",margin=c(3,4),trace="none",vline=NULL,hline=NULL,
    col=mycol,cexCol=1.5,cexRow=1.2,srtCol=0, xlab="Subclones", ylab="Samples",main="Estimated w",
    lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(0.8,6,1.5),lwid = c(0.1,4))


dev.new()
par(cex.main=1)
heatmap.2(data$w[,-1],Rowv=NA,Colv=NA,dendrogram="none", 
    key=TRUE,keysize=0.5,density.info="none",margin=c(3,4),trace="none",vline=NULL,hline=NULL,
    col=mycol,cexCol=1.5,cexRow=1.2,srtCol=0, xlab="Subclones", ylab="Samples",main="True w",
    lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(0.8,6,1.5),lwid = c(0.1,4))

