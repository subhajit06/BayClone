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
  
### This is the function to generate 
### random beta-dirichlet variables.
### It returns a n-by-length(b) matrix.

rdirichlet <- function(n, a){
    l = length(a)
    x = matrix(rgamma(n*l, a), ncol=l, byrow=TRUE)
    return(x/rowSums(x))
}

rbetaDir <- function(n,a1,a2,b){	
	y = rbeta(n,a1,a2);

    if(length(b)==2){
        h = rbeta(n,b[1],b[2]);
        x = cbind(h,1-h);
    } else x = rdirichlet(n,b);
    
    if(n>1) t = apply(x, 2, function(tmp) tmp*y)
    else t = y*x
    
    return(t);
}