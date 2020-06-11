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


##############################################################################

### This is the function to compute CPO.
### mchain is a list of MCMC samples,
### the first of which is initialization.

compute_CPO <- function(mchain,data){
    L = length(mchain)
    S = dim(data$n)[1]; T = dim(data$n)[2]
    cpo = matrix(0, S, T)
    for(s in 1:S){
        for(t in 1:T){
            n = data$n[s,t];
            N = data$N[s,t];
            cpo[s,t] = (L-1) / sum(sapply(2:L, function(l) 
                1/dbinom(n,N,mchain[[l]]$p[s,t])))
        }
    }
    return(cpo)
}