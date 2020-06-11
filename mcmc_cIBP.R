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
  
### MCMC sampler for finite cIBP,
### return a list of MCMC samples,
### based on "BayClone" [Sengupta, 2014] model.

source("rbetaDir.R")

mcmc_cIBP <- function(data, mcmcVar, hparam, seed){
    if(!missing(seed)) set.seed(seed);
   
    ## load the data; fixed across simulation
    N = data$N;
    n = data$n;

    ## hyperparameters
    alpha0 = hparam$alpha0; # scaler
    beta0 = hparam$beta0; # 2 by 1 vector
    a00 = hparam$a00; # scaler
    b00 = hparam$b00; # scaler
    alpha = hparam$alpha; # C by 1 vector
    sigma_theta = hparam$sigma_theta
    sigma_p0 = hparam$sigma_p0

    ## dimensions
    S = dim(n)[1]; T = dim(n)[2]
    C = length(alpha);

    ## MCMC variables
    n_burnin = mcmcVar$n_burnin;
    n_samples = mcmcVar$n_samples;
    n_thin = mcmcVar$n_thin;

    ## randomly initialize Z
    Z = matrix(0,S,C);
    Z = matrix(sample(c(0,0.5,1),S*C,replace=TRUE),S,C);
    Z[,1] = rep(1,S);

    ## initialize pi (require rbetaDir func)
    piZ = rbetaDir(C, alpha0/(C-1), 1, beta0);  
    piZ = cbind(1-rowSums(piZ), piZ)  
    
    ## initialize theta and w (same as rdirichlet func)
    theta = matrix(rgamma(T*C, alpha), ncol=C, byrow=TRUE)
    w = theta / rowSums(theta)

    ## initialize p0 and p; p is a S-by-T matrix
    p0 = rbeta(1, a00, b00)
    p = compute_p(Z, w, p0)
   
    ## record the initialize to result list mchain
    mchain = list(); # a list of length n_samples+1
    ii=1;kk=1; 

    mchain[[kk]] = list('Z'=Z, 'pi'=piZ, 'w'=w, 'p0'=p0, 'p'=p)
    kk = kk+1;

    cat('initialize done\n');

    ## loop started                           
    for(ii in 2:(n_burnin + n_samples * n_thin)){
       
        Z = update_Z(Z, piZ, w, p0, N, n);
                                                                                
        piZ = update_pi(Z, piZ, alpha0, beta0);
              
        theta = update_theta(Z, piZ, w, theta, p0, N, n, alpha, sigma_theta)
        w = theta / rowSums(theta)
   
        p0 = update_p0(Z, piZ, w, p0, N, n, a00, b00, sigma_p0)

        p = compute_p(Z, w, p0)

        ## save the result after burn in
        if((ii > n_burnin) && (ii%%n_thin == 0)){
        	mchain[[kk]] = list('Z'=Z, 'pi'=piZ, 'w'=w, 'p0'=p0, 'p'=p)
        	kk=kk+1
        }
        
        if(ii%%1000 == 0)                      
            cat('iter = ',ii,' update done\n');        
    }
    return(mchain)
}


update_Z <- function(Z, piZ, w, p0, N, n){

	S = dim(Z)[1]; C = dim(Z)[2];
	q_st = matrix(0, dim(w)[1], 3)   
    for(s in 1:S){
	    # c=1 is background, stays unchanged
	    for(c in 2:C){               
	        # q_st is a T-by-3 matrix;       
            if(C<=3) q_st[,1] = p0 * w[,1] + w[,-c(1,c)] * Z[s,-c(1,c)] else
	        q_st[,1] = p0 * w[,1] + w[,-c(1,c)] %*% Z[s,-c(1,c)];
	        q_st[,2] = q_st[,1] + 0.5 * w[,c];
	        q_st[,3] = q_st[,1] + w[,c];
	        
	        # logprob is a 3-by-1 vector
	        logprob = log(piZ[c,]) + n[s,]%*%log(q_st) + (N[s,]-n[s,])%*%log(1-q_st);
	        prob = exp(logprob-max(logprob))
	        prob = prob / sum(prob)
	        
	        # sample Z from a multinomial distribution    
	        k = which(rmultinom(1,1,prob)==1);
	        Z[s,c] = (k-1)/2;
	    }
	}
    return(Z)
}



update_pi <- function(Z, piZ, alpha0, beta0){
	C = dim(piZ)[1]
    # due to conjugacy we have a direct formula
    for(c in 2:C){
    	z = Z[,c]
    	m = c(sum(z==0.5), sum(z==1))
    	m0 = sum(z==0)
        h = rbetaDir(1, alpha0/(C-1)+sum(m), 1+m0, beta0+m)
        piZ[c,] = c(1-sum(h), h);
    }
    return(piZ)
}



update_theta <- function(Z, piZ, w, theta, p0, N, n, alpha, sigma_theta){
	p = compute_p(Z, w, p0)
	
    for(t in 1:T){
    	for(c in 1:C){
		    u = rnorm(1,0,sigma_theta);
		    prop_theta = theta;
		    prop_theta[t,c] = exp(log(theta[t,c])+u);
		    
		    tmp1 = alpha[c]*log(prop_theta[t,c])-prop_theta[t,c];
		    tmp2 = alpha[c]*log(theta[t,c])-theta[t,c];

		    # find new p_{st} for all s and this t
		    # any theta_{t,c} update lead to a p_{st} update 
		    prop_w = prop_theta[t,] / sum(prop_theta[t,])
		    prop_p = p0 * prop_w[1] + Z[,-1] %*% prop_w[-1]
		    curr_p = p[,t];
		    
		    tmp1 = tmp1 + n[,t]%*%log(prop_p) + (N[,t] - n[,t])%*%log(1-prop_p);
		    tmp2 = tmp2 + n[,t]%*%log(curr_p) + (N[,t] - n[,t])%*%log(1-curr_p);           
		 
		    v = log(runif(1))
			if(v < (tmp1 - tmp2)){
		        theta[t,c] = prop_theta[t,c];
		        p[,t] = prop_p
		    }    
    	}
    }   
    return(theta)
}


update_p0 <- function(Z, piZ, w, p0, N, n, a00, b00, sigma_p0){    

    u = rnorm(1,0,sigma_p0);
    
    prop_p0 = exp(log(p0)+u);
    if(prop_p0 >=1) prop_p0 = exp(log(p0)-u);

    prop_p = compute_p(Z, w, prop_p0);
    curr_p = compute_p(Z, w, p0);
    
    tmp1 = a00*log(prop_p0)+(b00-1)*log(1-prop_p0);
    tmp2 = a00*log(p0)+(b00-1)*log(1-p0);
    
    tmp1 = tmp1 + sum(n*log(prop_p)+(N-n)*log(1-prop_p));
    tmp2 = tmp2 + sum(n*log(curr_p)+(N-n)*log(1-curr_p));

    v = log(runif(1))
	if(v < (tmp1 - tmp2))
        p0 = prop_p0;

    return(p0)
}


compute_p <- function(Z, w, p0){
	p = t(p0 * w[,1] + t(Z[,-1] %*% t(w[,-1])))
	return(p)
}

