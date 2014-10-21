#
#     Description of this R script:
#     R tests for linear multiple output sparse group lasso routines.
#
#     Intended for use with R.
#     Copyright (C) 2014 Martin Vincent
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>
#

library(lsgl)

set.seed(100) #  ensures consistency of tests

##########################
###
### Dual test
###

# Simulate from Y=XB+E, the dimension of Y is N x K, X is N x p, B is p x K 

N <- 10 #number of samples
p <- 10 #number of features
K <- 15  #number of groups

B<-matrix(sample(c(rep(1,p^2*K*0.1),rep(0, p^2*K-as.integer(p^2*K*0.1)))),nrow=p^2,ncol=K) 

X<-matrix(rnorm(N*p,1,1),nrow=N,ncol=p)

Y<-kronecker(X,X)%*%B+matrix(rnorm(N*K,0,1),N^2,K)	


lambda<-lsgl.lambda(kron(X, X), Y, alpha=1, lambda.min=2, intercept=FALSE)

fit <-lsgl(kron(X, X), Y, alpha=1, lambda = lambda, intercept=FALSE)

## ||B - \beta||_F
if(min(sapply(fit$beta, function(beta) sum((B - beta)^2))) > 100) stop()

### Test predict
res <- predict(fit, kron(X, X))
if(min(sapply(res$Yhat, function(y) sum((Y - y)^2))) > 1e4) stop()

##########################
###
### Trple test
###

#FIXME

# Simulate from Y=XB+E, the dimension of Y is N x K, X is N x p, B is p x K 

N <- 50 #number of samples
p <- 10 #number of features
K <- 10  #number of groups

B<-matrix(sample(c(rep(1,p^3*K*0.1),rep(0, p^2*K-as.integer(p^3*K*0.1)))),nrow=p^3,ncol=K) 

X<-matrix(rnorm(N*p,1,1),nrow=N,ncol=p)

Y<-kronecker(X,kronecker(X,X))%*%B+matrix(rnorm(N*K,0,1),N^3,K)	

lambda<-lsgl.lambda(kron(X, X, X), Y, alpha=1, lambda.min=9000, intercept=FALSE)

fit <-lsgl(kron(X, X, X), Y, alpha=1, lambda = lambda, intercept=FALSE)

## ||B - \beta||_F
#if(min(sapply(fit$beta, function(beta) sum((B - beta)^2))) > 100) stop()

### Test predict
res <- predict(fit, kron(X, X))
#if(min(sapply(res$Yhat, function(y) sum((Y - y)^2))) > 1e4) stop()
