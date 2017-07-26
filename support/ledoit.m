function sigma=ledoit(ret)
% This is a support function for CovSmooth, please see help CovSmooth.
% function sigma=ledoit(ret)
% ret (t*n): t iid observations on n random variables
% sigma (n*n): invertible covariance matrix estimator
%
% This function computes the covariance matrix estimator 
% introduced by Olivier Ledoit in 
% "Portfolio Selection: Improved Covariance Matrix Estimation"
% (Job market paper, November 1994, reproduced by the UCLA
% Finance Department as Working Paper #5-96) and also
% 1st chapter of my Finance PhD thesis at MIT's Sloan
% School of Management called "Essays on Risk and Return
% in the Stock Market" (June 1995)
%
% This estimator is a weighted average of the sample
% covariance matrix and a "prior" or "shrinkage target".
% Here, the prior is given by a one-factor model.
% The factor is equal to the cross-sectional average
% of all the random variables.
% The weight, or "shrinkage intensity" is chosen to
% minimize quadratic loss measured by the Frobenius norm.
% The estimator is valid as the number of variables and/or the 
% number of observations go to infinity, but Monte-Carlo 
% simulations show that it works well for values as low as 10.
% The main advantage is that this estimator is guaranteed to
% be INVERTIBLE and well-conditioned even if variables 
% outnumber observations.
%
% Written by Olivier Ledoit (ledoit@ucla.edu) on 1/29/1996.
% Feel free to use it, modify it and/or distribute it, 
% as long as you keep this whole header intact.
% (c) Olivier Ledoit 1996

% de-mean returns
t=size(ret,1);
n=size(ret,2);
x=ret;
meanx=mean(x);
x=x-meanx(ones(t,1),:);
xmkt=mean(x')';

%Compute sample covariance matrix and prior
sample=cov([x xmkt]);
covmkt=sample(1:n,n+1);
varmkt=sample(n+1,n+1);
sample(:,n+1)=[];
sample(n+1,:)=[];
prior=covmkt*covmkt'./varmkt;
prior(logical(eye(n)))=diag(sample);

%Compute shrinkage parameters (as per Theorem 10)
d=1/n*norm(sample-prior,'fro')^2;
y=x.^2;
r2=1/n/t^2*sum(sum(y'*y))...
	-1/n/t*sum(sum(sample.^2));

% phi from section B.4 is divided into diagonal
% and off-diagonal terms, and the off-diagonal term
% is itself divided into smaller terms 
phidiag=1/n/t^2*sum(sum(y.^2))...
	-1/n/t*sum(diag(sample).^2);
z=x.*xmkt(:,ones(1,n));
v1=1/t^2*y'*z-1/t*covmkt(:,ones(1,n)).*sample;
phioff1=1/n*sum(sum(v1.*covmkt(:,ones(1,n))'))/varmkt...
	-1/n*sum(diag(v1).*covmkt)/varmkt;
v3=1/t^2*z'*z-1/t*varmkt*sample;
phioff3=1/n*sum(sum(v3.*(covmkt*covmkt')))/varmkt^2 ...
	-1/n*sum(diag(v3).*covmkt.^2)/varmkt^2;
phioff=2*phioff1-phioff3;
phi=phidiag+phioff;

%Compute the shrinkage covariance estimator

sigma=(r2-phi)/d*prior+(1-(r2-phi)/d)*sample;