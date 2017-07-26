function [y,x,trueM] = SimReg_Mixt(n,numred,mv)
%{ [responses, predictors, true model] = SimReg_Mixt(number observations, number junk variables, multivariate flag)
  Simulate regression data where the first five (after constant)
  predictors are correlated.  The responses are created from from the first
  four (including constant) predictors with mixtures of normals error
  terms.  The true model is [x1,x2,x3,x4].

  Where
  number observations --- (1 x 1) number observations required
  number redundant vars --- (1 x 1) number 'junk' variables added
  multivariate flag --- 1 (default) = multivariate Y, 0 = univariate Y
  responses --- (n x p) matrix of responses
  predictors --- (n x q) matrix of predictors
  true model --- structure having 3 elements:
     .preds: (1 x q*) vector indicating predictors in true model
     .coefs: (q* x p) matrix of true model coefficients
     .sigma: (p x p) matrix of true error covariance matrix
     .gen: string name of the generator

  [y,x,trueM] = SimReg_Mixt(100,15,1);

  See Also SimReg_B075, SimReg_B1, SimReg_B15, SimReg_B2, SimReg_Skew
  
  Copyright (C) 2007 J. Andrew Howe; see below
%}

if not(isscalar(n)) || (nargin < 2) || (nargin > 3) || not(isscalar(numred))
    % inputs not scalar, and not 2/3 of them
    fprintf('SimReg_Mixt: INVALID USAGE-Please read the following instructions!\n'), help SimReg_Mixt, return
end

if (nargin == 2); mv = 1; end; % default is multivariate regression

% params
cor = 0.3;                                      % multicollinearity
error_terms_variance = 1;                       % variance for x errors
% noise mixture parameters
Mu1 = [-2,0]; Sigma1 = [2,0;0,5]; pi1 = 0.33;
Mu2 = [8.3,8.1]; Sigma2 = [4.44,2.17;2.17,6.18]; pi2 = 0.33;
Mu3 = [5,8]; Sigma3 = [3.1,2.0;2.0,4.45]; pi3 = 0.33;
% normalize pis
sumpis = sum([pi1,pi2,pi3]);
pi1 = pi1/sumpis; pi2 = pi2/sumpis; pi3 = pi3/sumpis;

% generate the predictor variables
alpha = sqrt(1 - cor^2);                        % correlation factor
error_terms = randn(n,5)*sqrt(error_terms_variance);
x(:,1) = 10 + error_terms(:,1);
x(:,2) = 10 + cor*error_terms(:,1) + alpha*error_terms(:,2);
x(:,3) = 10 + cor*error_terms(:,1) + 0.5604*alpha*error_terms(:,2) + 0.8282*alpha*error_terms(:,3);
x(:,4) = -8 + x(:,1) + 0.5*x(:,2) + cor*x(:,3) + 0.5*error_terms(:,4);
x(:,5) = -5 + 0.5*x(:,1) + x(:,2) + 0.5*error_terms(:,5);

% generate the junk
jnk = rand(n,numred).*repmat([6:(numred + 5)], n,1);

% generate the noise from a mixture of three normals
if mv == 1
    X1 = mvnrnd(Mu1,Sigma1,ceil(n*pi1));
    X2 = mvnrnd(Mu2,Sigma2,ceil(n*pi2));
    X3 = mvnrnd(Mu3,Sigma3,ceil(n*pi3));    
    % simultaneously take only n obs and randomly permute the ordering
    E = [X1;X2;X3]; E = E(randperm(n),:);
else
    X1 = mvnrnd(Mu1(1),Sigma1(1,1),ceil(n*pi1));
    X2 = mvnrnd(Mu2(1),Sigma2(1,1),ceil(n*pi2));
    X3 = mvnrnd(Mu3(1),Sigma3(1,1),ceil(n*pi3));
    % simultaneously take only n obs and randomly permute the ordering
    E = [X1;X2;X3]; E = E(randperm(n));    
end

% generate the response variables
B = [-8,-5;1,0.5;0.5,0;0.3,0.3];                % real coefficient matrix
if (mv == 0); B = B(:,1); end;
x = [ones(n,1),x,jnk];                          % add constant & junk
y = x(:,[1:4])*B + E;                           % construct the responses

% true model
trueM.preds = [1,2,3,4];
trueM.coefs = B;
trueM.sigma = nan;
trueM.gen = 'SimReg_Mixt';

%{ JAH 20070502, 20120218

Copyright (C) 2007 J. Andrew Howe
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.%}
