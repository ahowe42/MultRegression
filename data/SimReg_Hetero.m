function [y,x,trueM] = SimReg_Hetero(n,numred,mv)
%{ [responses, predictors, true model] = SimReg_Hetero(number observations,number redundant variables,multivariate flag)
  Simulate regression data such that selection of the correct model should
  not violate the assumption that error terms have constant variance.  However,
  if an incorrect model is fit, the error terms be heteroskedastic.  Also
  present in a "wrong" model is multicolinearity. The true model is [x1,x2,x3,x4].

  Where
  number observations --- (1 x 1) number observations required
  number redundant vars --- (1 x 1) number 'junk' variables added
  multivariate flag --- 1 (default) = multivariate Y, 0 = univariate Y
  responses --- (n x p) matrix of responses
  predictors --- (n x q) matrix of predictors
  true model --- structure having 4 elements:
     .preds: (1 x q*) vector indicating predictors in true model
     .coefs: (q* x p) matrix of true model coefficients
     .sigma: (p x p) matrix of true error covariance matrix
     .gen: string name of the generator

  [y,x,trueM] = SimReg_Hetero(100,5,1);

  See Also SimReg_B075, SimReg_B15, SimReg_B2, SimReg_Mixt, SimReg_Skew
  
  Copyright (C) 2015 J. Andrew Howe; see below
%}

if not(isscalar(n)) || (nargin < 2) || (nargin > 3) || not(isscalar(numred))
    % inputs not scalar, and not 2/3 of them
    fprintf('SimReg_Hetero: INVALID USAGE-Please read the following instructions!\n'), help SimReg_Hetero, return
end

if (nargin == 2); mv = 1; end; % default is multivariate regression

% random noise error terms
if mv == 1
	mu = [0;0];                                 % errors mean
	sigma = [1,0.5;0.5,1];                      % errors variance matrix
else
	mu = 0; sigma = 1;
end

% make x1
l = linspace(1,n,n)';
x1 = randn(n,1).*exp(ceil(l/(n/3)));
% make the rest of X
x2 = sqrt(l)+randn(n,1);                        % linear term for y
x3 = 0.6*x2 + randn(n,1)*sqrt(3);               % something with mild correlation
x45 = rand(n,2)*chol([1,0.90;0.90,1]);
x67 = rand(n,2)*chol([1,0.90;0.90,1]) + x45.*repmat([0.9,0.8],n,1);
% generate the junk
jnk = rand(n,numred).*repmat([8:(numred + 7)], n,1);

% generate the response variables
B = [10,-10;3,4;4,5;2,3];               % real coefficient matrix
if (mv == 0); B = B(:,1); end; 
E = mvnrnd(mu,sigma,n);                 % Gaussian errors
x = [ones(n,1),x1,x2,x3,x45,x67,jnk];	% add constant & junk
y = x(:,[1:4])*B + E;                   % construct the responses

% true model
trueM.preds = [1,2,3,4];
trueM.coefs = B;
trueM.sigma = sigma;
trueM.gen = 'SimReg_Hetero';

%{ JAH 20150703

Copyright (C) 2015 J. Andrew Howe
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
