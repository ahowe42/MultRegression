function [ICScore,lackoffit,penalty,betas,conststatus] = MultRegGASub_IC(IC,Y,x,bin,smth,constlim)
%{
  [ICScore,lack-of-fit,penalty,coefs,constraints status] = MultRegGASub_IC(information criteria, ...
  response matrix, predictor matrix, subset vector, robust covariance, constraint test)

  Compute the information criteria score for a subset multiple
  regression model with the assumption of gaussian errors.

  Where
  information criteria --- code for IC to compute, choices are:
     AIC, SBC, CAIC, ICOMP_IFIM, ICOMP_IFIM_PEU, ICOMP_IFIM_PEULN; if 'ALL' is 
     passed, a vector of scores for all criteria in this order is given
  response vector --- (n x 1) vector of response variales
  predictor matrix --- (n x q) matrix of q predictor variables
  subset vector --- (1 x q) row vector of logicals indicating subsets
  robust covariance --- covariance smoother code for CovSmooth
  constraint test --- optional struct of constraint test inputs for violationchecker function
  ICScore --- information criteria score for model
  lack-of-fit --- 2*log(L(theta))
  penalty -- penalty term
  coefs --- regression coeficients
  constraints status --- struct of constraint statuses from violationchecker function
  
  Copyright (C) 2009 J. Andrew Howe; see below
%}  

% if no parameters, pass out list of choices
if (nargin == 0)
    ICScore = {'NEG_ADJR2';'AIC';'SBC';'CAIC';'ICOMP_IFIM';'ICOMP_IFIM_PEU';'ICOMP_IFIM_PEULN'};
    return;
end

[n, p] = size(Y); [nX, q] = size(x); IC = upper(IC);

if (n ~= nX) || (length(bin) ~= q) || (nargin < 5) || (nargin >6)
    % dimensional mismatches or wrong # arguments input
    fprintf('MultRegGASub_IC: INVALID USAGE-Please read the following instructions!\n'), help MultRegGASub_IC, return
end

% compute the ordinary least squares estimates
q = sum(bin); numparmsest = q + 1;
X = x(:,logical(bin));                          % subset selection
Xt = X';                                        % only transpose this once
des = Xt*X;                                     % model design matrix
invdes = inv(des);                              % inverse model design matrix
betas = invdes*(Xt*Y);                          % regression coefficients
preds = X*betas;                                % predicted values
resids = Y - preds;                             % residuals
mse = sum(resids.^2)/n;                         % mean squared error
Hat = X*invdes*Xt;								% hat matrix
SSR_dof = [Y'*(Hat-ones(n)/n)*Y;q-1];			% SSR & dof
SSE_dof = [Y'*(eye(n)-Hat)*Y,n-q];				% SSE & dof

% compute the information criteria
lackoffit = n*log(2*pi) + n*log(mse) + n;
switch IC
	case 'NEG_ADJR2'							% adj r^2: 1 - (SSE/(SSE+SSR))*((n-1)/(n-preds))
		penalty = nan;
        ICScore = -1*(1 - (SSE_dof(1)/(SSE_dof(1)+SSR_dof(1)))*((n-1)/(n-q)));
    case 'AIC'                                  % 2*(# parms)
        penalty = 2*numparmsest;
        ICScore = lackoffit + penalty;
    case 'SBC'                                  % log(n)*(# parms)
        penalty = log(n)*numparmsest;
        ICScore = lackoffit + penalty;
    case 'CAIC'                                 % (log(n)+1)*(# parms)
        penalty = (log(n) + 1)*numparmsest;
        ICScore = lackoffit + penalty;
    case {'ICOMP_IFIM','ICOMP_IFIM_PEU','ICOMP_IFIM_PEULN'}
        % inverse fisher information matrix
        zmat = zeros(q,1);
        finv = [mse*invdes,zmat;zmat',2*(mse^2)/n];
        penalty = 2*EntComp(finv,1);            % 2*C1(Finv)
        if isequal(IC([(end-3):end]),'_PEU')% m + 2*C1(Finv)
            penalty = numparmsest + penalty;
        elseif isequal(IC([(end-1):end]),'LN')  % m + log(n)*C1(Finv)
            penalty = numparmsest + log(n)*penalty/2;
        end
        ICScore = lackoffit + penalty;
end
% JAH 20140330 do ALL of them
if strcmp(IC, 'ALL')
    % stuff for ICOMPS first
    zmat = zeros(q,1);
    penalty = 2*EntComp([mse*invdes,zmat;zmat',2*(mse^2)/n],1);
	ICScore = zeros(1,7);
	ICScore(1) = -1*(1 - (SSE_dof(1)/(SSE_dof(1)+SSR_dof(1)))*((n-1)/(n-q)));
	ICScore(2) = lackoffit + 2*numparmsest;
	ICScore(3) = lackoffit + log(n)*numparmsest;
	ICScore(4) = lackoffit + (log(n) + 1)*numparmsest;
	ICScore(5) = lackoffit + penalty;
	ICScore(6) = lackoffit + (numparmsest + penalty);
	ICScore(7) = lackoffit + (numparmsest + log(n)*penalty/2);
    penalty = nan;
end

% JAH 20131110 only continue with the constraints checking if limits passed in, or if not empty
if (nargin == 5); return; end;

% compute some stats for the constraints check function
% NOTE: calcs for DoF consider the number betas estimated INCLUDING the constant if it is in X;
% unsure if this should include sigma or not (it currently doesn't).  MATLAB's regress function
% seems consistent with this. My LSRegM function counts betas EXCLUDING constant + sigma; 
stats = zeros(1,3);
stats(1,1) = (SSR_dof(1)/SSR_dof(2))/(SSE_dof(1)/SSE_dof(2));			% F stat for model; SSR/dof / SSE/dof
stats(1,2) = 1-fcdf(stats(1,1),SSR_dof(2),SSE_dof(2));					% p-value of F stat
stats(1,3) = SSE_dof(1)/SSE_dof(2);										% estimate of error variance

% JAH 20131028 added code to check satisfaction of regression constraints
conststatus = violationchecker(X, betas, resids, stats, constlim);

%{
JAH 20090228, checked for octave 3.4.3 20120309

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
along with this program.  If not, see <http://www.gnu.org/licenses/>.
%}