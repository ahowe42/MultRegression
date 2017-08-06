function [ICScore,lackoffit,penalty,betas,conststatus] = MultVarRegGASub_IC(IC,Y,x,bin,smth,trueM)
%{
  ICScore,lack-of-fit,penalty,coefs,constraints status] = MultVarRegGASub_IC(information criteria, ...
  response matrix, predictor matrix, subset vector, robust covariance, true model)

  Compute the information criteria score for a subset multivariate
  regression model with the assumption of gaussian errors.

  Where
  information criteria --- code for IC to compute, choices are:
     AIC, SBC, GAIC, ICOMP_IFIM, ICOMP_IFIM_PEU, ICOMP_IFIM_PEULN,
     ICOMP_MISP, ICOMP_MISP_PEU, ICOMP_MISP_PEULN, KLDist
  response matrix --- (n x p) matrix of p response variables
  predictor matrix --- (n x q) matrix of q predictor variables
  subset vector --- (1 x q) row vector of logicals indicating subsets
  robust covariance --- covariance smoother code for CovSmooth
  true model --- structure having 3 elements (only used for KL Distance):
     .preds: (1 x q*) vector indicating predictors in true model
     .coefs: (q* x p) matrix of true model coefficients
     .sigma: (p x p) matrix of true error covariance matrix
  ICScore --- information criteria score for model
  lack-of-fit --- optional -2*log(L(theta))
  penalty -- optional penalty term
  coefs --- regression coeficients

  Copyright (C) 2007 J. Andrew Howe; see below
%}

% if no parameters, pass out list of choices
% if 1 parameter, include debug choices
if (nargin == 0)
    ICScore = {'AIC';'SBC';'GAIC';'ICOMP_IFIM';'ICOMP_IFIM_PEU';...
        'ICOMP_IFIM_PEULN';'ICOMP_MISP';'ICOMP_MISP_PEU';...
        'ICOMP_MISP_PEULN';'KLDist'};
    return;
%elseif(nargin == 1)
%    ICScore = {'AIC';'SBC';'GAIC';'ICOMP_IFIM';'ICOMP_IFIM_PEU';...
%       'ICOMP_IFIM_MISP';'ICOMP_IFIM_MISP_PEU';'KLDist';'ICOMP_IFIMC';...
%        'ICOMP_IFIMC_PEU';'ICOMP_IFIMC_MISP';'ICOMP_IFIMC_MISP_PEU';...
%        'ICOMP_IFIMO_MISP';'ICOMP_IFIMO_MISP_PEU'};
%    return;    
end

[n, p] = size(Y); [nX, q] = size(x); IC = upper(IC);

if (n ~= nX) || (length(bin) ~= q) || (nargin ~= 6)
    % dimensional mismatches or wrong # arguments input
    fprintf('MultVarRegGASub_IC: INVALID USAGE-Please read the following instructions!\n'), help MultVarRegGASub_IC, return
end

% compute the ordinary least squares estimates
X = x(:,logical(bin));                          % subset selection
Xt = X';                                        % only transpose this once
des = Xt*X;                                     % model design matrix
invdes = inv(des);                              % inverse model design matrix
betas = invdes*(Xt*Y);                          % regression coefficients
preds = X*betas;                                % predicted values
resids = Y - preds;                             % residuals
S = (resids'*resids)/n;                         % estimate the covariance matrix
S = CovSmooth(S,smth,0,1,n);                    % smooth it
% Can't do below because ICOMP_IFIM and ICOMP_IFIM_MISP don't work with R
% compute correlation form of sigma, to make all IC scale invariant
%d = inv(sqrt(diag(diag(sigmahat_sm)))); R = d*S*d;
dS = det(S);

% compute the information criteria
lackoffit = n*p*log(2*pi) + n*log(dS) + n*p;
q = sum(bin); numparmsest = p*q + p*(p+1)/2;
switch IC
    case 'AIC'
        penalty = 2*numparmsest;
    case 'SBC'
        penalty = log(n)*numparmsest;
    case 'GAIC'
        Ip = eye(p);
        Z = resids*(inv(real(sqrtm(S))));           % standardized residuals
        ZtZ = Z'*Z;
        gamma2 = ZtZ(:)*(ZtZ(:)');
        gamma2st = gamma2 + (n^2)*Ip(:)*(Ip(:)');   % kurtosis matrix
        gamma2st = gamma2st.*sign(gamma2st);        % force to be positive
        penalty = 2*p*q + trace(gamma2st)/n;
    case {'ICOMP_IFIM','ICOMP_IFIM_PEU','ICOMP_IFIM_PEULN'}
        trS = trace(S);
        part1 = trS*trace(invdes) + (1/(2*n))*(trace(S^2) + trS^2 + 2*sum(diag(S).^2));
        part1 = numparmsest*log(part1/numparmsest);
        part2 = -(p + q)*log(dS) + p*log(det(des)) + 0.5*p*(p + 1)*log(n) - p*log(2);
        penalty = part1 + part2;
%        % closed version
%        [Dp, DpPlus] = DupMatrix(p);
%        Q1 = kron(S,invdes);
%        Q4 = (2/n)*(DpPlus*kron(S,S)*DpPlus');
%        Q2 = zeros(size(Q1,1),size(Q4,2));
%        IFIM = [Q1,Q2;Q2',Q4];
%        detIFIM = det(IFIM); traIFIM= trace(IFIM); rnkIFIM = rank(IFIM);
%        penalty = rnkIFIM*log(traIFIM/rnkIFIM) - log(detIFIM);
        % do heavier penalty if needed
        if isequal(IC([(end-3):end]),'_PEU')
            penalty = numparmsest + penalty;
        elseif isequal(IC([(end-5):end]),'_PEULN')
            penalty = numparmsest + log(n)*penalty/2;
        end
    case {'ICOMP_MISP','ICOMP_MISP_PEU','ICOMP_MISP_PEULN'}
        % compute matrix generalizations of skewness and kurtosis
        Ip = eye(p);
        Z = resids*(inv(real(sqrtm(S))));           % standardized residuals
        ZtZ = Z'*Z;
        gamma1 = ZtZ - n*Ip;
        gamma1 = Z(:)*(gamma1(:)');                 % skewness matrix
        gamma2 = ZtZ(:)*(ZtZ(:)');
        gamma2st = gamma2 + (n^2)*Ip(:)*(Ip(:)');   % kurtosis matrix
        gamma2st = gamma2st.*sign(gamma2st);        % force to be positive
        % other useful stuff
        [Dp, DpPlus] = DupMatrix(p); Dpt = Dp';
        sqS = real(sqrtm(S)); KronsqS = kron(sqS,sqS);
%        % closed version
%        iD = inv(sqS); % just reusing variable
%        Del = Dpt*kron(iD,iD)*Dp; iD = inv(Del);
%        Q1 = kron(S,invdes);
%        Q2 = (1/n)*kron(sqS,invdes*Xt)*gamma1*Dp*iD;
%        Q3 = (1/n)*iD*Dpt*gamma1'*kron(sqS,X*invdes);
%        Q4 = (1/n^2)*iD*Dpt*gamma2st*Dp*iD;
%        CovM = [Q1,Q2;Q3,Q4];
%        CovM = CovSmooth(CovM,smth,0,1,n);      % Cov(misspec) is rank deficient, so fix
%        detCovM = det(CovM); traCovM = trace(CovM); rnkCovM = rank(CovM);
        % sandwich covariance matrix
        isqrtS = inv(sqS); DpPlust = DpPlus';
        Del = Dpt*kron(isqrtS,isqrtS)*Dp;
        % inner product form
        Q1 = kron(S,invdes);
        Q4 = (2/n)*(DpPlus*kron(S,S)*DpPlust);
        Q2 = zeros(size(Q1,1),size(Q4,2));
        IPF = [Q1,Q2;Q2',Q4];            
        % outer product form
        Q1R = kron(inv(S),des);
        Q2R = 0.5*kron(isqrtS,Xt)*gamma1*DpPlust*Del;
        Q3R = 0.5*Del*DpPlus*gamma1'*kron(isqrtS,X);
        Q4R = 0.25*Del*DpPlus*gamma2st*DpPlust*Del;
        OPF = [Q1R,Q2R;Q3R,Q4R];
        CovM = IPF*OPF*IPF;                     % Finv*R*Finv
        CovM = CovSmooth(CovM,smth,0,1,n);      % Cov(misspec) is rank deficient, so fix
        detCovM = det(CovM); traCovM = trace(CovM); rnkCovM = rank(CovM);
%        % opened up version
%        Big = trace(DpPlus*KronsqS*gamma2st*KronsqS*(DpPlus'));
%        % compute trace(cov(theta_hat)_misp)
%        traCovM = trace(S)*trace(invdes) + (n^-2)*Big;
%        Big = det(Dpt*(gamma2st - gamma1'*kron(Ip,X*invdes*Xt)*gamma1)*Dp);
%        % compute det(cov(theta_hat)_misp) - ->0 since hat matrix is singular
%        detCovM = (2^(-p*(p - 1)))*(n^(-p*(p + 1)))*(dS^(p + q))*(det(des)^(-p))*Big;
%        rnkCovM = numparmsest;
        if (traCovM == 0) || (detCovM == 0)
            disp('debug')
            [st,i] = dbstack; eval(['dbstop in ''MultVarRegGASub_IC.m'' at ',num2str(st(i).line+1)]);
        end        
        penalty = rnkCovM*log(traCovM/rnkCovM) - real(log(detCovM));
        if isequal(IC([(end-3):end]),'_PEU')
            penalty = numparmsest + penalty;
        elseif isequal(IC([(end-5):end]),'_PEULN')
            penalty = numparmsest + log(n)*penalty/2;            
        end
    case 'KLDIST'
        invS = inv(S);
        moddiff = (x(:,trueM.preds)*trueM.coefs - preds);
        penalty = log(dS/det(trueM.sigma)) + trace(invS*trueM.sigma) + ...
            trace(invS*(moddiff'*moddiff)) - p;
        % subtract lackoffit, since it will be added on below, making this IC just the KL distance
        penalty = penalty - lackoffit;
end
ICScore = lackoffit + penalty;

%{
JAH 20070304, checked for octave 3.4.3 20120211

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
