function viola = violationchecker(x, beta, residuals, stats, confviola)
%{ violation_flags = violationchecker(X, betas, residuals, regress_stats, test_options)
  Check the results from application of a multiple regression model against
  several modeling assumptions and requirements as follows:
  p: valid model: F-test for regression model must be < specified confidence
  m: multicollinearity: VIF for the model must be < specified bound
  e: Gaussian errors: use the Kolmogorov-Smirnov test for normality on the residuals
  h: Homoskedastic errors: use the specified test for heteroskedasticity on the residuals

  Where
  X --- (n x p) matrix of independent variables (includinig the constant term if used)
  betas --- (p x 1) vector of regression coefficients
  residuals --- (n x 1) vector of residual error terms
  regress_stats --- (1 x 4) vector of regression statistics:
     R^square, F-statistic, p-value for F-statistic, estimated error variance
  test_options --- struct holding three parameters for the various tests:
     confviola.alpha: confidence level for hypothesis tests
     confviola.VIFbound: upper bound for VIF to declare the presence of multicollinearity
     confviola.HeteroTest: test for homoskedasticity
     'BPK' = Breush-Pagan, Koenker modification (Breush-Pagan 1979; Koenker 1981)           
     'W' = White (White 1980b)
     'Ws' = White, Wooldridge special case (White 1980b; Wooldridge 2006, p.286)
     and also the weights for each of the tests in order:
     modelp, multicollin, errornormal, heteros
  violation_flags --- struct holding results from the various tests:
     modelp: boolean flag indicating violation of test p
     multicollin: boolean flag indicating violation of test m
     errornormal: boolean flag indicating violation of test e
     heteros: boolean flag indicating violation of test h
     STR: 4-character string indicating which tests were violated
     TOTAL: integer number of violations
     SCORE: weighted score of violations (using input weights)

  Dr. Oguz Akbilgic 20131028
  Dr. John Andrew Howe

  Copyright (C) 2013 Oguz Akbilgic & J. Andrew Howe
%}

viola=[];
% test p:
viola.modelp = (stats(1,2) > confviola.alpha);

% test m:
% JAH 20150630 first check if the first column is the intercept
% if so, and it's the only column, just say the m test is 0
% otherwise strip it off
if isequal(x(:,1), ones(size(x,1),1))
    if size(x,2) == 1
        viola.multicollin = 0;
    else
        viola.multicollin = (sum(diag(inv(corr(x(:,[2:end])))) > confviola.VIFbound)>0);
    end
else
	viola.multicollin = (sum(diag(inv(corr(x))) > confviola.VIFbound)>0);
end

% test h:
% JAH 20150703 first check if the first column is the intercept
% if so, and it's the only column, just say the m test is 0
% otherwise strip it off
if isequal(x(:,1), ones(size(x,1),1))
    if size(x,2) == 1
        pVal = 0; % JAH 20150770 not sure what to do here, since the design matrix is uninvertible
    else
        pVal = TestHet(residuals, x(:,[2:end]), confviola.HeteroTest, x*beta);
    end
else
    pVal = TestHet(residuals, x, confviola.HeteroTest, x*beta);
end
viola.heteros = (pVal < confviola.alpha);

% test e:
% JAH 20150703 this test was not identifying nonnormality with
% PE beta = 0.5 tails, so I changed to Mardia
%viola.errornormal = (kstest(residuals/std(residuals)) == 1);
[K,S] = MVNormalityTest(residuals,1);
viola.errornormal = or(K(end) <  confviola.alpha, S(end) <  confviola.alpha);

% JAH 20140330 2nd to last item is a word coding which are violated
tmp = ['****';'pmeh'];
viola.STR = [tmp(viola.modelp+1,1),tmp(viola.multicollin+1,2),tmp(viola.errornormal+1,3),tmp(viola.heteros+1,4)];
% JAH 20131110 next to last item is count of all previous
viola.TOTAL = viola.modelp + viola.multicollin + viola.errornormal + viola.heteros;
% JAH 20150630 last item is a weighted score from 0 to 1 of violations
viola.SCORE = viola.modelp*confviola.p_weight + viola.multicollin*confviola.m_weight + viola.errornormal*confviola.e_weight + viola.heteros*confviola.h_weight;

%{
Copyright (C) 2013 J. Andrew Howe
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
