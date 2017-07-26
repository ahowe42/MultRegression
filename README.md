# MultRegression
## MATLAB / Octave toolbox for Multiple &amp; Multivariate Regression with several advanced modeling techniques.

This is a "toolbox" of code for performing both Multiple and Multivariate regression with several advanced modeling techniques.  I have worked on this code, off and on, since 2007.  It was adopted to work with Octave 3.4.3 in 2012.  Specifically, it includes:
- Genetic algorithm to identify optimal subset regression model
- Computation of various information-theoretic model selection criteria
- Creation of complex regression simulation protocols
- Regularization of ill-conditioned covariance matrices with a variety of smoothed estimators
- Optimization in the GA for multiple regression may be further constrained to minimize violation of certain regression assumptions

There are four driver scripts which control everything.

## MultRegAllSub_DRV
Perform all subsets analysis for multiple regression (only the regressors are subset), and compute various information criteria scores.  If you tell this to use simulated data, this can run multiple simulations, and determine how often each criterion selected the "true" model.

## MultVarRegAllSub_DRV
Same as above, but for multivariate regression.

## MultRegGASub_DRV
Perform optimal subsets analysis for multiple regression (only the regressors are subset) driven by the genetic algorithm, and compute various information criteria scores.  If you tell it to use simulated data, this can run multiple simulations, and determine how often each criterion selected the "true" model.

## MultVarRegGASub_DRV
Same as above, but for multivariate regression.
