# MultRegression
## MATLAB / Octave toolbox for Multiple &amp; Multivariate Regression with several advanced modeling techniques.

This is a "toolbox" of code for performing both Multiple and Multivariate regression with several advanced modeling techniques.  I have worked on this code, off and on, since 2007.  It was adopted to work with Octave 3.4.3 in 2012.  Specifically, it includes:
- Genetic algorithm to identify optimal subset regression model
- Computation of various information-theoretic model selection criteria
- Creation of complex regression simulation protocols
- Regularization of ill-conditioned covariance matrices with a variety of smoothed estimators
- Optimization in the GA for multiple regression may be further constrained to minimize violation of certain regression assumptions

Parts of this code were used in these publications:
- Bozdogan, H., Howe, J. A., April 2012. [Misspecified Multivariate Regression Models Using the Genetic Algorithm and Information  Complexity as the Fitness Function](http://ejpam.com/index.php/ejpam/article/view/1597). European Journal of Pure and Applied Mathematics 5 (2), 211–249.
- Bozdogan, H., Howe, J. A., Katragadda, S., Liberati, C., 2013. [Misspecification Resistant Model Selection Using Information Complexity with Applications](http://link.springer.com/chapter/10.1007%2F978-3-642-28894-4_20). Springer-Verlag, Berlin, pp. 165–172.

There are four driver scripts which control everything.

## MultRegAllSub_DRV
Perform all subsets analysis for multiple regression (only the regressors are subset), and compute various information criteria scores.  If you tell this to use simulated data, this can run multiple simulations, and determine how often each criterion selected the "true" model.

## MultVarRegAllSub_DRV
Same as above, but for multivariate regression.

## MultRegGASub_DRV
Perform optimal subsets analysis for multiple regression (only the regressors are subset) driven by the genetic algorithm, and compute various information criteria scores.  If you tell it to use simulated data, this can run multiple simulations, and determine how often each criterion selected the "true" model.

## MultVarRegGASub_DRV
Same as above, but for multivariate regression.
