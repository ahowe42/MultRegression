%{
usage: MVRGASB_LoadData

This scripts loads the data (nonsimulated) from data_file if data_type is
2.  Variables created are: n = number observations, py = number
responses, px = number predictors, HasConstant = boolean, Y = responses,
X = covariates.  data_file must be a delimited text file with the first
row indicating py and px.  After this, the data must be arranged with all
columns of Y and then all columns of X.  The first column of X must be
the intercept, for HasConstant to be accurate, if it is included.
If KeepConstant is TRUE and HasConstant is FALSE, a constant column will
be added.

Copyright (C) 2007 J. Andrew Howe; see below
%}

if data_type == 1
    % simulated data = do nothing, the simulator will take care of it
elseif data_type == 2;
    % real data must be organized in tab-delimited file such that:
    % row 1 = num_response   num_predictor
    % the rest = num_response y columns then num_predictor x columns
    datinput = dlmread([mydir,filesep,'data',filesep,data_file]);    % load the datafile
    % first, get number of responses and predictors, and num observations
    py = datinput(1,1); px = datinput(1,2); n = size(datinput,1) - 1;
    % JAH 20140405 ensure py = 1 if only multiple regression
    if (py > 1) && (mvflag == 0)
		error(sprintf('Input data has multivariate Y (py = %d), but this is multiple regression!',py))
    end
    % then get the data
    Y = datinput([2:end],[1:py]); X = datinput([2:end],[(py+1):(py + px)]);
    % do the predictors have a constant in the first column?
    HasConstant = (sum(X(:,1) == ones(n,1)) == n);
    if (KeepConstant == 1) && (HasConstant == 0)
		px = px + 1;
		X = [ones(n,1),X];
		HasConstant = 1;
    end
    clear datinput
end

%{
JAH 20070404, adapted for octave 3.4.3 20120213

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