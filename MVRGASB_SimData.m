%{
usage: MVRGASB_SimData

If data_type is 1, data_file is the program used to simulate the data
it can only accept the number observations as an argument, and must
return an [n x py] matrix of responses, [n x px] matrix of predictors,
and a structure identifying the true predictors, coefficients, and
errors (stored in variable trueM).  The responses are stored in Y and the
predictors are in X.  HasConstant indicates whether or not the first
column of X is a constant (all 1s).

Copyright (C) 2007 J. Andrew Howe; see below
%}

% simulate the data
func = data_file([1:strfind(data_file,'.')-1]);
[Y,X,trueM] = feval(func,n,extra_vars,mvflag); % JAH added mvflag param 20120309
% get the number responses and predictors
py = size(Y,2); px = size(X,2);
% do the predictors have a constant in the first column?
HasConstant = (sum(X(:,1) == ones(n,1)) == n);
% JAH 20140330 if no constant is already included, add it, and shift the true vector by 1
if (KeepConstant == 1) && (HasConstant == 0)
	px = px + 1;
	X = [ones(n,1),X];
	HasConstant = 1;
	trueM.preds = trueM.preds + 1;
end

clear func

%{
JAH 20070404, checked for octave 3.4.3 20120219

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