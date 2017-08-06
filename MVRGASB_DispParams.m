%{
usage: MVRGASB_DispParams

The purpose of this script is to print a "table" of the currently
selected parameters for the robust misspecified regression with GA.

Copyright (C) 2007 J. Andrew Howe; see below
%}

disp(repmat('#',1,50))
disp(sprintf('Population Size: %0.0f\nGeneration Size: %0.0f\nMinimum # of Generations: %0.0f',popul_size,num_generns,nochange_terminate))
disp(sprintf('Crossover Rate: %0.2f\nMutation Rate: %0.2f',prob_xover,prob_mutate))
if elitism; disp('Elitism is: ON'); else; disp('Elitism is: OFF'); end;
if pltflg; disp('Update Plot-on-the-Fly is: ON'); else; disp('Update Plot-on-the-Fly is: OFF'); end;
if plt3d; disp('Create 3d Score Plot is: ON'); else; disp('Create 3d Score Plot is: OFF'); end;
if mvflag == 1; disp(['Covariance Smoothing Code: ',regul_func]); end;
disp(['Objective Function: ',objec_func])
switch data_type
    case 1     % simulated data
		disp(sprintf('Simulated %d observations with %d noise variables from: %s',n,extra_vars,data_file))
    case 2      % real dta
        disp(['Real Data: ',data_file])        
end
if constraint_check == 1
	disp('Constraint Enforcement is: ON');
elseif mvflag == 0
	disp('Constraint Enforcement is: OFF');
end
disp(sprintf('\tParameters: Model alpha: %0.2f, VIF Limit: %0.2f, Heteroskedasticity test: %s',constlim.alpha,constlim.VIFbound,constlim.HeteroTest))
disp(sprintf('\tWeights: Model alpha = %0.2f, Multicollinearity = %0.2f, Nonnormality = %0.2f, Heterskedasticity = %0.2f',constlim.p_weight, constlim.m_weight,constlim.e_weight,constlim.h_weight))

disp(sprintf('Random State: %10.0f',rnd_stat))
if (KeepConstant == 1)
	disp('Force Constant is: ON');
else
	disp('Force Constant is: OFF');
end
disp(repmat('#',1,50))

%{
JAH 20070403, checked for octave 3.4.3 20120219

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