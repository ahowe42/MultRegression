%{
usage: MVRGASB_MRSummary

After performing many simulations/replications of the GA with information
criteria for subset regression, this script will provide
summary results across all runs.  The workspace .mat file from the
multirun must already be loaded in the workspace.

Copyright (C) 2007 J. Andrew Howe; see below
%}

clc
if (exist([mydir,mult_run([1:(end-4)]),'_SMRY.out'],'file') == 2)
    delete([mydir,mult_run([1:(end-4)]),'_SMRY.out']);   % delete existing summary file
end
diary([mydir,mult_run([1:(end-4)]),'_SMRY.out'])

%[st,i] = dbstack; eval(['dbstop in ',st(i).name,' at ',num2str(st(i).line+1),';']);
showtopsubs = 10;
% get the unique sorted solutions from all replications
[unique_solns,uniind] = unique(MRchroms(:,[2:end]),'rows');
numuni = size(unique_solns,1);
MR_BEST = [zeros(numuni,6),unique_solns];

% now get the min IC score (1), weights (2, later) replication index (3), frequency (4), time (5), and # constraints violated (6)
for rcnt = 1:numuni
    % find all solutions that match in all loci
    mtchs = (sum(repmat(unique_solns(rcnt,:),num_reps,1) == MRchroms(:,[2:end]),2) == px);
    MR_BEST(rcnt,4) = sum(mtchs);                           % frequency
    mtchs = MRscores(mtchs,:);
    [MR_BEST(rcnt,1),MR_BEST(rcnt,3)] = min(mtchs(:,2));    % min IC score
    MR_BEST(rcnt,5) = mtchs(MR_BEST(rcnt,3),4);             % computation time
    MR_BEST(rcnt,6) = mtchs(MR_BEST(rcnt,3),3);             % constraint violated scores JAH 20140328
    MR_BEST(rcnt,3) = mtchs(MR_BEST(rcnt,3),1);             % rep index
end                 % replications loop
[MR_BEST,srtind] = sortrows(MR_BEST,1);

% JAH 20140330 best constraint violation words
if mvflag == 0
	MR_BESTWords = MRwords(MR_BEST(:,3),:);
end

% compute IC weights for all generation results
wghts = exp(-0.5*(MR_BEST(:,1) - MR_BEST(1,1)));
MR_BEST(:,2) = wghts/sum(wghts);

top_best = min(showtopsubs,numuni); chromloc = [7:(px + 6)];
% prepare row headers
vars = [1:px]; rwhds = cell(top_best,1);
for rcnt = 1:top_best
    vrs = find(vars.*MR_BEST(rcnt,chromloc)) - HasConstant;
    vrs = sprintf('%d,',vrs);
    rwhds{rcnt,1} = ['{',vrs([1:(end-1)]),'}'];
end                 % replications loop
% compute the relative frequency of true model hits
select_true = 0; includ_true = 0;
if (data_type == 1)
    for rcnt = 1:numuni
        mod = find(vars.*MR_BEST(rcnt,chromloc));
        % equals true model
        if isequal(trueM.preds,mod)
            select_true = select_true + MR_BEST(rcnt,4);
            includ_true = includ_true + MR_BEST(rcnt,4);
            continue;
        end
        % includes true model
        if (sum(ISIN(trueM.preds,mod,1)) == length(trueM.preds))
            includ_true = includ_true + MR_BEST(rcnt,4);
        end
    end             % replications loop
end

% display the best chromosomes and scores
if mvflag == 1
	tab = table2str({objec_func,'Weights','Replication','Frequency','Time(s)'},MR_BEST([1:top_best],[1,2,3,4,5]),{'%0.2f','%0.3f','%d','%d','%0.1f'},0,rwhds);
else
	tab = table2str({objec_func,'Const. Viol.','Weights','Replication','Frequency','Time(s)'},MR_BEST([1:top_best],[1,6,2,3,4,5]),...
            {'%0.2f','%0.2f','%0.3f','%d','%d','%0.1f'},0,rwhds);
	tab = [tab,[repmat(['-';' ';'-'],1,size(MR_BESTWords,2));MR_BESTWords([1:top_best],:);repmat('-',1,size(MR_BESTWords,2))]];
end
MVRGASB_DispParams
lin = repmat('=',1,60); disp(' '), disp(lin)
if (mvflag == 1) % JAH added 20120309 
	disp('SUMMARY: Multivariate Regression Subsetting with the Genetic Algorithm')
else 
	disp('SUMMARY: Multiple Regression Subsetting with the Genetic Algorithm')
end
disp([mult_run])
if (data_type == 1)
    disp(sprintf('Results from %0.0f Simulations',num_reps))
else
    disp(sprintf('Results from %0.0f Replications',num_reps))
end
disp(tab)
disp(sprintf('Average Solution Order: %0.0f',mean(sum(MRchroms(:,[2:end]),2))))
% if simulated, show how well it hit the true model
if (data_type == 1)
    disp(['True Model: ',strrep(strrep(strrep(mat2str(trueM.preds - HasConstant),' ',','),'[','{'),']','}')])
    disp(sprintf([objec_func,' selected the true model in %0.2f%% of the simulations.'],100*select_true/num_reps))
    disp(sprintf([objec_func,' included the true model in %0.2f%% of the simulations.'],100*includ_true/num_reps))
    if isequal(trueM.preds,find(vars.*MR_BEST(1,chromloc)))
        tmp = 'selected';
    else
        tmp = 'didn''t select';
    end
    disp([objec_func,' ',tmp,' the true model as the best.'])
	if (mvflag == 1) % JAH 20120309
    	trueIC = MultVarRegGASub_IC(objec_func,Y,X(:,trueM.preds),ones(1,length(trueM.preds)),regul_func,trueM);
    	obj = sprintf('\n%s for true model (last run): %0.2f',objec_func,trueIC);
	else
    	[trueIC,lof,pen,bet,trueconst] = MultRegGASub_IC(objec_func,Y,X(:,trueM.preds),ones(1,length(trueM.preds)),regul_func,constlim);
		obj = sprintf('\n%s for true model (last run): %0.2f, Constraints Violated: %s (score = %0.2f)',objec_func,trueIC,trueconst.STR,trueconst.SCORE);
	end
    disp(obj)
    if (mvflag == 1) % JAH 20120309
        satIC = MultVarRegGASub_IC(objec_func,Y,X,ones(1,px),regul_func,trueM);
        obj = sprintf('%s for saturated model (last run): %0.2f',objec_func,satIC);
    else
        [satIC,lof,pen,bet,satconst] = MultRegGASub_IC(objec_func,Y,X,ones(1,px),regul_func,constlim);
        obj = sprintf('%s for saturated model (last run): %0.2f, Constraints Violated: %s (score = %0.2f)',objec_func,satIC,satconst.STR,satconst.SCORE);
    end
    disp(obj)
end
disp(lin)
disp(sprintf('Modeling Completed in \n\t%1.4f Seconds\n\t%1.4f Minutes\n\t%1.4f Hours',tottim./[1,60,3600]));
diary off

%{
JAH 20070405, checked for octave 3.4.3 20120219

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