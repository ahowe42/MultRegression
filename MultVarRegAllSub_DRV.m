%{
usage: MultVarRegAllSub_DRV

Perform all subsets analysis for multivariate regression (only the
regressors are subset), and compute various information criteria scores.
If you tell this to use simulated data, this can run multiple simulations,
and determine how often each criterion selected the "true" model.

Copyright (C) 2007 J. Andrew Howe; see below
%}

clear
drvstt = clock;
stt = sprintf('%4.0f%02.0f%02.0f_%02.0f%02.0f%02.0f',drvstt);
mvflag = 1;	% JAH 20120309 this is for multivariate regression

% find and store where I am located
mydir = dbstack; mydir = mydir(end);
mydir = which(mydir.file);
tmp = strfind(mydir,filesep);
mydir = mydir([1:(tmp(end)-1)]);

% initialize variables
data_type = 1+isequal('Real',questdlg('Use Simulated or Real Data?',...
    'MultVarRegAllSub_DRV','Simulated','Real','Simulated'));    % 1=sim, 2 = real

if data_type == 1   % simulated
    [data_file,data_path] = uigetfile('*.m', 'Load Simulated Data File',[mydir,filesep,'data',filesep]);
    answer = inputdlg({'Number Simulations','Number Noise Variables',...
        'Number Observations','Covariance Smoother','Random State','Keep Constant'},'Experiment Parameters',...
        1,{'100','4','500','MLE/EB','0','1'},'off');
    % get the parameters
    mciter = str2num(char(answer(1)));
    extra_vars = str2num(char(answer(2)));
    n = str2num(char(answer(3)));
    regul_func = char(answer(4));
    % JAH 20140327 added rand state
    if char(answer(5)) == '0';
		rnd_stat = sum(drvstt*1000000);
	else
		rnd_stat = str2num(char(answer(5)));
    end
	rand('state',rnd_stat); randn('state',rnd_stat); %randg('state',rnd_stat);
    KeepConstant = str2double(answer{6});		% JAH added keep constant flag 20140330 
else                % real
    [data_file,data_path] = uigetfile('*.*', 'Load Real Data File',[mydir,filesep,'data',filesep]);
    answer = inputdlg({'Covariance Smoother','Keep Constant'},'Experiment Parameters',1,{'MLE/EB','1'},'off');
    regul_func = char(answer{1});
    KeepConstant = str2double(answer{2});		% JAH added keep constant flag 20140330 
    mciter = 1;
end
dispgran = ceil(0.2*mciter);

% fix path
addpath(data_path,'-end');

% Information criteria used
IC = MultVarRegGASub_IC; numIC = length(IC);
% if real data, remove KLDist
if data_type == 2
    newIC = cell(numIC-1);
    for iccnt = 1:(numIC-1)
        newIC{iccnt} = IC{iccnt};
    end                 % information criteria loop
    IC = newIC; numIC = length(IC);
end

% get the matrix of subsets just once
MVRGASB_LoadData
if data_type == 1
    MVRGASB_SimData
    trusize = length(trueM.preds);
end
q = px; vars = [1:q]; 
subsets = VarSubset(X(:,[1:q])); numsubs = size(subsets,1);
qs = subsets(:,1); subsets = subsets(:,[2:end]);    % split number predictors from binary matrix

% initialize storage table
subset_scores = ones(numsubs,numIC)*Inf;
modpicks = zeros(numsubs,numIC);
equl_truemod = zeros(numIC,1);
incl_truemod = equl_truemod;

% output JAH changed 20120211 in octave 'output' has to go on mydir
outdir = [data_file([1:strfind(data_file,'.')-1]),filesep];
if (exist([mydir,filesep,'output',filesep,outdir],'dir') ~= 7)
    mkdir([mydir,filesep,'output',filesep],outdir);
end

outdir = [mydir,filesep,'output',filesep,outdir];
outfil = ['AllSubs_',stt];
clc,diary([outdir,outfil,'.out']);

% display settings
disp(repmat('@',1,50))
disp(sprintf('MultVarRegAllSub_DRV run on %4.0f%02.0f%02.0f_%02.0f%02.0f%02.0f',clock))
disp(repmat('#',1,50))
if data_type == 1
    disp(['Using Simulated Data: Data File: ',data_file])
    disp(sprintf('Random State: %10.0f',rnd_stat))
else
    disp(['Using Real Data: Data File: ',data_file])
end
disp(sprintf('Number Simulations: %0.0f\nNumber Observations: %0.0f',mciter,n))
disp(['Variables Included: ',strrep(strrep(strrep(mat2str(vars - HasConstant),' ',','),'[','{'),']','}')])
disp(['Covariance Smoother: ',regul_func])
disp(repmat('#',1,50))

% Perform the experiments
for mccnt = 1:mciter    
    if data_type == 1
        % simulate the data & remove redundant variables, if any
        MVRGASB_SimData; X = X(:,[1:q]);
        if rem(mccnt,dispgran) == 0; disp(sprintf('\tSimulation %0.0f of %0.0f',mccnt,mciter)); end;
    end
    for subcnt = 1:numsubs
        for iccnt = 1:numIC
            subset_scores(subcnt,iccnt) = MultVarRegGASub_IC(IC{iccnt},Y,X,subsets(subcnt,:),regul_func,trueM);
        end         % information criteria loop
    end             % subsets loop
    % now that we've got all the IC values, increment count of how many
    % times each picked the correct model (and others)
    if data_type == 1
        for iccnt = 1:numIC
            [val,ind] = min(subset_scores(:,iccnt));    % get the min of the ic
            mod = find(subsets(ind,:).*vars);       %    get the corresponding subset
            % equals true model
            equl_truemod(iccnt) = equl_truemod(iccnt) + isequal(trueM.preds,mod);
            % includes true model
            incl_truemod(iccnt) = incl_truemod(iccnt) + (sum(ISIN(trueM.preds,mod,1)) == trusize);
            % which did it pick?
            modpicks(ind,iccnt) = modpicks(ind,iccnt) + 1;
        end         % information criteria loop
    end
end                 % simulations loop

% final display of all model hits
if data_type == 1
    % drop all completely 0 rows
    allzero = (sum(modpicks,2) == 0);
    modpicks = modpicks(not(allzero),:);
    % prepare row headers
    rowheads = cell(sum(not(allzero)),1); zcnt = 0;
    for subcnt = 1:numsubs
        vrs = find(vars.*subsets(subcnt,:));
        % find the true model
        if isequal(vrs,trueM.preds)
            truesub = zcnt + 1;
        end    
        % make the row header
        if (allzero(subcnt) == 0)
            zcnt = zcnt + 1;
            vrs = sprintf('%d,',vrs - HasConstant);
            rowheads{zcnt,1} = ['{',vrs([1:(end-1)]),'}'];
        end
    end                 % subsets loop
    modpicks = 100*modpicks/mciter;
    tab1 = table2str(IC,modpicks,{'%0.2f'},1,rowheads);
    if truesub ~= 0
        tab1(truesub+3,1) = '*';
    end
    disp(' ')
    disp(['Model Hit Rates from ',num2str(mciter),' Simulations'])
    disp(tab1);
    disp('* Indicates the True Data Generating Model.')
    
    % display of % choosing exactly true model
    equl_truemod = 100*equl_truemod/mciter;
    [val,ind] = max(equl_truemod);
    tab2 = table2str({''},equl_truemod,{'%0.2f'},0,IC);
    disp(' ')
    disp(['Percent of ',num2str(mciter),' Simulations that Selected the True Structure'])
    disp(tab2([3:end],:))
    
    % display of % choosing including true model
    incl_truemod = 100*incl_truemod/mciter;
    [val,ind] = max(incl_truemod);
    tab2 = table2str({''},incl_truemod,{'%0.2f'},0,IC);
    disp(' ')
    disp(['Percent of ',num2str(mciter),' Simulations that Selected a Model Including the True Structure'])
    disp(tab2([3:end],:))    
end

% display scores for all subsets from final run
% prepare row headers
rwhds = cell(numsubs,1);
for rcnt = 1:numsubs
    vrs = find(vars.*subsets(rcnt,:)) - HasConstant;
    vrs = sprintf('%d,',vrs);
    rwhds{rcnt} = ['{',vrs([1:(end-1)]),'}'];
end                 % replications loop
tab = table2str(IC,subset_scores,{'%0.2f'},0,rwhds);
% find the min IC values
[minvals,mininds] = min(subset_scores);
% insert lines into table, at end of each section of same model size
grps = tabulate(qs); grps = grps([end:-1:1],2);
% take off the top and bottom lines, so tab size equals numsubs
tabtop = tab([1,2,3],:); tab = tab([4:(end-1)],:);
tabnew = []; grpcnt = 1; currgrp = 1;
for subcnt = 1:numsubs
    tabnew = [tabnew;tab(subcnt,:)];
    % insert dividing line if last of this group
    if grps(currgrp) == grpcnt
        tabnew = [tabnew;tabtop(1,:)];
        currgrp = currgrp + 1; grpcnt = 0;
    end
    grpcnt = grpcnt + 1;
end                 % subsets count
% display model scores table
disp(' '),disp('All Subsets Information Criteria Scores'),disp([tabtop; tabnew])
for iccnt = 1:numIC
    disp([IC{iccnt},' best subset ',rwhds{mininds(iccnt)},sprintf(': %0.2f, row %0.0f.',minvals(iccnt),mininds(iccnt))])
end                 % information criteria loop

% reset path
rmpath(data_path)

% clear junk
clear dispgran *cnt ind mod num* qs val vrs rowheads answer allzero newIC
clear stats vars rwhds numIC tab* currgrp grps mydir tmp

disp(sprintf('\nMultVarRegAllSub_DRV Required \n\t%1.4f Seconds\n\t%1.4f Minutes\n\t%1.4f Hours',etime(clock,drvstt)./[1,60,3600]));
disp([outfil , ' .mat and .out saved to']),disp(outdir)
disp(repmat('@',1,50))
diary off
save([outdir,outfil,'.mat']);

%{
JAH 20070421, added code to do real data 20070826, adapted for octave 3.4.3 20120211

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