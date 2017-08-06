%{
usage: MultRegGaSub_DRV

Perform optimal subsets analysis for multiple regression (only the
regressors are subset) driven by the genetic algorithm, and compute
various information criteria scores.  If you tell it to use simulated
data, this can run multiple simulations, and determine how often
each criterion selected the "true" model.

Copyright (C) 2009 J. Andrew Howe; see below
%}

% clean up and prepare
close all; clc; clear;
nw = clock;
mvflag = 0;	% JAH 20120309 this is for multivariate regression

% find and store where I am located
mydir = dbstack; mydir = mydir(end);
mydir = which(mydir.file);
tmp = strfind(mydir,filesep);
mydir = mydir([1:(tmp(end)-1)]);

% first check for a multirun
[pn,fn] = uigetfile('*.mat', 'Existing Run to Augment',[mydir,filesep,'output',filesep]);
if pn == 0
	mult_run = {};
else
	mult_run = [fn,pn]; mult_run = mult_run([(length(mydir)+1):end]);
end

if isempty(mult_run)
	% setup parms data_type: 1=sim, 2 = real
	data_type = 1+isequal('Real',questdlg('Use Simulated or Real Data?','MultRegGA','Simulated','Real','Simulated')); 
	if data_type == 1   % simulated
		[data_file,data_path] = uigetfile('*.m', 'Load Simulated Data File',[mydir,filesep,'data',filesep]);
	else                % real
		[data_file,data_path] = uigetfile('*.*', 'Load Real Data File',[mydir,filesep,'data',filesep]);
	end

	% fix path
	addpath(data_path,'-end');

	% show objective functions & regul functions on the screen
	RFs = CovSmooth; ICs = MultRegGASub_IC;
	disp('Objective Functions'), disp(ICs)
	disp('Covariance Regularization Functions'), disp(RFs)

	% GA parameters
    prompt = {'Number Reps','Population Size','Number Generations','Early Termination','Elitism (1=on,0=off)',...
		'Crossover Prob.','Mutation Prob.','Progress Plot (1=during,0=end)','3d Plot (1=on,0=off)',...
		'Objective Func. (see terminal)','Pause Between Reps (s)','Random State (0=random)'};
	def ={'100','30','60','40','0','0.75','0.10','0','0','SBC','2','0'};
	answer = inputdlg(prompt,'MultRegGA Parameters',1,def,'off');
	num_reps = str2double(answer{1});
	popul_size = str2double(answer{2});
	num_generns = str2double(answer{3});
	nochange_terminate = str2double(answer{4});
	elitism = str2double(answer{5});
	prob_xover = str2double(answer{6});
	prob_mutate = str2double(answer{7});
	pltflg = str2double(answer{8});
	plt3d = str2double(answer{9});
	objec_func = char(answer{10});	
	pause_sec = str2double(answer{11});
    % JAH 20140327 added rand state
    if char(answer(12)) == '0';
		rnd_stat = sum(clock*1000000);
	else
		rnd_stat = str2num(char(answer(12)));
    end
	rand('state',rnd_stat); randn('state',rnd_stat);% randg('state',rnd_stat);
	if elitism == 1
		plt3d = 0;
    end
    
    % other parameters
    prompt = {'Regularization Func. (see terminal)','Simulated Observations (if sim)','Extra Noise Variables (if sim)','Keep Constant'};
	def ={'MLE/EB','500','10','1'};    
    answer = inputdlg(prompt,'Data / Model Parameters',1,def,'off');
    regul_func = char(answer{1});
	if data_type == 1; n = str2double(answer{2}); end;
	if data_type == 1; extra_vars = str2double(answer{3}); end;
	KeepConstant = str2double(answer{4});		% JAH added keep constant flag 20140330    
    
	% JAH added constraint checking inputs 20131214
	prompt = {'Enforce Constraints','Hypo Tests alpha','VIF Limit','Heteroskedasticity Test (-BPK,-W,-Ws)',...
		'Model F-test Weight','Multicollinearity Weight','Error Nonnormality Weight','Error Heteroskedasticity Weight'};
	def = {'1','0.05','7','-BPK','0.25','0.25','0.25','0.25'};
	answer = inputdlg(prompt,'Constraint Parameters',1,def,'off');
	constlim=[];
	constraint_check = str2double(answer{1});
	constlim.alpha = str2double(answer{2});
	constlim.VIFbound = str2double(answer{3});
	constlim.HeteroTest = answer{4};
	% get the constraint weights and normalize added JAH 20150630
	tmp = str2double(answer(5:end));
	tmp(tmp < 0) = 0;
	tmp = tmp / sum(tmp);
	constlim.p_weight = tmp(1);
	constlim.m_weight = tmp(2);
	constlim.e_weight = tmp(3);
	constlim.h_weight = tmp(4);	

	% make the data - range specific output directory, and the MULTI-RUN prefix
	df = data_file([1:strfind(data_file,'.')-1]);
	outdir = [filesep,df];
	if exist([mydir,filesep,'output',outdir],'dir') ~= 7
		mkdir([mydir,filesep,'output'],outdir);
	end
	outdir = [filesep,'output',outdir,filesep];
else
	answer = inputdlg({'Number Repetitions to Add'},'MultVarRegGA Parameters',1,{'100'},'off');
	num_reps = str2double(answer{1});
end

if isempty(mult_run)
    % BEGIN A NEW MULTI-RUN
    MRscores = []; MRfil = []; MRchroms = []; mrcnt = 1; pretottim = 0; MRwords = [];
    mult_run = [sprintf(['MR_',df,'_%4.0f%02.0f%02.0f_%02.0f%02.0f%02.0f.mat'],nw)];
    mult_run = [outdir,mult_run]; save([mydir,mult_run]);
else
    % AUGMENT A MULTI-RUN
    t = num_reps;% save number reps to add
    load([mydir,mult_run]);
    num_reps = num_reps + t; % augment number reps
    mrcnt = mrcnt + 1;
    clear t; save([mydir,mult_run]);
end

% DO THE WORK
MVRGASB_LoadData
while mrcnt <= num_reps
    % update progress
    clc;
	disp(sprintf('Starting Rep %d of %d; %0.2f%% completed!',mrcnt,num_reps,100*(mrcnt-1)/num_reps))
	disp(sprintf('Modeling Required \n\t%1.4f Seconds\n\t%1.4f Minutes\n\t%1.4f Hours',etime(clock,nw)./[1,60,3600]));
	pause(pause_sec);
    % do the work
    MVRGASB
    clear f t allscores gen_results stt
    % close all figures
    if num_reps > 1; close all; end;
    % save the stuff
    load([mydir,mult_run])
    MRscores = [MRscores;[mrcnt,GA_BEST(1,[1,2]),ga_toc]];	% JAH 20140328 added in # constraints violated
    MRwords = [MRwords;GA_BEST_constwords(1,:)];
    MRfil = [MRfil;StrPad(save_prefix([(length(mydir)+1):length(save_prefix)]),100,'R',' ')];
    MRchroms = [MRchroms;[mrcnt,GA_BEST(1,chromloc)]];
    clear GA_BEST save_prefix chromloc ga_toc GA_BEST_constwords
    mrcnt = mrcnt + 1; save([mydir,mult_run]);
end
tottim = etime(clock,nw);           % get total modeling time
tottim = pretottim + tottim; pretottim = tottim;
num_reps = size(MRchroms,1);        % ensure num_reps is accurate
mrcnt = num_reps;                   % same for mrcnt
save([mydir,mult_run]);                % final save

% reset path
rmpath(data_path)

% SUMMARY - only if there were multiple runs
if num_reps > 1; MVRGASB_MRSummary; end;

%{
JAH 20090228, adapted for octave 3.4.3 20120309

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