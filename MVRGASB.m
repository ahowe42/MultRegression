%{
usage: MUST BE CALLED FROM MultVarRegGASub_DRV OR MultRegGASub_DRV

Perform one simulation or replication of using the GA with information
criteria to perform robust subset multiple / multivariate regression.

Copyright (C) 2007 J. Andrew Howe; see below
%}

tic; drvstt = clock;
stt = sprintf('%4.0f%02.0f%02.0f_%02.0f%02.0f%02.0f',drvstt);
if (mvflag == 1)
	save_prefix = [mydir,outdir,'MVRGASB_',objec_func,'_',stt];
else
	save_prefix = [mydir,outdir,'MRGASB_',objec_func,'_',stt];
end
diary([save_prefix,'.out']);

% simulate data if required
if (data_type == 1); MVRGASB_SimData; end

% Display parameters
MVRGASB_DispParams
disp(save_prefix),disp([' started on ',stt])
%disp(sprintf('Random State: %10.0f',rnd_stat))
disp(repmat('#',1,50)), disp(' ')

% initialize population - start out with about half 1s
population = zeros(popul_size,px);
population(unidrnd(popul_size*px,1,ceil(popul_size*px/2))) = 1;
all_zero = find((sum(population == zeros(popul_size,px),2) == px).*[1:popul_size]');
population(all_zero,unidrnd(px,length(all_zero),1)) = 1;
population = logical(population);

% initialize more "things"
genscores = zeros(num_generns,3); termcount = 0; best_score = Inf;
genbest = zeros(num_generns,px); pltxyy = []; allscores = [];
allchroms = []; allscores = []; allconsts = cell(0); % JAH added allscores 20090803; JAH 20131214 allconsts added
best_chrom = []; % JAH 20090804
allscoresfor3d = []; % JAH 20100313
genconstwords = []; % JAH added 20140330

% Begin Genetic Algorithm
fhga = figure;
for gencnt = 1:num_generns
    %num_mut = ceil(prob_mutate*popul_size); % number of chromosomes in population to mutate    
	% JAH 20140330 add enforcer to keep the constant
	if (KeepConstant == 1) && (HasConstant == 1); population(population(:,1) == 0,1) = 1; end;

    pop_fitness = ones(popul_size,2)*Inf; % JAH 20140324 change pop_fitness to have 2nd column for total # constraints violated
    conststatus = cell(popul_size,1);	% JAH 20131028 will hold structs with all constraint test results
    % COMPUTE OBJECTIVE FUNCTION VALUES
    for popcnt = 1:popul_size
        % JAH check if this chromosome already evaluated 20090803
        % JAH 20131028 edited to accept all arguments from MultVarRegGASub_IC, and pass in the struct of constraint violation
        if gencnt > 1
            preveval = find(sum(allchroms == repmat(population(popcnt,:),size(allchroms,1),1),2) == px);
            if isempty(preveval)
				if (mvflag == 1)
                	pop_fitness(popcnt,1) = MultVarRegGASub_IC(objec_func,Y,X,population(popcnt,:),regul_func,trueM);
				else
					[pop_fitness(popcnt,1),lof,pen,bet,tmp] = MultRegGASub_IC(objec_func,Y,X,population(popcnt,:),regul_func,constlim);
                    conststatus(popcnt) = {tmp}; % JAH 20150630 can't put a struct directly into a cell without the {} in MATLAB, though it worked ok in octave
					pop_fitness(popcnt,2) = conststatus{popcnt}.SCORE;
				end
            else
                pop_fitness(popcnt,1) = allscores(preveval(1));
				if mvflag == 0
					conststatus(popcnt) = allconsts(preveval(1));
					pop_fitness(popcnt,2) = conststatus{popcnt}.SCORE;
				end
            end
        else
			if (mvflag == 1)
		    	pop_fitness(popcnt,1) = MultVarRegGASub_IC(objec_func,Y,X,population(popcnt,:),regul_func,trueM);
			else
				[pop_fitness(popcnt,1),lof,pen,bet,tmp] = MultRegGASub_IC(objec_func,Y,X,population(popcnt,:),regul_func,constlim);
                conststatus(popcnt) = {tmp}; % JAH 20150630 can't put a struct directly into a cell without the {} in MATLAB, though it worked ok in octave
				pop_fitness(popcnt,2) = conststatus{popcnt}.SCORE;
			end
        end
    end             % chromosomes loop
    % JAH 20140324 sort by first total # constraints violated then fitness score
    if constraint_check == 1
		[sortval,sortind] = sortrows(pop_fitness,[2,1]);
	else
		[sortval,sortind] = sortrows(pop_fitness,[1]);
	end
	minval = sortval(1,1); minind = sortind(1);
    
    % same some stuff before moving along JAH 20140324 only average 1st column, and add in constraints violated
    genscores(gencnt,:) = [minval,mean(pop_fitness(not(isnan(pop_fitness(:,1))),1)),sortval(1,2)];
    genbest(gencnt,:) = population(minind,:);
    if mvflag == 0; genconstwords = [genconstwords;conststatus{minind}.STR]; end;
    if plt3d == 1
		allscoresfor3d = [allscoresfor3d;sort(pop_fitness(:,1),'descend')'];
    end
    
    % UPDATE DATA FOR PLOT
    pltxyy = [pltxyy;[gencnt,genscores(gencnt,1),genscores(gencnt,2)]];
    if (pltflg == 1)
        figure(fhga)
        [AX,H1,H2] = plotyy(pltxyy(:,1),pltxyy(:,2),pltxyy(:,1),pltxyy(:,3)); xlabel('Generation');
        title(['GA Progress: Objective function ',objec_func],'interpreter','none');
        set(get(AX(1),'Ylabel'),'String','Minimum Value (o)','color','b');
        set(H1,'color','b','marker','o'); set(AX(2),'ycolor','b');
        set(get(AX(2),'Ylabel'),'String','Average Value (*)','color','r');
        set(AX(2),'ycolor','r'); set(H2,'color','r','marker','*');
        drawnow
    end

    if (gencnt == num_generns); break; end;       % don't bother with the next generation

    % EARLY TERMINATION ALLOWED?
    if genscores(gencnt,1) < best_score
        best_score = genscores(gencnt,1);
        best_chrom = population(minind,:);
        termcount = 1;
    elseif (genscores(gencnt,1) > best_score) && not(elitism)
        % if elitism is on, we can still do early termination with this
        termcount = termcount + 1;
    elseif genscores(gencnt,1) == best_score
        termcount = termcount + 1;
        if termcount >= nochange_terminate
            disp(['Early Termination On Generation ',num2str(gencnt),' of ',num2str(num_generns)]);
            genscores = genscores([1:gencnt],:);
            genbest = genbest([1:gencnt],:);
            break;
        end
    end

    % SELECTION OF NEXT GENERATION (JAH 20140324 why does this not use GAselect?)
    % JAH 20140324 above already sorted scores first by total # of constraints violated then scores
    %[val, stdindex] = sort(pop_fitness);
    % prepare bins for roulette - bigger bins at the beginning with lower scores
    bins =  cumsum([popul_size:-1:1]/(popul_size*(popul_size + 1)/2))';
    % roulette selection - find into which bin the random falls
    new_pop = sum(repmat(rand(1,popul_size),popul_size,1) >= repmat(bins,1,popul_size))+1;
    new_pop = population(sortind(new_pop),:);

    % CROSSOVER OPERATION ON NEW POPULATION (JAH 20140324 why does this not use GAcrossover?)
    new_pop = new_pop(randperm(popul_size),:);  % randomly permute rows
    for popcnt = 2:2:popul_size
        if rand <= prob_xover
            xoverpoint = unidrnd(px - 2) + 1;  % ensure xover point not on ends
            tmp = new_pop(popcnt - 1,:);
            % trade right-most portions
            new_pop(popcnt-1,[(xoverpoint + 1):px]) = new_pop(popcnt,[(xoverpoint + 1):px]);
            new_pop(popcnt,[(xoverpoint + 1):px]) = tmp([(xoverpoint + 1):px]);
        end
    end             % chromosomes loop
    
    % CHROMOSOME MUTATION JAH 20140324 changed to not looop
	mutation_chances = prob_mutate > rand(popul_size,px);
	new_pop(mutation_chances) = not(new_pop(mutation_chances));    
    
    % FIX ALL-ZERO CHROMOSOMES
    all_zero = find((sum(new_pop == zeros(popul_size,px),2) == px).*[1:popul_size]');
    new_pop(all_zero,unidrnd(px,length(all_zero),1)) = 1;    
    
    % SAVE ALL UNIQUE CHROMOSOMES & SCORES JAH 20090803
    [tab,tmp] = unique(population,'rows');
    allchroms = [allchroms;tab];
    allscores = [allscores;pop_fitness(tmp,1)];	% JAH 20140324 only get 1st column of pop_fitness
    % also save all constraint statuses JAH 20131214
	allconsts = [allconsts;conststatus(tmp)];

    % CONVEY BEST INDIVIDUAL INTO NEW POPULATION
    if elitism == 1
        population = [new_pop;population(minind(1),:)];
    else
        population = new_pop;
    end
    popul_size = size(population,1);
    if rem(gencnt,2) == 1
		if mvflag == 1
			disp(sprintf('Generation %0.0f of %0.0f: Best Score = %0.4f (p=%d), Early Termination = %0.0f',gencnt,num_generns,best_score,sum(best_chrom),termcount));
		else
			disp(sprintf('Generation %0.0f of %0.0f: Best Score = %0.4f (p=%d), Early Termination = %0.0f\n\tConstraints met = %d, violated = %d',gencnt,num_generns,best_score,sum(best_chrom),termcount,sum(pop_fitness(:,2)==0),sum(pop_fitness(:,2)>0)));
		end
    end
end                 % generations loop
if mvflag == 1
	disp(sprintf('Generation %0.0f of %0.0f: Best Score = %0.4f (p=%d), Early Termination = %0.0f',gencnt,num_generns,best_score,sum(best_chrom),termcount));
else
	disp(sprintf('Generation %0.0f of %0.0f: Best Score = %0.4f (p=%d), Early Termination = %0.0f\n\tConstraints met = %d, violated = %d',gencnt,num_generns,best_score,sum(best_chrom),termcount,sum(pop_fitness(:,2)==0),sum(pop_fitness(:,2)>0)));
end
% End Genetic Algorithm

figure(fhga)
[AX,H1,H2] = plotyy(pltxyy(:,1),pltxyy(:,2),pltxyy(:,1),pltxyy(:,3)); xlabel('Generation');
title(['GA Progress: Objective function ',objec_func],'interpreter','none');
set(get(AX(1),'Ylabel'),'String','Minimum Value (o)','color','b');
set(H1,'color','b','marker','o'); set(AX(2),'ycolor','b');
set(get(AX(2),'Ylabel'),'String','Average Value (*)','color','r');
set(AX(2),'ycolor','r'); set(H2,'color','r','marker','*');
% JAH 20170726 save both an .eps and .fig / .figo
drawnow; print(fhga,[save_prefix,'_GA.eps'],'-depsc');
hgsave(fhga,[save_prefix,'_GA']);   % save the figure JAH changed to print for octave 20120219

% do 3-dimensional IC score surface plot
if plt3d
    fhIC = figure;
	surf(allscoresfor3d),title(['3-dimensional Surface Plot of ',objec_func,' scores.'],'interpreter','none');
    xlabel('Population'),ylabel('Generations')
	% JAH 20170726 save both an .eps and .fig / .figo
    drawnow; print(fhIC,[save_prefix,'_3d.eps'],'-depsc');
	hgsave(fhIC,[save_prefix,'_3d']);   % save the figure 20120219 JAH changed to print for octave
end

% SUMMARY: GA_BEST = [scores,constraints violated,weights,frequencies,solutions] % JAH 20140328 added constraints violated and edited below accordingly
showtopsubs = 10;
% get the unique sorted solutions and scores from all generations
gen_results = [genscores(:,[1,3]),zeros(gencnt,2),genbest]; % combine min score and chromosome
[unique_gens,uniind] = unique(gen_results,'rows');			% get just the unique results
numuni = size(unique_gens,1);                           	% number unique results
[GA_BEST,srtind] = sortrows(unique_gens,1);					% sort in ascending order
top_best = min(showtopsubs,numuni); chromloc = [5:(px + 4)];
% JAH 20140330 get the constraint violation words; uniqued and sorted same as GA_BEST
if mvflag == 0
	GA_BEST_constwords = genconstwords(uniind,:);
	GA_BEST_constwords = GA_BEST_constwords(srtind,:);
end

% compute IC weights for all generation results
wghts = exp(-0.5*(GA_BEST(:,1) - GA_BEST(1,1)));
GA_BEST(:,3) = wghts/sum(wghts);

%[st,i] = dbstack; eval(['dbstop in ',st(i).name,' at ',num2str(st(i).line+1),';']);
% compute frequencies & prepare row headers
vars = [1:px];
rwhds = cell(top_best,1);
for pcnt = 1:numuni
    % compute frequency
    GA_BEST(pcnt,4) = sum(sum(repmat(GA_BEST(pcnt,chromloc),gencnt,1) == gen_results(:,chromloc),2) == px);
    % prepare row header
    if pcnt <= top_best
        vrs = find(vars.*GA_BEST(pcnt,chromloc)) - HasConstant;
        vrs = sprintf('%d,',vrs);
        rwhds{pcnt,1} = ['{',vrs([1:(end-1)]),'}'];
    end
end                 % dimensions loop

% display the best chromosomes and scores
tab = table2str({objec_func,'Const. Viol.','Weights','Frequency'},GA_BEST([1:top_best],[1,2,3,4]),{'%0.2f','%0.2f','%0.3f','%0.0f'},0,rwhds);
if mvflag == 0
	tab = [tab,[repmat(['-';' ';'-'],1,size(GA_BEST_constwords,2));GA_BEST_constwords([1:top_best],:);repmat('-',1,size(GA_BEST_constwords,2))]];
end
lin = repmat('=',1,60); disp(' '), disp(lin), disp('GA Complete')
disp(sprintf('\tTotal Solutions Evaluated - %0.0f\n\tUnique Solutions Evaluated - %0.0f\n\tTotal Solutions Possible - %0.0f',size(allchroms,1),size(unique(allchroms,'rows'),1),2^px-1))
disp(tab)
disp(sprintf('Average Solution Order: %0.0f',mean(sum(genbest,2))))
if (data_type == 1)
    disp(['True Model: ',strrep(strrep(strrep(mat2str(trueM.preds - HasConstant),' ',','),'[','{'),']','}')])
	if (mvflag == 1)
    	trueIC = MultVarRegGASub_IC(objec_func,Y,X(:,trueM.preds),ones(1,length(trueM.preds)),regul_func,trueM);
    	obj = sprintf('%s for true model: %0.2f',objec_func,trueIC);
	else
		[trueIC,lof,pen,bet,trueconst] = MultRegGASub_IC(objec_func,Y,X(:,trueM.preds),ones(1,length(trueM.preds)),regul_func,constlim);
		obj = sprintf('%s for true model: %0.2f\nConstraints Violated: %s (score = %0.2f)',objec_func,trueIC,trueconst.STR,trueconst.SCORE);
	end
    disp(obj)
end

% also show the score for the saturated model JAH 20140328 changed to test constraints and display the result
if (mvflag == 1)
	satIC = MultVarRegGASub_IC(objec_func,Y,X,ones(1,px),regul_func,trueM);
	obj = sprintf('%s for saturated model: %0.2f',objec_func,trueIC);
else
	[satIC,lof,pen,bet,satconst] = MultRegGASub_IC(objec_func,Y,X,ones(1,px),regul_func,constlim);
	obj = sprintf('%s for saturated model: %0.2f, Constraints Violated: %s (score = %0.2f)',objec_func,satIC,satconst.STR,satconst.SCORE);
end
disp(obj)
disp(lin)

% get rid of unneeded variables
clear all_zero AX H1 H2 best_chrom best_score bins fhga genbest gencnt xoverpoint
clear genscores ind lin minind minval mutcnt new_pop numuni pcnt wghts trueIC pen
clear pltxyy pop_fitness popcnt popmut population rwhds samp_covar showtopsubs
clear tab termcount tmp top_best unique_gens val vars vrs fhIC drvstt satIC bet
clear allscoresfor3d lof pen bet satconst sortval sortind obj srtind uniind lof

ga_toc = toc;
disp(sprintf('GA Completed in \n\t%1.4f Seconds\n\t%1.4f Minutes\n\t%1.4f Hours',ga_toc./[1,60,3600]));
diary off, save([save_prefix,'.mat']);

%{
JAH 20070403, copied from GAtemplate, adapted for octave 3.4.3 20120219

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