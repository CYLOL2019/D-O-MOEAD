function [pareto] = do_moead(mop)
%MOEAD run moea/d algorithms for the given mop.
    %global variable definition.
    global optimal_inds
    evolve(mop);
    pareto=[optimal_inds];
end

% The evoluation setp in MOEA/D
function evolve(mop)
    global params;
    selindex = 1:params.popsize;  
    selindex = random_shuffle(selindex);
    for i=1:params.popsize
        index = selindex(i);
        r = rand;
        useneighbour = r < params.neighbormating;
        ind = genetic_op_DEPM_swap(index,mop.domain,useneighbour);
        update(index, ind, useneighbour);

        clear ind obj useneighbour;
    end
end
% update the population.
% index is the subproblem's index in the main population.
% ind is the individual structure.
% useneighbour is a bool determine whether the neighbourhood of index, or the whole population should be updated.
% this procedure is also governed by a parameter from params: params.updatesize, which determine how many subproblem
% should be updated at most by this new individual in MOEAD-DE.

%%
function update(index, ind, updateneighbour)
global idealpoint params optimal_inds;

% collect the updation index
if (updateneighbour)
updateindex = params.neighbor(1:params.Tr,index);
else
updateindex = [1:params.popsize]';
end

updateindex = random_shuffle(updateindex);
neighborpoints = [optimal_inds(updateindex)]; 
% Incomplete ICDV for each neighbour
tmp_inds = repmat(ind,length(updateindex),1);
[v, tmp_inds] = evaluate_mocpo(tmp_inds,updateindex);
idealpoint = min(idealpoint, min(v,[],2)); 
newobj  = subobjective([params.W(:,updateindex)], v, idealpoint, params.dmethod);    %objective values of current solution
oldobj  = subobjective([params.W(:,updateindex)], [neighborpoints.objective], idealpoint, params.dmethod);        %previous objective values
C       = newobj < oldobj;    %new solution is better?
if (length(C(C==1)) <= params.updatesize)
    toupdate = updateindex(C);
else
    toupdate = randsample(updateindex(C), params.updatesize);
end  
for i = toupdate'
    [optimal_inds(i)]=deal(tmp_inds(updateindex==i));   %repace with the new one 
end
end



