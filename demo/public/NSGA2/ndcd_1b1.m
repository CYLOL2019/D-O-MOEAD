%% non-dominance and crowding distance, 1-by-1
function [archive] = ndcd_1b1(optimal_inds, archive)

%% Problem Definition





%% NSGA-II Parameters

nPop = length(optimal_inds);        % Population Size




%if params.fes == 0    
%% initialize the probability vector

%xpop      = (rand(nPop,nVar))';


%end

%if params.fes ~= 0
    
%% sample initial population

%xpop      = population.parameter';

%end
%% Initialization

% parameters for ind
empty_individual.parameter=[];
empty_individual.combination=[];
empty_individual.weights=[];
%
empty_individual.Cost=[];
empty_individual.Rank=[];
empty_individual.DominationSet=[];
empty_individual.DominatedCount=[];
empty_individual.CrowdingDistance=[];







pop    = repmat(empty_individual,nPop,1);

for i=1:nPop
    
    pop(i).parameter = optimal_inds(i).parameter;
    pop(i).combination = optimal_inds(i).combination;
    pop(i).weights = optimal_inds(i).weights;
    pop(i).Cost = optimal_inds(i).objective;

end



pop_Archive = repmat(empty_individual,nPop,1);

for i = 1 : nPop
    pop_Archive(i).parameter = archive(i).parameter;
    pop_Archive(i).combination = archive(i).combination;
    pop_Archive(i).weights = archive(i).weights;
    pop_Archive(i).Cost = archive(i).objective;
end


%oldpop = repmat(empty_individual,

    % Merging
    pop=[pop pop_Archive];  
 %%Non-Dominated Sorting  

    % Non-Dominated Sorting
    [pop, F]=NonDominatedSorting(pop);

    % Calculate Crowding Distance
    pop=CalcCrowdingDistance(pop,F);

    % Sort Population
    [pop, F]=SortPopulation(pop);
    

    % Truncate
size_combine = 0;
Ranks=[pop.Rank];
MaxRank=max(Ranks);
for r=1:MaxRank
    if size_combine + length(F{r}) > nPop
        break;
    end
    size_combine = size_combine + length(F{r});
end
% eliminate one by one
pop_eliminate = pop(F{r});
while (size_combine + length(pop_eliminate) > nPop)
    F_eliminate = 1:length(pop_eliminate);
    pop_eliminate = CalcCrowdingDistance(pop_eliminate,{F_eliminate});
    % Sort Based on Crowding Distance
    [~, CDSO]=sort([pop_eliminate.CrowdingDistance],'descend');
    pop_eliminate(CDSO(end))=[];
end

for i = 1 : size_combine
    archive(i).parameter = pop(i).parameter;
    archive(i).combination = pop(i).combination;
    archive(i).weights = pop(i).weights;
    archive(i).objective = pop(i).Cost;
end
% if the length of F{1} has exceeded popsize
if isempty(i)
    i = 0;
end
for j = 1 : length(pop_eliminate)
    archive(i+j).parameter = pop_eliminate(j).parameter;
    archive(i+j).combination = pop_eliminate(j).combination;
    archive(i+j).weights = pop_eliminate(j).weights;
    archive(i+j).objective = pop_eliminate(j).Cost;
end
    
    
end

