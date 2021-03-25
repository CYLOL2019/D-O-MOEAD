function [v, x] = evaluate_mocpo(x,index)
%EVALUATE function evaluate an individual structure of a vector point with
%the given multiobjective problem.

%   Detailed explanation goes here
%   prob: is the multiobjective problem.
%   x: is a vector point, or a individual structure.
%   v: is the result objectives evaluated by the mop.
%   x: if x is a individual structure, then x's objective field is modified
%   with the evaluated value and pass back.

    v = [];
    if isstruct(x)
        for i  = 1 : length(index)
            sel = x(i).combination;
            % Get the local optimal CDV with an ICDV
            x(i).weights  = GAMSeachWs(sel,index(i));
            % Get the F-value of the CDV, which consumes one fitness
            % evaluation
            x(i).objective = evaluate_FV(x(i));
            v = [v, x(i).objective];
        end
    else
        X = x;
        v = prob.func(X);
    end
    
    % set the evaluation counter.
    global evalCounter;
    % Assume the cost of each grid is the same as one evaluation
    evalCounter = evalCounter+length(index); 
end

%% augmecon2 gams driver
function xval = GAMSeachWs(combination,subpindex)

    m2wgdx(combination,subpindex);
%     gdxWhos 'MtoG'
    % GAMS Tchebycheff
%     system('gams AUGMECON2_MIQP2I_Portfolio_eachWs --gdxin=MtoG '); 
    system('gams AUGMECON2_MIQP2I_Portfolio_eachWs --gdxin=MtoG logOption=0'); 
    irgdx 'GtoM'
%     xval 

    
    
%     % return mat directly*************************
%     load('variables.csv');
%     variables = variables(:,end-1:end)';
%     pf = zeros(mop.od,mocpo_params.grid);    
%     pf(:,1:size(variables,2)) = variables;
%     % AUGMECON2 will jump the same Pareto optimal, making solutions less
%     % then the grid number
%     if size(variables,2)<mocpo_params.grid
%         pf(:,size(variables,2)+1:end) = NaN;
%     end  
end

%% wgdx for GAMS
function m2wgdx(combination, subpindex)
global gamsparams mocpo_params params
    
    mu.name = 'mu';
    mu.val = mocpo_params.u(combination);
    mu.form = 'full';
%     mu.val = zeros(mocpo_params.K,2);
%     mu.val(:,1) = 1:mocpo_params.K;
%     mu.val(:,2) = mocpo_params.u(combination);
    mu.type = 'parameter';
    mu.uels = gamsparams.i.uels;
    
    norg.name = 'nor';
    norg.val = mocpo_params.nor;
    norg.form = 'full';
%     norg.val = zeros(4,2);
%     norg.val(:,1) = 1:4;
%     norg.val(:,2) = mocpo_params.nor;
    norg.type = 'parameter';
    norg.uels = gamsparams.n.uels;
    
    cov.name = 'cov';
    cov.val =  mocpo_params.Covariance(combination,combination);
    cov.form = 'full';
    cov.type = 'parameter';
    cov.uels = {gamsparams.i.uels,gamsparams.i.uels};
    
    LambdaWs.name = 'LambdaWs';
    LambdaWs.val = params.W(:,subpindex);
    LambdaWs.form = 'full';
    LambdaWs.type = 'parameter';
    LambdaWs.uels = gamsparams.k.uels;
    
    
        
    wgdx('MtoG',gamsparams.i,gamsparams.n,gamsparams.k,...
        mu,norg,cov,LambdaWs,...
        gamsparams.SumOneg,gamsparams.FloorC,gamsparams.RoundLot);
%     gdxInfo MtoG
end