function [v] = initidealpoint(x)

v = [];
for i  = 1 : length(x)
    sel = x(i).combination;
    vi   = augmecon2(sel);
    x(i).objective = vi;
    v = [v, vi];
end
v = min(v,[],2); 
end

%% augmecon2 gams driver
function [pf]= augmecon2(combination)
    m2wgdx(combination);
%     system('gams AUGMECON2_MIQP2I_Portfolio_sub --gdxin=MtoG'); 
    % AUGMECON2&CPLEX for mocpo_params.grid
    system('gams AUGMECON2_MIQP2I_Portfolio_payoff --gdxin=MtoG logOption=0'); 
    load('variables.csv');
    variables = variables(:,end-1:end)';
    pf = variables;
    % AUGMECON2 will jump the same Pareto optimal, making solutions less
    % then the grid number 
end

%% wgdx for GAMS
function m2wgdx(combination)
global gamsparams mocpo_params mop
    
    mu.name = 'mu';
    mu.val = zeros(mocpo_params.K,2);
    mu.val(:,1) = 1:mocpo_params.K;
    mu.val(:,2) = mocpo_params.u(combination);
    mu.type = 'parameter';
    mu.uels = gamsparams.i.uels;
    
    norg.name = 'nor';
    norg.val = zeros(4,2);
    norg.val(:,1) = 1:4;
    norg.val(:,2) = mocpo_params.nor;
    norg.type = 'parameter';
    norg.uels = gamsparams.n.uels;
    
    cov.name = 'cov';
    cov.val =  mocpo_params.Covariance(combination,combination);
    cov.form = 'full';
    cov.type = 'parameter';
    cov.uels = {gamsparams.i.uels,gamsparams.i.uels};
    
    g.name = 'g';
    g.uels = [];
    for jm = 1 : mop.od
        gstr = sprintf('g%d',jm-1);
        g.uels = [g.uels {gstr}];
    end
    
    
    wgdx('MtoG',gamsparams.i,gamsparams.n,...
        mu,norg,cov,g,...
        gamsparams.SumOneg,gamsparams.FloorC,gamsparams.RoundLot,gamsparams.k);
%     gdxInfo MtoG
end