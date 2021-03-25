% Main function for ICDV-MOEA/D
function main(run_index)
SetPATH;
global mop params mocpo_params gamsparams evalCounter optimal_inds T1;
path('../../problem',path); 
path('../../problem/portfolio problem',path); 
path('../public',path);
path('../public/NSGA2',path);

problems = {'port1','port2','port3','port4','port5',...
    'port6','port7','port8','port9','port10',...
    'port11','port12','port13','port14','port15',...
    'port16','port17','port18','port19','port20'};
nproblem = length(problems);
run_index = 1; 
for pn = nproblem:nproblem%1 : nproblem
    for r = run_index : run_index
        T1 = clock;
        % paramters for ICDV-MOEA/D
        params = loadparams_do();    
        % parameters for the multiobjective constrained portfolio optimization
        [mocpo_params,gamsparams] = loadparams_mocpo(char(problems(pn)));
        evalCounter = 0;
        tic;
        %init
        init(mop)
        % main loop
        ps = [optimal_inds.combination;optimal_inds.weights];
        pf = [optimal_inds.objective];
        T2 = clock;
        fes = [evalCounter;etime(T2,T1)];
        while ~terminate()
            evalCounter         
            [pareto] = do_moead(mop);
            ps = [ps,[optimal_inds.combination;optimal_inds.weights]];
            pf = [pf,pareto.objective];
            T2 = clock;
            fes = [fes,[evalCounter;etime(T2,T1)]];
        end
    end
    pf(2,:) = -pf(2,:);
    pf = mapminmax.reverse(pf,mocpo_params.objref);
    sdir = sprintf("../data/%s/run%d",char(problems(pn)),r);
    if ~exist(sdir, 'dir')
       mkdir(sdir)
    end
    runtime = toc;
    sfile = sprintf('%s/data.mat',sdir);
    disp(sfile);
    save(sfile,'pf','ps','fes','runtime');        
end


end

function y =terminate()
    global params evalCounter T1;
    T2 = clock;
    y = etime(T2,T1)>= params.stoptime;
%     y = evalCounter>=params.evaluation;
end
