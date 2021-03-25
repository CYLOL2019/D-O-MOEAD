function v = evaluate_FV(ind)
global mocpo_params mop
v = zeros(mop.od,1);
portfolio = zeros(mop.pd,1);
portfolio(ind.combination) = ind.weights*mocpo_params.Rl;
% risk
v(1) = portfolio'*mocpo_params.Covariance*portfolio;
% return
v(2) = portfolio'*mocpo_params.u;
% normalization
v = mapminmax.apply(v,mocpo_params.objref);
% minimization ws
v(2) = -v(2);
end