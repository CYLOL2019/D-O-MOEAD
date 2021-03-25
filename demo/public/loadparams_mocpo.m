function [mocpo_params,gamsparams] = loadparams_mocpo(testname)
global mop
%LOADPARMS_ICDV init the parameters for D&O-MOEA/D

% Parameters for MOCPO (constraints)
mocpo_params.K = 10;
mocpo_params.Lb = 0.01;
mocpo_params.Ub = 1;
mocpo_params.Rl = 0.008;
mocpo_params.PreAss = [30];
% Since 1 is divisible for Rl, it can be regarded as mix integer
% programming
mocpo_params.SumOne = 1/mocpo_params.Rl; 
mocpo_params.nor = zeros(4,1);
% Parameters from files
pfile = sprintf('%s.txt',testname);
input = textread(pfile);
[mocpo_params.NoA, mocpo_params.u, mocpo_params.Covariance] = DataInput(input);
% Normalization for the AUGMECON2&CPLEX
norfile = sprintf('../../problem/gamsdata/%s_n.txt',testname);
tmp = load(norfile);
[mu,~] = sort(mocpo_params.u,'descend');
% Least multipliers that y*Rl>=Lb
least_multipliers = ceil(mocpo_params.Lb/mocpo_params.Rl);
maxweight = mocpo_params.SumOne-least_multipliers*(mocpo_params.K-1);
weight = [maxweight; least_multipliers*ones(mocpo_params.K-1,1)]*mocpo_params.Rl;
tmp(3,2) = mu(1:mocpo_params.K)'*weight;
% normalization scalar
portnor=[tmp(1,2),tmp(2,2);tmp(3,2),tmp(4,2)];
[NULL,mocpo_params.objref] = mapminmax(portnor,0,1);
mocpo_params.nor(3) = portnor(2,1);
mocpo_params.nor(4) = portnor(2,2);
mocpo_params.nor(1) = portnor(1,1);
mocpo_params.nor(2) = portnor(1,2);
% testmop
testmop(testname, mocpo_params.NoA);% MOP could be obtained by function 'testmop.m'.

  % gamsparams wgdx
gamsparams.i.name = 'i';
gamsparams.i.uels = [];
for jm = 1 : mocpo_params.K
    uelstr = sprintf('%d',jm);
    gamsparams.i.uels = [gamsparams.i.uels {uelstr}];
end
gamsparams.n.name = 'n';
gamsparams.n.uels = [];
for jm = 1 : mop.od*2
    uelstr = sprintf('%d',jm);
    gamsparams.n.uels = [gamsparams.n.uels {uelstr}];
end
gamsparams.k.name = 'k';% number of objectives
gamsparams.k.uels = [];
for jm = 1 : mop.od
    uelstr = sprintf('%d',jm);
    gamsparams.k.uels = [gamsparams.k.uels {uelstr}];
end

gamsparams.SumOneg.name = 'SumOne';
gamsparams.SumOneg.val = mocpo_params.SumOne;
gamsparams.SumOneg.type = 'parameter';

gamsparams.FloorC.name = 'FloorC';
gamsparams.FloorC.val = least_multipliers;
gamsparams.FloorC.type = 'parameter';

gamsparams.RoundLot.name = 'RoundLot';
gamsparams.RoundLot.val = mocpo_params.Rl;
gamsparams.RoundLot.type = 'parameter';

end
%% DataInput
function [NoA, u, Covariance] = DataInput(input)
PortSize1 = [31 85 89 98 225];
NoA = input(1);
if ismember(NoA, PortSize1)
diagData = ones(1,NoA);
u = input(2:NoA+1,1);
variance = input(2:NoA+1,2);
cov = input(NoA+2:end,3);
A = tril(ones(NoA),0);

A(logical(A)) = cov;

cov = A + A' - diag(diagData);

Covariance = variance*variance'.*cov;
else
u = input(2:NoA+1,1);
Covariance = input(NoA+2:end,3);
diagData = ones(1,NoA);
A = tril(ones(NoA),0);

A(logical(A)) = Covariance;
diagData = diag(diagData)+1;
Covariance = A + A';
Covariance = Covariance./diagData;

end

end