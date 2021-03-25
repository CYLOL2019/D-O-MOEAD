*TITLE eps-Constraint Method for Multiobjective Optimization (EPSCM,SEQ=319)
$ontext
The eps-Constraint Method

$offtext

$inlinecom [ ]
$eolcom //
$STitle MIQCP model definitions
Option MIQCP = cplex;
Option IntVarUp = 0;



set i,n,k;
Alias (i,j);

Parameter mu(i),nor(n),cov(i,i),SumOne,FloorC,RoundLot,LambdaWs(k);
*----------------------------gdxin
$if not set gdxin $set gdxin MtoG
$GDXIN %gdxin%
$LOAD i n k mu nor cov SumOne FloorC RoundLot LambdaWs
$GDXIN
*-----------------------
Parameter dir(k) direction of the objective functions 1 for max and -1 for min
   / 1  -1
     2  -1
   /;


Integer Variables
    x(i) 'weights of assets';
Variables
    z(k)     'objective function variable'
    Ws  'Weights approach';        
    



    

Equations
    SumOneCon  sum to one constraint
    Ws_equation Weights equation
    FloorCon(i) floor constraints
    CeilCon(i) ceil consrraints;
    

Ws_equation .. Ws =e= LambdaWs('1')*(sum(i, x(i)*sum(j,cov(i,j)*x(j)))*sqr(RoundLot)-nor('2'))/(nor('1')-nor('2'))-LambdaWs('2')*(sum(i, x(i)*mu(i)*RoundLot)-nor('4'))/(nor('3')-nor('4'));
SumOneCon .. sum(i, x(i)) =e= SumOne;
FloorCon(i) .. x(i) =g= FloorC;
CeilCon(i) .. x(i) =l= SumOne;

*---------------------------------------------------------------------------
Model mod_Ws    / Ws_equation, SumOneCon, FloorCon, CeilCon / ;
*---------------------------------------------------------------------------
solve mod_Ws using MIQCP minimizing Ws;


Parameter xval(i);
xval(i)=x.l(i);
Execute_UnloadIdx 'GtoM' xval;


