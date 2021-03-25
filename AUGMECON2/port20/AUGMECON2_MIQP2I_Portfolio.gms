*TITLE eps-Constraint Method for Multiobjective Optimization (EPSCM,SEQ=319)
$ontext
The eps-Constraint Method

$offtext

$inlinecom [ ]
$eolcom //
$STitle MIQCP model definitions
Option MIQCP = cplex;
Option IntVarUp = 0;
Sets
* Asset def
  i    assests
*  / 1*31 /
* k(1) risk, k(2) return
  K 'objective functions' /1*2/
* urisk, lrisk, ureturn, lreturn/
  n   normalization;
*  /1*4/;
Alias (i,j);


Parameter dir(k) direction of the objective functions 1 for max and -1 for min
   / 1  -1
     2  1
   /;

Parameter mu(i) expected return for each asset
*/
*$include 'port16_mu.txt' ;
*/;
Parameter nor(n) nomalization scalar
*/
*$include 'port16_n.txt';
*/;
Parameter cov(i,j) co-variance matrix for risk
*$offlisting
*$include 'port16_cov.txt' ;
*$onlisting
*----------------------------gdxin
$if not set gdxin $set gdxin MtoG
$GDXIN %gdxin%
$LOADIDX nor mu cov
$GDXIN
*-----------------------
Binary Variables
    Sel(i)   Selection signal;
Integer Variables
    x(i) 'weights of assets';
Variables
    z(k)     'objective function variable';


Parameter
    FloorP(i) floor contraints
    CeilP(i)    ceil contraints;
    FloorP(i) = 0.01;
    CeilP(i) = 1;
    
Scalar    
    PreA    pre-assignment constraint /30/
    
    Cardinality cardinality constraint /10/
       
    RoundLot round lot contraint /0.008/
    
    SumOne sum to one with respect to round lot /125/;
    

Equations
    ObjRisk    total risk
    ObjReturn  total profit
    SumOneCon  sum to one constraint
    CardiCon cardinality constraint
    FloorCon(i) floor constraints
    CeilCon(i) ceil consrraints
    PreACon pre-assianment constraint;
    
ObjRisk ..  (sum(i, x(i)*sum(j,cov(i,j)*x(j)))*sqr(RoundLot)-nor('2'))/(nor('1')-nor('2')) =e= z('1');
ObjReturn ..  (sum(i, x(i)*mu(i)*RoundLot)-nor('4'))/(nor('3')-nor('4')) =e= z('2');
*ObjRisk ..  sum(i, x(i)*sum(j,cov(i,j)*x(j)))*sqr(RoundLot) =e= z('1');
*ObjReturn ..  nor('3')-sum(i, x(i)*mu(i)*RoundLot) =e= z('2');
SumOneCon .. sum(i, x(i)) =e= SumOne;
CardiCon .. sum(i, Sel(i)) =e= Cardinality;
FloorCon(i) .. x(i) =g= 2*Sel(i);
CeilCon(i) .. x(i) =l= SumOne*Sel(i);
PreAcon .. Sel('30') =e= 1;

*---------------------------------------------------------------------------
$STitle eps-constraint method

Set k1(k) the first element of k, km1(k) all but the first elements of k;
k1(k)$(ord(k)=1) = yes; km1(k)=yes; km1(k1) = no;
Set kk(k)     active objective function in constraint allobj

Parameter
   rhs(k)     right hand side of the constrained obj functions in eps-constraint
   maxobj(k)  maximum value from the payoff table
   minobj(k)  minimum value from the payoff table
   numk(k) ordinal value of k starting with 1

Scalar
iter   total number of iterations
infeas total number of infeasibilities
elapsed_time elapsed time for payoff and e-sonstraint
start start time
finish finish time

Variables
   a_objval   auxiliary variable for the objective function
   obj        auxiliary variable during the construction of the payoff table
Positive Variables
   sl(k)      slack or surplus variables for the eps-constraints
Equations
   con_obj(k) constrained objective functions
   augm_obj   augmented objective function to avoid weakly efficient solutions
   allobj     all the objective functions in one expression;

con_obj(km1)..   z(km1) - dir(km1)*sl(km1) =e= rhs(km1);

* We optimize the first objective function and put the others as constraints
* the second term is for avoiding weakly efficient points

* objfun=max z1 + 0.001*(s1/r1+0.1 s2/r2+ 0.01*s3/r3+...)
augm_obj..
  sum(k1,dir(k1)*z(k1))+1.0e-3*sum(km1,power(10,-(numk(km1)-1))*sl(km1)/(maxobj(km1)-minobj(km1))) =e= a_objval;

allobj..  sum(kk, dir(kk)*z(kk)) =e= obj;

Model mod_payoff    / ObjRisk, ObjReturn, SumOneCon, CardiCon, FloorCon, CeilCon, PreACon, allobj / ;
Model mod_epsmethod / ObjRisk, ObjReturn, SumOneCon, CardiCon, FloorCon, CeilCon, PreACon, con_obj, augm_obj / ;
*--------------------------------------------------------------------------    
Parameter
  payoff(k,k)  payoff tables entries;
Alias(k,kp);

option optCr = 0, optca = 0;
option limrow=0, limcol=0, solprint=off ;
$offlisting;
$offsymxref;
$offsymlist;
$offuelxref;
$offuellist;

*file cplexopt /cplex.opt/;
*put cplexopt;
*put 'threads 4'/;
*put 'parallelmode 1'/;
*putclose cplexopt;
*mod_epsmethod.optfile=1;
*option optca=0.;
*mod_payoff.optfile=1;
*mod_epsmethod.optfile=1;

* Generate payoff table applying lexicographic optimization
*    z.fx(kk) = z.l(kk); // freeze the value of the last objective optimized 
*   kk(k++1) = kk(k);   // cycle through the objective functions
File fx  / 'solutions.txt' /;
File fx_  / 'variables.csv' /;
fx_.pw = 1e5;
*-------------------No Augment
*loop(kp,
*  kk(kp)=yes;
*  repeat
*    display kk;
*    solve mod_payoff using MIQCP maximizing obj;
*    display x.l, z.l;
*    display mod_payoff.modelstat;
*    if(    mod_payoff.modelStat <> %modelStat.optimal%
*        and mod_payoff.modelStat <> %modelStat.locallyOptimal%
*        and mod_payoff.modelStat <> %modelStat.feasibleSolution%
*        and mod_payoff.modelStat <> %modelStat.integerSolution%,
*        abort 'no optimal solution for mod_payoff';
*         );
*    payoff(kp,kk) = z.l(kk);
*    z.fx(kk) = z.l(kk);
*    kk(k++1) = kk(k);
*    
*
*  until kk(kp); kk(kp) = no;
** release the fixed values of the objective functions for the new iteration
*  put fx_ 0:5:0;
*  loop(i, put fx_ x.l(i):4:0);
*  loop(k, put fx_ z.l(k):16:8);
*  put fx_ /;
*  z.up(k) = inf; z.lo(k) =-inf;
*
*);
loop(kp,
  kk(kp)=yes;
    display kk;
    solve mod_payoff using MIQCP maximizing obj;
    display x.l, z.l;
    display mod_payoff.modelstat;
    if(    mod_payoff.modelStat <> %modelStat.optimal%
        and mod_payoff.modelStat <> %modelStat.locallyOptimal%
        and mod_payoff.modelStat <> %modelStat.feasibleSolution%
        and mod_payoff.modelStat <> %modelStat.integerSolution%,
        abort 'no optimal solution for mod_payoff';
         );
    payoff(kp,k) = z.l(k);
    kk(k++1) = kk(k);
    kk(kp) = no;
* release the fixed values of the objective functions for the new iteration
  z.up(k) = inf; z.lo(k) =-inf;

);
*---------------------------


*if (mod_payoff.modelstat<>1 and mod_payoff.modelstat<>8, abort 'no optimal solution for mod_payoff');


PUT fx ' PAYOFF TABLE'/   ;
loop (kp,
        loop(k, put payoff(kp,k):16:8);
        put /;
     );
put fx /;

*display payoff;
minobj(k)=smin(kp,payoff(kp,k));
*minobj(k)=0;
maxobj(k)=smax(kp,payoff(kp,k));

**$ontext
*$set fname h.%scrext.dat%

*gridpoints=max integer of km1 = 4149  or 1093
$if not set gridpoints $set gridpoints 2000
Set g grid points /g0*g%gridpoints%/
    grid(k,g) grid
Parameter
    gridrhs(k,g) rhs of eps-constraint at grid point
    maxg(k) maximum point in grid for objective
    posg(k) grid position of objective
    firstOffMax, lastZero some counters
    numk(k) ordinal value of k starting with 1
    numg(g) ordinal value of g starting with 0
    step(k) step of grid points in objective functions
    jump(k) jumps in the grid points' traversing
;
lastZero=1; loop(km1, numk(km1)=lastZero; lastZero=lastZero+1); numg(g) = ord(g)-1;
*grid(km1,g) = yes; // Here we could define different grid intervals for different objectives
grid(km1,g) = yes; 
maxg(km1) = smax(grid(km1,g), numg(g));
step(km1)=(maxobj(km1)- minobj(km1))/maxg(km1);
gridrhs(grid(km1,g))$(dir(km1)=-1) = maxobj(km1) - numg(g)/maxg(km1)*(maxobj(km1)- minobj(km1));
gridrhs(grid(km1,g))$(dir(km1)=1) = minobj(km1) + numg(g)/maxg(km1)*(maxobj(km1)- minobj(km1));
*display gridrhs;

PUT fx ' Grid points'/   ;
loop (g,
        loop(km1, put gridrhs(km1,g):16:8);
        put /;
     );
put fx /;
put fx 'Efficient solutions'/;

* Walk the grid points and take shortcuts if the model becomes infeasible
posg(km1) = 0;
iter=0;
infeas=0;

start=jnow;
*GtoM
*Parameter xval(i);
*Parameter zval(k);
*----------------------------
repeat
  rhs(km1) = sum(grid(km1,g)$(numg(g)=posg(km1)), gridrhs(km1,g));
  solve mod_epsmethod maximizing a_objval using MIQCP;
  iter=iter+1;
*if (mod_epsmethod.modelstat<>1 and mod_epsmethod.modelstat<>8,  // not optimal is in this case infeasible
  if (mod_payoff.modelStat <> %modelStat.optimal%
        and mod_payoff.modelStat <> %modelStat.locallyOptimal%
        and mod_payoff.modelStat <> %modelStat.feasibleSolution%
        and mod_payoff.modelStat <> %modelStat.integerSolution%,
    infeas=infeas+1;
    put fx iter:5:0, '  infeasible'/;
    lastZero = 0; loop(km1$(posg(km1)>0 and lastZero=0), lastZero=numk(km1));
*posg(km1)$(numk(km1)<=lastZero) = maxg(km1); // skip all solves for more demanding values of rhs(km1)
    posg(km1)$(numk(km1)<=lastZero) = maxg(km1);
  else
    put fx iter:5:0;
    loop(k, put fx z.l(k):16:8);
*    put fx_ iter:5:0;
    loop(i, put fx_ x.l(i):4:0);
    loop(k, put fx_ z.l(k):16:8);
    put fx_ /;
* GtoM code
*    xval(i)=x.l(i);
*    zval(k)=z.l(k);
*    Execute_UnloadIdx 'GtoM' xval zval;
*    put fx ' *** '; // put /;
    put fx ' *** '; 
    loop(km1, put fx sl.l(km1):16:8);
*    put fx ' *** '; // put /;
    put fx ' *** '; 
    jump(km1)=1;
*   find the first off max (obj function that hasn't reach the final grid point).
*   If this obj.fun is k then assign jump for the 1..k-th objective functions
    firstOffMax = 0;
    loop(km1$(posg(km1)<maxg(km1) and firstOffMax=0), firstOffMax=numk(km1));
    jump(km1)$(numk(km1)<=firstOffMax)=1+floor(abs(sl.L(km1))/step(km1));
*    put ' lastzero= ', lastzero:3:0 ;
    loop(km1, put fx jump(km1):5:0) ;
    loop(km1$(jump(km1)>1),put '   jump');
    put /;
    );
* Proceed forward in the grid
  firstOffMax = 0;
  loop(km1$(posg(km1)<maxg(km1) and firstOffMax=0), posg(km1)=min((posg(km1)+jump(km1)),maxg(km1)); firstOffMax=numk(km1));
  posg(km1)$(numk(km1)<firstOffMax) = 0;
*until sum(km1$(posg(km1)=maxg(km1)),1)= card(km1) and firstOffMax=0;
until sum(km1$(posg(km1)=maxg(km1)),1)= card(km1) and firstOffMax=0;

finish=jnow;
elapsed_time=(finish-start)*86400;

put /;
put 'Infeasibilities = ', infeas:5:0 /;
put 'Elapsed time: ',elapsed_time:10:2, ' seconds' / ;
*$offtext
*putclose fx; // close the point file
putclose fx;
putclose fx_;
**$offtext


*---------------------------------------------------
$stitle outpuf solutions
*File opt /MOPO_Solutions.opt/;
*
*put opt ' Solutions'/;
*loop(i, put x.l(i):10:5);
*put  /' Variance'/ z.l('1'):10:5;
*put /' fsum'/ sum(i,x.l(i)):10:5;
*put /'dmeam'/ z.l('2'):10:5;
*putclose opt;
*
**p1.optFile = 1;
**option limCol = 0, limRow = 0;
**
**solve p1 using minlp minimizing variance;
**
**if(    p1.modelStat <> %modelStat.optimal%
**   and p1.modelStat <> %modelStat.locallyOptimal%
**   and p1.modelStat <> %modelStat.feasibleSolution%
**   and p1.modelStat <> %modelStat.integerSolution%,
**   abort 'Could not solve p1 minimizing variance';
**);
*------------------------------------------------------
$stitle enumeration method
* just to be sure we also do complete enumeration and put results
* on a separate file.
*Set b / zero, one /;
*
*Parameter boole(b) / zero 0, one 1 /;
*
*Alias (b,b1,b2,b3,b4);
*
*File res / results.put /;
*put  res;
*
*Scalar min / 1.0e10 /;
*
*p1.solPrint = %solPrint.Summary%;
*p1.optFile  = 0;
*
*loop(i, put ' ',i.tl:4 );
*put '  variance';
*loop((b1,b2,b3,b4),
*   active.fx('hardware') = boole(b1);
*   active.fx('software') = boole(b2);
*   active.fx('show-biz') = boole(b3);
*   active.fx('t-bills')  = boole(b4);
*
*   solve p1 minimizing variance using rminlp;
*   put / boole(b1):5:0 boole(b2):5:0 boole(b3):5:0 boole(b4):5:0;
*
*   if(p1.solveStat <> %solveStat.normalCompletion%,
*      put '            *** failed solveStat=' p1.solveStat:0:0 ' modelStat=', p1.modelStat:0:0;
*      display 'Solver failed', p1.solveStat, p1.modelStat;
*   else
*      if(p1.modelStat <= %modelStat.locallyOptimal% or p1.modelStat = %modelStat.feasibleSolution%,
*         put variance.l:15:5;
*         if(variance.l < min,
*            put ' *';
*            min = variance.l;
*         );
*      else
*         put '    infeas, etc';
*      );
*   );
*);