$title "let's get started"
$onEolCom

SETS
    i       groups / i1*i4 /
    j       properties / j1*j12 /
    k       indices /1*2/
    l       indices /1*7/
    Bo(i)   oddly-bonded bonds /i1,i3/
    SDx(i)  /system.empty/
    SDy(i)  /i1*i4/
    SDz(i)  /system.empty/
    H(i)    groups with higher-order bonds /i2*i4/
    O(i)    groups with single bond /i1/ !! define relation with H
;

PARAMETERS
    pL(j)   property lower bounds
    pU(j)   property upper bounds
    nL(i)   groups lower bounds
    nU(i)   groups upper bounds
    b(i)    group valency
        / i1 1,
          i2 2,
          i3 3,
          i4 4/
    S(i)    # of single bonds on group i
        / i1 1,
          i2 2,
          i3 3,
          i4 4/
    D(i)  # of double bonds on group i
    T(i)  # of triple bonds on group i
    Nmax  max # of any group
    Te    evaporation temperature
    Tc    condensation temperature
    Tm    mean temperature
    Tb0
    Tc0
    Pc1
    Pc2
    Hv0
    c(i,l)  group constants
;

pL('j1') = 0;
pL('j2') = 0;
pL('j3') = 0;
pL('j4') = 0;
pL('j5') = 0;
pL('j6') = 1.1;
pL('j7') = 0;
pL('j8') = 0;
pL('j9') = -1000;
pL('j10') = -1000;
pL('j11') = 0;
pL('j12') = 0;

pU('j1') = 1000;
pU('j2') = 1000;
pU('j3') = 1000;
pU('j4') = 10000;
pU('j5') = 1000;
pU('j6') = 100;
pU('j7') = 14;
pU('j8') = 143.67;
pU('j9') = 1000;
pU('j10') = 1000;
pU('j11') = 1;
pU('j12') = 20.353;

D(i) = 0;
T(i) = 0;
nU(i) = 3;
nL(i) = 0;
Nmax = 3;
Te   = 272;
Tc   = 316;
Tm   = 294;

* the below constants are from M&G Table 2
Tb0  = 222.543; !! K
Tc0  = 231.239; !! K
Pc1  = 5.9827; !! bar
Pc2  = 0.108998; !!bar^{-0.5}
Hv0  = 11.733; !! kJ/mol

Table c(i,l)
         1       2       3         4      5       6       7
i1       0.8491  1.7506  0.018615  0.217  35.1152 39.5923 -9.9232
i2       0.7141  1.3327  0.013547  4.91   22.6346 45.0933 -15.7033
i3       0.2925  0.596   0.007259  7.962  8.9272  59.9786 -29.5143
i4       -0.0671 0.0306  0.001219  10.73  0.3456  74.0368 -45.7878
;


VARIABLES
    p(j)    property value
    f0(k)
    f1(k)
    f2(k)
    theta
    OF
;

INTEGER VARIABLES
    n(i)    group multiplicity
    ZBo     number of odd bonds
    ZS      number of single bonds
    ZD      number of double bonds
    ZT      number of triple bonds
;
ZBo.up = 10; ZS.up = 10; ZD.up = 10; ZT.up = 10;

BINARY VARIABLES
    YSDx    exist singly and doubly-bonded groups
    YSDy    exist singly not doubly-bonded groups
    YSDz    exist doubly not singly-bonded groups
    YH      exist higher-order groups
;


EQUATIONS
    eq1a    Normal boiling point (Tb) !! 1a-d are from M&G Table 1
    eq1b    Critical temperature (Tc)
    eq1c    Critical pressure (Pc)
    eq1d    Standard enthalpy of vaporization at 298K (Hv)
    eq2a,eq2b,eq3a,eq3b,eq4a1,eq4a2,eq4b1,eq4b2,eq4c1,eq4c2,eq5,eq6,eq7,eq8,eq9,eq10
* eq11, eq12, eq13, eq17a, eq17b, eq18a, eq18b, eq19a, eq19b, eq20, eq21, eq22, eq23, eq24, eq25(i), eq27a, eq27b, eq28, eq29, eq30, eq31,eq32(i)
* eq33(j),eq34(j),
ObjFun
;

eq1a..     p('j1') =e= Tb0*log(sum(i,n(i)*c(i,'1')));
eq1b..     p('j2') =e= Tc0*log(sum(i,n(i)*c(i,'2')));
eq1c..     p('j3') =e= power(sum(i,n(i)*c(i,'3'))+Pc2,-2)+Pc1;
eq1d..     p('j4') =e= sum(i,n(i)*c(i,'4'))+Hv0;
eq2a..     p('j5') =e= (sum(i,n(i)*c(i,'5'))-19.7779)+(sum(i,n(i)*c(i,'6'))+22.5981)*theta+(sum(i,n(i)*c(i,'7'))-10.7983)*power(theta,2);
eq2b..     theta =e= (Tm-298)/700;
eq3a..     p('j6') =e= exp(f0('1')+p('j11')*f1('1')+power(p('j11'),2)*f2('1'))*p('j3');
eq3b..     p('j7') =e= exp(f0('2')+p('j11')*f1('2')+power(p('j11'),2)*f2('2'))*p('j3');
eq4a1..    f0('1') =e= (-5.97616*(1-Te/p('j2'))+1.29874*(1-Te/p('j2'))**1.5-0.60394*(1-Te/p('j2'))**2.5-1.06841*(1-Te/p('j2'))**5)/(Te/p('j2'));
eq4b1..    f1('1') =e= (-5.03365*(1-Te/p('j2'))+1.11505*(1-Te/p('j2'))**1.5-5.41217*(1-Te/p('j2'))**2.5-7.46628*(1-Te/p('j2'))**5)/(Te/p('j2'));
eq4c1..    f2('1') =e= (-0.64771*(1-Te/p('j2'))+1.29874*(1-Te/p('j2'))**1.5-4.26979*(1-Te/p('j2'))**2.5-3.25259*(1-Te/p('j2'))**5)/(Te/p('j2'));
eq4a2..    f0('2') =e= (-5.97616*(1-Tc/p('j2'))+1.29874*(1-Tc/p('j2'))**1.5-0.60394*(1-Tc/p('j2'))**2.5-1.06841*(1-Tc/p('j2'))**5)/(Tc/p('j2'));
eq4b2..    f1('2') =e= (-5.03365*(1-Tc/p('j2'))+1.11505*(1-Tc/p('j2'))**1.5-5.41217*(1-Tc/p('j2'))**2.5-7.46628*(1-Tc/p('j2'))**5)/(Tc/p('j2'));
eq4c2..    f2('2') =e= (-0.64771*(1-Tc/p('j2'))+1.29874*(1-Tc/p('j2'))**1.5-4.26979*(1-Tc/p('j2'))**2.5-3.25259*(1-Tc/p('j2'))**5)/(Tc/p('j2'));
eq5..      p('j8') =e= 1/4.1868*(p('j5')+8.314*(1.45+0.45/(1-Tm/p('j2'))+0.25*p('j12')*(17.11+25.22*(1-Tm/p('j2'))**(1/3)/(Tm/p('j2'))+1.742/(1-Tm/p('j2')))));
eq6..      p('j9') =e= -5.97214-log(p('j3')/1.013)+6.09648*p('j2')/p('j1')+1.28862*log(p('j1')/p('j2'))-0.169347*(p('j1')/p('j2'))**6;
eq7..      p('j10') =e= 15.2518-15.6875*p('j2')/p('j1')-13.4721*log(p('j1')/p('j2'))+0.4357*(p('j1')/p('j2'))**6;
eq8..      p('j11') =e= p('j9')/p('j10');
eq9..      p('j12') =e= p('j4')*((1-Te/p('j2'))/(1-p('j1')/p('j2')))**0.38;
eq10..      sum(i, n(i)) =g= 2;
* eq11..      sum(Bo, n(Bo)) =e= 2*ZBo;
* eq12..      sum(i, n(i)*b(i)) =g= 2*sum(i, n(i)-1);
* eq13..      sum(i, n(i)*b(i)) =l= sum(i, n(i))*sum(i, n(i)-1);
* eq17a..     YSDx =l= sum(SDx, n(SDx));
* eq17b..     sum(SDx, n(SDx)) =l= Nmax*YSDx*card(SDx);
* eq18a..     YSDy =l= sum(SDx, n(SDx));
* eq18b..     sum(SDy, n(SDy)) =l= Nmax*YSDy*card(SDy);
* eq19a..     YSDz =l= sum(SDx, n(SDx));
* eq19b..     sum(SDz, n(SDz)) =l= Nmax*YSDz*card(SDz);
* eq20..      YSDy + YSDz -1 =l= YSDx;
* eq21..      sum(i$(b(i)=1), n(i)) - sum(i$(b(i)=3), n(i)) - 2*sum(i$(b(i)=4), n(i)) =e= 2;
* eq22..      sum(i$(S(i)>0), n(i)*S(i)) =e= 2*ZS;
* eq23..      sum(i$(D(i)>0), n(i)*D(i)) =e= 2*ZD;
* eq24..      sum(i$(T(i)>0), n(i)*T(i)) =e= 2*ZT;
* Alias (i, ii);
* eq25(i)..       sum(ii, n(ii)) =g= n(i)*(b(i) - 1)+2;
* eq27a..     YH =l= sum(H, n(H));
* eq27b..     sum(H, n(H)) =l= Nmax*YH*card(H);
* eq28..      sum(O$(b(O) = 1), n(O)) =l= sum(H, n(H)*S(H)) + Nmax*(1-YH)*card(SDy);
* eq29..      sum(O$(b(O) = 2), n(O)) =l= sum(H, n(H)*D(H)) + Nmax*(1-YH)*card(SDy);
* eq30..      sum(O$(b(O) = 3), n(O)) =l= sum(H, n(H)*T(H)) + Nmax*(1-YH)*card(SDy);
* eq31..      sum(H, n(H)*(S(H)+D(H)+T(H))) - sum(O, n(O)) =e= 2*sum(H, n(H)-1);
* eq32(i)..   n(i) =l= Nmax;
* eq33(j)..   p(j) =l= pU(j);
* eq34(j)..   p(j) =g= pL(j);

ObjFun..        OF   =e= p('j8')/p('j12');

n.lo(i)=nL(i);
n.up(i)=nU(i);

p.lo(j)=pL(j);
p.up(j)=pU(j);

* set intial values
p.l(j) = 1; !! completely arbitrary
p.l('j2') = max(Tc,Te); !! avoid negative bases for powers in eqs 4-5
n.l('i1') = 1; !! start with one of group 1

OPTION SYSOUT = ON;

option minlp = baron;

Model CAMP /all/;

CAMP.OPTFILE=1;
$onecho > baron.opt
# https://www.gams.com/latest/docs/S_BARON.html#BARONCompIIS
CompIIS = 5
IISOrder = 1
$offecho

Solve CAMP using MINLP minimise OF;