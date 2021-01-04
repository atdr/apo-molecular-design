$title "let's get started"
$onEolCom

SETS
    i       groups / i1*i36 / !! CH3, CH2, CH C Br
    j       properties / j1*j12 /
    k       indices /1*2/
    l       indices /1*7/
    Bo(i)   oddly-bonded bonds /i1,i3,i5,i8,i10,i12,i13,i15*i21,i23*i26,i28,i31*i33,i35/
    SDx(i)  /system.empty/
    SDy(i)  /i1*i34/
    SDz(i)  /system.empty/
    H(i)    groups with higher-order bonds /i2*i4,i6*i9,i11,i12,i14,i22,i23,i25,i27,i29*i31,i34,i36/
    O(i)    groups with single bonds /i1,i5,i10,i13,i15*i21,i24,i26,i28,i32,i33,i35/ !! define relation with H

    x       integer cuts /1*10/
    dyn(x)  dynamic set of c
;

dyn(x) = no;

PARAMETERS
    pL(j)   property lower bounds
        /j1 0,
        j2  0,
        j3  0,
        j4  0,
        j5  0,
        j6  1.1,
        j7  0,
        j8  0,
        j9  -10,
        j10 -10,
        j11 0,
        j12 20.353/
    pU(j)   property upper bounds
        /j1 1000,
        j2  1000,
        j3  1000,
        j4  10000,
        j5  1000,
        j6  100,
        j7  14,
        j8  143.67,
        j9  10,
        j10 10,
        j11 1,
        j12 10000/
    nL(i)   groups lower bounds
    nU(i)   groups upper bounds
    b(i)    group valency
        /i1 1, i2  2, i3  3, i4  4, i5  1, i6  2, i7  2, i8  3, i9  4, i10 1,
        i11 2, i12 3, i13 1, i14 2, i15 1, i16 1, i17 1, i18 1, i19 1, i20 1,
        i21 1, i22 2, i23 3, i24 1, i25 1, i26 1, i27 2, i28 1, i29 2, i30 2,
        i31 3, i32 1, i33 1, i34 2, i35 1, i36 2/
    S(i)    # of single bonds on group i
        /i1 1, i2  2, i3  3, i4  4, i5  1, i6  2, i7  2, i8  3, i9  4, i10 1,
        i11 2, i12 3, i13 1, i14 2, i15 1, i16 1, i17 1, i18 1, i19 1, i20 1,
        i21 1, i22 2, i23 3, i24 1, i25 1, i26 1, i27 2, i28 1, i29 2, i30 2,
        i31 3, i32 1, i33 1, i34 2, i35 1, i36 2/
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
    nv(i,x)
    OFv(x)
    pv(j,x)
;

D(i) = 0;
T(i) = 0;
Nmax = 3;
nU(i) = Nmax;
nL(i) = 0;
Te   = 272;
Tc   = 316;
Tm   = 294;

* the below constants are from M&G Table 2
Tb0  = 222.543;     !! K
Tc0  = 231.239;     !! K
Pc1  = 5.9827;      !! bar
Pc2  = 0.108998;    !! bar^{-0.5}
Hv0  = 11.733;      !! kJ/mol

Table c(i,l)
        1        2        3         4       5        6          7
i1    0.8491   1.7506   0.018615   0.217   35.1152  39.5923   -9.9232
i2    0.7141   1.3327   0.013547   4.91    22.6346  45.0933   -15.7033
i3    0.2925   0.596    0.007259   7.962   8.9272   59.9786   -29.5143
i4    -0.0671  0.0306   0.001219   10.73   0.3456   74.0368   -45.7878
i5    1.5596   3.2295   0.025745   4.031   49.2506  59.384    -21.7908
i6    1.5597   3.0741   0.023003   9.456   35.2248  62.1924   -24.8156
i7    1.3621   2.7717   0.021137   8.602   37.6299  62.1285   -26.0637
i8    1.2971   2.5666   0.019609   14.095  21.3528  66.3947   -29.3703
i9    1.2739   2.6391   0.014114   19.91   10.2797  65.5372   -30.6057
i10   2.6364   6.2561   6.2561     11.754  48.4648  37.237    -13.0635
i11   2.0246   4.3756   4.3756     12.048  36.5885  47.6004   -22.8148
i12   1.7049   3.7063   3.7063     16.597  29.1848  52.3817   -30.8526
i13   3.9093   8.8073   0.036746   20.55   56.1685  46.9337   -31.3325
i14   0.5149   0.8543   0.018572   1.621   44.3567  44.5875   -23.282
i15   1.1916   1.7737   0.048565   7.352   63.2024  51.9366   -28.6308
i16   2.4231   4.5036  -0.00146    9.888   28.026  -7.1651     2.4332
i17   3.342    7.8956   0.028236   17.251  60.8262  41.9908   -20.4091
i18   2.567    5.2188  -0.005401   24.214  27.2107  2.7609     1.306
i19   5.1108   14.6038  0.009885   17.002  46.5577  48.2322   -20.4868
i20   3.1178   7.0058   0.025227   15.195  59.3032  67.8149   -20.9948
i21   1.7703   3.4393   0.020084   5.783   50.5604  38.9681   -4.7799
i22   1.3368   2.4217   0.017954   9.997   39.5784  41.8177   -11.0837
i23   0.8924   0.7889   0.014487   14.62   25.675   24.7281    4.2419
i24   2.5972   5.6064   0.015249   15.422  51.5048  44.4133   -19.6155
i25   2.5388   5.8013   0.010204   12.37   66.8423  102.4553  -43.3306
i26   2.7987   8.1745   0.011413   15.432  57.6861  64.0768   -21.048
i27   2.0948   4.2847   0.013049   16.048  44.1122  77.2155   -33.5086
i28   1.6525   2.8546   0.01079    17.257  53.7012  71.7948   -22.9685
i29   2.2514   4.5529   0.015863   11.831  44.6388  68.5041   -26.7106
i30   1.3841   3.0106   0.021186   9.493   41.4064  85.0996   -35.6318
i31   1.1222   2.1673   0.027454   12.636  30.1561  81.6814   -36.1441
i32   4.5871   12.9827  0.036523   21.923  58.2837  49.6388   -15.6291
i33   4.5311   10.9507  0.021056   29.64   63.7851  83.4744   -35.1171
i34   3.8069   9.5487   0.014899   29.173  51.1442  94.2934   -45.2029
i35   1.7618   3.7897   0.01401    6.144   45.9768  20.6417    -8.3297
i36   1.6767   4.587    0.010888   12.54   26.7371  21.7676    -6.4481
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
ZBo.l = 1;

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
    eq2a    vapour heat capacity (Poling)
    eq2b    theta for ditto
    eq3a    Pvp (evap)
    eq3b    Pvp (cond)
    eq4a1,eq4a2,eq4b1,eq4b2,eq4c1,eq4c2 !! some parameter for Pvp (woo woo) - index 1/2 refer to evap/cond
    eq5     liquid heat capacity
    eq6     alpha (for eq8)
    eq7     beta (for eq8)
    eq8     acentricity factor
    eq9     enthalpy of vaporization at Te
    eq10
    eq11 !! makes it unhappy
    eq12
    eq13
* eq17a, eq17b
    eq18a, eq18b
* eq19a, eq19b
    eq20
    eq21 !! makes it unhappy
    eq22
* eq23, eq24
    eq25(i)
    eq27a, eq27b, eq28
* eq29, eq30
    eq31
    ObjFun
    IntCut(x)
* nUp(i), nLo(i), pUp(j), pLo(j)
;

eq1a..      p('j1') =e= Tb0*log(sum(i,n(i)*c(i,'1')));
eq1b..      p('j2') =e= Tc0*log(sum(i,n(i)*c(i,'2')));
eq1c..      p('j3') =e= power(sum(i,n(i)*c(i,'3'))+Pc2,-2)+Pc1;
eq1d..      p('j4') =e= sum(i,n(i)*c(i,'4'))+Hv0;
eq2a..      p('j5') =e= (sum(i,n(i)*c(i,'5'))-19.7779)+(sum(i,n(i)*c(i,'6'))+22.5981)*theta+(sum(i,n(i)*c(i,'7'))-10.7983)*power(theta,2);
eq2b..      theta =e= (Tm-298)/700;
eq3a..      p('j6') =e= exp(f0('1')+p('j11')*f1('1')+power(p('j11'),2)*f2('1'))*p('j3');
eq3b..      p('j7') =e= exp(f0('2')+p('j11')*f1('2')+power(p('j11'),2)*f2('2'))*p('j3');
eq4a1..     f0('1') =e= (-5.97616*(1-Te/p('j2'))+1.29874*(1-Te/p('j2'))**1.5-0.60394*(1-Te/p('j2'))**2.5-1.06841*(1-Te/p('j2'))**5)/(Te/p('j2'));
eq4b1..     f1('1') =e= (-5.03365*(1-Te/p('j2'))+1.11505*(1-Te/p('j2'))**1.5-5.41217*(1-Te/p('j2'))**2.5-7.46628*(1-Te/p('j2'))**5)/(Te/p('j2'));
eq4c1..     f2('1') =e= (-0.64771*(1-Te/p('j2'))+2.41539*(1-Te/p('j2'))**1.5-4.26979*(1-Te/p('j2'))**2.5-3.25259*(1-Te/p('j2'))**5)/(Te/p('j2'));
eq4a2..     f0('2') =e= (-5.97616*(1-Tc/p('j2'))+1.29874*(1-Tc/p('j2'))**1.5-0.60394*(1-Tc/p('j2'))**2.5-1.06841*(1-Tc/p('j2'))**5)/(Tc/p('j2'));
eq4b2..     f1('2') =e= (-5.03365*(1-Tc/p('j2'))+1.11505*(1-Tc/p('j2'))**1.5-5.41217*(1-Tc/p('j2'))**2.5-7.46628*(1-Tc/p('j2'))**5)/(Tc/p('j2'));
eq4c2..     f2('2') =e= (-0.64771*(1-Tc/p('j2'))+2.41539*(1-Tc/p('j2'))**1.5-4.26979*(1-Tc/p('j2'))**2.5-3.25259*(1-Tc/p('j2'))**5)/(Tc/p('j2'));
eq5..       p('j8') =e= p('j5')+8.314*(1.45+0.45/(1-Tm/p('j2'))+0.25*p('j11')*(17.11+25.22*(1-Tm/p('j2'))**(1/3)/(Tm/p('j2'))+1.742/(1-Tm/p('j2'))));
eq6..       p('j9') =e= -5.97214-log(p('j3')/1.013)+6.09648*p('j2')/p('j1')+1.28862*log(p('j1')/p('j2'))-0.169347*(p('j1')/p('j2'))**6;
eq7..       p('j10') =e= 15.2518-15.6875*p('j2')/p('j1')-13.4721*log(p('j1')/p('j2'))+0.4357*(p('j1')/p('j2'))**6;
eq8..       p('j11') =e= p('j9')/p('j10');
eq9..       p('j12') =e= p('j4')*((1-Te/p('j2'))/(1-p('j1')/p('j2')))**0.38;
eq10..      sum(i, n(i)) =g= 2;
eq11..      sum(Bo, n(Bo)) =e= 2*ZBo;
eq12..      sum(i, n(i)*b(i)) =g= 2*(sum(i, n(i))-1);
eq13..      sum(i, n(i)*b(i)) =l= sum(i, n(i))*(sum(i, n(i))-1);
* eq17a..     sum(SDx, n(SDx)) =g= YSDx;
* eq17b..     sum(SDx, n(SDx)) =l= Nmax*YSDx*card(SDx);
eq18a..     sum(SDy, n(SDy)) =g= YSDy;
eq18b..     sum(SDy, n(SDy)) =l= Nmax*YSDy*card(SDy);
* eq19a..     sum(SDz, n(SDz)) =g= YSDz;
* eq19b..     sum(SDz, n(SDz)) =l= Nmax*YSDz*card(SDz);
eq20..      YSDy + YSDz -1 =l= YSDx;
* eq21..      sum(i$(b(i)=1), n(i)) - sum(i$(b(i)=3), n(i)) - 2*sum(i$(b(i)=4), n(i)) =e= 2;
eq21..      sum(i,n(i)*(2-b(i))) =e= 2; !! Odele-Macchietto version of above (simpler)
eq22..      sum(i$(S(i)>0), n(i)*S(i)) =e= 2*ZS;
* eq23..      sum(i$(D(i)>0), n(i)*D(i)) =e= 2*ZD;
* eq24..      sum(i$(T(i)>0), n(i)*T(i)) =e= 2*ZT;
Alias (i, ii);
eq25(i)..   sum(ii, n(ii)) =g= n(i)*(b(i) - 1)+2;
eq27a..     YH =l= sum(H, n(H));
eq27b..     sum(H, n(H)) =l= Nmax*YH*card(H);
eq28..      sum(O$(b(O) = 1), n(O)) =l= sum(H, n(H)*S(H)) + Nmax*(1-YH)*card(SDy);
* eq29..      sum(O$(b(O) = 2), n(O)) =l= sum(H, n(H)*D(H)) + Nmax*(1-YH)*card(SDy);
* eq30..      sum(O$(b(O) = 3), n(O)) =l= sum(H, n(H)*T(H)) + Nmax*(1-YH)*card(SDy);
eq31..      sum(H, n(H)*(S(H)+D(H)+T(H))) - sum(O, n(O)) =e= 2*(sum(H, n(H))-1);

ObjFun..    OF   =e= p('j8')/p('j12');

IntCut(x)$(dyn(x)).. sum(i,abs(nv(i,x)-n(i)))=g=1;

* nUp(i)..   n(i) =l= nU(i);
* nLo(i)..   n(i) =g= nL(i);
* pUp(j)..   p(j) =l= pU(j);
* pLo(j)..   p(j) =g= pL(j);
n.lo(i)=nL(i);
n.up(i)=nU(i);
p.lo(j)=pL(j);
p.up(j)=pU(j);

* set intial values
p.l('j1') = 250; !! completely arbitrary
p.l('j2') = 400;
p.l('j3') = 20;
p.l('j4') = 20;
p.l('j5') = 100;
p.l('j6') = 1.5;
p.l('j7') = 4.5;
p.l('j8') = 130;
p.l('j9') = -0.1;
p.l('j10') = -4;
p.l('j11') = 0.03;
p.l('j12') = 15;
n.l('i1') = 2;
* n.l('i16') = 1; !! start with some of group 1

OPTION SYSOUT = ON;

option minlp = baron;
* options optcr=0;
* options optca=0;

Model CAMD /all/;

CAMD.OPTFILE=1;
$onecho > baron.opt
# https://www.gams.com/latest/docs/S_BARON.html#BARONCompIIS
CompIIS = 1
$offecho

* Solve CAMD using MINLP minimise OF;
*Integer cut
nv(i,x)=0;
alias(x,xx);
loop(xx,
         solve CAMD minimising OF using MINLP;
*         solve modelPO minimising OF using MIP;
         nv(i,xx)=n.l(i);
         OFv(xx)=OF.l;
         pv(j,xx) = p.l(j);
* n.l(i) = 0;
* n.l('i1') = 1;
         dyn(xx)=yes;
);

display nv, pv, OFv;
