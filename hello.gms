$title "let's get started"

SETS
    i       groups / i1*i20 /
    j       properties / j1*j8 /
    Bo(i)   oddly-bonded bonds
    SDx(i)
    SDy(i)
    SDz(i)
    H(i)    groups with higher-order bonds
    O(i)    groups with single bond !! define relation with H
;

PARAMETERS
    po(j)   property target
    ps(j)   property corresponding scale
    pL(j)   property lower bounds
    pU(j)   property upper bounds
    b(i)    group valency
    S(i)    # of single bonds on group i
    D(i)    # of double bonds on group i
    T(i)    # of triple bonds on group i
    Nmax    max # of any group
;

INTEGER VARIABLES
    n(i)    group multiplicity
    ZBo     number of odd bonds
    ZS
    ZD
    ZT
;

BINARY VARIABLES
    YSDx    exist singly and doubly-bonded groups
    YSDy    exist singly but not doubly-bonded groups
    YSDz    exist doubly but not singly-bonded groups
    YH      exist higher-order groups
;
SCALAR theta

EQUATIONS
* eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9,
eq1a,eq1b,eq1c,eq1d,eq2a,eq2b,eq3a,eq3b,eq4a1,eq4a2,eq4b1,eq4b2,eq4c1,eq4c2,eq5,eq6,eq7,eq8,eq9,eq10, eq11, eq12, eq13, eq14, eq15, eq16, eq17a, eq17b, eq18a, eq18b, eq19a, eq19b, eq20, eq21, eq22, eq23, eq24, eq25(i), eq26, eq27a, eq27b, eq28, eq29, eq30, eq31
;

eq1a..     po(1) =e= Tb0*log(sum(i,n(i)*c(i,1)));
eq1b..     po(2) =e= Tc0*log(sum(i,n(i)*c(i,2)));
eq1c..     po(3) =e= (sum(i,n(i)*c(i,3))+Pc2)^-2+Pc1;
eq1d..     po(4) =e= sum(i,n(i)*c(i,4))+Hv0;
eq2a..     po(5) =e= (sum(i,n(i)*c(i,5))-19.7779)+(sum(i,n(i)*c(i,6))+22.5981)*theta+(sum(i,n(i)*c(i,7))-10.7983)&theta^2;
eq2b..     theta =e= (Tm-298)/700;
eq3a..     po(6) =e= exp(f0(1)+po(11)*f1(1)+po(11)^2*f2(1))*po(3);
eq3b..     po(7) =e= exp(f0(2)+po(11)*f1(2)+po(11)^2*f2(2))*po(3);
eq4a1..    f0(1) =e= (-5.97616*(1-Te/po(2))+1.29874*(1-Te/po(2))^1.5-0.60394*(1-Te/po(2))^2.5-1.06841*(1-Te/po(2))^5)/(Te/po(2));
eq4b1..    f1(1) =e= (-5.03365*(1-Te/po(2))+1.11505*(1-Te/po(2))^1.5-5.41217*(1-Te/po(2))^2.5-7.46628*(1-Te/po(2))^5)/(Te/po(2));
eq4c1..    f2(1) =e= (-0.64771*(1-Te/po(2))+1.29874*(1-Te/po(2))^1.5-4.26979*(1-Te/po(2))^2.5-3.25259*(1-Te/po(2))^5)/(Te/po(2));
eq4a2..    f0(2) =e= (-5.97616*(1-Tc/po(2))+1.29874*(1-Tc/po(2))^1.5-0.60394*(1-Tc/po(2))^2.5-1.06841*(1-Tc/po(2))^5)/(Tc/po(2));
eq4b2..    f1(2) =e= (-5.03365*(1-Tc/po(2))+1.11505*(1-Tc/po(2))^1.5-5.41217*(1-Tc/po(2))^2.5-7.46628*(1-Tc/po(2))^5)/(Tc/po(2));
eq4c2..    f2(2) =e= (-0.64771*(1-Tc/po(2))+1.29874*(1-Tc/po(2))^1.5-4.26979*(1-Tc/po(2))^2.5-3.25259*(1-Tc/po(2))^5)/(Tc/po(2));
eq5..      po(8) =e= 1/4.1868*(po(5)+8.314*(1.45+0.45/(1-Tm/po(2))+0.25*po(12)*(17.11+25.22*(1-Tm/po(2))^(1/3)/(Tm/po(2))+1.742/(1-Tm/po(2)))));
eq6..      po(9) =e= -5.97214-log(po(3)/1.013)+6.09648*po(2)/po(1)+1.28862*log(po(1)/po(2))+-0.169347*(po(1)/po(2))^6;
eq7..      po(10) =e= 15.2518-15.6875*po(2)/po(1)-13.4721*log(po(1)/po(2))+0.4357*(po(1)/po(2))^6;
eq8..      po(11) =e= po(9)/po(10);
eq9..      po(12) =e= po(4)*((1-Te/po(2))/(1-po(1)/po(2)))^0.38;
eq10..      sum(i, n(i)) =g= 2;
eq11..      sum(Bo, n(Bo)) =e= 2*ZBo;
eq12..      sum(i, n(i)*b(i)) =g= 2*sum(i, n(i)-1);
eq13..      sum(i, n(i)*b(i)) =l= sum(i, n(i))*sum(i, n(i)-1);
eq14..      YSDx$(sum(SDx, n(SDx)) >= 1) =e= 1;
eq15..      YSDy$(sum(SDy, n(SDy)) >= 1) =e= 1;
eq16..      YSDz$(sum(SDz, n(SDz)) >= 1) =e= 1;
eq17a..     YSDx =l= sum(SDx, n(SDx));
eq17b..     sum(SDx, n(SDx)) =l= Nmax*YSDx*card(SDx);
eq18a..     YSDy =l= sum(SDx, n(SDx));
eq18b..     sum(SDy, n(SDy)) =l= Nmax*YSDy*card(SDy);
eq19a..     YSDz =l= sum(SDx, n(SDx));
eq19b..     sum(SDz, n(SDz)) =l= Nmax*YSDz*card(SDz);
eq20..      YSDy + YSDz -1 =l= YSDx;
eq21..      sum(i$(b(i)=1), n(i)) - sum(i$(b(i)=3), n(i)) - 2*sum(i$(b(i)=4), n(i)) =e= 2;
eq22..      sum(i$(S(i)>0), n(i)*S(i)) =e= 2*ZS;
eq23..      sum(i$(D(i)>0), n(i)*D(i)) =e= 2*ZD;
eq24..      sum(i$(T(i)>0), n(i)*T(i)) =e= 2*ZT;
Alias (i, ii);
eq25(i)..       sum(ii, n(ii)) =g= n(i)*(b(i) - 1)+2;
eq26..      YH$(sum(H, n(H) >= 1)) =e= 1;
eq27a..     YH =l= sum(H, n(H));
eq27b..     sum(H, n(H)) =l= Nmax*YH*card(H);
eq28..      sum(O$(b(O) = 1), n(O)) =l= sum(H, n(H)*S(H)) + Nmax*(1-YH)*card(SDy);
eq29..      sum(O$(b(O) = 2), n(O)) =l= sum(H, n(H)*D(H)) + Nmax*(1-YH)*card(SDy);
eq30..      sum(O$(b(O) = 3), n(O)) =l= sum(H, n(H)*T(H)) + Nmax*(1-YH)*card(SDy);
eq31..      sum(H, n(H)*(S(H)+D(H)+T(H))) - sum(O, n(O)) =e= 2*sum(H, n(H)-1);

