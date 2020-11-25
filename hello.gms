$title let's get started

SETS
    i       groups / i1*i20 /
    j       properties / j1*j8 /
    B(i)    oddly-bonded bonds
    SDx(i)
    SDy(i)
    SDz(i)
    H(i)    groups with higher-order bonds
    O(i)    groups with single bond /*define relation with H*/
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
;

INTEGER VARIABLES
    n(i)    group multiplicity
    ZB      number of odd bonds
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

EQUATIONS
eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq11, eq12, eq13, eq14, eq15, eq16, eq17, eq18, eq19, eq20, eq21, eq22, eq23, eq24, eq25(i), eq26, eq27, eq28, eq29, eq30, eq31, eq32, eq33, eq34, eq35, eq36, eq37, eq38, eq39, eq40, eq41, eq42, eq43, eq44, eq45, eq46, eq47
;

eq1..	;
eq2..	;
eq3..	;
eq4..	;
eq5..	;
eq6..	;
eq7..	;
eq8..	;
eq9..	;
eq10..	sum(i, n(i)) =g= 2;
eq11..	sum(B, n(B)) =e= 2*ZB;
eq12..	sum(i, n(i)*b(i)) =g= 2*sum(i, n(i)-1);
eq13..	sum(i, n(i)*b(i)) =l= sum(i, n(i))*sum(i, n(i)-1);
eq14..  YSDx$(sum(SDx, n(SDx)) >= 1) =e= 1;
eq15..  YSDy$(sum(SDy, n(SDy)) >= 1) =e= 1;
eq16..  YSDz$(sum(SDz, n(SDz)) >= 1) =e= 1;
eq17a..	YSDx =l= sum(SDx, n(i));
eq17b..	sum(SDx, n(SDx)) =l= Nmax*YSDx*card(SDx);
eq18a..	YSDy =l= sum(SDx, n(SDx));
eq18b..	sum(SDy, n(SDy)) =l= Nmax*YSDy*card(SDy);
eq19a..	YSDz =l= sum(SDx, n(SDx));
eq19b..	sum(SDz, n(SDz)) =l= Nmax*YSDz*card(SDz);
eq20..	YSDy + YSDz -1 =l= YSDx;
eq21..	sum(i$(b(i)=1), n(i)) - sum(i$(b(i)=3), n(i)) - 2*sum(i$(b(i)=4), n(i)) =e= 2;
eq22..	sum(i$(S(i)>0), n(i)*S(i)) = 2*ZS;
eq23..	sum(i$(D(i)>0), n(i)*D(i)) = 2*ZD;
eq24..	sum(i$(T(i)>0), n(i)*T(i)) = 2*ZT;
eq25(i)..	sum(i, n(i)) =g= n(i)*(b(i) - 1)+2; !! this will probably break due to indexing control
eq26..	YH$(sum(H, n(i) >= 1)) =e= 1;
eq27a..	YH =l= sum(H, n(i));
eq27b..	sum(H, n(H)) =l= Nmax*YH*card(H);
eq28..	sum(O$(b(O) = 1), n(i)) =l= sum(H, n(H)*S(H)) + Nmax*(1-YH)*card(SDy);
eq29..	sum(O$(b(O) = 2), n(i)) =l= sum(H, n(H)*D(H)) + Nmax*(1-YH)*card(SDy);
eq30..	sum(O$(b(O) = 3), n(i)) =l= sum(H, n(H)*T(H)) + Nmax*(1-YH)*card(SDy);
eq31..	sum(H, n(H)*(S(H)+D(H)+T(H)) - sum(O, n(O)) =e= 2*sum(H, n(H)-1);
eq32..	;

