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

EQUATIONS
* eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9,
eq10, eq11, eq12, eq13, eq14, eq15, eq16, eq17a, eq17b, eq18a, eq18b, eq19a, eq19b, eq20, eq21, eq22, eq23, eq24, eq25(i), eq26, eq27a, eq27b, eq28, eq29, eq30, eq31
;

* eq1..	    ;
* eq2..	    ;
* eq3..	    ;
* eq4..	    ;
* eq5..	    ;
* eq6..	    ;
* eq7..	    ;
* eq8..	    ;
* eq9..	    ;
eq10..	    sum(i, n(i)) =g= 2;
eq11..	    sum(Bo, n(Bo)) =e= 2*ZBo;
eq12..	    sum(i, n(i)*b(i)) =g= 2*sum(i, n(i)-1);
eq13..	    sum(i, n(i)*b(i)) =l= sum(i, n(i))*sum(i, n(i)-1);
eq14..      YSDx$(sum(SDx, n(SDx)) >= 1) =e= 1;
eq15..      YSDy$(sum(SDy, n(SDy)) >= 1) =e= 1;
eq16..      YSDz$(sum(SDz, n(SDz)) >= 1) =e= 1;
eq17a..	    YSDx =l= sum(SDx, n(SDx));
eq17b..	    sum(SDx, n(SDx)) =l= Nmax*YSDx*card(SDx);
eq18a..	    YSDy =l= sum(SDx, n(SDx));
eq18b..	    sum(SDy, n(SDy)) =l= Nmax*YSDy*card(SDy);
eq19a..	    YSDz =l= sum(SDx, n(SDx));
eq19b..	    sum(SDz, n(SDz)) =l= Nmax*YSDz*card(SDz);
eq20..	    YSDy + YSDz -1 =l= YSDx;
eq21..	    sum(i$(b(i)=1), n(i)) - sum(i$(b(i)=3), n(i)) - 2*sum(i$(b(i)=4), n(i)) =e= 2;
eq22..	    sum(i$(S(i)>0), n(i)*S(i)) =e= 2*ZS;
eq23..	    sum(i$(D(i)>0), n(i)*D(i)) =e= 2*ZD;
eq24..	    sum(i$(T(i)>0), n(i)*T(i)) =e= 2*ZT;
Alias (i, ii);
eq25(i)..	sum(ii, n(ii)) =g= n(i)*(b(i) - 1)+2;
eq26..	    YH$(sum(H, n(H) >= 1)) =e= 1;
eq27a..	    YH =l= sum(H, n(H));
eq27b..	    sum(H, n(H)) =l= Nmax*YH*card(H);
eq28..	    sum(O$(b(O) = 1), n(O)) =l= sum(H, n(H)*S(H)) + Nmax*(1-YH)*card(SDy);
eq29..	    sum(O$(b(O) = 2), n(O)) =l= sum(H, n(H)*D(H)) + Nmax*(1-YH)*card(SDy);
eq30..	    sum(O$(b(O) = 3), n(O)) =l= sum(H, n(H)*T(H)) + Nmax*(1-YH)*card(SDy);
eq31..	    sum(H, n(H)*(S(H)+D(H)+T(H))) - sum(O, n(O)) =e= 2*sum(H, n(H)-1);

