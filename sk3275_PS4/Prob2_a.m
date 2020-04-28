    syms A B C v1 v2 v3 v4;
    [Vm1,Vm2] = deal(5.0);
    Vm3 = 1.0;
    Vm4 = 1.0;
    [Ks1,Ks2,Ks3,Ks4] = deal(5.0);
    [Ki1,Ki2] = deal(1.0);
    S_tot = 100.0;
    [I1,I2] = deal(0);
    v1 = (Vm1*A)/((1+(I1/Ki1))*(Ks1+A));
    v2 = (Vm2*A)/((1+(I2/Ki2))*(Ks2+A));
    v3 = (Vm3*B)/(Ks3+B);
    v4 = (Vm4*C)/(Ks3+C);

    eq1 = v1 - v3 == 0;
    eq2 = v2 - v4 == 0;
    eq3 = A + B + C == S_tot;
    eq4 = A >= 0;
    eq5 = B >= 0;
    eq6 = C >= 0;

    eqns = [eq1,eq2,eq3,eq4,eq5,eq6];

    S = solve(eqns,[A,B,C]);

    (double(S.A)),(double(S.B)),(double(S.C))
