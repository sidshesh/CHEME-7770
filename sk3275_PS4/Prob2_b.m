clf
syms A B C v1 v2 v3 v4;

% I1_values = 0.01:10:1010.01;
% I2_values = 0.01:10:1010.01;
I1_values = logspace(-2,3,10);
I2_values = logspace(-2,3,10);
A_values = zeros(numel(I1_values),numel(I2_values));
tic
I1_count = 1;
for i= I1_values
    I2_count=1;
    for j = I2_values
        if (A_values(I1_count,I2_count)~=0)
            A_values(I2_count,I1_count) = A_values(I1_count,I2_count);
            I2_count=I2_count+1;
            continue
        end
            [Vm1,Vm2] = deal(5.0);
            Vm3 = 1.0;
            Vm4 = 1.0;
            [Ks1,Ks2,Ks3,Ks4] = deal(5.0);
            [Ki1,Ki2] = deal(1.0);
            S_tot = 100.0;
            I1 = i;
            I2 = j;
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

            ans_A = real(double(S.A));

            A_values(I2_count,I1_count) = ans_A;

            I2_count=I2_count+1;
    end
    I1_count = I1_count+1;
    toc
    disp(I1_count);
end

[X,Y]=meshgrid(I1_values,I2_values);
%mesh(X,Y,A_values)
h = axes;
surf(X,Y,A_values)
colormap winter;
title('Variation of [A] with I_1 and I_2 concentrations')
xlabel('[I_1]')
ylabel('[I_2]')
zlabel('[A]')
set(h,'xscale','log')
set(h,'yscale','log')