clf
syms x y;

%[k1,k2,k3,k4] = deal(10.0);
[k1,k2,k3,k4] = deal(0.1);

%gamma_1 * R_T/ V2 = a1 
%gamma_3 * X_t/ V4 = a2
a1 = 5.0;
a2 = 10.0;

%I've used logspace() instead of constructing a normally separated array.
%This is because I want more points at lower range of kd_inv 
%where there is more variation in the output values
%rather than in the higher range where the outputs plateau
%This gives us a much better fit
kd_inv_values = logspace(-2,1,100);
count=1;
tic
x_values = zeros(1,numel(kd_inv_values));
y_values = zeros(1,numel(kd_inv_values));
theta_b_values = zeros(1,numel(kd_inv_values));
for i=kd_inv_values
    kd_inv = i;
    theta_b = kd_inv/(1+kd_inv);

    eqn1 = a1*theta_b == ((k1+1-x)/(k2+x))*(x/(1-x)) ;
    eqn4 = x>=0 ;
    eqn7 = x<=1 ;

    X = solve([eqn1,eqn4,eqn7],x);

    eqn2 = a2*X == ((k3+1-y)/(k4+y))*(y/(1-y)) ;
    eqn5 = y>=0 ;
    eqn8 = y<=1 ;
    Y = solve([eqn2,eqn5,eqn8],y);

    double(X);
    double(Y);
    disp("Iteration number :")
    disp(count)
    x_values(count)=double(X);
    y_values(count)=double(Y);
    theta_b_values(count)=double(theta_b);
    toc
    count=count+1;
end

hill_fit = @(A,x)  (A(3).*(x.^A(1)))./((A(2)^A(1))+(x.^A(1)));

x0 = [0.1 0.1 0.1];
hill_parameters_x = lsqcurvefit(hill_fit,x0,kd_inv_values,x_values)
hill_coefficient_x = real(hill_parameters_x(1));
new_x = hill_fit(hill_parameters_x,kd_inv_values);
figure(1)
hold on
plot(kd_inv_values,new_x,'b')
plot(kd_inv_values,x_values,'p')
legend("x^* calculated from Hill function","Actual x^*",'Location','southeast')
xlabel("1/k_D")
hold off

hill_parameters_y = lsqcurvefit(hill_fit,x0,kd_inv_values,y_values)
hill_coefficient_y = real(hill_parameters_y(1));
new_y = hill_fit(hill_parameters_y,kd_inv_values);
figure(2)
hold on
plot(kd_inv_values,new_y,'b')
plot(kd_inv_values,y_values,'p')
legend("y^* calculated from hill function","Actual y^*",'Location','southeast')
xlabel("1/k_D")
hold off

hill_parameters_theta_b = lsqcurvefit(hill_fit,x0,kd_inv_values,theta_b_values)
hill_coefficient_theta_b = real(hill_parameters_theta_b(1));
new_theta_b = hill_fit(hill_parameters_theta_b,kd_inv_values);
figure(3)
hold on
plot(kd_inv_values,new_theta_b,'b')
plot(kd_inv_values,theta_b_values,'p')
legend("\theta_B calculated from hill function","Actual \theta_B",'Location','southeast')
xlabel("1/k_D")
hold off

disp(hill_coefficient_x);
disp(hill_coefficient_y);
disp(hill_coefficient_theta_b);