clf
syms x y;

[k1,k2,k3,k4] = deal(10.0);
%[k1,k2,k3,k4] = deal(0.1);

%gamma_1 * R_T/ V2 = a1 
%gamma_3 * X_t/ V4 = a2

a1 = 5.0;
a2 = 10.0;
kd_inv_values = 0:0.5:2
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

hold on
plot(kd_inv_values,x_values)
plot(kd_inv_values,y_values)
plot(kd_inv_values,theta_b_values)
legend("x^*","y^*","\theta_B")
xlabel('1/k_D')
ylabel('Response values')
hold off