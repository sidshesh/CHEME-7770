using DifferentialEquations
using Plots
using Statistics

function prob2d(X,Y,Z)
    alpha_x = 3.9e-2
    alpha_y = 4.3e-3
    beta_x = 6.1
    beta_y = 5.7
    x_y = 7.9e-4
    x_z = 12e-2
    y_z = 11e-3
    z_x = 1.3e-5
    n_zx = 2.32
    n_xy = 2
    n_xz = 2
    n_yz = 2
    del_z = 1.04
    del_y = 1.05

    S = 0.02
    #S = 10
    #S = 1e5
    #S = 0.2 #For Problem 2D
    #S = 50000 #For Problem 2D

    dX = ((alpha_x + (beta_x*S))/(1+S+((Z/z_x)^n_zx))) - X
    dY = ((alpha_y + (beta_y*S))/(1+S+((X/x_y)^n_xy))) - (del_y*Y)
    dZ = (1/(1+((X/x_z)^n_xz)+((Y/y_z)^n_yz))) - (del_z*Z)
    return [dX,dY,dZ]
end

function ODE_model!(du,u,model,t)
    du[1], du[2], du[3] = model(u[1], u[2], u[3])
end

init_conc = [0.0,0.0,0.0]
t_span = [0.0,150.0]

prob = ODEProblem(ODE_model!,init_conc,t_span,prob2d)
sol = solve(prob)
t=sol.t

X_values = sol[1,:]
Y_values = sol[2,:]
Z_values = sol[3,:]



#Small code modification for Problem 2D

#X_last = sol[1,floor(length(sol))]
#Y_last = sol[2,floor(length(sol))]
#Z_last = sol[3,floor(length(sol))]
#println("X Steady State value ", X_last)
#println("Y Steady State value ", Y_last)
#println("Z Steady State value ", Z_last)

#End of code modification for Problem 2D

plot(t,X_values,xlabel="Time",ylabel="Concentration",label="X vs t")
