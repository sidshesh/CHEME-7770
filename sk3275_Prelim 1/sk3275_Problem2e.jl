using DifferentialEquations
using Plots
function prob2e1(X,Y,Z)
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

    S = 100
    #S = 5e5

    dX = ((alpha_x + (beta_x*S))/(1+S+((Z/z_x)^n_zx))) - X
    dY = ((alpha_y + (beta_y*S))/(1+S+((X/x_y)^n_xy))) - (del_y*Y)
    dZ = (1/(1+((X/x_z)^n_xz)+((Y/y_z)^n_yz))) - (del_z*Z)
    return [dX,dY,dZ]
end

function ODE_model!(du,u,model,t)
    du[1], du[2], du[3] = model(u[1], u[2], u[3])
end

#For Cell 1
#Initial concentration calculated using code for Problem 2d.
#init_conc = [0.0006443075432132409  0.5842951962116366  0.00034066995015958647] #S=0.2
init_conc = [5.734114476572251  0.005147127512792916  0.0004208860487737176] #S=50000
t_span = [0,100.0]
prob1 = ODEProblem(ODE_model!,init_conc,t_span,prob2e1)
sol = solve(prob1)
t = sol.t
Z = sol[3,:]

plot(t,Z,xlabel="t",ylabel="Concentration of Z",label="Cell 1")


#For Cell2
#init_conc = [0.0008053844290165512  0.7303689952645458  0.0004258374376994831] #S=0.2
init_conc = [7.167643095715314  0.0064339093909911455  0.0005261075609671469] #S=50000
t_span = [0,100.0]
prob2 = ODEProblem(ODE_model!,init_conc,t_span,prob2e1)
sol = solve(prob2)
t = sol.t
Z = sol[3,:]

plot!(t,Z,xlabel="t",ylabel="Concentration of Z",label="Cell 2")


#For Cell3
#init_conc = [0.00048323065740993066  0.43822139715872743  0.00025550246261968985] #S=0.2
init_conc = [4.300585857429188  0.003860345634594687  0.0003156645365802882] #S=50000
t_span = [0,100.0]
prob3 = ODEProblem(ODE_model!,init_conc,t_span,prob2e1)
sol = solve(prob3)
t = sol.t
Z = sol[3,:]

plot!(t,Z,xlabel="t",ylabel="Concentration of Z",label="Cell 3")
