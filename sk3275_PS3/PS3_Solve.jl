include("Include.jl")

s_matrix = DataFrame(CSV.File("S_Matrix.csv",header=false))
s_matrix = convert(Matrix,s_matrix)

E = 0.01E-3 # steady-state enzyme concentration (units [mmol/gDW])

# Metabolic kcats for v1, v2, v3, v4, v5f, v5r
kcat_v1 = (203)*(3600)   # v1
kcat_v2 = (34.5)*(3600)  # v2
kcat_v3 = (249)*(3600)   # v3
kcat_v4 = (88.1)*(3600)  # v4
kcat_v5f = (13.7)*(3600) # v5f
kcat_v5r = (13.7)*(3600) # v5r

# Metabolic Vmax for v1, v2, v3, v4, v5f, v5r
v1_max = (kcat_v1)*(E)
v2_max = (kcat_v2)*(E)
v3_max = (kcat_v3)*(E)
v4_max = (kcat_v4)*(E)
v5f_max = (kcat_v5f)*(E)
v5r_max = (kcat_v5r)*(E)

# Using Park et al as reference
v1 = v1_max*((4.67E-3)/(3.92E-4 + 4.67E-3))*((1.49E-2)/(1.54E-4 + 1.49E-2))
v2 = v2_max*(1)
v3 = v3_max*((2.55E-4)/(1.55E-3 + 2.55E-4))
v4 = v4_max*(1)
v5f = v5f_max*((2.55E-4)/(3.50E-6 + 2.55E-4))
v5r = v5r_max*(1)

metabolic_flux_bounds_array = [
    0.0     v1;       # v1 ATP + L-Citrulline + L-Aspartate --> AMP + Diphosphate + 2-(Nomega-L-Arginino)succinate (units [mmol/gDW-hr])
    0.0     v2;       # v2 2-(Nomega-L-Arginino)succinate --> Fumarate + L-Arginine (units [mmol/gDW-hr])
    0.0     v3;       # v3 L-Arginine + H2O --> L-Ornithine + Urea (units [mmol/gDW-hr])
    0.0     v4;       # v4 Carbamoyl phosphate + L-Ornithine --> phosphate + L-Citrulline (units [mmol/gDW-hr])
    0.0     v5f;      # v5f 2.0*L-Arginine + 4.0*O2 + 3.0*NADPH + 3.0*H --> 2.0*Nitric oxide + 2.0*L-Citrulline + 3.0*NADP + 4.0*H2O (units [mmol/gDW-hr])
    -v5r   0.0 ;      # v5r 2.0*Nitric oxide + 2.0*L-Citrulline + 3.0*NADP + 4.0*H2O --> 2.0*L-Arginine + 4.0*O2 + 3.0*NADPH + 3.0*H (units [mmol/gDW-hr])
    0.0     10.0;     # b1 [] -> Carbamoyl phosphate (units [mmol/gDW-hr])
    0.0     10.0;     # b2 [] -> L-Aspartate (units [mmol/gDW-hr])
    0.0     10.0;     # b3 Fumarate -> [] (units [mmol/gDW-hr])
    0.0     10.0;     # b4 Urea -> [] (units [mmol/gDW-hr])
    0.0     10.0;     # b5 [] -> ATP (units [mmol/gDW-hr])
    0.0     10.0;     # b6 AMP -> [] (units [mmol/gDW-hr])
    0.0     10.0;     # b7 Diphosphate -> [] (units [mmol/gDW-hr])
    0.0     10.0;     # b8 Phosphate -> [] (units [mmol/gDW-hr])
    -10.0     10.0;   # b9 [] -> NADPH (units [mmol/gDW-hr])
    -10.0     10.0;   # b10 [] -> H (units [mmol/gDW-hr])
    -10.0     10.0;   # b11 [] -> O2 (units [mmol/gDW-hr])
    -10.0     10.0;   # b12 Nitric oxide -> [] (units [mmol/gDW-hr])
    -10.0     10.0;   # b13 NADP -> [] (units [mmol/gDW-hr])
    -10.0     10.0;   # b14 [] -> H20 (units [mmol/gDW-hr])
    -10.0     10.0;   # b15 H20 -> [] (units [mmol/gDW-hr])
    ];

flux_array=Variable(21)

b = s_matrix*flux_array

#Objective:  Maximize Urea flux (b4) (flux_array[10])
problem = maximize(flux_array[10])

#Constraints :
#Fluxes must be within bounds
#stoichiometric_matrix * Flux array = b = 0

for i in 1:length(flux_array)
    problem.constraints += flux_array[i] <= metabolic_flux_bounds_array[i,2]
    problem.constraints += flux_array[i] >= metabolic_flux_bounds_array[i,1]
end

for i in 1:length(b)
    problem.constraints += b[i] == 0
end

solve!(problem,SCSSolver())

max_urea_flux = problem.optval

println("Maximum Urea Flux is ",round(max_urea_flux,digits=8)," mmol/gDW-hr")
println("Flux Array is", flux_array)
