#Submitted by Siddharth Krishnan - sk3275
#CHEME7770 - HW5
#PRoblem 2 Phase Plots

include("PhasePortraitV2.jl")

# Function for a dual repression system without cooperativity
# x1: range of x1 values (i.e. A values)
# x2: range of x2 values (i.e. R values)
# We use `@.` to apply the calculations across all rows.
# Note that model parameters are specified within the function
# Returns computed (dx1/dt, dx2/dt) over the range of (x1, x2)
function toggleMono(x1, x2)
    d1 = 30            #degradation rate const. for repressor 1
                       
    r10 = 100             #max generation rate of repressor 1
    r20= 1             #max generation rate of repressor 2
    r2 = 100            #Association constant for repressor 1 - promot. 2
    r1 = 5000             #Association constant for repressor 2 - promot. 1

    u = @. -d1*x1 + (r10+r1*(x1^2))/(1+x1^2+(x2^2)) #eqn 11
    v = @. -x2 + (r20+r2*(x1^2))/(1+(x1^2))   #eqn 12

    return (u,v)
end

#Range of x1, x2 values
x1range = (0,200.0,30)          #Has the form (min, max, num points)
x2range = (0,100.0,30)          #Has the form (min, max, num points)
x₀ = ([1.0,0.0],[0.0,10.0])  #initial state vectors; a common must be included after the first array
tspan=(0.0,30.0)             #time span

#Call the phaseplot functon to construct the phase portrait
phaseplot(toggleMono, x1range, x2range, xinit=x₀, t=tspan, clines=true,
        norm=true, scale=0.5)
