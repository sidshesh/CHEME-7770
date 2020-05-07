#Submitted by Siddharth Krishnan - sk3275
#CHEME 7770 - HW5
#Problem 3 Phase Plots

include("PhasePortraitV2.jl")

# Function for a dual repression system without cooperativity
# x1: range of x1 values (i.e. A values)
# x2: range of x2 values (i.e. R values)
# We use `@.` to apply the calculations across all rows.
# Note that model parameters are specified within the function
# Returns computed (dx1/dt, dx2/dt) over the range of (x1, x2)
function toggleMono(x1, x2)
    a=10.0
    n=1
    u = @. -x1 + a/(1+x2^n) #eqn 11
    v = @. -x2 + a/(1+x1^n)   #eqn 12

    return (u,v)
end

#Range of x1, x2 values
x1range = (0,20,15)          #Has the form (min, max, num points)
x2range = (0,20,15)          #Has the form (min, max, num points)
x₀ = ([0.9,0.35],[0.1, 0.0])  #initial state vectors; a common must be included after the first array
tspan=(0.0,30.0)             #time span

#Call the phaseplot functon to construct the phase portrait
phaseplot(toggleMono, x1range, x2range, xinit=(), t=tspan, clines=true,
        norm=true, scale=0.5)
