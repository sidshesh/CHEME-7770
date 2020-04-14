include("Include.jl")
s_matrix = DataFrame(CSV.File("S_Matrix.csv",header=false))
s_matrix = convert(Matrix,s_matrix)

atom_matrix = DataFrame(CSV.File("Atom_Matrix.csv",header=false))
atom_matrix = convert(Matrix,atom_matrix)

checkbalance = transpose(atom_matrix)*s_matrix

print(checkbalance)
#The first six columns of the checkbalance matrix are all zeros
#that is, the elements are all balanced within the cycle
