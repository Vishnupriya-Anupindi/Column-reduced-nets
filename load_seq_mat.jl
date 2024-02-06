using CSV, DataFrames

function row_to_mat(row,b,m)
    return stack(first(digits(x, base = b, pad = m),m) for x in first(row,m))
end

function load_seq_mat(filename,b,m,s)
    df = CSV.read(filename, DataFrame, header = false, delim = ' ')
    K = Matrix(df[1:s,1:m])
    C = [row_to_mat(K[i,:],b,m) for i in 1:s]
    return C
end

# C = load_seq_mat("Data/sobol_Cs.txt",2,4,2)


###################### SOME ATTEMPTS ############################

#readdlm("Data/sobol_Cs.txt", ' ', Int128, '\n')

# df = CSV.read("Data/sobol_Cs.txt", DataFrame, header = false, delim = ' ')
# C = Matrix(df[1:3,1:5])

#row_to_mat(C[2,:],2,4)

# function load_seq_mat(filename,b,m,s)
#     df = CSV.read(filename, DataFrame, header = false, delim = ' ')
#     K = Matrix(df[1:s,1:m])
#     C = []
#     for i in 1:s
#         push!(C,row_to_mat(K[i,:],b,m))
#     end
#     return C
# end
