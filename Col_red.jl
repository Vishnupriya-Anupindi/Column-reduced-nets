using Random, DataFrames, CSV, LinearAlgebra, ProgressMeter
include("Col_red_utils.jl")

mkpath("data")
fn_postfix = "i1_14408200_less_check"

#@profview let 
begin 
    b = 3
    m = 4
    s = 3
    ρ = m 
    t = m - ρ 
    N = b^m

    bf = float(b)
    c_z = 500
    pts = [zeros(s) for i in 1:N ]
    badic = get_badic(b,m)

    println("b = $b, m = $m, s = $s, # matrices = $(b^(m*m)), # zero row matrices = $(b^(m*(m-1)))" )
end

begin
    matrix_range = 0:b^(m*m)-1
    @assert length(matrix_range) <= b^(m*m) "matrix range large"
    @assert last(matrix_range) >= b^(m*(m-1)) "zero rows"

    if s==2
        C = (zeros(Int,m,m),zeros(Int,m,m))
    
        df = DataFrame(i1 = Int64[], i2 = Int64[])
        CSV.write("data/output_b$(b)_m$(m)_s$(s)$(fn_postfix).csv ", df)
    
        prog = Progress(Int(length(matrix_range)*(length(matrix_range)-1)/2))

        for i1 in matrix_range

            # write in the beginning, to ensure that we don't accidently shift the writing part
            if i1 % 100 == 0 
                CSV.write("data/output_b$(b)_m$(m)_s$(s)$(fn_postfix).csv ", df, append = true)
                empty!(df)
            end

            C[1] .= reshape(int_to_matrix(i1, b, m), m, m)

            if det(C[1]) % b == 0
                next!(prog,step=i1)
                continue
            end

            for i2 in 0:i1-1
                next!(prog)
                
                C[2] .= reshape(int_to_matrix(i2, b, m), m, m)

                if det(C[2]) % b == 0
                    continue
                end

                if compare_rows_matrix(C...)
                    continue
                end
                    
                get_points!(pts, C, badic, b, m, bf)  
                
                if is_NNLD(c_z,s,pts)
                    push!(df, (i1, i2))
                end
            end
        end
        CSV.write("data/output_b$(b)_m$(m)_s$(s)$(fn_postfix).csv ", df, append = true)


    elseif s==3
        C = (zeros(Int,m,m),zeros(Int,m,m),zeros(Int,m,m))
    
        df = DataFrame(i1 = Int64[], i2 = Int64[], i3 = Int64[])
        CSV.write("data/output_b$(b)_m$(m)_s$(s)$(fn_postfix).csv ", df)
        
        

        prog = Progress(binomial(length(matrix_range),2)) #Change this to 3 later.

        for i1 in [14408200] #matrix_range

            # write in the beginning, to ensure that we don't accidently shift the writing part
            #Shifting to i2 since we are keeping i1 fixed.
            #if i1 % 100 == 0 
            #    CSV.write("data/output_b$(b)_m$(m)_s$(s)$(fn_postfix).csv ", df, append = true)
            #    println("write dataframe with i1 = $i1")
            #    empty!(df)
            #end

            C[1] .= reshape(int_to_matrix(i1, b, m), m, m)

            if det(C[1]) % b == 0
                next!(prog,step= binomial(i1,2))
                continue
            end

            for i2 in matrix_range

                #Comment this out when computing for all i1.
                if i2 % 100 == 0 
                    CSV.write("data/output_b$(b)_m$(m)_s$(s)$(fn_postfix).csv ", df, append = true)
                    println("write dataframe with i2 = $i2")
                    empty!(df)
                end
                #End of commenting out
                
                C[2] .= reshape(int_to_matrix(i2, b, m), m, m)

                if det(C[2]) % b == 0
                    next!(prog,step= binomial(i2,1))
                    continue
                end

                for i3 in 0:i2-1
                    next!(prog)
                
                    C[3] .= reshape(int_to_matrix(i3, b, m), m, m)

                    if det(C[3]) % b == 0
                        continue
                    end
                

                    if compare_rows_matrix(C...)
                        continue
                    end
                    
                    get_points!(pts, C, badic, b, m, bf)  
                
                    if is_NNLD(c_z,s,pts)
                        push!(df, (i1, i2, i3))
                    end
                end
            end
        end
        CSV.write("data/output_b$(b)_m$(m)_s$(s)$(fn_postfix).csv ", df, append = true)

    else 
        @warn "not implemented"
    end

end


df_result = CSV.read("data/output_b$(b)_m$(m)_s$(s)$(fn_postfix).csv ", DataFrame)

filter!(row -> validate_NNLD(row,c_z*100,b,m,s), df_result)
filter!(row -> validate_NNLD(row,c_z*70,b,m,s), df_result)
filter!(row -> validate_NNLD(row,c_z*50,b,m,s), df_result)
#CSV.write("filter.csv",df_result)


CSV.write("data/filter_b$(b)_m$(m)_s$(s)$(fn_postfix).csv ",df_result)
df_result = CSV.read("data/filter_b$(b)_m$(m)_s$(s)$(fn_postfix).csv ", DataFrame)

begin
    idxs = values(df_result[11040,:]) # Use Vector() if number of row entries are large
    C = [get_matrix(i,b,m) for i in idxs]
    J = hcat(C...)
    display(J)
    @show idxs
    
    pts = get_points(C, get_badic(b,m), b, m,s)

    @show is_NNLD(10000,s,pts)



    if s == 2
        include("NNLD_plots.jl")

        plot_points(pts)
    end
end

#idxs = values(df_result[1762,:]) 
#C = [get_matrix(i,b,m) for i in idxs]
#K = hcat(C[1])
#for i in 2:s
#    K = hcat(K,C[i])
#end
    