using LinearAlgebra
using Test
using BenchmarkTools
using Statistics
using CairoMakie
using GLM, StatsBase, DataFrames, CSV

@inline function norm_coord(v, b, bf = float(b))
    v_1 = 0.0
    for i in eachindex(v)
        v_1 += v[i] * bf^(-i)
    end
    return v_1
end

@test norm_coord([1 1 0 1],2) == 13/16

function red_mat(C,w)
    CT = deepcopy(C)
    for i in eachindex(w)
        for j in 1:w[i]
            CT[i][:,end-j+1] .= 0  
        end  
    end
    return CT
end

# CT[i][:,end-w[i]+1 : end] .= 0 replace this instead of the j for loop

function compute_P_j(P, A, w, st, τ)
(;C,b,m,s) = P
    P_j = zeros(1,τ)
    Cn = zeros(Int64, m)
    #X_j = Float64[]; sizehint!(X_j,b^(m - minimum(w))) #zeros(Float64,b^(m - minimum(w))) this gave error
    for j in st:-1:1
        # computing X_j
        X_j = zeros(Float64,b^(m-w[j])) # Array{Float64}(undef,b^(m-w[j]))
        for i in 1:b^(m-w[j])
            Cn .= (C[j]*(digits(i-1, base=b, pad=m))) .% b  # Can improve this
            X_j[i] = norm_coord(Cn,b)
        end
        @views q_j = X_j*A[j:j,:]

        # computing P_j
        n_w = min( get(w,j+1,m) , m) - w[j]
        P_j = repeat(P_j,b^n_w) + q_j      #Can improve this
        # @show size(P_j) n_w q_j
        
    end
    return P_j
end

# pts needs to be of size (s, N) where N is length(badic) = b^m
function get_points!(pts, P) 
    (;C,b,m,s) = P  
    Cn = zeros(Int64, m) 
    @inbounds for j in eachindex(pts)  # 1:N
        n = digits(j-1, base=b, pad=m)

        for i in eachindex(C)
            mul!(Cn, C[i], n)
            for k in eachindex(Cn)
                Cn[k] = Cn[k] % b  # @. Cn = Cn %b
            end 
            pts[j][i] = norm_coord(Cn, b)
        end
    end
end

function get_points(P) 
    (;C,b,m,s) = P 
    pts = [zeros(s) for i in 1:b^m]
    get_points!(pts, P)
    return pts
end

function mat_mul_Pj(P, A, w)
    CT = red_mat(P.C,w)
    PT = (;P...,C = CT)
    X = stack(get_points(PT))'
    return X*A
end

@testset "column reduced computation" begin 
    P = (C = ( [1 0; 0 1], [0 1; 1 0] ),
        b = 2, m = 2,s = 2)
    A = [1 2 3; 2 5 1]
    w = [0,1]
    st = 2
    τ = 3
    @time prod_alg = compute_P_j(P, A, w, st, τ)
    CT = red_mat(P.C,w)
    PT = (;P...,C = CT)
    X = stack(get_points(PT))'

    @test prod_alg == X*A
end

@testset "get points" begin
    C = ( [1 0; 0 1], [0 1; 1 0] )
    b = 2
    m = 2
    s = 2
    P = (;C,b,m,s)
    pts = get_points(P)
    @test pts == [[0.0, 0.0], [0.5, 0.25],[0.25; 0.5], [0.75, 0.75]]
    
end







   # C = ([2 1 1; 1 1 0; 1 0 0], [0 1 1; 1 1 0; 1 0 0])
   # b = 3
   # m = 3
   # s = 2


function runtime_theory(τ, b, m, s, w_s = @. min(floor(Int64,log2(1:s)),m))
    st = findlast(w_s.< m)
    run_theory = 0.0
    for j in 1:st
        run_theory += τ*b^(m-w_s[j])
    end
    return run_theory
end 

function runtime_theory_rr(τ, b, m, s, w_s = @. min(floor(Int64,log2(1:s)),m))
    st = findlast(w_s.< m)
    run_theory_rr = 0.0
    for j in 1:st
        run_theory_rr += τ*b^(m-w_s[j]) + b^m
    end
    return run_theory_rr
end 

function regres_comp(s_val,T_val)
    df = DataFrame(x = log.(s_val),y = log.(T_val))
    ols = lm(@formula(y ~ x), df)
    return exp(coef(ols)[1]), coef(ols)[2]
end

function regres_theory(T_val, T_theory)
    df = DataFrame(x = T_val,y = T_theory)
    ols = lm(@formula(x ~ y), df)
    return coef(ols)[1], coef(ols)[2]
end

function plot_lines!(s_val,T_val,label)
    lines!(s_val,T_val; linestyle = :solid, label , linewidth = 2)
    c_1,c_2 = regres_comp(s_val,T_val)
    lines!(s_val,c_1.*(s_val.^c_2), color = "light gray", linestyle = :dot)
end
    

#################################################################################

function reduced_mv_product_qmc(b, m, z, w, A)
    N = b^m;         #total number of cubature points N
    M = @. b^(m-w);    #number of reduced points per dimension
    s = length(z);   #number of dimensions s
    #push!(w,m) # extend reduction indices artificially to compute differences

    P = zeros(1, size(A,2));
    for j = s:-1:1
         Xj = @. mod((0:M[j]-1)*z[j], N)/N;
         #@show size(Xj') size(A[j:j,:])
         P = repeat(P, Int64(b^(get(w,j+1,m)-w[j]))) + Xj * A[j:j,:];
    end
    return P
end


####################################  ROW REDUCED MATRIX PRODUCT #############################################

function row_red_mat(C,w)
    CT = deepcopy(C)
    for i in eachindex(w)
        for j in 1:w[i]
            CT[i][end-j+1,:] .= 0  
        end  
    end
    return CT
end


function row_red_prod(P, A, w, st, τ, pts)
    (;C,b,m,s) = P
    P_j = zeros(b^m,τ)
    c = [zeros(1,τ) for i in 1:b^(m- minimum(w))]
    #c = Vector{Matrix{Float64}}(undef,b^(m- minimum(w))) 
    @inbounds for j in st:-1:1
        #Computing row vectors
        for k in 0:b^(m-w[j])-1
           @views @. c[k+1] = k/(b^(m-w[j]))*A[j:j,:] # Julia has a problem with index
        end
        #Compute P_j
        for i in 1:b^m
            k_i = floor(Int,pts[i][j]*b^(m-w[j]))
            @views @. P_j[i:i,:] = P_j[i:i,:] + c[k_i + 1]
        end
    end
    return P_j
end

@testset "row reduced computation" begin 
    P = (C = ( [1 0; 0 1], [0 1; 1 0] ),
        b = 2, m = 2,s = 2)
    pts = get_points(P)
    A = [1 2 3; 2 5 1]
    w = [0,1]
    st = 2
    τ = 3
    @time prod_alg = row_red_prod(P, A, w, st, τ, pts)
    CT = row_red_mat(P.C,w)
    PT = (;P...,C = CT)
    X = stack(get_points(PT))'

    @test prod_alg == X*A
end

# begin
#     P = (C = ( [1 0; 0 1], [0 1; 1 0] ),
#         b = 2, m = 12,s = 2)
#     τ = 20
#     s = 2000
#     m = 12
#     C = Matrix{Int64}[]
#     for i in 1:s
#         C_i = rand(0:1,m,m)
#         push!(C,C_i)
#     end
#     P = (;P...,C,s,m)
#     A = rand(s,τ) 
#     w = @. min(floor(Int64,log2(1:s)),m)
#     st = findlast(w.< m)
#     pts = get_points(P)
#     pts = SizedVector{s}.(pts)
#     @time prod_alg = row_red_prod(P, A, w, st, τ, pts);
#     @profview prod_alg = row_red_prod(P, A, w, st, τ, pts)
#     #CT = row_red_mat(P.C,w)
#     #PT = (;P...,C = CT)
#     #X = stack(get_points(PT))'
# end

# begin 
#     P = (C = ( [1 0; 0 1], [0 1; 1 0] ),
#         b = 2, m = 2,s = 2)
#     τ = 20
#     s = 2000
#     m = 12
#     C = Matrix{Int64}[]
#     for i in 1:s
#         C_i = rand(0:1,m,m)
#         push!(C,C_i)
#     end
#     P = (;P...,C,s,m)
#     A = rand(s,τ) 
#     w = @. min(floor(Int64,log2(1:s)),m)
#     st = findlast(w.< m)
#     @time prod_alg = compute_P_j(P, A, w, st, τ)
#     @profview prod_alg = compute_P_j(P, A, w, st, τ)
#     CT = red_mat(P.C,w)
#     PT = (;P...,C = CT)
#     X = stack(get_points(PT))'

#     @test prod_alg ≈ X*A
# end


# begin
#     s = 2000
#     m = 12
#     C = Matrix{Int64}[]
#     for i in 1:s
#         C_i = rand(0:1,m,m)
#         push!(C,C_i)
#     end
#     P = (;P...,C,s,m)
#     A = rand(s,τ) 
#     w_s = @. min(floor(Int64,log2(1:s)),m)
#     st = findlast(w_s.< m)
# end
# @profview compute_P_j(P, A, w_s, st, τ)

# begin
#     s = 1
#     Ps = (;P...,s)
#     A_s = rand(s,τ) 
#     w_s = @. min(floor(Int64,log2(1:s)),m)
#     st = findlast(w_s.< m)
#     C_sobol = load_seq_mat("Data/sobol_Cs.txt",P.b,m,s)
#     Ps_sobol = (;Ps...,C=C_sobol,s)
#     #T_sobol_test = @belapsed compute_P_j($Ps_sobol, $A_s, $w_s, $st, $τ)
#     T_2_sobol_test = @belapsed mat_mul_Pj($Ps_sobol, $A_s, $w_s)

# end