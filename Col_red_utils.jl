using LinearAlgebra
using Test

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
    for j in st:-1:1
        # computing X_j
        X_j = zeros(Float64,b^(m-w[j])) # Array{Float64}(undef,b^(m-w[j]))
        for i in 1:b^(m-w[j])
            Cn = (C[j]*(digits(i-1, base=b, pad=m))) .% b
            X_j[i] = norm_coord(Cn,b)
        end
        q_j = X_j*A[j:j,:]

        # computing P_j
        n_w = min( get(w,j+1,m) , m) - w[j]
        P_j = repeat(P_j,b^n_w) + q_j
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

@testset "column reduced computation" begin 
    P = (C = ( [1 0; 0 1], [0 1; 1 0] ),
        b = 2, m = 2,s = 2)
    A = [1 2 3; 2 5 1]
    w = [0,1]
    st = 2
    τ = 3
    prod_alg = compute_P_j(P, A, w, st, τ)
    
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
    
