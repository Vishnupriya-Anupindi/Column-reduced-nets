using Test

@inline function norm_coord(v, b, bf = float(b))
    v_1 = 0.0
    for i in eachindex(v)
        v_1 += v[i] * bf^(-i)
    end
    return v_1
end

@test norm_coord([1 1 0 1],2) == 13/16

function compute_P_j(pts, A, w, b, m, s, st, τ)
    P_j = zeros(1,τ)
    for j in st:-1:1
        n_w = min((w[j+1]) , m) - w[j]
        P_j = repeat(P_j,b^n_w) 
        @show(size(P_j))
        @show n_w
    end
    return P_j
end

compute_P_j(0,0,[2,2,3,4],2,5,4,3,5)


# pts needs to be of size (s, N) where N is length(badic) = b^m
function get_points!(pts, C, badic, b, m,bf = float(b))  
    Cn = zeros(Int64, m) 
    @inbounds for j in eachindex(pts)  # 1:N
        n = badic[j]

        for i in eachindex(C)
            mul!(Cn, C[i], n)
            for k in eachindex(Cn)
                Cn[k] = Cn[k] % b
            end 
            pts[j][i] = norm_coord(Cn, b, bf)
        end
    end
end

function get_points(C, badic, b, m, s, bf = float(b))  
    pts = [zeros(s) for i in 1:length(badic)]
    get_points!(pts, C, badic, b, m, bf)
    return pts
end

get_badic(b, m) = collect.(Iterators.product(fill(0:b-1, m)...))[:]

@testset "get points" begin
    C = ( [1 0; 0 1], [0 1; 1 0] )
    b = 2
    m = 2
    s = 2
    badic = get_badic(b,m)
    bf = float(b)
    pts = [zeros(s) for i in 1:length(badic)]
    get_points!(pts, C, badic, b, m, bf)
    @test pts == [[0.0, 0.0], [0.5, 0.25],[0.25; 0.5], [0.75, 0.75]]
    
end

@testset "get points" begin
    C = ([2 1 1; 1 1 0; 1 0 0], [0 1 1; 1 1 0; 1 0 0])
    b = 3
    m = 3
    s = 2
    badic = get_badic(b,m)
    bf = float(b)
    pts = [zeros(s) for i in 1:length(badic)]
    get_points!(pts, C, badic, b, m, bf)
    @test pts == [0 0.25 0.5 0.75; 0 0.5 0.25 0.75]
end


# generates the ith matrix (returned as a vector) 
# in base b with m rows and m columns 
function int_to_matrix!(C, i, b, m) 
    C .= 0
    if i > 0
        digits!(C, i, base=b)
    end
    return C
end 

function int_to_matrix(i, b, m) 
    C = zeros(Int, m*m)
    return int_to_matrix!(C, i, b, m)
end 

@test int_to_matrix(0, 3, 2) == [0, 0, 0, 0]
@test int_to_matrix(5, 3, 2) == [2, 1, 0, 0]
@test int_to_matrix(3^(2*2)-1, 3, 2) == [2, 2, 2, 2]



get_matrix(i, b, m) = reshape(int_to_matrix(i, b, m), m, m)
#get_matrices(i1,i2,b,m) = [ get_matrix(i1,b,m) get_matrix(i2,b,m) ]



function compare_rows_matrix(C1,C2)
    C = (C1,C2)
    if C[1][1,:] == C[2][1,:]
        return true
    end

    if C[1][1,:] == C[2][2,:]
        return true
    end

    if C[1][2,:] == C[2][1,:]
        return true
    end
    return false
end

function compare_rows_matrix(C1,C2,C3)
    C = (C1,C2,C3)
    for d1 in 1:3
        for d2 in 1:3
            if d1<d2 && C[d1][1,:] == C[d2][1,:]
                return true
            end

            #if d1!=d2 && C[d1][1,:] == C[d2][2,:]
            #    return true
            #end
        end
    end
    return false
end



#hshojhsfjkhsfjdlk