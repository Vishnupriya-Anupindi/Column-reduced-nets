include("Col_red_utils.jl")



function compute_P_j_row_col(P, A, w, st, τ)
    (;C,b,m,s) = P
    P_j = zeros(1,τ)
    Cn = zeros(Int64, m)
    base_n = zeros(Int64, m)
    X_j = zeros(Float64,b^(m - minimum(w))) #Float64[]; sizehint!(X_j,b^(m - minimum(w)))this gave error
    for j in st:-1:1
        # computing X_j
        #X_j = zeros(Float64,b^(m-w[j])) # Array{Float64}(undef,b^(m-w[j]))
        resize!(Cn,m-w[j])
        resize!(base_n, m-w[j])
        C_j_red = C[j][1:m-w[j],1:m-w[j]]
        for i in 1:b^(m-w[j])
            digits!(base_n, i-1, base=b)
            mul!(Cn,C_j_red,base_n)
            @. Cn = Cn % b  # Can improve this
            X_j[i] = norm_coord(Cn,b)
        end
        @views q_j = X_j[1:b^(m-w[j])]*A[j:j,:] #affected by row reduction

        # computing P_j
        n_w = min( get(w,j+1,m) , m) - w[j]
        P_j = repeat(P_j,b^n_w) + q_j      #Can improve this
        # @show size(P_j) n_w q_j
        
    end
    return P_j
end

red_mat(C, u, w) = map( i -> C[i][1:end-u[i], 1:end-w[i]], eachindex(w))


@testset "column reduced computation" begin 
    P = (C = ( [1 0; 0 1], [0 1; 1 0] ),
        b = 2, m = 2,s = 2)
    A = [1 2 3; 2 5 1]
    w = [0,1]
    st = 2
    τ = 3
    @time prod_alg = compute_P_j_row_col(P, A, w, st, τ)
    # CT = red_mat(P.C,w,w)
    # PT = (;P...,C = CT)
    # X = stack(get_points(PT))'

    # @test prod_alg == X*A 
end

begin
    b = 2;  m = 12;  s = 500 ; τ = 20
    C = [rand(0:1,m,m) for i in 1:s]
    P = (;C , b, m , s )
    A = rand(s,τ)
    w = @. min(floor(Int64,log2(1:s)),m)
    st = findlast(w.< m)
    compute_P_j_row_col(P, A, w, st, τ)
    @time prod_alg = compute_P_j_row_col(P, A, w, st, τ)
    # test_1  0.010893 seconds (1.70 k allocations: 16.995 MiB)
   # test_2 naive row_col: 0.024204 seconds (38.40 k allocations: 21.356 MiB, 36.77% gc time)
   #test 3 0.019666 seconds (38.90 k allocations: 21.489 MiB, 44.39% gc time)
   # test 4: with mul!: 0.016291 seconds (2.20 k allocations: 17.129 MiB, 53.02% gc time)
    #@profview prod_alg = compute_P_j_row_col(P, A, w, st, τ)
end

#mkpath("Output")
case = 3   # 1 means plot with random matrices, 2 means plot with sobol and niederreiter
s_range = 1600
step_size = 200
m = 12
fn_postfix = "case$(case)_m$(m)_s$(s_range)"

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 0.1
#BenchmarkTools.DEFAULT_PARAMETERS.samples = 2

τ = 20

P = (C = ( [1 0; 0 1], [0 1; 1 0] ),
        b = 2, m, s = 2)


begin 
    s_val = collect(1:step_size:s_range)
    T_val = Float64[]
    T_val2 = Float64[]
    T_val_r = Float64[]
    T_val_gp = Float64[]
    T_val3 = Float64[]

    T_val_sobol = Float64[]
    T_val2_sobol = Float64[]

    T_val_nied = Float64[]
    T_val2_nied = Float64[]

    (;m) = P
    for s in s_val

        C = Matrix{Int64}[]
        for i in 1:s
            C_i = rand(0:1,m,m)
            push!(C,C_i)
        end

        Ps = (;P...,C,s)
        A_s = rand(s,τ) 
        
        #w_s = @. min(floor(Int64,log2(1:s)),m)

        if case == 11
            w_s = @. min(floor(Int64,log2(1:s)),m)
        elseif case == 12
            w_s = @. min(floor(Int64,(log2(1:s))^(1/2)),m)
        elseif case == 13
            w_s = @. min(floor(Int64,(log2(1:s))^(1/4)),m)
        end

        st = findlast(w_s.< m)

        # @time compute_P_j(Ps, A_s, w_s, st, τ)
        T = @belapsed compute_P_j($Ps, $A_s, $w_s, $st, $τ)
        push!(T_val,T)

        T_2 = @belapsed mat_mul_Pj($Ps, $A_s, $w_s)
        push!(T_val2,T_2)

        local z_r = rand(1:1000,s);
        T_3 = @belapsed reduced_mv_product_qmc($P.b, $m, $z_r, $w_s, $A_s);
        push!(T_val3,T_3)

        pts = get_points(Ps)
        T_gp = @belapsed get_points($Ps)
        T_r = @belapsed row_red_prod($Ps, $A_s, $w_s, $st, $τ, $pts)
        push!(T_val_gp, T_gp)
        push!(T_val_r, T_r)
    

        println("Finished s=", s)

    end
end

T_theory = [runtime_theory(τ, P.b, P.m, s) for s in s_val]


df = DataFrame(s = s_val, row_red = T_val_r, col_red = T_val, std_mat = T_val2, lat_red_row = T_val3, theo_col = T_theory)
CSV.write("runtime_$(fn_postfix)_b$(P.b).csv", df)