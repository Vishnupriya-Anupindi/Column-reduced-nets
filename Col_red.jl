include("Col_red_utils.jl")

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 0.1





τ = 20

P = (C = ( [1 0; 0 1], [0 1; 1 0] ),
        b = 2, m = 12,s = 2)


begin
    s = 5000
    m = 15
    C = Matrix{Int64}[]
    for i in 1:s
        C_i = rand(0:1,m,m)
        push!(C,C_i)
    end
    P = (;P...,C,s,m)
    A = rand(s,τ) 
    w_s = @. min(floor(Int64,log2(1:s)),m)
    st = findlast(w_s.< m)
end
@profview compute_P_j(P, A, w_s, st, τ)

begin 
    s_val = collect(1:200:1500)
    T_val = Float64[]
    T_val2 = Float64[]
    (;m) = P
    for s in s_val

        C = Matrix{Int64}[]
        for i in 1:s
            C_i = rand(0:1,m,m)
            push!(C,C_i)
        end
        Ps = (;P...,C,s)
        A_s = rand(s,τ) 
        w_s = @. min(floor(Int64,log2(1:s)),m)
        st = findlast(w_s.< m)

        # @time compute_P_j(Ps, A_s, w_s, st, τ)
        T = @belapsed compute_P_j($Ps, $A_s, $w_s, $st, $τ)
        push!(T_val,T)

        T_2 = @belapsed mat_mul_Pj($Ps, $A_s, $w_s)
        push!(T_val2,T_2)
    end
end

begin
    fig = Figure()
    ax = Axis(fig[1,1])
    lines!(s_val,T_val)
    lines!(s_val,T_val2,linestyle = :dash)
    fig
end