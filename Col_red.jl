include("Col_red_utils.jl")
using CairoMakie



b=2
m=12
τ = 20

P = (C = ( [1 0; 0 1], [0 1; 1 0] ),
        b = 2, m = 12,s = 2)
A = [1 2 3; 2 5 1]
w = [0,1]
st = 2
τ = 3
@time prod_alg = compute_P_j(P, A, w, st, τ)
CT = red_mat(P.C,w)
PT = (;P...,C = CT)
X = stack(get_points(PT))'

s_val = collect(1:200:1500)
T_val = Float64[]
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

    @time compute_P_j(Ps, A_s, w_s, st, τ)
    T = @elapsed compute_P_j(Ps, A_s, w_s, st, τ)
    push!(T_val,T)
end


fig = Figure()
ax = Axis(fig[1,1])
lines!(s_val,T_val)
fig