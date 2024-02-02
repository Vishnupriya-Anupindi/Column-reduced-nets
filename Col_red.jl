include("Col_red_utils.jl")
include("row_latt_red.jl")

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 0.1

τ = 20

P = (C = ( [1 0; 0 1], [0 1; 1 0] ),
        b = 2, m = 12,s = 2)


begin
    s = 2000
    m = 12
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
# @profview compute_P_j(P, A, w_s, st, τ)

begin 
    s_val = collect(1:200:1500)
    T_val = Float64[]
    T_val2 = Float64[]
    T_val_r = Float64[]
    T_val_gp = Float64[]
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

        local z_r = rand(1:1000,s);
        T_3 = @belapsed reduced_mv_product_qmc($b, $m, $z_r, $w_s, $A_s);
        push!(T_val3,T_3)

        pts = get_points(Ps)
        T_gp = @belapsed get_points($Ps)
        T_r = @belapsed row_red_prod($Ps, $A_s, $w_s, $st, $τ, $pts)
        push!(T_val_gp, T_gp)
        push!(T_val_r, T_r)
        println("Finished s=", s)
    end
end

df = DataFrame(s = s_val, row_red = T_val_r, col_red = T_val, std_mat = T_val2, lat_red_row = T_val3 )
CSV.write("runtime_b$(b)_m$(m)_s$(s).csv", df)

function runtime_theory(τ, b, m, s, w_s = @. min(floor(Int64,log2(1:s)),m))
    st = findlast(w_s.< m)
    run_theory = 0.0
    for j in 1:st
        run_theory += τ*b^(m-w_s[j])
    end
    return run_theory
end 

T_theory = [runtime_theory(τ, P.b, P.m, s) for s in s_val]


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

begin
    fig = Figure()
    ax = Axis(fig[1,1], title = "log_plot", xlabel = "log s", ylabel = "Runtime (log seconds)",xscale = log10, yscale = log10)
    lines!(s_val,T_val, label="col_red", linewidth = 2)
    c_1,c_2 = regres_comp(s_val,T_val)
    lines!(s_val,c_1.*(s_val.^c_2), color = "light gray", linestyle = :dot)

    lines!(s_val,T_val2,linestyle = :solid, label="std_mul",linewidth = 2)
    c_3,c_4 = regres_comp(s_val,T_val2)
    lines!(s_val,c_3.*(s_val.^c_4), color = "light gray", linestyle = :dot)

    lines!(s_val,T_val_r, label="row_red", linewidth = 2)
    c_7,c_8 = regres_comp(s_val,T_val_r)
    lines!(s_val,c_7.*(s_val.^c_8), color = "light gray", linestyle = :dot)

    #lines!(s_val,T_val3,linestyle = :solid, label="row_latt_red",linewidth = 2)
    #c_5,c_6 = regres_comp(s_val,T_val3)
    #lines!(s_val,c_5.*(s_val.^c_6), color = "light gray", linestyle = :dot)

    d_1,d_2 =  regres_theory(T_val, T_theory)
    lines!(s_val, d_2.*T_theory,linestyle = :dash, label="theoretical",linewidth = 2)
    
    Legend(fig[1,2], ax)
    save("log_plot_comparision_v1.png", fig)
    fig
end

begin
    fig = Figure()
    ax = Axis(fig[1,1], title = "plot", xlabel = "s", ylabel = "Runtime (seconds)")
    lines!(s_val,T_val, label="col_red", linewidth = 2)
    c_1,c_2 = regres_comp(s_val,T_val)
    lines!(s_val,c_1.*(s_val.^c_2), color = "light gray", linestyle = :dot)

    lines!(s_val,T_val2,linestyle = :solid, label="std_mul",linewidth = 2)
    c_3,c_4 = regres_comp(s_val,T_val2)
    lines!(s_val,c_3.*(s_val.^c_4), color = "light gray", linestyle = :dot)

    lines!(s_val,T_val_r, label="row_red", linewidth = 2)
    c_7,c_8 = regres_comp(s_val,T_val_r)
    lines!(s_val,c_7.*(s_val.^c_8), color = "light gray", linestyle = :dot)

    lines!(s_val,T_val3,linestyle = :solid, label="row_latt_red",linewidth = 2)
    c_5,c_6 = regres_comp(s_val,T_val3)
    lines!(s_val,c_5.*(s_val.^c_6), color = "light gray", linestyle = :dot)

    d_1,d_2 =  regres_theory(T_val, T_theory)
    lines!(s_val, d_2.*T_theory,linestyle = :dash, label="theoretical",linewidth = 2)
    
    Legend(fig[1,2], ax)
    save("theory_plot_comparision_v1.png", fig)
    fig
end
