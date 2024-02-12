include("Col_red_utils.jl")
#include("load_seq_mat.jl")
#include("row_latt_red.jl") #already included

#mkpath("Output")
case = 31   # 1 means plot with random matrices, 2 means plot with sobol and niederreiter
tau_range = 1000
step_size = 60
m = 12
s = 800
fn_postfix = "case$(case)_m$(m)_s$(s)_tau_$(tau_range)"

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 0.1
#BenchmarkTools.DEFAULT_PARAMETERS.samples = 2

τ = 20

P = (C = ( [1 0; 0 1], [0 1; 1 0] ),
        b = 2, m, s)


begin 
    τ_val = collect(1:step_size:tau_range)
    T_val = Float64[]
    T_val2 = Float64[]
    T_val_r = Float64[]
    T_val_gp = Float64[]
    T_val3 = Float64[]

    T_val_sobol = Float64[]
    T_val2_sobol = Float64[]

    T_val_nied = Float64[]
    T_val2_nied = Float64[]

    C = Matrix{Int64}[]
    for i in 1:s
        C_i = rand(0:1,m,m)
        push!(C,C_i)
    end
    Ps = (;P...,C)
    w_s = @. min(floor(Int64,log2(1:s)),m)
    st = findlast(w_s.< m)

    #(;m) = P
    for τ in τ_val

        A_s = rand(s,τ) 
        
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

        println("Finished tau=", τ)

    end
end

#T_theory = [runtime_theory(τ, P.b, P.m, s) for s in τ_val]


df = DataFrame(τ = τ_val, row_red = T_val_r, col_red = T_val, std_mat = T_val2, lat_red_row = T_val3)
CSV.write("runtime_$(fn_postfix)_b$(P.b).csv", df)

begin
    fig = Figure()
    ax = Axis(fig[1,1], title = "log_plot", xlabel = "log τ", ylabel = "Runtime (log seconds)",xscale = log10, yscale = log10)
    
    plot_lines!(τ_val,T_val,"col_red")
    plot_lines!(τ_val,T_val2,"std_mul")
    plot_lines!(τ_val,T_val_r,"row_red")
    plot_lines!(τ_val,T_val3,"row_latt_red")

    Legend(fig[1,2], ax)
    save("Output/logplot_$(fn_postfix).png", fig)
    fig
end


begin
    fig = Figure()
    ax = Axis(fig[1,1], title = "plot", xlabel = "τ", ylabel = "Runtime (seconds)")

    plot_lines!(τ_val,T_val,"col_red")
    plot_lines!(τ_val,T_val2,"std_mul")
    plot_lines!(τ_val,T_val_r,"row_red")
    plot_lines!(τ_val,T_val3,"row_latt_red")
    
    
    Legend(fig[1,2], ax)
    save("Output/plot_$(fn_postfix).png", fig)
    fig
end
