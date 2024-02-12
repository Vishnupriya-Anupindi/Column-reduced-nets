include("Col_red_utils.jl")
#include("load_seq_mat.jl")
#include("row_latt_red.jl") #already included

#mkpath("Output")
case = 21   # 1 means plot with random matrices, 2 means plot with sobol and niederreiter
s = 800
step_size = 2
m_range = 20
fn_postfix = "case$(case)_m$(m_range)_s$(s)"

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 0.1
#BenchmarkTools.DEFAULT_PARAMETERS.samples = 2

τ = 20

P = (C = ( [1 0; 0 1], [0 1; 1 0] ),
        b = 2, m=2, s)

begin 
    m_val = collect(10:step_size:m_range)
    T_val = Float64[]
    T_val2 = Float64[]
    T_val_r = Float64[]
    T_val_gp = Float64[]
    T_val3 = Float64[]

    T_val_sobol = Float64[]
    T_val2_sobol = Float64[]

    T_val_nied = Float64[]
    T_val2_nied = Float64[]

    (;s) = P
    A_s = rand(s,τ) 
    
    for m in m_val

        C = Matrix{Int64}[]
        for i in 1:s
            C_i = rand(0:1,m,m)
            push!(C,C_i)
        end

        Pm = (;P...,C,m)

        w_s = @. min(floor(Int64,log2(1:s)),m)

        st = findlast(w_s.< m)

        # @time compute_P_j(Pm, A_s, w_s, st, τ)
        T = @belapsed compute_P_j($Pm, $A_s, $w_s, $st, $τ)
        push!(T_val,T)

        T_2 = @belapsed mat_mul_Pj($Pm, $A_s, $w_s)
        push!(T_val2,T_2)

        local z_r = rand(1:1000,s);
        T_3 = @belapsed reduced_mv_product_qmc($P.b, $m, $z_r, $w_s, $A_s);
        push!(T_val3,T_3)

        pts = get_points(Pm)
        T_gp = @belapsed get_points($Pm)
        T_r = @belapsed row_red_prod($Pm, $A_s, $w_s, $st, $τ, $pts)
        push!(T_val_gp, T_gp)
        push!(T_val_r, T_r)
        

        if case == 2
            
            C_sobol = load_seq_mat("Data/sobol_Cs.txt",P.b,m,s)
            Pm_sobol = (;Pm...,C=C_sobol,s)
            T_sobol = @belapsed compute_P_j($Pm_sobol, $A_s, $w_s, $st, $τ)
            push!(T_val_sobol,T_sobol)
            @show s w_s
            #@show Pm_sobol.C
            #@show s

            T_2_sobol = @belapsed mat_mul_Pj($Pm_sobol, $A_s, $w_s)
            push!(T_val2_sobol,T_2_sobol)

            C_nied = load_seq_mat("Data/niederreiter_Cs.txt",P.b,m,s)
            Pm_nied = (;Pm...,C=C_nied,s)
            T_nied = @belapsed compute_P_j($Pm_nied, $A_s, $w_s, $st, $τ)
            push!(T_val_nied,T_nied)
            @show s w_s
            #@show Pm_sobol.C
            #@show s

            T_2_nied = @belapsed mat_mul_Pj($Pm_nied, $A_s, $w_s)
            push!(T_val2_nied,T_2_nied)

        end


        println("Finished m=", m)

    end
end

#T_theory = [runtime_theory(τ, P.b, P.m, s) for s in s_val]


df = DataFrame(m = m_val, row_red = T_val_r, col_red = T_val, std_mat = T_val2, lat_red_row = T_val3)
CSV.write("runtime_$(fn_postfix)_b$(P.b).csv", df)

begin
    fig = Figure()
    ax = Axis(fig[1,1], title = "log_plot", xlabel = "log m", ylabel = "Runtime (log seconds)",xscale = log10, yscale = log10)
    
    plot_lines!(m_val,T_val,"col_red")
    plot_lines!(m_val,T_val2,"std_mul")
    plot_lines!(m_val,T_val_r,"row_red")
    plot_lines!(m_val,T_val3,"row_latt_red")

    #d_1,d_2 =  regres_theory(T_val, T_theory)
    #lines!(m_val, d_2.*T_theory,linestyle = :dash, label="theoretical",linewidth = 2)
    
    if case ==2
        plot_lines!(m_val,T_val_sobol,"Sobol_col_red")
        plot_lines!(m_val,T_val2_sobol,"Sobol_std_mul")

        plot_lines!(m_val,T_val_nied,"Nied_col_red")
        plot_lines!(m_val,T_val2_nied,"Nied_std_mul")
    end

    Legend(fig[1,2], ax)
    save("Output/logplot_$(fn_postfix).png", fig)
    fig
end


begin
    fig = Figure()
    ax = Axis(fig[1,1], title = "plot", xlabel = "m", ylabel = "Runtime (log seconds)", yscale = log10)

    plot_lines!(m_val,T_val,"col_red")
    plot_lines!(m_val,T_val2,"std_mul")
    plot_lines!(m_val,T_val_r,"row_red")
    plot_lines!(m_val,T_val3,"row_latt_red")
    
    #d_1,d_2 =  regres_theory(T_val, T_theory)
    #lines!(m_val, d_2.*T_theory,linestyle = :dash, label="theoretical",linewidth = 2)
    
    Legend(fig[1,2], ax)
    save("Output/plot_$(fn_postfix).png", fig)
    fig
end
