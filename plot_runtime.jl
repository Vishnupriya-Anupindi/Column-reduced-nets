using BenchmarkTools
using Statistics
using CairoMakie
using GLM, StatsBase, DataFrames, CSV

b = 2
case = 12   # 1 means plot with random matrices, 2 means plot with sobol and niederreiter
s_range = 1600
step_size = 200
m = 12
fn_postfix = "case$(case)_m$(m)_s$(s_range)"

#df = DataFrame(s = df.s, row_red = df.row_red, col_red = df.col_red, std_mat = df.std_mat, lat_red_row = df.lat_red_row, theo_col = df.theo_col)
df = CSV.read("runtime_$(fn_postfix)_b$(b).csv", DataFrame)

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
    scatter!(s_val,T_val; marker = :dot, label , markersize = 10)
    #c_1,c_2 = regres_comp(s_val,T_val)
    #lines!(s_val,c_1.*(s_val.^c_2), color = "light gray", linestyle = :dot)
end



begin
    fig = Figure()
    ax = Axis(fig[1,1], title = "log_plot", xlabel = "log s", ylabel = "Runtime (log seconds)",xscale = log10, yscale = log10, xminorticksvisible = true, xminorgridvisible = true,
    xminorticks = IntervalsBetween(5))
    
    plot_lines!(df.s,df.col_red,"col_red")
    plot_lines!(df.s,df.std_mat,"std_mul")
    plot_lines!(df.s,df.row_red,"row_red")
    plot_lines!(df.s,df.lat_red_row,"row_latt_red")

    d_1,d_2 =  regres_theory(df.col_red, df.theo_col)
    lines!(df.s, d_2.*df.theo_col,linestyle = :dash, label="theoretical",linewidth = 2)
    
    if case ==2
        plot_lines!(df.s,df.col_red_sobol,"Sobol_col_red")
        plot_lines!(df.s,df.std_mat_sobol,"Sobol_std_mul")

        plot_lines!(df.s,df.col_red_nied,"Nied_col_red")
        plot_lines!(df.s,df.std_mat_nied,"Nied_std_mul")
    end

    axislegend(ax, merge = true, position = :lt)
    save("Output/logplot_$(fn_postfix).png", fig)
    fig
end


begin
    fig = Figure()
    ax = Axis(fig[1,1], title = "plot", xlabel = "s", ylabel = "Runtime (seconds)" , xminorticksvisible = true, xminorgridvisible = true,
    xminorticks = IntervalsBetween(5))


    plot_lines!(df.s,df.col_red,"col_red")
    plot_lines!(df.s,df.std_mat,"std_mul")
    plot_lines!(df.s,df.row_red,"row_red")
    plot_lines!(df.s,df.lat_red_row,"row_latt_red")
    
    d_1,d_2 =  regres_theory(df.col_red, df.theo_col)
    lines!(df.s, d_2.*df.theo_col,linestyle = :dash, label="theoretical",linewidth = 2)
    
    Legend(fig[1,2], ax)
    save("Output/plot_$(fn_postfix).png", fig)
    fig
end