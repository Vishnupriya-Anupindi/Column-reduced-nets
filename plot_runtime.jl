using BenchmarkTools
using Statistics
using CairoMakie
using GLM, StatsBase, DataFrames, CSV

b=2
case = 11   # 1 means plot with random matrices, 2 means plot with sobol and niederreiter
s_range = 1600
step_size = 200
m = 12
τ = 20
s_val = collect(1:step_size:s_range)
fn_postfix = "case$(case)_m$(m)_s$(s_range)"

#df = DataFrame(s = df.s, row_red = df.row_red, col_red = df.col_red, std_mat = df.std_mat, lat_red_row = df.lat_red_row, theo_col = df.theo_col)
df = CSV.read("runtime_$(fn_postfix)_b$(b).csv", DataFrame)

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

function plot_lines!(s_val,T_val,label,ptstyle, colour)
    lines!(s_val,T_val; linestyle = :solid, color = colour, label , linewidth = 1.5)
    scatter!(s_val,T_val; marker = ptstyle, color = colour, label , markersize = 12)
    #c_1,c_2 = regres_comp(s_val,T_val)
    #lines!(s_val,c_1.*(s_val.^c_2), color = "light gray", linestyle = :dot)
end


T_theory_rr = [runtime_theory_rr(τ, b, m, s) for s in s_val]

begin
    fig = Figure()
    ax = Axis(fig[1,1], title = "", xlabel = "log s", ylabel = "Runtime (log seconds)",xscale = log10, yscale = log10, xminorticksvisible = true, xminorgridvisible = true,
    xminorticks = IntervalsBetween(5),yminorticksvisible = true, yminorgridvisible = true,
    yminorticks = IntervalsBetween(5))
    
    plot_lines!(df.s,df.col_red,"Column reduced mm product",:circle, :blue)
    plot_lines!(df.s,df.std_mat,"Standard mm product",:rect, :orange)
    plot_lines!(df.s,df.row_red,"Row reduced mm product", :xcross, :darkred)
    #plot_lines!(df.s,df.lat_red_row,"Lattice row reduced mm product",:utriangle)

    d_1,d_2 =  regres_theory(df.col_red, df.theo_col)
    lines!(df.s, d_2.*df.theo_col,linestyle = :dash, label="Theoretical estimate column reduced",linewidth = 1.5, color = :black)

    d_3,d_4 =  regres_theory(df.row_red, T_theory_rr)
    lines!(df.s, d_4.*T_theory_rr,linestyle = :dash, label="Theoretical estimate row reduced",linewidth = 1.5, color = :gray)
    
    if case ==2
        plot_lines!(df.s,df.col_red_sobol,"Sobol_col_red")
        plot_lines!(df.s,df.std_mat_sobol,"Sobol_std_mul")

        plot_lines!(df.s,df.col_red_nied,"Nied_col_red")
        plot_lines!(df.s,df.std_mat_nied,"Nied_std_mul")
    end

    axislegend(ax, merge = true, position = :lt)
    save("Output/logplot_$(fn_postfix)_rr.png", fig)
    fig
end


begin
    fig = Figure()
    ax = Axis(fig[1,1], title = "", xlabel = "s", ylabel = "Runtime (log seconds)" , yscale = log10, yminorticksvisible = true, yminorgridvisible = true,
    yminorticks = IntervalsBetween(5))


    plot_lines!(df.s,df.col_red,"Column reduced mm product",:circle, :blue)
    plot_lines!(df.s,df.std_mat,"Standard mm product", :rect, :orange)
    plot_lines!(df.s,df.row_red,"Row reduced mm product", :xcross, :darkred)
    #plot_lines!(df.s,df.lat_red_row,"row_latt_red", :utriangle)
    
    d_1,d_2 =  regres_theory(df.col_red, df.theo_col)
    lines!(df.s, d_2.*df.theo_col,linestyle = :dash, label="Theoretical estimate column reduced",linewidth = 1.5, color = :black)

    d_3,d_4 =  regres_theory(df.row_red, T_theory_rr)
    lines!(df.s, d_4.*T_theory_rr,linestyle = :dash, label="Theoretical estimate row reduced",linewidth = 1.5, color = :gray)
    
    axislegend(ax, merge = true, position = :rb)
    save("Output/semilog_plot_$(fn_postfix)_rr.png", fig)
    fig
end