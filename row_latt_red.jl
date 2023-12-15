include("Col_red_utils.jl")


BenchmarkTools.DEFAULT_PARAMETERS.seconds = 0.1


b = 2; m = 12; s = 10;
τ = 20
local z_r = rand(1:1000,s);
local w_r = @. min(floor(2*log2(1:s)), m); local A_r = rand(s,τ);

#@profview reduced_mv_product_qmc(b, m, z_r, w_r, A_r)
     

#P = reduced_mv_product_qmc(b, m, z_r, w_r, A_r);

begin
    s_val = collect(1:200:1500)
    T_val3 = Float64[]
    for s in s_val

        local z_r = rand(1:1000,s)
        local w_r = @. min(floor(2*log2(1:s)), m)
        local A_r = rand(s,τ)
        # @btime reduced_mv_product_qmc(b, m, z_r, w_r, A_r);
        T = @belapsed reduced_mv_product_qmc($b, $m, $z_r, $w_r, $A_r);
        push!(T_val3,T)
    end

end


#fig = Figure()
#ax = Axis(fig[1,1], title = "lat_row_red", xlabel = "s", ylabel = "t")
#lines!(s_val,T_val3)
#save("plot_lat_row_red_1.png", fig)
#fig