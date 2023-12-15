include("Col_red_utils.jl")


BenchmarkTools.DEFAULT_PARAMETERS.seconds = 0.1


b = 2; m = 22; s = 10000;
τ = 6 
z = rand(1:1000,s);
w = @. min(floor(2*log2(1:s)), m); A = rand(s,τ);

@profview reduced_mv_product_qmc(b, m, z, w, A)
     

P = reduced_mv_product_qmc(b, m, z, w, A);

begin
    s_val = collect(1:200:1500)
    T_val = Float64[]
    for s in s_val

        z = rand(1:1000,s)
        w = @. min(floor(2*log2(1:s)), m)
        A = rand(s,τ)
        # @btime reduced_mv_product_qmc(b, m, z, w, A);
        T = @belapsed reduced_mv_product_qmc($b, $m, $z, $w, $A);
        push!(T_val,T)
    end

end


fig = Figure()
ax = Axis(fig[1,1])
lines!(s_val,T_val)
fig