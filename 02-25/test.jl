##
using Plots
include("LagrangeInterp.jl")

Xdata = [1, 2,3, 4, 5.5];
Xdata = LinRange(1,5.5, 30);
Ydata = sin.(Xdata)
f = LagrangeInterp.interp(Xdata, Ydata)

plot( t->sin(t), 1,5.5, lw=5)
plot!(f, 1, 5.5, lw=5)
scatter!(Xdata, Ydata)