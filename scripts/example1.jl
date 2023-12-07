using Pkg
Pkg.activate("..")
using MFDecoupling

fp1 = ARGS[1]
fp2 = ARGS[2]

LL = parse(Int, ARGS[3])
Uin = parse(Float64, ARGS[4]);
Vin = parse(Float64, ARGS[5]);
UU = parse(Float64, ARGS[6])
VV = parse(Float64, ARGS[7]);
tmin = parse(Float64, ARGS[8]);
tmax = parse(Float64, ARGS[9]);
outputf = ARGS[10]


tspan = (tmin,tmax)

read_inputs(fp1, fp2, LL)
# solve
