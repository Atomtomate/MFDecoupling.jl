data_dir = "./data"

const Uin_4::Float64 = 3.0
const Vin_4::Float64 = 0.9
const Uin_3::Float64 = 0.3
const Vin_3::Float64 = 0.3
const Uin_2::Float64 = 0.1
const Vin_2::Float64 = 0.5
const Uin::Float64 = 0.5
const Vin::Float64 = 0.75
const tmin::Float64 = 0.0
const tmax::Float64 = 950.0
tspan = (tmin,tmax)
tsave = LinRange(0.0,950.0,5000)


UList_1 = union(LinRange(0.0,5.0,60), [Uin])
VList_1 = LinRange(0.0,1.0,60)

# coarse
UList_2 = union(LinRange(0.0,12.0,50), [Uin])#union(LinRange(0.0,0.2,10),LinRange(3.0,5.0,10),LinRange(10.0,20.0,10))
VList_2 = union(LinRange(0.0,1.0,50), [Vin])#union(LinRange(0.0,2.0,25),[0.4,0.5,0.6])

# Upt
UList_3 = LinRange(3.0,5.0,40)
VList_3 = union(LinRange(0.0,1.0,30), [Vin])


UList = union(UList_1, UList_2, UList_3)
VList = union(VList_1, VList_2, VList_3)
UV_list = Base.product(UList,VList)
println("#of jobs:", length(UV_list))

fl = readdir(abspath(data_dir))
fi = sort(map(x->parse(Int,first(split(last(split(x,"_")),"."))), fl))

println("There are $(length(fl)) files in the data dir")
println("Done check: ", length(fi) == last(fi))
