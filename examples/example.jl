using Gallium
fileLoc,_ = splitdir(@__FILE__)
push!(LOAD_PATH,"$(fileLoc)/../src/")
using CCBlade
using PyPlot
PyPlot.close("all")
# -------- propeller example ----------------


# geometry
scale = 5 #for verification purposes: increase Re and M
Rhub =.0254*.5*scale
Rtip = .0254*3.0*scale


r = .0254*[0.7526, 0.7928, 0.8329, 0.8731, 0.9132, 0.9586, 1.0332,
     1.1128, 1.1925, 1.2722, 1.3519, 1.4316, 1.5114, 1.5911,
     1.6708, 1.7505, 1.8302, 1.9099, 1.9896, 2.0693, 2.1490, 2.2287,
     2.3084, 2.3881, 2.4678, 2.5475, 2.6273, 2.7070, 2.7867, 2.8661, 2.9410]*scale
chord = .0254*[0.6270, 0.6255, 0.6231, 0.6199, 0.6165, 0.6125, 0.6054, 0.5973, 0.5887,
          0.5794, 0.5695, 0.5590, 0.5479, 0.5362, 0.5240, 0.5111, 0.4977,
          0.4836, 0.4689, 0.4537, 0.4379, 0.4214, 0.4044, 0.3867, 0.3685,
          0.3497, 0.3303, 0.3103, 0.2897, 0.2618, 0.1920]*scale

theta = pi/180.0*[40.2273, 38.7657, 37.3913, 36.0981, 34.8803, 33.5899, 31.6400,
                   29.7730, 28.0952, 26.5833, 25.2155, 23.9736, 22.8421, 21.8075,
                   20.8586, 19.9855, 19.1800, 18.4347, 17.7434, 17.1005, 16.5013,
                   15.9417, 15.4179, 14.9266, 14.4650, 14.0306, 13.6210, 13.2343,
                   12.8685, 12.5233, 12.2138]
B = 2  # number of blades

# aftype = CCBlade.af_from_aerodynfile("airfoils/xf-clarky-il-1000000.csv")
#
# n = length(r)
# af = Array(CCBlade.AirfoilData, n)
# for i = 1:n
#     af[i] = aftype
# end

#Load the prop airfoil ND table
ClarkY_non_extrap = JLD.load("$(fileLoc)/airfoils/af_prop_ClarkY.jld")
ClarkY_non_extrap = ClarkY_non_extrap["NDtable"]
TSR = .231

grid_alphas=[i for i in -180:1.0:180]

n = length(r)
# af = Array{Any}(n)
# for i = 1:n
#
#     r_over_R = r[i]
#     c_over_r = chord[i]/r[i]
#     NDextrap3D_3Dtable = AirfoilPrep.NDTable_correction3D_extrap(ClarkY_non_extrap,r_over_R,c_over_r,TSR;grid_alphas=grid_alphas)
#     afspl = AirfoilPrep.SplineND_from_tableND(NDextrap3D_3Dtable)
#     af[i] = afspl
# end
# JLD.save("$(fileLoc)/airfoils/af.jld", "af", af)
af = JLD.load("$(fileLoc)/airfoils/af.jld")
af = af["af"]


precone = 0.0

rho = 1.225
mu = 1.789E-5
a = 1225 #km/h

Vinf = 10.0
Omega = 8000.0*pi/30.0*scale/4.19

inflow = CCBlade.simpleinflow(Vinf, Omega, r, precone, rho, mu, a)
rotor = CCBlade.Rotor(r, chord, theta, af, Rhub, Rtip, B, precone)

turbine = false

# @enter CCBlade.distributedloads(rotor, inflow, turbine)
Np, Tp = CCBlade.distributedloads(rotor, inflow, turbine)

PyPlot.figure("loads")
PyPlot.plot(r/Rtip, Np,label = "Normal")
PyPlot.plot(r/Rtip, Tp, label = "Tangential")
PyPlot.legend(loc = "best")


J = linspace(0.1, 0.9, 20)

# Omega = 8000.0*pi/30
n = Omega/(2*pi)
D = 2*Rtip*cos(precone)

eff = zeros(20)
CT = zeros(20)
CQ = zeros(20)

N = 20

Re = zeros(N)
M = zeros(N)
TSR_arry = zeros(N)

for i = 1:N
    Vinf = J[i] * D * n
    Vinfeff = sqrt((Omega*Rtip*0.7)^2+Vinf^2)
    Re[i] = rho/mu*Vinfeff.*chord[22]
    M[i] = Vinfeff/a
    TSR_arry[i] = TSR = 1/(4*pi*J[i])

    inflow = CCBlade.simpleinflow(Vinf, Omega, r, precone, rho, mu, a)

    T, Q = CCBlade.thrusttorque(rotor, [inflow], turbine)
    eff[i], CT[i], CQ[i] = CCBlade.nondim(T, Q, Vinf, Omega, rho, Rtip, precone, turbine)

end

Re_ave = mean(Re)
M_ave = mean(M)
TSR_ave = mean(TSR_arry)

println("Re_ave $Re_ave \n M_ave $M_ave")

PyPlot.figure("CT")
PyPlot.plot(J, CT)
PyPlot.xlabel("J")
PyPlot.xlabel("CT")

PyPlot.figure("CQ")
PyPlot.plot(J, CQ)
PyPlot.xlabel("J")
PyPlot.xlabel("CQ")

PyPlot.figure("eta")
PyPlot.plot(J, eff)
PyPlot.xlabel("J")
PyPlot.xlabel(L"\eta_prop")

println("J $J")
println("eta $eff")

#M =
J_old = 0.1:0.042105263157894736:0.9
eta_old = [0.254429, 0.343448, 0.422759, 0.493193, 0.555543, 0.610559, 0.658931, 0.701294, 0.738223, 0.770222, 0.797684, 0.820705, 0.838741, 0.850184, 0.850937, 0.828907, 0.731107, 0.0, 0.0, 0.0]
