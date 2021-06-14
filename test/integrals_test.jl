using PyPlot


## -- Problem definition --
# definition of the configuration
R = 1
χ = 30*pi/180

# sampling position
ψ = 90 #Caution: not defined as the standard one
r = .5
z = 0

# wake intensity
γt = -1
uz0 = γt/2

eps = 1e-3 #small tolerance to avoid evaluating r=R




## -- Perform the integral --

# init, span over a radius, fore-aft diameter if psi = 0
rr = range(-R-eps,R+eps,length=101)

ut = eval_ut(rr, ψ, z, χ, R) .* γt

#--  approximate formulation
Kzt_approx = rr./R .* tan(χ/2) #another approx exist, only valid on psi=0,z=0
uzt_approx = uz0 .* ( 1 .+ Kzt_approx .* cos(ψ))
Ft = Kzt_approx ./2 ./ tan(χ/2)

Kξt_approx = Kzt_approx ./ sin(χ)
Kxt_approx = Kξt_approx .* cos(χ)
uxt_approx = uz0 * (tan(χ/2) .- Kxt_approx .* cos(ψ))

# uyt_approx = -uz0 .* Ft .* sec(χ/2)^2 .* sin(ψ)

urt_approx = tan(χ/2) .* cos(ψ) .* uzt_approx - uz0 .* Ft .* sec(χ/2)^2
uψt_approx = -tan(χ/2) .* sin(ψ) .* uzt_approx



## -- Plots --

plt.figure(1)
plot(rr,ut[1,:])
plot(rr,urt_approx)

plt.figure(2)
plot(rr,ut[2,:])
plot(rr,uψt_approx)
plt.show()

plt.figure(3)
plot(rr,ut[3,:])
plot(rr,uzt_approx)

## -- Perform the integral --
#polar plots
rr = range(0,R-eps,length=17)
ψψ = range(0,2*pi,length=13)'

ut = eval_ut(rr, ψψ, z, χ, R) .* γt
ur = eval_ur(rr, ψψ, z, χ, R) .* γt

##

f = plt.figure(4)
ax = f.add_subplot(111, polar=true)

u1 = ut[1,:,:] .* cos.(ψψ) - ut[2,:,:] .* sin.(ψψ)
u2 = ut[1,:,:] .* sin.(ψψ) + ut[2,:,:] .* cos.(ψψ)
ax.quiver(ψψ, rr, u1, u2)


f = plt.figure(5)
ax = f.add_subplot(111, polar=true)
ax.contour(ψψ, rr, ut[3,:,:],[-.7,-.6,-.5,-.4,-.3])


f = plt.figure(6)
ax = f.add_subplot(111, polar=true)
ax.contour(ψψ, rr, ur[3,:,:],[-.2,-.1,-.05,0,.05,.1,.2])



##
using Test

no, we = gausslegendre( 10000 );
no .= 2 .* pi .* .5 * (no .+ 1)
k_u = zeros(3, length(no))
u = zeros(3)

r = .33
ψ = 1.
z = .05

#1
eval_u!(u, r, ψ, z, χ, R, k_u, no, we)

#2
ut = eval_ut(r, ψ, z, χ, R; n=10000 )
ur = eval_ur(r, ψ, z, χ, R)

for i = 1:3
    @test isapprox(u[i], ut[i]+ur[i], atol=1e-12)
end