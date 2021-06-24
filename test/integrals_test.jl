using PyPlot
using CCBlade
using FastGaussQuadrature
using Test

## -- Problem definition --
# definition of the configuration
R = 1
χ = 30*pi/180
Θ = 5*pi/180
λ = 8

# sampling position
x = 0
ψ = 90 #Caution: not defined as the standard one
r = .5

# wake intensity
γt = -1
uz0 = γt/2

eps = 1e-3 #small tolerance to avoid evaluating r=R




## -- Perform the integral --

# init, span over a radius, fore-aft diameter if psi = 0
rr = range(-R-eps,R+eps,length=101)

ut = CCBlade.eval_ut(x, ψ, rr, R, χ) .* γt

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

ut_approx2 = zeros(3,length(rr))
CCBlade.eval_ut_approx!(ut_approx2, x, ψ, rr, R, χ, CCBlade.BranlardApprox())
ut_approx2.*=γt

for j = 1:4:101
    @test isapprox(uzt_approx[j], ut_approx2[1,j], atol=1e-12)
end

## -- Linear plots --

plt.figure(1)
plot(rr,ut[1,1,1,:])
plot(rr,uzt_approx)
plot(rr,ut_approx2[1,:])

plt.figure(2)
plot(rr,ut[2,1,1,:])
plot(rr,uψt_approx)
plot(rr,ut_approx2[2,:])
plt.show()

plt.figure(3)
plot(rr,ut[3,1,1,:])
plot(rr,urt_approx)
plot(rr,ut_approx2[3,:])



## -- Polar plots --
rr = range(0,R-eps,length=17)
ψψ = range(0,2*pi,length=13)'

ut = CCBlade.eval_ut(x, ψψ, rr, R, χ) .* γt
ur = CCBlade.eval_ur_0(ψψ, rr, R, χ) .* γt
ul = CCBlade.eval_ul(x, ψψ, rr, R, χ) .* 2. ./tan(χ/2) # =ul/uy0

## ut

f = plt.figure(4)
ax = f.add_subplot(111, polar=true)

u1 = ut[3,1,:,:]' .* cos.(ψψ) - ut[2,1,:,:]' .* sin.(ψψ)
u2 = ut[3,1,:,:]' .* sin.(ψψ) + ut[2,1,:,:]' .* cos.(ψψ)
# ax.quiver(ψψ, rr, u1, u2)
ax.quiver(ψψ, rr, -u1, -u2)
ax.invert_xaxis()
ax.set_theta_zero_location("W")  # theta=0 at the left
# ax.set_theta_direction(-1)  # theta increasing clockwise

f = plt.figure(5)
ax = f.add_subplot(111, polar=true)
ax.contour(ψψ, rr, ut[1,1,:,:]',[-.7,-.6,-.5,-.4,-.3])
ax.invert_xaxis()
ax.set_theta_zero_location("W")  # theta=0 at the left
ax.set_theta_direction(-1)  # theta increasing clockwise

## ur

f = plt.figure(6)
ax = f.add_subplot(111, polar=true)
ax.contour(ψψ, rr, ur[1,1,:,:]',[-.2,-.1,-.05,0,.05,.1,.2])
ax.invert_xaxis()
ax.set_theta_zero_location("W")  # theta=0 at the left
ax.set_theta_direction(-1)  # theta increasing clockwise

## ul

f = plt.figure(7)
ax = f.add_subplot(111, polar=true)

u1 = ul[3,1,:,:]' .* cos.(ψψ) - ul[2,1,:,:]' .* sin.(ψψ)
u2 = ul[3,1,:,:]' .* sin.(ψψ) + ul[2,1,:,:]' .* cos.(ψψ)
# ax.quiver(ψψ, rr, u1, u2)
ax.quiver(ψψ, rr, -u1, -u2)
ax.invert_xaxis()
ax.set_theta_zero_location("W")  # theta=0 at the left
# ax.set_theta_direction(-1)  # theta increasing clockwise

f = plt.figure(8)
ax = f.add_subplot(111, polar=true)
ax.contour(ψψ, rr, ul[1,1,:,:]'.*100.,[-10.,-5.,-2.,-.5,0.,.5,2.,5.,10.])
ax.invert_xaxis()
ax.set_theta_zero_location("W")  # theta=0 at the left
ax.set_theta_direction(-1)  # theta increasing clockwise



## -- verif at a given point that in-place and out-of-place agree --

no, we = gausslegendre( 10000 );
no .= 2 .* pi .* .5 * (no .+ 1)
k_u = zeros(3, length(no))
u = zeros(3)

x = .0
ψ = 1.
r = .33

#1
CCBlade.eval_u!(u, ψ, r, R, χ, k_u, no, we)

#2
ut = CCBlade.eval_ut(x, ψ, r, R, χ; n=10000 )
ur = CCBlade.eval_ur_0(ψ, r, R, χ)
ul = CCBlade.eval_ul(x, ψ, r, R, χ; n=10000 )


for i = 1:3
    @test isapprox(u[i], ut[i]+ur[i]+ul[i], atol=1e-12)
end


## -- verif in-place and out-of-place epsilons --

k_u = zeros(3, length(no))
I = zeros(3)
Iff = zeros(3)

ψ = 1.
r = .33
CT = 0.2

ϵx,ϵψ,ϵr = CCBlade.epsilons(ψ, r, R, χ, Θ, λ, CT)
ϵx2,ϵψ2,ϵr2 = CCBlade.epsilons!(ψ, r, R, χ, Θ, λ, CT, no, we, k_u, I, Iff )

@test isapprox(ϵx, ϵx2, atol=1e-12)
@test isapprox(ϵψ, ϵψ2, atol=1e-12)
@test isapprox(ϵr, ϵr2, atol=1e-12)
