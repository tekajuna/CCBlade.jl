using Base: _linspace
using PyPlot
using FastGaussQuadrature
using LinearAlgebra

## -- Problem definition --
# definition of the configuration
R = 1
χ = 30*pi/180

# sampling position
ψ = 0
r = .5
z = 0

# wake intensity
γt = 1
uz0 = γt/2


#θ ; vect - integration variable
#ψ,r,z ; scalars
#m,R ; scalars
#k_u ; [size(θ) x 3] - output
function integ_ut!(k_u,θ,r,ψ,z,m,R)

    apz = R .* (R .- r .* cos.(θ .- ψ))
    bpz = R .* m .* cos.(θ)

    apr = R .* z .* cos.(θ .- ψ)
    bpr =-R .* cos.(θ .- ψ)

    apψ = R .* z .* sin.(θ .- ψ)
    bpψ =-R .* sin.(θ .- ψ)

    a = R^2 .+ r.^2 .+ z.^2 .- 2 .*r.*R .* cos.(θ .- ψ)
    b = 2*m*R .* cos.(θ) .- 2 .*m.*r .* cos.(ψ) .- 2 .* z
    c = 1 + m^2

    den = (sqrt.(a) .* (2*sqrt.(a.*c) .+ b) )

    k_u[1,:] = 2 .* (apz .* sqrt.(c) .+ bpz .* sqrt.(a)) ./ den
    k_u[2,:] = 2 .* (apr .* sqrt.(c) .+ bpr .* sqrt.(a)) ./ den
    k_u[3,:] = 2 .* (apψ .* sqrt.(c) .+ bpψ .* sqrt.(a)) ./ den
end

function eval_ut(rr, ψψ, zz, χ, R)

    # integration stuff
    nodes, weights = gausslegendre( 10000 );
    nodes .= 2 .* pi .* .5 * (nodes .+ 1)

    # params
    m = tan(χ)  # x = m*z
    nr = length(rr)
    nψ = length(ψψ)
    nz = length(zz)

    # prealloc
    k_u = zeros(3, length(nodes))
    ut = zeros(3, nr*nψ*nz)

    for i = 0 : nr*nψ*nz - 1
        ir = mod(i,nr) +1
        iψ = mod(floor(Int64,i/nr),nψ) +1
        iz = floor(Int64,i/(nr*nψ)) +1
        integ_ut!(k_u,nodes,rr[ir],ψψ[iψ],zz[iz],m,R)
    
        for j = 1:3
            ut[j,i+1] = dot( weights, k_u[j,:] ) * pi # quadrature with interval [0,2pi]
        end
    end
    ut .*= 1 / (4*pi) 

    ut = reshape(ut,(3,nr,nψ,nz))
    
    return ut
end

## -- Perform the integral --

# init, span over a radius, fore-aft diameter if psi = 0
rr = range(-R,R,length=101)

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
plot(rr,uzt_approx)

plt.figure(2)
plot(rr,ut[2,:])
plot(rr,urt_approx)

plt.figure(3)
plot(rr,ut[3,:])
plot(rr,uψt_approx)
plt.show()

