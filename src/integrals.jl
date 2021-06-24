using FastGaussQuadrature
using LinearAlgebra


# Most of the formulas in this file come from Branlard2016
# - Be very careful here: he defines the psi angle with respect to the horizontal (i.e. not the standard definition of psi)
# - all induced velocities here are for a unit gamma_t = uz0*2
# - all inplace functions are additive

## ----------------- velocity induced by the tangential vorticity in the outer region -----------------

"""
    integ_ut!(k_u, θ, r, ψ, z, m, R)

Integrand for the velocity induced by the tangantial vorticity.

**Arguments**
- `k_u :: Array{Float,2}`: preallocated integrand of size [3*`Ni`], for the computation of [ux,uψ,ur]
- `θ :: Array{Float,1}`: integration variable of size `Ni`
- `x,ψ,r :: Float`: coordinates in the rotor c.s.
- `m,R :: Float`: skew angle param and scalars and rotor radius
"""
function integ_ut!(k_u, θ, x, ψ, r, m, R)

    apz = R .* (R .- r .* cos.(θ .- ψ))
    bpz = R .* m .* cos.(θ)

    apr = R .* x .* cos.(θ .- ψ)
    bpr =-R .* cos.(θ .- ψ)

    apψ = R .* x .* sin.(θ .- ψ)
    bpψ =-R .* sin.(θ .- ψ)

    a = R^2 .+ r.^2 .+ x.^2 .- 2 .*r.*R .* cos.(θ .- ψ)
    b = 2*m*R .* cos.(θ) .- 2 .*m.*r .* cos.(ψ) .- 2 .* x
    c = 1 + m^2

    den = (sqrt.(a) .* (2*sqrt.(a.*c) .+ b) )

    k_u[1,:] = 2 .* (apz .* sqrt.(c) .+ bpz .* sqrt.(a)) ./ den
    k_u[2,:] = 2 .* (apψ .* sqrt.(c) .+ bpψ .* sqrt.(a)) ./ den
    k_u[3,:] = 2 .* (apr .* sqrt.(c) .+ bpr .* sqrt.(a)) ./ den
end


"""
    eval_ut!(u, x, ψ, r, χ, R, k_u, no, we; γt=1.0)

Add the velocity induced by the tangantial component of vorticity in the wake to `u`.

**Arguments**
- `u :: Array{Float,1}`: preallocated velocity vector [ux,uψ,ur]
- `x,ψ,r :: Float`: coordinates where to evaluate the velocity
- `χ :: Float`: skew angle
- `R :: Float`: rotor radius
- `k_u :: Array{Float,2}`: preallocated integrand of size [3*`Ni`], for the computation of [ux,uψ,ur]
- `no :: Array{Float,1}`: `Ni` nodes for the Gauss-Legendre quadrature between [0,2pi]
- `we :: Array{Float,1}`: `Ni` weights for the Gauss-Legendre quadrature 
- `γt :: Float`: tangential vorticity in the outer sheet
"""
function eval_ut!(u, x, ψ, r, χ, R, k_u, no, we; γt=1.0)

    m = tan(χ)  # x = m*z

    integ_ut!(k_u,no,x,ψ,r,m,R)
    
    for j = 1:3 #is this faster than a matmul?
        u[j] += γt * dot( we, k_u[j,:] ) * 0.25 # quadrature with interval [0,2pi]: *pi, factor: 1/(4pi)
    end
    # u = γt .* k_u * we .* 0.25
end


"""
    eval_ut(xx, ψψ, rr, χ, R; n=10000)

Convenience function to compute the velocity induced by the tangential vorticity of unit intensity, mainly for plotting

**Arguments**
- `xx,ψψ,rr  :: Array{Float,1}`: list of coordinates where to evaluate the velocity, respectively of size `nx,nψ,nr`
- `χ :: Float`: skew angle
- `R :: Float`: rotor radius
- `n :: Int`: number of integration points in the quadrature

**Returns**
- `u::Array{Float,4}`: induced velocity, of size [3 x `nx` x `nψ` x `nr`] (meshgrid)
"""
function eval_ut(xx, ψψ, rr, χ, R; n=10000)

    # integration stuff
    nodes, weights = gausslegendre( n );
    nodes .= 2 .* pi .* .5 * (nodes .+ 1)
    # weights .*= 2 * pi / 2 # the change of variable [-1:1]->[0:2pi] can be included in the weights. We don't do it though, since the pi will cancel out withint the integrals.

    # params
    nx = length(xx)
    nψ = length(ψψ)
    nr = length(rr)

    # prealloc
    k_u = zeros(3, length(nodes))
    ut = zeros(3, nr*nψ*nx)

    for i = 0 : nr*nψ*nx - 1
        ix = mod(i,nx) +1
        iψ = mod(floor(Int64,i/nx),nψ) +1
        ir = floor(Int64,i/(nx*nψ)) +1

        eval_ut!( view(ut,:,i+1), xx[ix],ψψ[iψ],rr[ir], χ, R, k_u, nodes, weights)
    end
    
    ut = reshape(ut,(3,nx,nψ,nr))
    
    return ut
end


# -- approximated versions --
# TODO: regroup with the not approximated version

abstract type AbstractApproxUt end

struct ColemanApprox <: AbstractApproxUt end

struct PittPetersApprox <: AbstractApproxUt end

struct BranlardApprox <: AbstractApproxUt end

"""
    eval_ut_approx!(ut_approx, rr, ψψ, zz, χ, R ; γt=1.0)

Add the velocity induced by the tangantial component of vorticity in the wake to `u`.

**Arguments**
- `u_approx :: Array{Float,1}`: preallocated velocity vector [ux,uψ,ur]
- `x,ψ,r :: Float`: coordinates where to evaluate the velocity
- `χ :: Float`: skew angle
- `R :: Float`: rotor radius
- `approx<:AbstractApproxUt`: approximation type
- `γt :: Float`: tangential vorticity in the outer sheet
"""
eval_ut_approx!

function eval_ut_approx!(u_approx, x, ψ, r, χ, R, approx::AbstractApproxUt ; γt=1.0)
    @error "need to choose an ApproxUt"
end

function eval_ut_approx(u_approx, x, ψ, r, χ, R, approx::ColemanApprox ; γt=1.0)
    Kzt_approx = r./R .* tan(χ/2)
    u_approx[1,:] .+= ( 1 .+ Kzt_approx .* cos(ψ)) *.5 *γt
end

function eval_ut_approx!(u_approx, x, ψ, r, χ, R, approx::PittPetersApprox; γt=1.0)
    Kzt_approx = r./R .* tan(χ/2) .* (15.0*pi/32.)
    u_approx[1,:] .+= ( 1 .+ Kzt_approx .* cos(ψ)) *.5 *γt
end

function eval_ut_approx!(u_approx, x, ψ, r, χ, R, approx::BranlardApprox; γt=1.0)

    Kzt_approx = r./R .* tan(χ/2) #another approx exists, only valid on psi=0,z=0
    u_approx[1,:] .+= ( 1 .+ Kzt_approx .* cos(ψ)) *.5 *γt
    Ft = Kzt_approx ./2 ./ tan(χ/2)

    # Kξt_approx = Kzt_approx ./ sin(χ)
    # Kxt_approx = Kξt_approx .* cos(χ)
    # uxt_approx =  (tan(χ/2) .- Kxt_approx .* cos(ψ))
    # uyt_approx = -uz0 .* Ft .* sec(χ/2)^2 .* sin(ψ)

    u_approx[3,:] .+= (tan(χ/2) .* cos(ψ) .* u_approx[1,:] .- Ft .* sec(χ/2)^2 *.5) *γt
    u_approx[2,:] .+= (-tan(χ/2) .* sin(ψ) .* u_approx[1,:]) *γt
end



## ----------------- velocity induced by the root vorticity -----------------


"""
    eval_ur_0!(u, ψ, r, χ, R; Γr = 1.0)

Add the velocity induced by the root component of vorticity in the wake to `u`, in the rotor plane (x=0).

**Arguments**
- `u :: Array{Float,1}`: preallocated velocity vector [ux,uψ,ur]
- `ψ,r :: Float`: coordinates where to evaluate the velocity
- `χ :: Float`: skew angle
- `R :: Float`: rotor radius
- `Γr :: Float`: circulation of the root vortex
"""
function eval_ur_0!(u, ψ, r, χ, R; Γr = 1.0)

    u[1] += Γr .* sin(χ) .* sin.(ψ) ./ (4 .*pi .* r .* (1.0 .- cos.(ψ) .* sin(χ)) )
    u[3] += Γr .* cos(χ) ./ (4*pi .* r .* (1.0 .- cos.(ψ) .* sin(χ)) )

end

"""
    eval_ur_ff!(u, ψ, r, χ, R; Γr = 1.0)

Same as `eval_ur!` but for the far field ``(x \\rightarrow \\infty)``.
"""
function eval_ur_ff!(u, ψ, r, χ, R, Θ; Γr = 1.0)
    
    u[2] += Γr ./ (2*pi .* r .* cos.(χ) .* cos(Θ) )

end

"""
    eval_ur(xx, ψψ, rr, χ, R)

Convenience function to compute the velocity induced by the root vortex of unit circulation in the rotor plane, mainly for plotting

**Arguments**
- `xx,ψψ,rr  :: Array{Float,1}`: list of coordinates where to evaluate the velocity, respectively of size `nx,nψ,nr`
- `χ :: Float`: skew angle
- `R :: Float`: rotor radius

**Returns**
- `u::Array{Float,4}`: induced velocity, of size [3 x `nx` x `nψ` x `nr`] (meshgrid)
"""
function eval_ur_0(ψψ, rr, χ, R)
    
    nx = 1
    nψ = length(ψψ)
    nr = length(rr)

    ur = zeros(3, nr*nψ*nx)

    for i = 0 : nr*nψ*nx - 1
        # ix = mod(i,nx) +1
        iψ = mod(floor(Int64,i/nx),nψ) +1
        ir = floor(Int64,i/(nx*nψ)) +1

        eval_ur_0!( view(ur,:,i+1) , ψψ[iψ], rr[ir], χ, R)
    end
    ur = reshape(ur,(3,nx,nψ,nr))

    return ur
end


## ----------------- velocity induced by the longitudinal vorticity in the outer sheet -----------------

#TODO: add vel induced by longitundinal vorticity, u_l 



## ----------------- total induced velocity -----------------

"""
    eval_u!(u, x, ψ, r, χ, R, k_u, no, we)

Add the velocity induced by unit tangantial and root components of vorticity in the wake to `u`.

**Arguments**
- `u :: Array{Float,1}`: preallocated velocity vector [ux,uψ,ur]
- `x,ψ,r :: Float`: coordinates where to evaluate the velocity
- `χ :: Float`: skew angle
- `R :: Float`: rotor radius
- `k_u :: Array{Float,2}`: preallocated integrand of size [3*`Ni`], for the computation of [ux,uψ,ur]
- `no :: Array{Float,1}`: `Ni` nodes for the Gauss-Legendre quadrature between [0,2pi]
- `we :: Array{Float,1}`: `Ni` weights for the Gauss-Legendre quadrature 
"""
function eval_u!(u, x, ψ, r, χ, R, k_u, no, we)

    eval_ur_0!( u, ψ, r, χ, R)
    eval_ut!( u, x, ψ, r, χ, R, k_u, no, we)

end


## ----------------- epsilon functions for the BEM -----------------
