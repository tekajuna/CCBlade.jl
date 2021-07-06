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
    eval_ut!(u, x, ψ, r, R, χ, k_u, no, we; γt=1.0)

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
function eval_ut!(u, x, ψ, r, R, χ, k_u, no, we; γt=1.0)

    m = tan(χ)  # x = m*z

    integ_ut!(k_u,no,x,ψ,r,m,R)
    
    for j = 1:3 #is this faster than a matmul?
        u[j] += γt * dot( we, k_u[j,:] ) * 0.25 # quadrature with interval [0,2pi]: *pi, factor: 1/(4pi)
    end
    # u = γt .* k_u * we .* 0.25
end


"""
    eval_ut(xx, ψψ, rr, R, χ; n=10000)

Convenience function to compute the velocity induced by the tangential vorticity of unit intensity, mainly for plotting

**Arguments**
- `xx,ψψ,rr  :: Array{Float,1}`: list of coordinates where to evaluate the velocity, respectively of size `nx,nψ,nr`
- `χ :: Float`: skew angle
- `R :: Float`: rotor radius
- `n :: Int`: number of integration points in the quadrature

**Returns**
- `u::Array{Float,4}`: induced velocity, of size [3 x `nx` x `nψ` x `nr`] (meshgrid)
"""
function eval_ut(xx, ψψ, rr, R, χ; n=10000)

    # integration stuff
    nodes, weights = gausslegendre_0_2π(n)
    
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

        eval_ut!( view(ut,:,i+1), xx[ix],ψψ[iψ],rr[ir], R, χ, k_u, nodes, weights)
    end
    
    ut = reshape(ut,(3,nx,nψ,nr))
    
    return ut
end


### -- approximated versions --
#TODO: regroup with the not approximated version

abstract type AbstractApproxUt end

struct ColemanApprox <: AbstractApproxUt end

struct PittPetersApprox <: AbstractApproxUt end

struct BranlardApprox <: AbstractApproxUt end

"""
    eval_ut_approx!(ut_approx, rr, ψψ, zz, R, χ ; γt=1.0)

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

function eval_ut_approx!(u_approx, x, ψ, r, R, χ, approx::AbstractApproxUt ; γt=1.0)
    @error "need to choose an ApproxUt"
end

function eval_ut_approx(u_approx, x, ψ, r, R, χ, approx::ColemanApprox ; γt=1.0)
    Kzt_approx = r./R .* tan(χ/2)
    u_approx[1,:] .+= ( 1 .+ Kzt_approx .* cos(ψ)) *.5 *γt
end

function eval_ut_approx!(u_approx, x, ψ, r, R, χ, approx::PittPetersApprox; γt=1.0)
    Kzt_approx = r./R .* tan(χ/2) .* (15.0*pi/32.)
    u_approx[1,:] .+= ( 1 .+ Kzt_approx .* cos(ψ)) *.5 *γt
end

function eval_ut_approx!(u_approx, x, ψ, r, R, χ, approx::BranlardApprox; γt=1.0)

    Kzt_approx = r./R .* tan(χ/2) #another approx exists, only valid on psi=0,z=0
    u_approx[1,:] .+= ( 1 .+ Kzt_approx .* cos(ψ)) *.5 *γt
    Ft = Kzt_approx ./2 ./ tan(χ/2)

    u_approx[3,:] .+= (tan(χ/2) .* cos(ψ) .* u_approx[1,:] .- Ft .* sec(χ/2)^2 *.5) *γt
    u_approx[2,:] .+= (-tan(χ/2) .* sin(ψ) .* u_approx[1,:]) *γt
end



## ----------------- velocity induced by the root vorticity -----------------


"""
    eval_ur_0!(u, ψ, r, R, χ; Γr = 1.0)

Add the velocity induced by the root component of vorticity in the wake to `u`, in the rotor plane (x=0).

**Arguments**
- `u :: Array{Float,1}`: preallocated velocity vector [ux,uψ,ur]
- `ψ,r :: Float`: coordinates where to evaluate the velocity
- `χ :: Float`: skew angle
- `R :: Float`: rotor radius
- `Γr :: Float`: circulation of the root vortex
"""
function eval_ur_0!(u, ψ, r, R, χ; Γr = 1.0)

    u[1] += Γr .* sin(χ) .* sin.(ψ) ./ (4 .*pi .* r .* (1.0 .- cos.(ψ) .* sin(χ)) )
    u[2] += Γr .* cos(χ) ./ (4*pi .* r .* (1.0 .- cos.(ψ) .* sin(χ)) )

end

"""
    eval_ur_ff!(u, ψ, r, R, χ; Γr = 1.0)

Same as `eval_ur!` but for the far field ``(x \\rightarrow \\infty)``.
"""
function eval_ur_ff!(u, ψ, r, R, χ, Θ; Γr = 1.0)
    
    u[2] += Γr ./ (2*pi .* r .* cos.(χ) .* cos(Θ) )

end

"""
    eval_ur(xx, ψψ, rr, R, χ)

Convenience function to compute the velocity induced by the root vortex of unit circulation in the rotor plane, mainly for plotting

**Arguments**
- `xx,ψψ,rr  :: Array{Float,1}`: list of coordinates where to evaluate the velocity, respectively of size `nx,nψ,nr`
- `χ :: Float`: skew angle
- `R :: Float`: rotor radius

**Returns**
- `u::Array{Float,4}`: induced velocity, of size [3 x `nx` x `nψ` x `nr`] (meshgrid)
"""
function eval_ur_0(ψψ, rr, R, χ)
    
    nx = 1
    nψ = length(ψψ)
    nr = length(rr)

    ur = zeros(3, nr*nψ*nx)

    for i = 0 : nr*nψ*nx - 1
        # ix = mod(i,nx) +1
        iψ = mod(floor(Int64,i/nx),nψ) +1
        ir = floor(Int64,i/(nx*nψ)) +1

        eval_ur_0!( view(ur,:,i+1) , ψψ[iψ], rr[ir], R, χ)
    end
    ur = reshape(ur,(3,nx,nψ,nr))

    return ur
end


## ----------------- velocity induced by the longitudinal vorticity in the outer sheet -----------------

"""
    integ_ul!(k_u, θ, r, ψ, z, m, R)

Integrand for the velocity induced by the longitudinal vorticity.

**Arguments**
- `k_u :: Array{Float,2}`: preallocated integrand of size [3*`Ni`], for the computation of [ux,uψ,ur]
- `θ :: Array{Float,1}`: integration variable of size `Ni`
- `x,ψ,r :: Float`: coordinates in the rotor c.s.
- `m,R :: Float`: skew angle param and scalars and rotor radius
"""
function integ_ul!(k_u, θ, x, ψ, r, m, R)
    rt = r./R
    xt = x./R

    D1 = sqrt.(1 .+ rt.^2 .+ xt.^2 .- 2.0.*rt.* cos.(θ .- ψ))
    D2 = m .* (cos.(θ) .- rt .* cos.(ψ)) .- xt 
    D  = D1 .* (D2 .+ sqrt.(1+m^2) .* D1 )

    k_u[1,:] = m .* ( -sin.(θ) .+ rt .* sin.(ψ) ) ./ D
    k_u[2,:] = ( -m .* xt .* cos.(ψ) .- cos.(θ .- ψ) .+ rt ) ./ D
    k_u[3,:] = ( -m .* xt .* sin.(ψ) .+ sin.(θ .- ψ) ) ./ D
end


"""
    eval_ul!(u, x, ψ, r, R, χ, k_u, no, we; γl=1.0)

Add the velocity induced by the longitudinal component of vorticity in the wake to `u`.

**Arguments**
- `u :: Array{Float,1}`: preallocated velocity vector [ux,uψ,ur]
- `x,ψ,r :: Float`: coordinates where to evaluate the velocity
- `χ :: Float`: skew angle
- `R :: Float`: rotor radius
- `k_u :: Array{Float,2}`: preallocated integrand of size [3*`Ni`], for the computation of [ux,uψ,ur]
- `no :: Array{Float,1}`: `Ni` nodes for the Gauss-Legendre quadrature between [0,2pi]
- `we :: Array{Float,1}`: `Ni` weights for the Gauss-Legendre quadrature 
- `γl :: Float`: tangential vorticity in the outer sheet
"""
function eval_ul!(u, x, ψ, r, R, χ, k_u, no, we; γl=1.0)

    m = tan(χ)  # x = m*z

    integ_ul!(k_u,no,x,ψ,r,m,R)
    
    for j = 1:3 #is this faster than a matmul?
        u[j] += γl * dot( we, k_u[j,:] ) * 0.25 # quadrature with interval [0,2pi]: *pi, factor: 1/(4pi)
    end
    # u = γt .* k_u * we .* 0.25
end


"""
    eval_ul(xx, ψψ, rr, R, χ; n=10000)

Convenience function to compute the velocity induced by the longitudinal vorticity of unit intensity, mainly for plotting

**Arguments**
- `xx,ψψ,rr  :: Array{Float,1}`: list of coordinates where to evaluate the velocity, respectively of size `nx,nψ,nr`
- `χ :: Float`: skew angle
- `R :: Float`: rotor radius
- `n :: Int`: number of integration points in the quadrature

**Returns**
- `u::Array{Float,4}`: induced velocity, of size [3 x `nx` x `nψ` x `nr`] (meshgrid)
"""
function eval_ul(xx, ψψ, rr, R, χ; n=10000)

    # integration stuff
    nodes, weights = gausslegendre_0_2π(n)

    # params
    nx = length(xx)
    nψ = length(ψψ)
    nr = length(rr)

    # prealloc
    k_u = zeros(3, length(nodes))
    ul = zeros(3, nr*nψ*nx)

    for i = 0 : nr*nψ*nx - 1
        ix = mod(i,nx) +1
        iψ = mod(floor(Int64,i/nx),nψ) +1
        ir = floor(Int64,i/(nx*nψ)) +1

        eval_ul!( view(ul,:,i+1), xx[ix],ψψ[iψ],rr[ir], R, χ, k_u, nodes, weights)
    end
    
    ul = reshape(ul,(3,nx,nψ,nr))
    
    return ul
end

### --- approximations ---

# See Branlard2016
# function eval_ul_approx!(u_approx, x, ψ, r, R, χ, approx::BranlardApprox; γl=1.0)
#     # u_approx[1,:] .+= ...
#     # u_approx[3,:] .+= ...
#     # u_approx[2,:] .+= ...
# end


## ----------------- total induced velocity -----------------

"""
    eval_u!(u, ψ, r, R, χ, k_u, no, we; γt=1.0, γl=1.0, Γr=1.0 )

Add the velocity induced by unit tangantial, longitudinal and root components of vorticity in the wake to `u`, evaluated at the rotor disk (x=0).

**Arguments**
- `u :: Array{Float,1}`: preallocated velocity vector [ux,uψ,ur]
- `ψ,r :: Float`: coordinates where to evaluate the velocity
- `χ :: Float`: skew angle
- `R :: Float`: rotor radius
- `k_u :: Array{Float,2}`: preallocated integrand of size [3*`Ni`], for the computation of [ux,uψ,ur]
- `no :: Array{Float,1}`: `Ni` nodes for the Gauss-Legendre quadrature between [0,2pi]
- `we :: Array{Float,1}`: `Ni` weights for the Gauss-Legendre quadrature 
"""
function eval_u!(u, ψ, r, R, χ, k_u, no, we; γt=1.0, γl=1.0, Γr=1.0 )

    eval_ur_0!( u, ψ, r, R, χ; Γr)
    eval_ut!( u, 0.0, ψ, r, R, χ, k_u, no, we; γt)
    eval_ul!( u, 0.0, ψ, r, R, χ, k_u, no, we; γl)

end


## ----------------- convenience functions -------------------------

"""
(private)    gammas(r, R, χ, λ, CT ; γt=1.0)

Compute the vorticity/circulation ratio between various wake components.

**Arguments**
- `r :: Float`: radial location
- `χ :: Float`: skew angle
- `R :: Float`: rotor radius
- `λ :: Float`: tip speed ratio
- `CT :: Float`: expected/approximated rotor thrust coefficient
- `γt :: Float`: tangential vorticity in wake cylinder sheet (keep 1 to obtain ratios as an output)

**Returns**
- `Γr::Float`: circulation of root vortex
- `γl::Float`: longitudinal vorticity wake cylinder sheet
- `γb::Float`: radial bound vorticity on rotor disk
"""
@inline function gammas(r, R, χ, λ, CT ; γt=1.0)
    h = pi * R / λ * (1 + sqrt(1 - CT))

    Γr = -γt * h / cos(χ)
    γl = γt / (2 * pi * R) * h / cos(χ)
    γb = γt / (2 * pi * r) * h / cos(χ)

    return Γr, γl, γb
end

"""
(private)   gausslegendre_0_2π(n)

    Compute nodes and weights for a Gauss-Legendre quadrature between 0 and 2π such that
    ```math
     \\int_0^{2\\pi} f(x) dx = \\pi \\sum_n \\textrm{weights}_n f(\\textrm{nodes}_n)
    ```
"""
@inline function gausslegendre_0_2π(n)
    nodes, weights = gausslegendre( n );
    nodes .= 2 .* pi .* .5 * (nodes .+ 1)
    # weights .*= 2 * pi / 2 # the change of variable [-1:1]->[0:2pi] can be included in the weights. We don't do it though, since the pi will cancel out withint the integrals.
    return nodes, weights
end


## ----------------- epsilon functions for the BEM -----------------

"""
    epsilons(ψ, r, R, χ, Θ, λ, CT )

Compute the epsilon factors at the rotor disk (x=0). See Theory in the documentation.

**Arguments**
- `ψ,r :: Float`: azimuthal, radial location
- `R :: Float`: rotor radius
- `χ :: Float`: skew angle
- `Θ :: Float`: rotor tilt
- `λ :: Float`: tip speed ratio
- `CT :: Float`: expected/approximated rotor thrust coefficient
- `n :: Int`: number of quadrature points

**Returns**
- `ϵx,ϵψ,ϵr::Float`: epsilon factors
"""
function epsilons(ψ, r, R, χ, Θ, λ, CT ; n=10000)

    # integration stuff
    nodes, weights = gausslegendre_0_2π(n)

    # prealloc
    k_u = zeros(3, n)
    I = zeros(3)
    Iff = zeros(3)

    # compute
    ϵx,ϵψ,ϵr = epsilons!(ψ, r, R, χ, Θ, λ, CT, nodes, weights, k_u, I, Iff )

    return ϵx,ϵψ,ϵr
end


"""
    epsilons!(ψ, r, R, χ, Θ, λ, CT, no, we, k_u, I, Iff ; n=10000)

    See `epsilons`.

**Arguments**
- `no,we :: Array{Float,1}`: nodes and weights for the Gauss-Legendre quadrature
- `k_u, I, Iff :: Array{Float,1}`: buffer of size 3
"""
function epsilons!(ψ, r, R, χ, Θ, λ, CT, no, we, k_u, I, Iff )

    Γr, γl, _ = gammas(r, R, χ, λ, CT )

    I .= 0
    Iff .= 0
    eval_u!(I, ψ, r, R, χ, k_u, no, we; γt=1.0, γl, Γr )
    eval_ur_ff!(Iff, ψ, r, R, χ, Θ; Γr )
    
    ϵx = .5 / I[1]
    ϵψ = .5 * Iff[2] / I[2]
    ϵr = I[3] / I[1]

    return ϵx,ϵψ,ϵr
end