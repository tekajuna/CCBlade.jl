using FastGaussQuadrature
using LinearAlgebra


# Most of the formulas in this file come from Branlard2016
# - Be very careful here: he defines the psi angle with respect to the horizontal (i.e. not the standard definition of psi)
# - all induced velocities here are for a unit gamma_t = uz0*2
# - all inplace functions are additive

## ----------------- velocity induced by the tangential vorticity in the outer region -----------------

#θ ; vect - integration variable
#ψ,r,z ; scalars
#m,R ; scalars
#k_u ; [3 x size(θ)] - output
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

    k_u[1,:] = 2 .* (apr .* sqrt.(c) .+ bpr .* sqrt.(a)) ./ den
    k_u[2,:] = 2 .* (apψ .* sqrt.(c) .+ bpψ .* sqrt.(a)) ./ den
    k_u[3,:] = 2 .* (apz .* sqrt.(c) .+ bpz .* sqrt.(a)) ./ den
end


#additive to ut
function eval_ut!(ut, rr, ψψ, zz, χ, R, k_u, no, we)

    m = tan(χ)  # x = m*z

    integ_ut!(k_u,no,rr,ψψ,zz,m,R)
    
    for j = 1:3
        ut[j] += dot( we, k_u[j,:] ) * 0.25 # quadrature with interval [0,2pi]: *pi, factor: 1/(4pi)
    end
end

# convenient function, mainly for plotting
#ψ,r,z ; vectors/scalars
#χ,R ; scalars
#ur ; [2 x size(r) x size(ψ) x size(z)] - output
#n ; number of quadrature points
function eval_ut(rr, ψψ, zz, χ, R; n=10000)

    # integration stuff
    nodes, weights = gausslegendre( n );
    nodes .= 2 .* pi .* .5 * (nodes .+ 1)

    # params
    # m = tan(χ)  # x = m*z
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

        # integ_ut!(k_u,nodes,rr[ir],ψψ[iψ],zz[iz],m,R)
    
        # for j = 1:3
        #     ut[j,i+1] = dot( weights, k_u[j,:] ) * pi # quadrature with interval [0,2pi]
        # end
        eval_ut!( view(ut,:,i+1), rr[ir],ψψ[iψ],zz[iz], χ, R, k_u, nodes, weights)
    end
    # ut .*= 1 / (4*pi) 

    ut = reshape(ut,(3,nr,nψ,nz))
    
    return ut
end

# -- approximated versions --

# Mind the 1/2 factor so that the restuls has to be multiplied by gamma_t (and uz0 = gamma_t/2)

#neglects radial and tangential components
function eval_ut_Coleman!(ut_approx, rr, ψψ, zz, χ, R)
    Kzt_approx = rr./R .* tan(χ/2)
    ut_approx[3,:] .+= ( 1 .+ Kzt_approx .* cos(ψψ))
end

#neglects radial and tangential components
function eval_ut_PittPeters!(ut_approx, rr, ψψ, zz, χ, R)
    Kzt_approx = rr./R .* tan(χ/2) .* (15.0*pi/32.)
    ut_approx[3,:] .+= ( 1 .+ Kzt_approx .* cos(ψψ)) *.5
end

function eval_ut_Branlard!(ut_approx, rr, ψψ, zz, χ, R)

    Kzt_approx = rr./R .* tan(χ/2) #another approx exists, only valid on psi=0,z=0
    ut_approx[3,:] .+= ( 1 .+ Kzt_approx .* cos(ψψ)) *.5
    Ft = Kzt_approx ./2 ./ tan(χ/2)

    # Kξt_approx = Kzt_approx ./ sin(χ)
    # Kxt_approx = Kξt_approx .* cos(χ)
    # uxt_approx =  (tan(χ/2) .- Kxt_approx .* cos(ψ))
    # uyt_approx = -uz0 .* Ft .* sec(χ/2)^2 .* sin(ψ)

    ut_approx[1,:] .+= tan(χ/2) .* cos(ψψ) .* ut_approx[3,:] .- Ft .* sec(χ/2)^2 *.5
    ut_approx[2,:] .+= -tan(χ/2) .* sin(ψψ) .* ut_approx[3,:]
end



## ----------------- velocity induced by the root vorticity -----------------

#additive to ut
function eval_ur!(ur, rr, ψψ, zz, χ, R)

    ur[1] += cos(χ) ./ (4*pi .* rr .* (1.0 .- cos.(ψψ) .* sin(χ)) )
    ur[3] += sin(χ) .* sin.(ψψ) ./ (4 .*pi .* rr .* (1.0 .- cos.(ψψ) .* sin(χ)) )

end

#ψ,r,z ; vectors/scalars
#χ,R ; scalars
#ur ; [2 x size(r) x size(ψ) x size(z)] -  tangential component, and axial component (no induced radial vel)
### ==== CAUTION =====
#  This is obviously wrong since it misses the dependence in z!
#  Plus, it seems like this considers an infinite filament, not a semi-infinite filament !! 
#     For the semi-infinite filament with chi = 0, the tg vel should be 1/2 Gamma/4/pi/r
#     Intuition: must be a (1+tan(r/z))/2 with a cos(chi)*(1-cos(psi)) somewhere
# ======================
#isn't it r/R that we need to use?
function eval_ur(rr, ψψ, zz, χ, R)
    nr = length(rr)
    nψ = length(ψψ)
    nz = length(zz)

    ur = zeros(3, nr*nψ*nz)

    for i = 0 : nr*nψ*nz - 1
        ir = mod(i,nr) +1
        iψ = mod(floor(Int64,i/nr),nψ) +1
        iz = floor(Int64,i/(nr*nψ)) +1

        # ur[1,i+1] = cos(χ) ./ (4*pi .* rr[ir] .* (1.0 .- cos.(ψψ[iψ]) .* sin(χ)) )
        # ur[3,i+1] = sin(χ) .* sin.(ψψ[iψ]) ./ (4 .*pi .* rr[ir] .* (1.0 .- cos.(ψψ[iψ]) .* sin(χ)) )
        eval_ur!( view(ur,:,i+1) , rr[ir], ψψ[iψ], zz[iz], χ, R)
    end
    ur = reshape(ur,(3,nr,nψ,nz))

    return ur
end


#TODO: add vel induced by longitundinal vorticity, u_l ?


## ----------------- total induced velocity -----------------

#all induced velocities, sum of tangential component in the tip wake and axial root vortex
# We neglext the axial component of the tip wake.
function eval_u!(u, rr, ψψ, zz, χ, R, k_u, no, we)

    eval_ur!( u, rr, ψψ, zz, χ, R)
    eval_ut!( u, rr, ψψ, zz, χ, R, k_u, no, we)

end


## ----------------- epsilon functions for the BEM -----------------
