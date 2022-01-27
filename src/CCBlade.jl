#=
Author: Andrew Ning

A general blade element momentum (BEM) method for propellers/fans and turbines.

Some unique features:
- a simple yet very robust solution method ideal for use with optimization
- designed for compatibility with algorithmic differentiation tools
- allows for arbitrary inflow conditions, including reversed flow, hover, etc.
- convenience methods for common wind turbine inflow scenarios

=#

module CCBlade

import FLOWMath

export Rotor, Section, OperatingPoint, Outputs
export simple_op, windturbine_op, flexturbine_op
export solve, thrusttorque, nondim, localLoadsToBladeFrame


include("airfoils.jl")  # all the code related to airfoil data
include("integrals.jl") # the code related to wake models and integrals to compute induced velocities

# --------- structs -------------

"""
    Rotor(Rhub, Rtip, B; precone=0.0, tilt=0.0, turbine=false, 
        mach=nothing, re=nothing, rotation=nothing, tip=PrandtlTipHub(), wakeCyl=false)

Parameters defining the rotor (apply to all sections).  

**Arguments**
- `Rhub::Float64`: hub radius (along blade length)
- `Rtip::Float64`: tip radius (along blade length)
- `B::Int64`: number of blades
- `precone::Float64`: precone angle
- 
- `turbine::Bool`: true if using wind turbine conventions
- `mach::MachCorrection`: correction method for Mach number
- `re::ReCorrection`: correction method for Reynolds number
- `rotation::RotationCorrection`: correction method for blade rotation
- `tip::TipCorrection`: correction method for hub/tip loss
- `wakeCyl::Bool`: true to use the cylinder wake model
"""
struct Rotor{TF, TI, TB, 
        T1 <: Union{Nothing, MachCorrection}, T2 <: Union{Nothing, ReCorrection}, 
        T3 <: Union{Nothing, RotationCorrection}, T4 <: Union{Nothing, TipCorrection},
        TV1, TV2}
    Rhub::TF
    Rtip::TF
    B::TI
    precone::TF
    tilt::TF
    turbine::TB
    mach::T1
    re::T2
    rotation::T3
    tip::T4
    wakeCyl::TB
    no::TV1
    we::TV1
    k_u::TV2
    I::TV1
    Iff::TV1
end

# convenience constructor with keyword parameters
function Rotor(Rhub, Rtip, B; precone=0.0, tilt=0.0, turbine=false, mach=nothing, re=nothing, rotation=nothing, tip=PrandtlTipHub(), wakeCyl=false
    )
    
    #wke integration nodes
    n = 10000
    no, we = gausslegendre_0_2π(n)

    # prealloc
    k_u = zeros(3, n)
    I = zeros(3)
    Iff = zeros(3)

    return Rotor(Rhub, Rtip, B, precone, tilt, turbine, mach, re, rotation, tip, wakeCyl, no, we, k_u, I, Iff)
end


"""
    Section(r, chord, theta, af)
    Section(chord, theta, af, x_az, y_az, z_az, coning, sweep)
    Section(chord, theta, af, precone, x_def, y_def, z_def, coning, sweep)

Define sectional properties for one station along rotor
    
**Arguments**
- `r::Float64`: radial location along blade, measured in the blade root frame
- `chord::Float64`: corresponding local chord length
- `theta::Float64`: corresponding twist angle (radians)
- `af::Function or AFType`: if function form is: `cl, cd = af(alpha, Re, Mach)`
- `x_az::Float64`: x location when blade is deflected in the azimuthal frame (m)
- `y_az::Float64`: y location when blade is deflected in the azimuthal frame (m)
- `z_az::Float64`: z=radial location when blade is deflected (m)
or 
- `precone::Float`: precone angle (rad)
- `x_def::Float64`: x location when blade is deflected, in the blade root frame (m)
- `y_def::Float64`: y location when blade is deflected, in the blade root frame (m)
- `z_def::Float64`: z location when blade is deflected, in the blade root frame (m)
- `coning::Float64`: local coning angle w.r.t. blade root at that location, in blade root frame (rad), i.e. lcon>0 for a blade deflected downwind
- `sweep::Float64`: local sweep angle w.r.t. blade root at that location (rad)
"""
Section

struct Section{TF1, TF2, TF3, TAF, TF4, TF5}
    r::TF1  # different types b.c. of dual numbers.  often r is fixed, while chord/theta vary.
    chord::TF2
    theta::TF3
    af::TAF
    x_az::TF4 #could discuss if it is better to give this in blade root...: definitely more convenient from structure point of view
    y_az::TF4
    z_az::TF4
    coning::TF5 # coning contains the local cone in blade root frame (see thrusttorque routine)
    sweep::TF5
    #replace with curv object?? TFc<:Curvature
end  

# convenience constructor for undeflected blade. Avoiding kwargs so that it can still be broadcasted.
Section(r, chord, theta, af, precone=zero(r)
    ) = Section(r, chord, theta, af,-r*sin(precone), zero(r), r*cos(precone), zero(r), zero(r))

Section(chord, theta, af, precone, x_def, y_def, z_def, coning, sweep, coordsInAzimuthalFrame=false) = begin
    if coordsInAzimuthalFrame
        x_az = x_def
        y_az = y_def
        z_az = z_def
    else
        x_az = x_def * cos(precone) - z_def * sin(precone)
        y_az = y_def
        z_az = x_def * sin(precone) + z_def * cos(precone) 
    end
    r = z_def
    Section(r, chord, theta, af, x_az, y_az, z_az, coning, sweep)
end

# convenience function to access fields within an array of structs
function Base.getproperty(obj::Vector{Section{TF1, TF2, TF3, TAF, TF4}}, sym::Symbol) where {TF1, TF2, TF3, TAF, TF4}
    return getfield.(obj, sym)
end # This is not always type stable b/c we don't know if the return type will be float or af function.


"""
    OperatingPoint(Vx, Vy, Omega, rho)
    OperatingPoint(Vx, Vy, Vz, Omega, rho; pitch=0.0, mu=1.0, asound=1.0)

Operation point for a rotor.  
The x direction is the axial direction, and y direction is the tangential direction in the rotor plane.  
See Documentation for more detail on coordinate systems. #TODO: expand/rewrite
`Vx` and `Vy` vary radially at same locations as `r` in the rotor definition.

**Arguments**
- `Vx::Float64`: velocity in x-direction along blade
- `Vy::Float64`: velocity in y-direction along blade
- `Vz::Float64`: velocity in z-direction along blade
- `Vhub::Float64`: upstream velocity
- `Omega::Float64`: rotation rate
- `rho::Float64`: fluid density
- `pitch::Float64`: pitch angle (radians)
- `mu::Float64`: fluid dynamic viscosity (unused if Re not included in airfoil data)
- `asound::Float64`: fluid speed of sound (unused if Mach not included in airfoil data)
"""
struct OperatingPoint{TF1, TF2, TF3, TF4, TF5, TF6, TF7}
    Vx::TF1
    Vy::TF1
    Vz::TF1
    Vhub::TF1
    Omega::TF1
    rho::TF2  # different type to accomodate ReverseDiff
    pitch::TF3  
    yaw::TF4
    azimuth::TF5
    mu::TF6
    asound::TF7
    #extend this with gauss legendre stuff? or put a wake object in rotor
end

# convenience constructor when Re and Mach are not used.
OperatingPoint(Vx, Vy, Vz, Omega, rho; pitch=zero(rho), yaw=zero(rho), azimuth=zero(rho), mu=one(rho), asound=one(rho)) = OperatingPoint(Vx, Vy, Vz, Vx, Omega, rho, pitch, yaw, azimuth, mu, asound ) 

OperatingPoint(Vx, Vy, Omega, rho; kwargs...) = OperatingPoint(Vx, Vy, zero(Vy), Omega, rho; kwargs...) 

# convenience function to access fields within an array of structs
function Base.getproperty(obj::Vector{OperatingPoint{TF1, TF2, TF3, TF4, TF5}}, sym::Symbol) where {TF1, TF2, TF3, TF4, TF5}
    return getfield.(obj, sym)
end


"""
    Outputs(Np, Tp, a, ap, u, v, phi, alpha, W, cl, cd, cn, ct, F, G)

Outputs from the BEM solver along the radius.

**Arguments**
- `Np::Float64`: normal force per unit length
- `Tp::Float64`: tangential force per unit length
- `a::Float64`: axial induction factor
- `ap::Float64`: tangential induction factor
- `u::Float64`: axial induced velocity
- `v::Float64`: tangential induced velocity
- `phi::Float64`: inflow angle
- `alpha::Float64`: angle of attack
- `W::Float64`: inflow velocity
- `cl::Float64`: lift coefficient
- `cd::Float64`: drag coefficient
- `cn::Float64`: normal force coefficient
- `ct::Float64`: tangential force coefficient
- `F::Float64`: hub/tip loss correction
- `G::Float64`: effective hub/tip loss correction for induced velocities: `u = Vx * a * G, v = Vy * ap * G`
"""
struct Outputs{TF}
    Np::TF
    Tp::TF
    a::TF
    ap::TF
    u::TF
    v::TF
    phi::TF
    alpha::TF
    W::TF
    cl::TF
    cd::TF
    cn::TF
    ct::TF
    F::TF
    G::TF
end

# convenience constructor to initialize
Outputs() = Outputs(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

# convenience function to access fields within an array of structs
function Base.getproperty(obj::Vector{Outputs{TF}}, sym::Symbol) where TF
    return getfield.(obj, sym)
end

# -------------------------------




# ------------ BEM core ------------------

"""
(private) residual function
"""
function residual(phi, rotor, section, op; trial::Int64=1)

    # unpack inputs
    r = section.r
    # x_az = section.x_az
    y_az = section.y_az
    z_az = section.z_az

    chord = section.chord
    theta = section.theta
    af = section.af

    Rhub = rotor.Rhub
    Rtip = rotor.Rtip
    B = rotor.B
    
    # external velocities, expressed in the rotor plane c.s.
    Vx = op.Vx
    Vy = op.Vy
    Vz = op.Vz
    #TODO:  consider adding displacement velocity here

    #unpack
    rho = op.rho
    yaw = op.yaw
    tilt = rotor.tilt
    azimuth = op.azimuth
    precone = rotor.precone
    pitch = op.pitch
    cone = section.coning
    # swp = section.sweep
    
    # constants
    sigma_p = B*chord/(2.0*pi*z_az) #assuming all blades have the same loading for the computation of the induction
    sphi = sin(phi)
    cphi = cos(phi)
    ca = sin(azimuth)
    sa = cos(azimuth)
    sc0 = sin(-precone)
    cc0 = cos(-precone)
    sc = sin(cone)
    cc = cos(cone)
    sp = sin(pitch)
    cp = cos(pitch)
    cY = cos(yaw)
    cT = cos(rotor.tilt)
    λ  =  Rtip * op.Omega / op.Vhub 
    
    # blade coordinates in the hub frame ≡ cylinder wake frame, in cylindrical coordinates
    #  because we neglect the tilt angle in the computation of the epsilon factors
    y_hu = y_az*ca - z_az*sa
    z_hu = y_az*sa + z_az*ca

    ψ_Br = atan(z_hu,y_hu) #CAUTION: this is the psi as defined in Branlart2016 (not the same as the psi of the WT)
    r_Br   = sqrt(y_hu^2+z_hu^2)

    # trigonometric factors
    # p1f = (cc*cp*cc0 - sc*sc0)
    # p2f = (cc*cp*sc0 + sc*cc0)
    # p3 = cp 
    # p4 = cc0*sp
    # p5 = cc*sp
    # p6 = sc0*sp
    #APPROX PITCH:
    p1 = cos(-precone+cone)
    p2 = sin(-precone+cone)
    p3 = 1. #->sweep through shearing assumption...?
    p4 = 0.
    p5 = 0.
    p6 = 0.

    # epsilon ratios
    if rotor.wakeCyl
        a_tmp = (rotor.turbine  ? -.33 : +.1)     #arbitrary choice
        CT =  4*a_tmp*(1+a_tmp) #TENTATIVE APPROXIMATE VALUE (Branlart2016)
        ϵx,ϵψ,ϵr = epsilons!(ψ_Br, r_Br, Rtip, yaw, 0.0, λ, CT, rotor.no, rotor.we, rotor.k_u, rotor.I, rotor.Iff )
        #TODO: double check yaw convention: For positive chi, Branlart's x+ goes in the same direction as the wake, that is in my hub's y+. But I might use the wrong yaw convention.

    else   
        #reset values to match previous CCblade version
        ϵx = 1.
        ϵψ = 1.
        ϵr = 0.

        cT = 1.0
        # cY = 0. #it was so 

        sigma_p = B*chord/(2.0*pi*r)  #in the original CCBlade, r does not account for precone. It does in the new version, which leads to a small difference on sigma_p.
        #---> this is questionable, theoretically z_az should always be used in this sigma_p.

        Vx = (p1 * Vx - p2 * Vz) #velocity in the blade root frame
        p1 = 1. #forcing the forces to remain in the blade root frame
        p2 = 0.
    end

    # angle of attack
    # alpha = theta - phi 
    alpha = theta + pitch - phi #APPROX PITCH

    # Reynolds/Mach number
    W0 = sqrt(Vx^2 + Vy^2)  # ignoring induction, which is generally a very minor difference and only affects Reynolds/Mach number
    Re = rho * W0 * chord / op.mu
    Mach = W0/op.asound  # also ignoring induction

    # airfoil cl/cd
    if rotor.turbine
        cl, cd = afeval(af, -alpha, Re, Mach)
        cl *= -1
    else
        cl, cd = afeval(af, alpha, Re, Mach)
    end

    # airfoil corrections
    if !isnothing(rotor.re)
        cl, cd = re_correction(rotor.re, cl, cd, Re)
    end
    if !isnothing(rotor.mach)
        cl, cd = mach_correction(rotor.mach, cl, cd, Mach)
    end
    if !isnothing(rotor.rotation)
        cl, cd = rotation_correction(rotor.rotation, cl, cd, chord/r, r/Rtip, Vy/Vx*Rtip/r, alpha, phi)
    end

    # resolve into normal and tangential forces
    cn = cl*cphi - cd*sphi
    ct = cl*sphi + cd*cphi


    # changing frame from local airfoil to azm
    # TODO: use change of frame functions from supplemental?
    #   otherwise, A' in my notes, with sweep = 0 assumed
    cx = p1 * cn - p4 * ct
    cy = p5 * cn + p3 * ct #note: could neglect the first term as we do for the momentum
    cz =-p2 * cn + p6 * ct  #note: could also neglect the last term, for the same reason

    # hub/tip loss
    F = 1.0
    if !isnothing(rotor.tip)
        F = tip_correction(rotor.tip, r, Rhub, Rtip, phi, B)   
    end

    # sec parameters
    k  = cx*sigma_p/(4.0*F*sphi*sphi) /(ϵx*cY*cT)
    kp = cy*sigma_p/(4.0*F*sphi*cphi) /(ϵψ)

    # --- solve for induced velocities ------
    if isapprox(Vx, 0.0, atol=1e-6)

        #no yaw model for this yet:
        k  = cn*sigma_p/(4.0*F*sphi*sphi)
        kp = ct*sigma_p/(4.0*F*sphi*cphi)

        u = sign(phi)*kp*cn/ct*Vy
        v = zero(phi)
        Un = Vx + u
        Ut = Vy - v 
        a = zero(phi)
        ap = zero(phi)
        R = sign(phi) - k

    elseif isapprox(Vy, 0.0, atol=1e-6)

        #no yaw model for this yet:
        k  = cn*sigma_p/(4.0*F*sphi*sphi)
        kp = ct*sigma_p/(4.0*F*sphi*cphi)
        
        u = zero(phi)
        v = k*ct/cn*abs(Vx)
        Un = Vx + u
        Ut = Vy - v
        a = zero(phi)
        ap = zero(phi)
        R = sign(Vx) + kp
    
    else

        if phi < 0
            k *= -1
        end

        if isapprox(k, 1.0, atol=1e-6)  # state corresopnds to Vx=0, return any nonzero residual
            return 1.0, Outputs()
        end

        b1 = (p1 * Vx - p2 * Vz)
        b2 = (p1 - p2 * ϵr ) * Vx
        b3 = k / Vx^2

        if k >= -2.0/3  # momentum region
            # a = k/(1 - k)

            radical = 4* (b1^2 - b1*b2) *b3 +1
            if radical < 0.
                if radical > 1e-6
                    print("Warning: radical <0 :" )
                    println(radical)
                end
                radical = 0
            end
            a_  = ((- 2*b1*b2*b3 + 1) + sqrt(radical)) / (2 * (b2^2*b3 - 1) ) #always discard this one?
            a   = ((- 2*b1*b2*b3 + 1) - sqrt(radical)) / (2 * (b2^2*b3 - 1) )
            

        else  # empirical region
            g1 = F*(2*k - 1) + 10.0/9
            g2 = F*(F - 2*k - 4.0/3)
            g3 = 2*F*(1 - k) - 25.0/9

            if isapprox(g3, 0.0, atol=1e-6)  # avoid singularity
                a = 1.0/(2.0*sqrt(g2)) - 1
            else
                a = (g1 + sqrt(g2)) / g3
            end
        end

        u = a * Vx
        Un = b1 + b2 * a

        # -------- tangential induction ----------
        if Vx < 0
            kp *= -1
        end

        if isapprox(kp, -1.0, atol=1e-6)  # state corresopnds to Vy=0, return any nonzero residual
            return 1.0, Outputs()
        end

        kplhs = Un * p3 / ( Vx * (1+a) ) * kp

        # ap = kp/(1 + kp)
        ap = kplhs/(1 + kplhs)

        # println("ap")
        # # println(kplhs)
        # println(ap)
        # println(kp/(1 + kp))
        # println("-")

        v = ap * Vy
        Ut = p3 * Vy * (1 - ap)


        # ------- residual function -------------
        # R = sin(phi)/(1 + a) - Vx/Vy*cos(phi)/(1 - ap)
        #R = sin(phi)/Un - cos(phi)/Ut

        if trial == 1
            R = sin(phi)/Un - cos(phi)/Ut # Fails later portions --with corrected TSR it now is opposite
        end 
        if trial == 2
            R = cos(phi)/Ut - sin(phi)/Un
        end
        if trial == 3
            R = sin(phi)- cos(phi) * Un/Ut #Matches later portions --with corrected TSR it now is opposite
        end
        if trial == 4
            R = sin(phi)*(p3 * Vy )*(1 - ap)  - cos(phi) * (b1 + b2 * a)
        end
        # if trial == 5
        sin(phi)/(Vx + u) == cos(phi)/(Vy - v )
        # end

        # println(R)
        # # println(sin(phi)/(1 + a) - Vx/Vy*cos(phi)/(1 - ap))
    end


    # ------- loads ---------
    # W = sqrt((Vx + u)^2 + (Vy - v)^2)
    W = sqrt(Un^2 + Ut^2)

    #these are the forces in azm c.s.
    Np = cx*0.5*rho*W^2*chord
    Tp = cy*0.5*rho*W^2*chord
    #Zp = ... #Should add the spanwise force for completeness

    # The BEM methodology applies hub/tip losses to the loads rather than to the velocities.  
    # This is the most common way to implement a BEM, but it means that the raw velocities are misleading 
    # as they do not contain any hub/tip loss corrections.
    # To fix this we compute the effective hub/tip losses that would produce the same thrust/torque.
    # In other words:
    # CT = 4 a (1 + a) F = 4 a G (1 + a G)\n
    # This is solved for G, then multiplied against the wake velocities.
    
    if isapprox(Vx, 0.0, atol=1e-6)
        G = sqrt(F)
    elseif isapprox(Vy, 0.0, atol=1e-6)
        G = F
    else
        G = (-1.0 + sqrt(1.0 + 4*a*(1.0 + a)*F))/(2*a)
    end
    u *= G
    v *= G

    """
    if rotor.turbine
        return R, Outputs(-Np, -Tp, -a, -ap, -u, -v, phi, -alpha, W, -cl, cd, -cn, -ct, F, G)
    else
        return R, Outputs(Np, Tp, a, ap, u, v, phi, alpha, W, cl, cd, cn, ct, F, G)
    end
    """
    if rotor.turbine
        return R, Outputs(-Np, -Tp, -a, -ap, -u, -v, phi, -alpha, W, -cl, cd, -cn, -ct, F, G), Un, Ut, ap, Vy, b1,b2,p3,kplhs, [Un,Ut,ap,Vy,b1,b2,p3,kplhs,kp,ϵx,ϵψ,ϵr,u,a,Vx,Vy]
    else                                                                                                                        #  2o 2o       opp   hmm            d      off                       
        return R, Outputs(Np, Tp, a, ap, u, v, phi, alpha, W, cl, cd, cn, ct, F, G)#, Un, Ut, ap, Vy, b1,b2,p3,kplhs
    end

end




"""
(private) Find a bracket for the root closest to xmin by subdividing
interval (xmin, xmax) into n intervals.

Returns found, xl, xu.
If found = true a bracket was found between (xl, xu)
"""
function firstbracket(f, xmin, xmax, n, backwardsearch=false)

    xvec = range(xmin, xmax, length=n)
    if backwardsearch  # start from xmax and work backwards
        xvec = reverse(xvec)
    end

    fprev = f(xvec[1])
    for i = 2:n
        fnext = f(xvec[i])
        if fprev*fnext < 0  # bracket found
            if backwardsearch
                return true, xvec[i], xvec[i-1]
            else
                return true, xvec[i-1], xvec[i]
            end
        end
        fprev = fnext
    end

    return false, 0.0, 0.0
end


"""
    solve(rotor, section, op)

Solve the BEM equations for given rotor geometry and operating point.

**Arguments**
- `rotor::Rotor`: rotor properties
- `section::Section`: section properties
- `op::OperatingPoint`: operating point

**Returns**
- `outputs::Outputs`: BEM output data including loads, induction factors, etc.
"""
function solve(rotor, section, op; trial=1)

    # error handling
    if typeof(section) <: Vector
        error("You passed in an vector for section, but this funciton does not accept an vector.\nProbably you intended to use broadcasting (notice the dot): solve.(Ref(rotor), sections, ops)")
    end

    # check if we are at hub/tip
    if isapprox(section.r, rotor.Rhub, atol=1e-6) || isapprox(section.r, rotor.Rtip, atol=1e-6)
        return Outputs()  # no loads at hub/tip
    end

    # parameters
    npts = 10  # number of discretization points to find bracket in residual solve

    # unpack
    Vx = op.Vx
    Vy = op.Vy
    theta = section.theta + op.pitch

    # ---- determine quadrants based on case -----
    Vx_is_zero = isapprox(Vx, 0.0, atol=1e-6)
    Vy_is_zero = isapprox(Vy, 0.0, atol=1e-6)

    # quadrants
    epsilon = 1e-6
    q1 = [epsilon, pi/2]
    q2 = [-pi/2, -epsilon]
    q3 = [pi/2, pi-epsilon]
    q4 = [-pi+epsilon, -pi/2]

    if Vx_is_zero && Vy_is_zero
        return Outputs()

    elseif Vx_is_zero

        startfrom90 = false  # start bracket at 0 deg.

        if Vy > 0 && theta > 0
            order = (q1, q2)
        elseif Vy > 0 && theta < 0
            order = (q2, q1)
        elseif Vy < 0 && theta > 0
            order = (q3, q4)
        else  # Vy < 0 && theta < 0
            order = (q4, q3)
        end

    elseif Vy_is_zero

        startfrom90 = true  # start bracket search from 90 deg

        if Vx > 0 && abs(theta) < pi/2
            order = (q1, q3)
        elseif Vx < 0 && abs(theta) < pi/2
            order = (q2, q4)
        elseif Vx > 0 && abs(theta) > pi/2
            order = (q3, q1)
        else  # Vx < 0 && abs(theta) > pi/2
            order = (q4, q2)
        end

    else  # normal case

        startfrom90 = false

        if Vx > 0 && Vy > 0
            order = (q1, q2, q3, q4)
        elseif Vx < 0 && Vy > 0
            order = (q2, q1, q4, q3)
        elseif Vx > 0 && Vy < 0
            order = (q3, q4, q1, q2)
        else  # Vx[i] < 0 && Vy[i] < 0
            order = (q4, q3, q2, q1)
        end

    end

        

    # ----- solve residual function ------

    # # wrapper to residual function to accomodate format required by fzero
    R(phi) = residual(phi, rotor, section, op)[1]

    success = false
    for j = 1:length(order)  # quadrant orders.  In most cases it should find root in first quadrant searched.
        phimin, phimax = order[j]

        # check to see if it would be faster to reverse the bracket search direction
        backwardsearch = false
        if !startfrom90
            if phimin == -pi/2 || phimax == -pi/2  # q2 or q4
                backwardsearch = true
            end
        else
            if phimax == pi/2  # q1
                backwardsearch = true
            end
        end
        
        # force to dual numbers if necessary
        phimin = phimin*one(section.chord)
        phimax = phimax*one(section.chord)

        # find bracket
        success, phiL, phiU = firstbracket(R, phimin, phimax, npts, backwardsearch)

        # once bracket is found, solve root finding problem and compute loads
        if success
            phistar, _ = FLOWMath.brent(R, phiL, phiU)
            _, outputs = residual(phistar, rotor, section, op)
            return outputs
        end    
    end    

    # it shouldn't get to this point.  if it does it means no solution was found
    # it will return empty outputs
    # alternatively, one could increase npts and try again
    
    @warn "Invalid data (likely) for this section.  Zero loading assumed."
    return Outputs()
end



# ------------ inflow ------------------



"""
    simple_op(Vinf, Omega, r, rho; pitch=0.0, mu=1.0, asound=1.0, precone=0.0)

Uniform inflow through rotor.  Returns an OperatingPoint object.

**Arguments**
- `Vinf::Float`: freestream speed (m/s)
- `Omega::Float`: rotation speed (rad/s)
- `r::Float`: radial location where inflow is computed (m)
- `rho::Float`: air density (kg/m^3)
- `pitch::Float`: pitch angle (rad)
- `mu::Float`: air viscosity (Pa * s)
- `asounnd::Float`: air speed of sound (m/s)
- `precone::Float`: precone angle (rad)
"""
function simple_op(Vinf, Omega, r, rho; pitch=zero(rho), mu=one(rho), asound=one(rho), precone=zero(Vinf))

    # error handling
    if typeof(r) <: Vector
        error("You passed in an vector for r, but this function does not accept an vector.\nProbably you intended to use broadcasting")
    end

    #result in the AZM c.s.
    Vx = Vinf
    Vy = Omega * r * cos(precone)

    return OperatingPoint(Vx, Vy, Omega, rho; pitch=pitch, mu=mu, asound=asound)

end


"""
    windturbine_op(Vhub, Omega, pitch, r, precone, yaw, tilt, azimuth, hubHt, shearExp, rho, mu=1.0, asound=1.0)

Compute relative wind velocity components along a straight blade accounting for inflow conditions
and orientation of turbine.  See Documentation for angle definitions.

**Arguments**
- `Vhub::Float64`: freestream speed at hub (m/s)
- `Omega::Float64`: rotation speed (rad/s)
- `pitch::Float64`: pitch angle (rad)
- `r::Float64`: radial location along the blade axis (in the blade root frame) where inflow is computed (m)
- `precone::Float64`: precone angle (rad)
- `yaw::Float64`: yaw angle (rad)
- `tilt::Float64`: tilt angle (rad)
- `azimuth::Float64`: azimuth angle to evaluate at (rad)
- `hubHt::Float64`: hub height (m) - used for shear
- `shearExp::Float64`: power law shear exponent
- `rho::Float64`: air density (kg/m^3)
- `mu::Float64`: air viscosity (Pa * s)
- `asound::Float64`: air speed of sound (m/s)
"""
function windturbine_op(Vhub, Omega, pitch, r, precone, yaw, tilt, azimuth, hubHt, shearExp, rho, mu=one(rho), asound=one(rho))

    xaz = r * -sin(precone)  #AZM c.s.
    yaz = zero(r)
    zaz = r * cos(precone)
    lcon = zero(r)
    lswp = zero(r)
    return flexturbine_op(Vhub, Omega, pitch, xaz, yaz, zaz, yaw, tilt, azimuth, hubHt, shearExp, rho, mu, asound)
end

"""
    flexturbine_op(Vhub, Omega, pitch, r, precone, yaw, tilt, azimuth, hubHt, shearExp, rho, mu=1.0, asound=1.0)

Compute relative wind velocity components along blade accounting for inflow conditions,
 orientation of turbine and deflection of the blades.  See Documentation for angle definitions.

**Arguments**
- `Vhub::Float64`: freestream speed at hub (m/s)
- `Omega::Float64`: rotation speed (rad/s)
- `pitch::Float64`: pitch angle (rad)
- `x_az::Float64`: x location where inflow is computed in the azimuthal frame (m)
- `y_az::Float64`: y location where inflow is computed in the azimuthal frame (m)
- `z_az::Float64`: z location where inflow is computed in the azimuthal frame (m)
- `yaw::Float64`: turbine yaw angle (rad)
- `tilt::Float64`: turbine tilt angle (rad)
- `azimuth::Float64`: blade azimuth angle to evaluate at (rad)
- `hubHt::Float64`: hub height (m) - used for shear
- `shearExp::Float64`: power law shear exponent
- `rho::Float64`: air density (kg/m^3)
- `mu::Float64`: air viscosity (Pa * s)
- `asound::Float64`: air speed of sound (m/s)
or
- `section::Section`: a section struct
"""
flexturbine_op

function flexturbine_op(Vhub, Omega, pitch, x_az, y_az, z_az, yaw, tilt, azimuth, hubHt, shearExp, rho, mu=one(rho), asound=one(rho))
    
    sy = sin(yaw)
    cy = cos(yaw)
    st = sin(tilt)
    ct = cos(tilt)
    sa = sin(azimuth)
    ca = cos(azimuth)
    
    # get section heights in wind-aligned coordinate system
    heightFromHub = (y_az*sa + z_az*ca)*ct - x_az*st

    # velocity with shear
    V = Vhub*(1 + heightFromHub/hubHt)^shearExp

    # transform wind from inertial to AZM c.s.
    # (V,0,0) -> yaw -> tilt -> azimuth 
    # TODO: use supplemental?
    Vwind_xr = V * ( cy*ct )
    Vwind_yr = V * (cy*st*sa - sy*ca)
    Vwind_zr = V * (cy*st*ca + sy*sa) 

    # relative wind from rotation in blade root c.s. (i.e., Omega cross r in the azimuthal c.s., then tilted by precone angle)
    Vrot_xr = 0.0
    Vrot_yr = Omega*z_az
    Vrot_zr =-Omega*y_az
    
    # total velocity
    Vxr = Vwind_xr + Vrot_xr
    Vyr = Vwind_yr + Vrot_yr
    Vzr = Vwind_zr + Vrot_zr
    
    # operating point
    return OperatingPoint(Vxr, Vyr, Vzr, Vhub, Omega, rho, pitch, yaw, azimuth, mu, asound)

end

function flexturbine_op(Vhub, Omega, pitch, section::Section{TF1, TF2, TF3, TAF, TF4}, 
                        yaw, tilt, azimuth, hubHt, shearExp, rho, mu=one(rho), asound=one(rho)) where {TF1, TF2, TF3, TAF, TF4}
    x_az = section.x_az
    y_az = section.y_az
    z_az = section.z_az
    # lcon = section.coning #WARNING: depending on what we choose to store in Sections, this contains precone or not.
    # lswp = section.sweep
    return flexturbine_op(Vhub, Omega, pitch, x_az, y_az, z_az, yaw, tilt, azimuth, hubHt, shearExp, rho, mu, asound)
end
# -------------------------------------


# -------- convenience methods ------------

"""
    thrusttorque(rotor, sections, outputs::Vector{TO}) where TO

integrate the thrust/torque across the blade, 
including 0 loads at hub/tip, using a trapezoidal rule.

**Arguments**
- `rotor::Rotor`: rotor object
- `sections::Vector{Section}`: rotor object
- `outputs::Vector{Outputs}`: output data along blade

**Returns**
- `T::Float64`: thrust (along x-dir see Documentation).
- `Q::Float64`: torque (along x-dir see Documentation).
"""
#TODO: use buffers to avoid allocations?
function thrusttorque(rotor, sections, outputs::Vector{TO}) where TO

    # add hub/tip for complete integration.  loads go to zero at hub/tip.
    rvec = [s.r for s in sections]
    rfull = [rotor.Rhub; rvec; rotor.Rtip]

    # angles between the local deflected blade orientation and the rotor plane
    # confull = [0.0; [s.coning for s in sections]; 0.0] .+ rotor.precone #since conicity and precone are two successive transforms, I can just add up the angles
    confull = [0.0; [s.coning for s in sections]; 0.0] 
    swpfull = [0.0; [s.sweep for s in sections]; 0.0]

    # extract loading as force components in the rotor plane c.s.,
    fx, fy, fz = outputs.Np, -outputs.Tp, outputs.Np.*0.  #CAUTION: minus sign, see explanation in localLoadsToBladeFrame
    if !rotor.wakeCyl 
        #In the original CCBlade, the forces are computed in the blade root frame. Switching to Azm:
        fx .*= cos(rotor.precone)
    end

    # extend to the root and tip
    fxfull = [0.0; fx; 0.0] 
    fyfull = [0.0; fy; 0.0] 
    fzfull = [0.0; fz; 0.0]  #TODO: add Fz !!!

    # position of the blade sections in the azimuthal frame, as per the definition of Section. 
    x_az = [s.x_az for s in sections]
    y_az = [s.y_az for s in sections]
    z_az = [s.z_az for s in sections]

    # Section stores the local coning angle. From the AZM frame perspective, we need to add the precone (minus sign because sign convention are opposed)
    confull .-= rotor.precone

    # guessing the position of the root and the tip, in the azimuthal coordinates, after preconing and deflection
    x_az_tip = x_az[end] + sin(confull[end]) * (rotor.Rtip-rvec[end])
    y_az_tip = y_az[end] - sin(swpfull[end]) * cos(confull[end]) * (rotor.Rtip-rvec[end])
    z_az_tip = z_az[end] + cos(swpfull[end]) * cos(confull[end]) * (rotor.Rtip-rvec[end])
    x_az_root = x_az[1] + sin(confull[1]) * (rotor.Rhub-rvec[1])
    y_az_root = y_az[1] - sin(swpfull[1]) * cos(confull[1]) * (rotor.Rhub-rvec[1])
    z_az_root = z_az[1] + cos(swpfull[1]) * cos(confull[1]) * (rotor.Rhub-rvec[1])
    # assembling it:
    # xdeffull = [x_az_root ; x_az ; x_az_tip] #unused
    y_az_full = [y_az_root ; y_az ; y_az_tip]
    z_az_full = [z_az_root ; z_az ; z_az_tip]

    # due to the deflection, Np and Tp may contribute to the torque = r cross (fx fy fz)
    thrust = fxfull
    torque = y_az_full .* fzfull - z_az_full .* fyfull

    # integrate Thrust and Torque (trapezoidal)
    # Note: we do neglect the blade axial extension
    T = rotor.B * FLOWMath.trapz(rfull, thrust)
    Q = rotor.B * FLOWMath.trapz(rfull, torque)

    # println(rfull)
    # println(thrust)
    # println(torque)

    return T, Q
end


"""
    thrusttorque(rotor, sections, outputs::Matrix{TO}) where TO

Integrate the thrust/torque across the blade given an array of output data.
Generally used for azimuthal averaging of thrust/torque.
`outputs[i, j]` corresponds to `sections[i], azimuth[j]`.  Integrates across azimuth
"""
function thrusttorque(rotor, sections, outputs::Matrix{TO}) where TO

    T = 0.0
    Q = 0.0
    nr, naz = size(outputs)

    for j = 1:naz
        Tsub, Qsub = thrusttorque(rotor, sections, outputs[:, j])
        T += Tsub / naz
        Q += Qsub / naz
    end

    return T, Q
end

"""
    localLoadsToBladeFrame(rotor, sections, outputs::Vector{TO}) where TO

transform the loading along the blade to express it in force components (x,y,z)
in the blade root coordinate system (x towards suction side at root, y towards root leading edge, z towards the blade tip at root).
If `toRotorPlane=true`, the precone angle is also accounted for such that the output is expressed in the rotor plane.

**Arguments**
- `rotor::Rotor`: rotor object
- `sections::Vector{Section}`: rotor object
- `outputs::Vector{Outputs}`: output data along blade
- `toRotorPlane::Bool`: if true, the output are expressed in the rotor plane

**Returns**
- `fx::Vector{Float64}`: force along x-dir (see Documentation).
- `fy::Vector{Float64}`: force along y-dir (see Documentation).
- `fz::Vector{Float64}`: force along z-dir (see Documentation).
"""
#DEPRECATED. Use Andrew's utilisty functions instead.
function localLoadsToBladeFrame(rotor, sections, outputs::Vector{TO}; toRotorPlane=false) where TO

    # Caution: For props Ct is defined positive towards the trailing edge, 
    #           for turbines, it is positive towards the LE.
    #          Cn is always positive towards the suction side.
    #  Hence, we need to change the sign of Tp to match the definition of the blade/rotor frame.
    sig = -1.

    # loadings in their local deflected c.s.
    Np = outputs.Np
    Tp = outputs.Tp .* sig
    
    # prepare to add the precone to the coning=curvature if needed
    pcn = toRotorPlane ? rotor.precone : zero(rotor.precone)
    
    # angles between the local deflected blade orientation and the blade root c.s.
    con = [s.coning - pcn for s in sections] #coning positive bckwd
    swp = [s.sweep for s in sections]

    # switch loads from local deflected blade c.s. to azimuthal c.s. (that is, the blade root c.s. or rotor plane)
    fx = cos.(con).*Np .+ sin.(swp).*sin.(con).*Tp
    fy = cos.(swp).*Tp
    fz = -sin.(con).*Np .+ sin.(swp).*cos.(con).*Tp    

    return fx, fy, fz
end


"""
    nondim(T, Q, Vhub, Omega, rho, rotor, rotortype)

Nondimensionalize the outputs.

**Arguments**
- `T::Float64`: thrust (N)
- `Q::Float64`: torque (N-m)
- `Vhub::Float64`: hub speed used in turbine normalization (m/s)
- `Omega::Float64`: rotation speed used in propeller normalization (rad/s)
- `rho::Float64`: air density (kg/m^3)
- `rotor::Rotor`: rotor object
- `rotortype::String`: normalization type

**Returns**

`if rotortype == "windturbine"`
- `CP::Float64`: power coefficient
- `CT::Float64`: thrust coefficient
- `CQ::Float64`: torque coefficient

`if rotortype == "propeller"`
- `eff::Float64`: efficiency
- `CT::Float64`: thrust coefficient
- `CQ::Float64`: torque coefficient

`if rotortype == "helicopter"`
- `FM::Float64`: figure of merit
- `CT::Float64`: thrust coefficient
- `CQ or CP::Float64`: torque/power coefficient (they are identical)
"""
function nondim(T, Q, Vhub, Omega, rho, rotor, rotortype)

    P = Q * Omega
    Rp = rotor.Rtip*cos(rotor.precone)

    if rotortype == "windturbine"  # wind turbine normalizations

        q = 0.5 * rho * Vhub^2
        A = pi * Rp^2

        CP = P / (q * A * Vhub)
        CT = T / (q * A)
        CQ = Q / (q * Rp * A)

        return CP, CT, CQ

    elseif rotortype == "propeller"

        n = Omega/(2*pi)
        Dp = 2*Rp

        if T < 0
            eff = 0.0  # creating drag not thrust
        else
            eff = T*Vhub/P
        end
        CT = T / (rho * n^2 * Dp^4)
        CQ = Q / (rho * n^2 * Dp^5)

        return eff, CT, CQ

    elseif rotortype == "helicopter"

        A = pi * Rp^2

        CT = T / (rho * A * (Omega*Rp)^2)
        CP = P / (rho * A * (Omega*Rp)^3)  # note that CQ = CP
        FM = CT<0. ? NaN : CT^(3.0/2)/(sqrt(2)*CP)

        return FM, CT, CP
    end

end


end  # module
