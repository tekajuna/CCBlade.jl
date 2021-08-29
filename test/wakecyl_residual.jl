# --------- qualitative --------------
# The following tests are only qualitative.  They are not comparisons to known
# empirical data.  Rather they are from the figures used in the documentaiton.  Qualitatively
# the output is about right.  Minor changes in, for example, the airfoil interpolation method
# coudl slightly change the outputs.  The main purpose of these tests is to alert us if something
# significant changes.  

using CCBlade

# @testset "residual_behavior" begin


    Rhub = 1.5
    Rtip = 63.0
    B = 3
    precone = 2.5*pi/180
    tilt = 5.0*pi/180
    yaw = 0.0*pi/180
    
    # ============ PLAY WITH THIS : ==========================
    #THESE 3 SHOULD BE THE SAME AND THEY ARE NOT -> bug in precone???
    precone = 0.
    tilt = 2.5*pi/180 # -> works

    # precone = 2.5*pi/180
    # tilt = 5.0*pi/180

    # precone = -2.5*pi/180
    # tilt = 0.

    #wakecyl=true
    # tilt 0, precone 0 : some difference of course
    # tilt 5, precone 0 : very close
    # tilt 0, precone 2.5 : deviates at high tsr. CAUTION: must average over azm (because of shear)
    # tilt 5, precone 2.5 : deviates at high tsr. CAUTION: must average over azm (because of shear+tilt)

    #wakecyl=false
    # tilt 0, precone 0 : some difference of course, same as wakecyl
    # tilt 5, precone 0 : very close, same
    # tilt 0, precone 2.5 : very close, CT/CP slightly above
    # tilt 5, precone 2.5 : wrong
    
    Rtip_def = Rtip*cos(precone)
    Rhub_def = Rhub*cos(precone)

    rotor = Rotor(Rhub, Rtip, B, precone=precone, turbine=true)
    # rotor = Rotor(Rhub_def, Rtip_def, B; precone=precone, turbine=true, wakeCyl=true)
    
    r = [2.8667, 5.6000, 8.3333, 11.7500, 15.8500, 19.9500, 24.0500,
        28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500,
        56.1667, 58.9000, 61.6333]
    chord = [3.542, 3.854, 4.167, 4.557, 4.652, 4.458, 4.249, 4.007, 3.748,
        3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419]
    theta = pi/180*[13.308, 13.308, 13.308, 13.308, 11.480, 10.162, 9.011, 7.795,
        6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106]
    
    # Define airfoils.  In this case we have 8 different airfoils that we load into an array.
    # These airfoils are defined in files.
    aftypes = Array{AlphaAF}(undef, 8)
    aftypes[1] = AlphaAF("airfoils/Cylinder1.dat", radians=false) 
    aftypes[2] = AlphaAF("airfoils/Cylinder2.dat", radians=false)
    aftypes[3] = AlphaAF("airfoils/DU40_A17.dat", radians=false)
    aftypes[4] = AlphaAF("airfoils/DU35_A17.dat", radians=false)
    aftypes[5] = AlphaAF("airfoils/DU30_A17.dat", radians=false)
    aftypes[6] = AlphaAF("airfoils/DU25_A17.dat", radians=false)
    aftypes[7] = AlphaAF("airfoils/DU21_A17.dat", radians=false)
    aftypes[8] = AlphaAF("airfoils/NACA64_A17.dat", radians=false)
    
    # indices correspond to which airfoil is used at which station
    af_idx = [1, 1, 2, 3, 4, 4, 5, 6, 6, 7, 7, 8, 8, 8, 8, 8, 8]
    
    # create airfoil array 
    airfoils = aftypes[af_idx]
    
    sections = Section.(r, chord, theta, airfoils, precone)
    
    # operating point for the turbine
    

    hubHt = 90.0
    shearExp = 0.2
    
    Vinf = 10.0
    tsr = 7.55
    # tsr = 12. #...........................
    Omega = Vinf*tsr/Rtip_def
    azimuth = 0.0*pi/180
    rho = 1.225
    pitch = 0.0
    

# #### CHECK residual

# # TODO: COMPARE precone,tilt combinations
# precone = 0.
# tilt = 2.5*pi/180 # -> works
# Rtip_def = Rtip*cos(precone)
# # Rhub_def = Rhub*cos(precone)
# # sections = Section.(r, chord, theta, airfoils, precone)
# sections = Section.(r, chord, theta, airfoils, zero(r), zero(r), cos(precone)*r, zero(r), zero(r) )
# # sections = Section.(r, chord, theta, airfoils, zero(r), zero(r), r, zero(r), zero(r) )
# Omega = Vinf*tsr/Rtip_def
# op = windturbine_op.(Vinf, Omega, pitch, r, precone, yaw, tilt, azimuth, hubHt, shearExp, rho)
#     # rotor1 = Rotor(Rhub_def, Rtip_def, B; precone=precone, turbine=true, wakeCyl=false)
#     rotor1 = Rotor(Rhub, Rtip, B; precone=precone, turbine=true, wakeCyl=false)
#     resi1,_ = CCBlade.residual(10*pi/180, rotor1, sections[1], op[1]);

# println(op[1])

# precone = 2.5*pi/180
# tilt = 5.0*pi/180
# Rtip_def = Rtip*cos(precone)
# # Rhub_def = Rhub*cos(precone)
# # sections = Section.(r, chord, theta, airfoils, precone)
# sections = Section.(r, chord, theta, airfoils, zero(r), zero(r), cos(precone)*r, zero(r), zero(r) )
# # sections = Section.(r, chord, theta, airfoils, zero(r), zero(r), r, zero(r), zero(r) )
# Omega = Vinf*tsr/Rtip_def
# op = windturbine_op.(Vinf, Omega, pitch, r, precone, yaw, tilt, azimuth, hubHt, shearExp, rho)
#     # rotor2 = Rotor(Rhub_def, Rtip_def, B; precone=precone, turbine=true, wakeCyl=false)
#     rotor2 = Rotor(Rhub, Rtip, B; precone=precone, turbine=true, wakeCyl=false)
#     resi2,_ = CCBlade.residual(10*pi/180, rotor2, sections[1], op[1]);

# println(op[1])
# # OperatingPoint{Float64, Float64, Float64, Float64, Float64, Float64, Float64}(10.053271233693533, 3.438762622588977, 0.43893530137807996, 10.0, 1.1995544084100105, 1.225, 0.0, 0.0, 0.0, 1.0, 1.0)
# # OperatingPoint{Float64, Float64, Float64, Float64, Float64, Float64, Float64}(10.024556658162442, 3.4354896825396826, 0.877035064462535, 10.0, 1.1995544084100105, 1.225, 0.0, 0.0, 0.0, 1.0, 1.0)

#     println(resi1)
#     println(resi2)

#     error("stop")

    # ====================== SINGLE OP ========================9==================================

    
    op = windturbine_op.(Vinf, Omega, pitch, r, precone, yaw, tilt, azimuth, hubHt, shearExp, rho)
    
    out = solve.(Ref(rotor), sections, op)
    
    
    nr = length(r)
    Npnorm = [0.09718339956327719, 0.13149361678303992, 0.12220741751313423, 1.1634860128761517, 1.7001259694801125, 2.0716257635881257, 2.5120015027019678, 3.1336171133495045, 3.6916824972696465, 4.388661772599469, 5.068896486058948, 5.465165634664408, 6.035059239683594, 6.539134070994739, 6.831387531628286, 6.692665814418597, 4.851568452578296]
    Tpnorm = [-0.03321034106737285, -0.08727189682145081, -0.12001897678344217, 0.4696423333976085 , 0.6226256283641799 , 0.6322961942049257 , 0.6474145670774534 , 0.6825021056687035 , 0.6999861694557595 , 0.7218774840801262 , 0.7365515905555542 , 0.7493905698765747 , 0.7529143446199785 , 0.7392483947274653, 0.6981206044592225, 0.614524256128813, 0.40353047553570615]
    for i = 1:nr
        # @test isapprox.(out.Np[i]/1e3, Npnorm[i], atol=1e-3)
        # @test isapprox.(out.Tp[i]/1e3, Tpnorm[i], atol=1e-3)
    end
    

    # ---------------- plots loads -----------------------
    using PyPlot

    # plot distributed loads
    figure()
    plot(r/Rtip, out.Np/1e3)
    plot(r/Rtip, out.Tp/1e3)
    plot(r/Rtip, Npnorm, "--")
    plot(r/Rtip, Tpnorm, "--")
    xlabel("r/Rtip")
    ylabel("distributed loads (kN/m)")
    legend(["flapwise", "lead-lag"])
    # savefig("loads-turbine.svg"); nothing # hide

    # figure()
    # plot(r/Rtip, op.Vx)
    # plot(r/Rtip, op.Vy)


    # ---------------- plot residuals -----------------------
    figure()
    # plot(r/Rtip, out.a)
    # plot(r/Rtip, out.ap)
    plot(r/Rtip, out.phi.*180. /pi)

    #CHOOSE A SPANWISE LOCATION index
    isp = length(r)-1

    #creating phi vector refined near 0
    nn = 1000
    p1 = range(-180.,-1.,length=100)
    p2 = range(-1,1.,length=nn)
    phi = vcat(p1[1:end-1],p2,-p1[end-1:-1:1])
    
    Res1 = zero(phi)

    plot(r[isp]/Rtip, 0.,"x")

    # compute and plot the residual 
    for i =1:length(phi)
        Res1[i] = CCBlade.residual(phi[i].* pi/180, rotor, sections[isp], op[isp])[1]
    end

    figure()
    plot(out[isp].phi.* 180/pi, 0.,"x")
    plot(phi,Res1)
    ylim([-1,1])


    

    # rr = range(0,R-eps,length=17)
    # ψψ = range(0,2*pi,length=13)'
    
    # ut = CCBlade.eval_ut(x, ψψ, rr, R, χ) .* γt
    # ur = CCBlade.eval_ur_0(ψψ, rr, R, χ) .* γt
    # ul = CCBlade.eval_ul(x, ψψ, rr, R, χ) .* 2. ./tan(χ/2) # =ul/uy0
    
    # ## ut
    
    # f = plt.figure(4)
    # ax = f.add_subplot(111, polar=true)
    
    # u1 = ut[3,1,:,:]' .* cos.(ψψ) - ut[2,1,:,:]' .* sin.(ψψ)
    # u2 = ut[3,1,:,:]' .* sin.(ψψ) + ut[2,1,:,:]' .* cos.(ψψ)
    # # ax.quiver(ψψ, rr, u1, u2)
    # ax.quiver(ψψ, rr, -u1, -u2)
    # ax.invert_xaxis()
    # ax.set_theta_zero_location("W")  # theta=0 at the left
    # # ax.set_theta_direction(-1)  # theta increasing clockwise
    
    # f = plt.figure(5)
    # ax = f.add_subplot(111, polar=true)
    # ax.contour(ψψ, rr, ut[1,1,:,:]',[-.7,-.6,-.5,-.4,-.3])
    # ax.invert_xaxis()
    # ax.set_theta_zero_location("W")  # theta=0 at the left
    # ax.set_theta_direction(-1)  # theta increasing clockwise
    


    # error()




    
    # T, Q = thrusttorque(rotor, sections, out)
    
    # azangles = pi/180*[0.0, 90.0, 180.0, 270.0]
    # ops = windturbine_op.(Vinf, Omega, r, precone, yaw, tilt, azangles', hubHt, shearExp, rho)
    # outs = solve.(Ref(rotor), sections, ops)
    
    # T, Q = thrusttorque(rotor, sections, outs)
    
    # CP, CT, CQ = nondim(T, Q, Vinf, Omega, rho, rotor)
    

    # ====================== ACROSS MULTIPLE OP ==========================================================


    ntsr = 20  # number of tip-speed ratios
    tsrvec = range(2, stop=15, length=ntsr)
    cpvec = zeros(ntsr)  # initialize arrays
    ctvec = zeros(ntsr)
    
    azangles = pi/180*[0.0, 90.0, 180.0, 270.0]
    # azangles = [0.0,]

    # figure()
    for i = 1:ntsr
        Omega = Vinf*tsrvec[i]/rotorR
    
        ops = windturbine_op.(Vinf, Omega, pitch, r, precone, yaw, tilt, azangles', hubHt, shearExp, rho)     #tilt!!!
        outs = solve.(Ref(rotor), sections, ops)
    
        # for i= 1:4
        #     plot(r/Rtip, outs[:,i].Np/1e3,"-o")
        # end
    

        T, Q = thrusttorque(rotor, sections, outs)

    
        cpvec[i], ctvec[i], _ = nondim(T, Q, Vinf, Omega, rho, rotor, "windturbine")
    end
    
    cpvec_test = [0.02350364982213745, 0.07009444382725231, 0.13891408294580307, 0.21795999362154536, 0.30793657066245234, 0.39220453169574504, 0.43584242088313013, 0.4590001212905916, 0.4695346630018769, 0.4680559467091491, 0.45853627726342805, 0.44346560401350865, 0.4247660950655427, 0.4031517121734798, 0.37850875685858926, 0.3506421939620393, 0.31947670368255254, 0.28492212475211376, 0.24684802187254656, 0.2053516821754716]
    ctvec_test = [0.1234606206655405, 0.19135957414241966, 0.2753772801035205, 0.3637320188268945, 0.46087282130376256, 0.5702751788685602, 0.6510892321304749, 0.7117745993566226, 0.7600809643657047, 0.7994008741637466, 0.8332719078069747, 0.8638031244208144, 0.8924359419919573, 0.9199569340127343, 0.9465496984524919, 0.9723679911393518, 0.9975143642860824, 1.0220510627832238, 1.0458478359218666, 1.0685364647098319]
    
    for i = 1:ntsr
        # @test isapprox(cpvec[i], cpvec_test[i], atol=1e-3)
        # @test isapprox(ctvec[i], ctvec_test[i], atol=1e-3)
    end
    

    
    figure()
    plot(tsrvec,cpvec)
    plot(tsrvec,cpvec_test)
    plot(tsr,0.0,"x")
    xlabel(L"\lambda")
    ylabel(L"C_P")
    
    figure()
    plot(tsrvec,ctvec)
    plot(tsrvec,ctvec_test)
    xlabel(L"\lambda")
    ylabel(L"C_T")

    
    
# end