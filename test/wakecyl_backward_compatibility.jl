using CCBlade
using Test


@testset "pitch_and_twist" begin
#no wake cylinder at this point, just checking that the new formulation  still permits to exchange pitch and twist.


# -------- verification: turbine ------
# same turbine as in runtest

# --- rotor definition ---
Rhub = 0.01
Rtip = 5.0
Rtip_eff = Rtip*100  # to eliminate tip effects as consistent with their study.
B = 3  # number of blades

#let's start with no yaw, tilt... nothing
rotor = Rotor(Rhub, Rtip_eff, B, turbine=true)

# --- section definitions ---
#straight blades

r = [0.2, 1, 2, 3, 4, 5]
gamma = [61.0, 74.31002131, 84.89805553, 89.07195504, 91.25038415, 92.58003871]
theta = (90.0 .- gamma)*pi/180
chord = [0.7, 0.706025153, 0.436187551, 0.304517933, 0.232257636, 0.187279622]

function affunc(alpha, Re, M)

    cl = 0.084*alpha*180/pi

    return cl, 0.0
end 

sections = Section.(r, chord, theta, Ref(affunc))


# --- inflow definitions ---

Vinf = 7.0
tsr = 8
Omega = tsr*Vinf/Rtip
rho = 1.0

ops = simple_op.(Vinf, Omega, r, rho)

# --- preliminary test of the residuals ---

# phi = 30. *pi/180.
# CCBlade.residual(phi, rotor, sections[2], ops[2])



# -------- verification: turbine PITCH VS TWIST ------

# START WITH PITCH = 0, SHOULD GIVE EXACT SAME RES
# WITH PITCH =/= 0, SHOULD ALSO GIVE SAME RES even though some intermediate res may be different?

ptch = 3. *pi/180.

sections_pitch = Section.(r, chord, theta .+ ptch, Ref(affunc))
ops_nopitch = simple_op.(Vinf, Omega, r, rho)

sections_nopitch = Section.(r, chord, theta, Ref(affunc))
ops_pitch = simple_op.(Vinf, Omega, r, rho, pitch=ptch)


# - first simply on the residual -
phi = 30. *pi/180.

R1,_ = CCBlade.residual(phi, rotor, sections_pitch[2], ops_nopitch[2])
R2,_ = CCBlade.residual(phi, rotor, sections_nopitch[2], ops_pitch[2])

@test isapprox(R1, R2, atol=1e-8)

# - the the full solve -
#1
outputs_1 = solve.(Ref(rotor), sections_pitch, ops_nopitch)
T1, Q1 = thrusttorque(rotor, sections_pitch, outputs_1)

#2
outputs_2 = solve.(Ref(rotor), sections_nopitch, ops_pitch)
T2, Q2 = thrusttorque(rotor, sections_nopitch, outputs_2)

@test isapprox(T1, T2, atol=1e-8)

end


# #################################################################################################
# #################################################################################################
# #################################################################################################


@testset "other_tmp_test" begin



# -------- verification: propellers.  using script at http://www.aerodynamics4students.com/propulsion/blade-element-propeller-theory.php ------


# ************************************
# I increased their tolerance to 1e-6
# ************************************

# --- rotor definition ---
D = 1.6
Rhub = 0.0
Rtip = D/2
Rhub_eff = 1e-6  # something small to eliminate hub effects
Rtip_eff = 100.0  # something large to eliminate tip effects
B = 2  # number of blades

rotor_no_F = Rotor(Rhub_eff, Rtip_eff, B)
rotor = Rotor(Rhub, Rtip, B)


# --- section definitions ---

R = D/2.0
r = range(R/10, stop=R, length=11)
pitch = 1.0  # pitch distance in meters.
theta = atan.(pitch./(2*pi*r))
chord = 0.10


function affunc(alpha, Re, M)

    cl = 6.2*alpha
    cd = 0.008 - 0.003*cl + 0.01*cl*cl

    return cl, cd
end 

sections = Section.(r, chord, theta, Ref(affunc))


# --- inflow definitions ---
RPM = 2100
rho = 1.225



# --- preliminary test of the residuals ---
Vinf = 5.0
Omega = RPM * pi/30 

ops = simple_op.(Vinf, Omega, r, rho)

phi = 30. *pi/180.
CCBlade.residual(phi, rotor, sections[6], ops[6])




# # --- evaluate ---
tsim = 1e3*[1.045361193032356, 1.025630300048415, 1.005234466788998, 0.984163367036026, 0.962407923825140, 0.939960208707079, 0.916813564966455, 0.892962691000145, 0.868403981825492, 0.843134981103815, 0.817154838249790, 0.790463442573673, 0.763063053839278, 0.734956576558370, 0.706148261507327, 0.676643975451150, 0.646450304160057, 0.615575090105131, 0.584027074365864, 0.551815917391907, 0.518952127358381, 0.485446691671386, 0.451311288662196, 0.416557935286392, 0.381199277009438, 0.345247916141561, 0.308716772800348, 0.271618894441869, 0.233967425339051, 0.195775319296371, 0.157055230270717, 0.117820154495231, 0.078082266879117, 0.037854059197644, -0.002852754149850, -0.044026182837742, -0.085655305814570, -0.127728999394140, -0.170237722799272, -0.213169213043848, -0.256515079286031, -0.300266519551194, -0.344414094748869, -0.388949215983616, -0.433863576642539, -0.479150401337354, -0.524801553114807, -0.570810405128802, -0.617169893200684, -0.663873474163182, -0.710915862524620, -0.758291877949762, -0.805995685105502, -0.854022273120508, -0.902366919041604, -0.951025170820984, -0.999992624287163, -1.049265666456123, -1.098840222937414, -1.148712509929845]
qsim = 1e2*[0.803638686218187, 0.806984572453978, 0.809709290183008, 0.811743686838315, 0.813015017103876, 0.813446921530685, 0.812959654049620, 0.811470393912576, 0.808893852696513, 0.805141916379142, 0.800124489784850, 0.793748780791057, 0.785921727832179, 0.776548246109426, 0.765532528164390, 0.752778882688809, 0.738190986274448, 0.721673076180745, 0.703129918771009, 0.682467282681955, 0.659592296506578, 0.634413303042323, 0.606840565246423, 0.576786093366321, 0.544164450503912, 0.508891967461804, 0.470887571011192, 0.430072787279711, 0.386371788290446, 0.339711042057184, 0.290019539402947, 0.237229503458026, 0.181274942660876, 0.122093307308376, 0.059623821454727, -0.006190834182631, -0.075406684829235, -0.148076528546541, -0.224253047813501, -0.303980950928302, -0.387309291734422, -0.474283793689904, -0.564946107631716, -0.659336973911858, -0.757495165410553, -0.859460291551374, -0.965266648683888, -1.074949504731187, -1.188540970723477, -1.306072104649531, -1.427575034895290, -1.553080300508925, -1.682614871422754, -1.816205997296014, -1.953879956474228, -2.095662107769925, -2.241576439746701, -2.391647474158875, -2.545897099743367, -2.704346566395035]

for i = 1:60

    Vinf = float(i)
    Omega = RPM * pi/30 

    ops = simple_op.(Vinf, Omega, r, rho)

    # --- evaluate ---

    out = solve.(Ref(rotor_no_F), sections, ops)

    # Np, Tp = loads(outputs)
    # T, Q = thrusttorque(r[1], r, r[end], Np, Tp, B)
    # their spreadsheet did not use trapzezoidal rule, so just do a rect sum.
    T = sum(out.Np*(r[2]-r[1]))*B
    Q = sum(r.*out.Tp*(r[2]-r[1]))*B

    @test isapprox(T, tsim[i], atol=1e-2)  # 2 decimal places
    @test isapprox(Q, qsim[i], atol=2e-3)

end


Vinf = 20.0
Omega = RPM * pi/30

op = simple_op.(Vinf, Omega, r, rho)

out = solve.(Ref(rotor_no_F), sections, op)

T = sum(out.Np*(r[2]-r[1]))*B
Q = sum(r.*out.Tp*(r[2]-r[1]))*B
eff, CT, CQ = nondim(T, Q, Vinf, Omega, rho, rotor, "propeller")


@test isapprox(CT, 0.056110238632657, atol=1e-7)
@test isapprox(CQ, 0.004337202960642, atol=1e-8)
@test isapprox(eff, 0.735350632777002, atol=1e-6)


end


# #################################################################################################
# #################################################################################################
# #################################################################################################



