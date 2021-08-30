# Examples using the cylinder wake model

We provide examples using the wake model, in addition with supplemental material to understand the behavior of CCBlade when the wake model is on.

Useful link: [theory on the cylinder wake model](wake_model.md).



## Example 1: coned rotor

We first consider a simple rotor with straight blades under various coning angles. As the coning increases, we expect that the radial flow component (which is neglected in the standard BEM) plays an increasingly important role.

```@setup coned-rotor
using PyPlot
using CCBlade
nothing #hide
```


Let us use the same turbine as in the other [examples](howto.md).

```@example coned-rotor
    Rhub = 1.5
    Rtip = 63.0
    B = 3
    hubHt = 90.0

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
    aftypes[1] = AlphaAF("data/Cylinder1.dat", radians=false) 
    aftypes[2] = AlphaAF("data/Cylinder2.dat", radians=false)
    aftypes[3] = AlphaAF("data/DU40_A17.dat", radians=false)
    aftypes[4] = AlphaAF("data/DU35_A17.dat", radians=false)
    aftypes[5] = AlphaAF("data/DU30_A17.dat", radians=false)
    aftypes[6] = AlphaAF("data/DU25_A17.dat", radians=false)
    aftypes[7] = AlphaAF("data/DU21_A17.dat", radians=false)
    aftypes[8] = AlphaAF("data/NACA64_A17.dat", radians=false)

    # indices correspond to which airfoil is used at which station
    af_idx = [1, 1, 2, 3, 4, 4, 5, 6, 6, 7, 7, 8, 8, 8, 8, 8, 8]
    
    # create airfoil array 
    airfoils = aftypes[af_idx]
    
    sections = Section.(r, chord, theta, airfoils)
    nothing # hide
```


The main parameters that we want to explore the influence of are:

```@example coned-rotor
    precone = 10.0*pi/180  #the original value is -2.5 (negtive means forward)
    tilt = 0.0*pi/180     #the original value is 5.0 (positive means upwards)
    yaw = 0.0*pi/180
    shearExp = 0.0

    nothing # hide
```

The response to shear will also differ between the standard BEM and the one working with the wake model.


We also define common operating conditions.

```@example coned-rotor
    Vinf = 10.0
    tsr = 7.55
    rotorR = Rtip*cos(precone)
    Omega = Vinf*tsr/rotorR
    azimuth = 0.0*pi/180
    rho = 1.225
    pitch = 0.0
    nothing # hide
```




**Rotor 1**

As a reference, we will use the unconed rotor. Note the keyword argument `wakeCyl` to specify wheher or not you want to activate the wake model. If `false`, the standard BEM routines are used.

```@example coned-rotor

    rotor_nocone_nowake = Rotor(Rhub, Rtip, B; precone=0.0, turbine=true, wakeCyl=false)

    op = windturbine_op.(Vinf, Omega, pitch, r, 0.0, yaw, tilt, azimuth, hubHt, shearExp, rho)
    
    out_nocone_nowake = solve.(Ref(rotor_nocone_nowake), sections, op)

    nothing # hide
```


**Rotor 2**

The second rotor to compare with has some precone but does not use the wake model.

```@example coned-rotor

    rotor_cone_nowake = Rotor(Rhub, Rtip, B; precone=precone, turbine=true, wakeCyl=false)

    op = windturbine_op.(Vinf, Omega, pitch, r, precone, yaw, tilt, azimuth, hubHt, shearExp, rho)
    
    out_cone_nowake = solve.(Ref(rotor_cone_nowake), sections, op)
    
    nothing # hide
```

**Rotor 3**

The third rotor has precone and uses the wake model.

```@example coned-rotor
    
    rotor_cone_wake = Rotor(Rhub, Rtip, B; precone=precone, turbine=true, wakeCyl=true)

    op = windturbine_op.(Vinf, Omega, pitch, r, precone, yaw, tilt, azimuth, hubHt, shearExp, rho)
    
    out_cone_wake = solve.(Ref(rotor_cone_wake), sections, op)

    nothing # hide
```


**Comparison**

``` @example coned-rotor
    figure()
    plot(r/Rtip, out_nocone_nowake.Np/1e3)
    plot(r/Rtip, out_cone_nowake.Np/1e3,"--")
    plot(r/Rtip, out_cone_wake.Np/1e3)
    xlabel("r/Rtip")
    ylabel("normal loads (kN/m)")

    savefig("wakeEx1_Np.svg") # hide    

    figure()
    plot(r/Rtip, out_nocone_nowake.Tp/1e3)
    plot(r/Rtip, out_cone_nowake.Tp/1e3,"--")
    plot(r/Rtip, out_cone_wake.Tp/1e3)
    xlabel("r/Rtip")
    ylabel("tangential loads (kN/m)")

    legend(["out_nocone_nowake", "out_cone_nowake", "out_cone_wake"])

    savefig("wakeEx1_Tp.svg") # hide

    nothing # hide
```

![](wakeEx1_Np.svg)
![](wakeEx1_Tp.svg)





## Example 2: coned rotor with shear and tilt

TODO

## Example 3: influence of yaw

TODO

## Example 4: rotor with curved blades

TODO