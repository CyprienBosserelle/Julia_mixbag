module CoastWaves


    export hs2hrms hrms2hs qkhf OffshoreWaveheight

    g=9.81;

    """
        Hs2Hrms(Hs)

    Convert from Significant wave height (Hs) to RMS wave height (Hrms).
    See also Hrms2Hs
    CB 30/03/2020

    # Examples
    ```julia-repl
    julia> Hrms=Hs2Hrms(1.0)
    1
    ```
    """
    hs2hrms(Hs)=Hs./sqrt(2);

    hrms2hs(Hrms)=Hrms.*sqrt(2);

    """
        kh = qkhf( w, h )
    QKHFS - Quick iterative calculation of kh in dispersion relationship
    Input:
    w Angular wave frequency = 2*pi/T where T = wave period [1/s]
    h Water depth [m]
    Returns:
    kh = wavenumber * depth [ ]
    Hard-wired for MKS units.
    Orbital velocities from kh are accurate to 3e-12 !
    RL Soulsby (2006) "Simplified calculation of wave orbital velocities"
    HR Wallingford Report TR 155, February 2006
    Eqns. 12a - 14
    Originally coded in Matlab by:
    csherwood@usgs.gov
    Sept 10, 2006

    Julia version from CB: 30/03/2020
    """
    function qkhfs( w, h )
        x = w.^2*h./g;
        y = sqrt.(x) .* (x.<1) + x.* (x.>=1);
        #this appalling bit of code is about 25% faster than a for loop in Matlab
        t = tanh.( y );
        y = y-( (y.*t -x)./(t+y.*(1 .- t.^2)));
        t = tanh.( y );
        y = y-( (y.*t -x)./(t+y.*(1 .- t.^2)));
        t = tanh.( y );
        y = y-( (y.*t -x)./(t+y.*(1 .- t.^2)));
        kh=y;
        return kh;
    end

    """
        Hoff=OffshoreWaveheight(Hin,Tp,h)
    Calculate offshore wave height [Hoff] based on wave height [Hin] at depth [h]
    for waves of period [Tp]
    Hin, Hoff, h are in [m] and Tp in [s].
    This uses the linear wave theory.
    """
    function OffshoreWaveheight(Hin,Tp,h)
        #Calculate wave length and wave celerity
        sig=2.0*pi./Tp
        kh = qkhfs( sig, h );
        L=2*pi./(kh./h);
        C=L./Tp;
        # offshore group velocity
        Cgo=g .* Tp ./ (4 * pi);
        # Group velocity at depth h
        Cg=0.5 .* (1 + (2 .* kh)./(sinh.(2 .* kh))) .* C;
        # shoaling coeff
        Ks=sqrt.(Cgo ./ Cg)

        return Hin./Ks
    end
