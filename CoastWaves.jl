module CoastWaves


    export hs2hrms,hrms2hs,qkhfs,OffshoreWaveheight

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

    function Shoalcoeff(H,Tp,h)
        #Calculate wave length and wave celerity
        sig=2.0*pi./Tp
        kh = qkhfs(sig, h);
        L=2*pi./(kh./h);
        C=L./Tp;
        # offshore group velocity
        Cgo=g .* Tp ./ (4 * pi);
        # Group velocity at depth h
        Cg=0.5 .* (1 + (2 .* kh)./(sinh.(2 .* kh))) .* C;
        # shoaling coeff
        Ks=sqrt.(Cgo ./ Cg)

        return Ks
    end

    """
        Hoff=OffshoreWaveheight(Hin,Tp,h)
    Calculate offshore wave height [Hoff] based on wave height [Hin] at depth [h]
    for waves of period [Tp]
    Hin, Hoff, h are in [m] and Tp in [s].
    This uses the linear wave theory.
    """
    function OffshoreWaveheight(Hin,Tp,h)
       
        # shoaling coeff
        Ks=Shoalcoeff(Hin,Tp,h)

        return Hin./Ks
    end

    """
    Sd,Sp,Rp,Sm=Wsetup(Hs,T,D,slp)
    This function calculate the wave set down(Sd) wave set up (Sp) wave run up
    (Rp) and the distance of the inundation (Sm) for a given wave height(Hs),
    Period (T) at a given depth (D) and a given beach slope(slp(e.g.: slp=1/30)).
    Wave set-up and Run-up calculated using CEM (Weggel, longuet-Higgins)
    The program calculate the beaking height and depth for each wave using
    linear wave theory.
    Matlab transcripted by    Cyprien Bosserelle 22/05/2008
    Julia translation CB 24/03/2023
    """
    function Wsetup(Hs,Tp,h,slp)
        Hrms=hs2hrms(Hs)

        Ho=OffshoreWaveheight(Hrms,Tp,h)

        hh=20:-0.01:0.01;
        KSall=Shoalcoeff.(Ho,Tp,hh)

        Hin=Ho.*KSall;
        gamo=Hin./hh;

        a=43.8*(1-exp(-19*slp));
        b=1.56/(1+exp(-19.5*slp));

        gam1=b.-a.*Hin./(g*Tp^2);

        indmin=argmin(abs.(gamo.-gam1));

        Hb=Hin[indmin];
        hb=hh[indmin];

        kb=qkhfs( 2.0*pi./Tp, hb )

        gam=Hb/hb;

        Sd=(-1/8)*Hb^2*(kb)/(sinh(2*kb*hb));
        Sp=Sd+(1/(1+8/(3*gam^2)))*(hb-Sd);
        nn=1/(1+8/(3*gam^2))*slp;
        Rp=Sp/(slp-nn);
        Sm=Sp+nn*Rp;

        return Sd,Sp,Rp,Sm
    end


    function Stockdon2006(Hs,Tp,h,slp)
        sig=2.0*pi./Tp
        kh = qkhfs(sig, h);
        Lp=2*pi./(kh./h);
        Sp=0.35*slp*sqrt(Hs*Lp)
        return Sp
    end



end
