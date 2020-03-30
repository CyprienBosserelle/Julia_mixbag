module ReefRunup

    import CoastWaves
    export Becker2014 simplifiedGourlay

    # g=9.81 #already imported from CoastWaves
    """
        n=simplifiedGourlay(Hs, Tp, h, h3)
    Calculate wave setup on a reef using a simplified version of Gourlay eq..
    This uses a simplified version of Eq 11 asssuming Kp=0.44 KR=0 and γr=0.4
    T√[(η+h)/g]=15

    See also:
    Gourlay()

    CB converts a matlab script but I can't remember who wrote it in the first place
    """
    function simplifiedGourlay(Hs, Tp, h, hr)
        # Hs is teh significant wave height at a depth h
        # Tp is the peak period
        # h3 is the depth of the reef flat






        Hrms=hs2hrms.(Hs)

        #Offshore significant wave height
        Ho=OffshoreWaveheight(Hrms,Tp,h);

        err=Ho;
        nmax=0.0;
        while abs(err)>0.00001
            global nmax,err
            nmax=nmax+err;

            y=0.006565*(1-0.134*((nmax+hr)/Ho)^2)*(Ho/(nmax+hr))^(3/2);

            x=y*Tp*sqrt(g*Ho);
            if (x<=0.0) #if x is negative or null the solution cannot converge
                nmax=0.0; #It will happen when the wave will not break so no setup occurs
                break;
            end


            err=x-nmax;

        end

        return nmax;
    end





    function Gourlay()

        hr=0.5;
        Ho=5.1
        Tp=12
        Kp=0.44
        γ=0.4
        KR=0



        err=Ho;
        nmax=0.0;
        while abs(err)>0.00001
            global nmax,err
            nmax=nmax+err;



            y=3.0/(64.0*pi)*Kp*(1-KR*KR-4.0*pi*γ*γ*(((nmax+hr)/Ho)^2)/Tp*sqrt(abs(nmax+hr)/g))*(Ho/(nmax+hr))^(3/2);

            x=y*Tp*sqrt(g*Ho);
            if (x<=0.0) #if x is negative or null the solution cannot converge
                nmax=0.0; #It will happen when the wave will not break so no setup occurs
                break;
            end


            err=x-nmax;

            println(nmax,err)

        end





    end
    # Function to write md files (next gen?)
    function Becker2014setup(γb,Hb,Hi)
        # γb is breaking gamma (0.4 to 1.3)
        # Hb is breaker height as RMS wave height
        # Hi is residual rms wave height
        # n=5/16*γb(Hb-1.2Hi)
        return 5.0/16*γb*(Hb-1.2*Hi)
    end
end
