module ReefRunup

    import Printf
    export Becker2014

    # Function to write md files (next gen?)
    function Becker2014setup(γb,Hb,Hi)
        # γb is breaking gamma (0.4 to 1.3)
        # Hb is breaker height as RMS wave height
        # Hi is residual rms wave height
        # n=5/16*γb(Hb-1.2Hi)
        return 5.0/16*γb*(Hb-1.2*Hi)
    end
end
