module TEOS10_mixbag

    export o2sat

    function o2sat(SP, CT)

        #x = gsw_SP_from_SA(SA,p);
        x = SP;
        pt = CT;#gsw_pt_from_CT(SA,CT); # pt is potential temperature referenced to
                                    # the sea surface.
        pt68 = pt .* 1.00024; # pt68 is the potential temperature in degress C on
                      # the 1968 International Practical Temperature Scale IPTS-68.
        gsw_T0 = 273.15;
        y = log((298.15 - pt68) ./ (gsw_T0 + pt68));

        # The coefficents below are from the second column of Table 1 of Garcia and
        # Gordon (1992)
        a0 =  5.80871;
        a1 =  3.20291;
        a2 =  4.17887;
        a3 =  5.10006;
        a4 = -9.86643e-2;
        a5 =  3.80369;
        b0 = -7.01577e-3;
        b1 = -7.70028e-3;
        b2 = -1.13864e-2;
        b3 = -9.51519e-3;
        c0 = -2.75915e-7;

        O2sol = exp(a0 + y .* (a1 + y .* (a2 + y .* (a3 + y .* (a4 + a5 * y)))) + x .* (b0 + y .* (b1 + y .* (b2 + b3 * y)) + c0 * x));

        return O2sol;
    end
end
# function gsw_SP_from_SA(SA,p,long,lat)
    #gsw_SP_from_SA                  Practical Salinity from Absolute Salinity
    #Converted from the Matlab version of TEOS 10
    # Original AUTHOR:
    #  David Jackett, Trevor McDougall and Paul Barker    [ help_gsw@csiro.au ]
    # Converted to Julia by:
    # Cyprien Bosserelle
    #
    #DESCRIPTION:
    #  Calculates Practical Salinity from Absolute Salinity.
    #
    # INPUT:
    #  SA    =  Absolute Salinity                                      [ g/kg ]
    #  p     =  sea pressure                                           [ dbar ]
    #           ( i.e. absolute pressure - 10.1325 dbar )
    #  long  =  longitude in decimal degrees                     [ 0 ... +360 ]
    #                                                     or  [ -180 ... +180 ]
    #  lat   =  latitude in decimal degrees north               [ -90 ... +90 ]
    #
    #  p, lat and long may have dimensions 1x1 or Mx1 or 1xN or MxN,
    #  where SA is MxN.
    #
    # OUTPUT:
    #  SP        =  Practical Salinity  (PSS-78)
    #  Original software is available from http://www.TEOS-10.org
