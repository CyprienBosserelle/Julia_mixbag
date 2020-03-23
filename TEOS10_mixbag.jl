module TEOS10_mixbag

    export o2sat,gsw_specvol,gsw_rho

    function gsw_specvol(sa,ct,p)
        #  Calculates specific volume from Absolute Salinity, Conservative
		#  Temperature and pressure, using the computationally-efficient
		#  polynomial expression for specific volume (Roquet et al., 2014).
		#
		# sa     : Absolute Salinity                               [g/kg]
		# ct     : Conservative Temperature (ITS-90)               [deg C]
		# p      : sea pressure                                    [dbar]
		#          ( i.e. absolute pressure - 10.1325 dbar )
		#
		# specvol: specific volume                                 [m^3/kg]

		#Define "some" constant
		gsw_sfac = 0.0248826675584615;
		offset = 5.971840214030754e-1;
		v000 =  1.0769995862e-3;
		v001 = -6.0799143809e-5;	v002 =  9.9856169219e-6;	v003 = -1.1309361437e-6;
		v004 =  1.0531153080e-7;	v005 = -1.2647261286e-8;	v006 =  1.9613503930e-9;
		v010 = -3.1038981976e-4;	v011 =  2.4262468747e-5;	v012 = -5.8484432984e-7;
		v013 =  3.6310188515e-7;	v014 = -1.1147125423e-7;	v020 =  6.6928067038e-4;
		v021 = -3.4792460974e-5;	v022 = -4.8122251597e-6;	v023 =  1.6746303780e-8;
		v030 = -8.5047933937e-4;	v031 =  3.7470777305e-5;	v032 =  4.9263106998e-6;
		v040 =  5.8086069943e-4;	v041 = -1.7322218612e-5;	v042 = -1.7811974727e-6;
		v050 = -2.1092370507e-4;	v051 =  3.0927427253e-6;	v060 =  3.1932457305e-5;
		v100 = -1.5649734675e-5;	v101 =  1.8505765429e-5;	v102 = -1.1736386731e-6;
		v103 = -3.6527006553e-7;	v104 =  3.1454099902e-7;	v110 =  3.5009599764e-5;
		v111 = -9.5677088156e-6;	v112 = -5.5699154557e-6;	v113 = -2.7295696237e-7;
		v120 = -4.3592678561e-5;	v121 =  1.1100834765e-5;	v122 =  5.4620748834e-6;
		v130 =  3.4532461828e-5;	v131 = -9.8447117844e-6;	v132 = -1.3544185627e-6;
		v140 = -1.1959409788e-5;	v141 =  2.5909225260e-6;	v150 =  1.3864594581e-6;
		v200 =  2.7762106484e-5;	v201 = -1.1716606853e-5;	v202 =  2.1305028740e-6;
		v203 =  2.8695905159e-7;	v210 = -3.7435842344e-5;	v211 = -2.3678308361e-7;
		v212 =  3.9137387080e-7;	v220 =  3.5907822760e-5;	v221 =  2.9283346295e-6;
		v222 = -6.5731104067e-7;	v230 = -1.8698584187e-5;	v231 = -4.8826139200e-7;
		v240 =  3.8595339244e-6;	v300 = -1.6521159259e-5;	v301 =  7.9279656173e-6;
		v302 = -4.6132540037e-7;	v310 =  2.4141479483e-5;	v311 = -3.4558773655e-6;
		v312 =  7.7618888092e-9;	v320 = -1.4353633048e-5;	v321 =  3.1655306078e-7;
		v330 =  2.2863324556e-6;	v400 =  6.9111322702e-6;	v401 = -3.4102187482e-6;
		v402 = -6.3352916514e-8;	v410 = -8.7595873154e-6;	v411 =  1.2956717783e-6;
		v420 =  4.3703680598e-6;	v500 = -8.0539615540e-7;	v501 =  5.0736766814e-7;
		v510 = -3.3052758900e-7;	v600 =  2.0543094268e-7;


        xs	= sqrt(gsw_sfac*sa + offset);
		ys	= ct*0.025;
		z	= p*1e-4;



		value = v000
	    + xs*(v010 + xs*(v020 + xs*(v030 + xs*(v040 + xs*(v050
	    + v060*xs))))) + ys*(v100 + xs*(v110 + xs*(v120 + xs*(v130 + xs*(v140
	    + v150*xs)))) + ys*(v200 + xs*(v210 + xs*(v220 + xs*(v230 + v240*xs)))
	    + ys*(v300 + xs*(v310 + xs*(v320 + v330*xs)) + ys*(v400 + xs*(v410
	    + v420*xs) + ys*(v500 + v510*xs + v600*ys))))) + z*(v001 + xs*(v011
	    + xs*(v021 + xs*(v031 + xs*(v041 + v051*xs)))) + ys*(v101 + xs*(v111
	    + xs*(v121 + xs*(v131 + v141*xs))) + ys*(v201 + xs*(v211 + xs*(v221
	    + v231*xs)) + ys*(v301 + xs*(v311 + v321*xs) + ys*(v401 + v411*xs
	    + v501*ys)))) + z*(v002 + xs*(v012 + xs*(v022 + xs*(v032 + v042*xs)))
	    + ys*(v102 + xs*(v112 + xs*(v122 + v132*xs)) + ys*(v202 + xs*(v212
	    + v222*xs) + ys*(v302 + v312*xs + v402*ys))) + z*(v003 + xs*(v013
	    + v023*xs) + ys*(v103 + v113*xs + v203*ys) + z*(v004 + v014*xs + v104*ys
	    + z*(v005 + v006*z)))));

		return (value);
    end
    function gsw_rho(sa, ct, p)

    	return (1.0/gsw_specvol(sa,ct,p));
    end
    function o2sat(SP, CT)

        #x = gsw_SP_from_SA(SA,p);
        x = SP;#practical salinity
        pt = CT;#gsw_pt_from_CT(SA,CT); # pt is potential temperature referenced to
                                    # the sea surface.
        pt68 = pt .* 1.00024; # pt68 is the potential temperature in degress C on
                      # the 1968 International Practical Temperature Scale IPTS-68.
        gsw_T0 = 273.15;
        #   Calculate Ts from T (deg C)
        y = log((298.15 - pt68) ./ (gsw_T0 + pt68));

        # The coefficents below are from the second column of Table 1 of Garcia and
        # Gordon (1992) in umol/kg
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
        # in umol/kg
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
