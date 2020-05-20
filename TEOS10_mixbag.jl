module TEOS10_mixbag

    export o2sat,gsw_specvol,gsw_rho,gsw_CT_from_t

	function gsw_entropy_part_zerop(SA,pt0)
		# gsw_entropy_part_zerop          entropy_part evaluated at the sea surface
		#
		# This function calculates entropy at a sea pressure of zero, except that
		# it does not evaluate any terms that are functions of Absolute Salinity
		# alone.  By not calculating these terms, which are a function only of
		# Absolute Salinity, several unnecessary computations are avoided
		# (including saving the computation of a natural logarithm). These terms
		# are a necessary part of entropy, but are not needed when calculating
		# potential temperature from in-situ temperature.
		# The inputs to "gsw_entropy_part_zerop(SA,pt0)" are Absolute Salinity
		# and potential temperature with reference sea pressure of zero dbar.
		#
		# VERSION NUMBER: 3.05 (27th January 2015)
		#
		#
		sfac = 0.0248826675584615;               # sfac = 1/(40*(35.16504/35));

		x2 = sfac.*SA;
		x = sqrt(x2);
		y = pt0.*0.025;

		g03 =  y.*(-24715.571866078 + y.*(2210.2236124548363 +
		    y.*(-592.743745734632 + y.*(290.12956292128547 +
		    y.*(-113.90630790850321 + y.*21.35571525415769)))));

		g08 = x2.*(x.*( x.*(y.*(-137.1145018408982 + y.*(148.10030845687618 +
		    y.*(-68.5590309679152 + 12.4848504784754.*y)))) +
		    y.*(-86.1329351956084 + y.*(-30.0682112585625 + y.*3.50240264723578))) +
		    y.*(1760.062705994408 + y.*(-675.802947790203 +
		    y.*(365.7041791005036 + y.*(-108.30162043765552 + 12.78101825083098.*y)))));

		entropy_part_zerop = -(g03 + g08).*0.025;

		return entropy_part_zerop

	end
	function gsw_entropy_part(SA,t,p)
		#gsw_entropy_part    entropy minus the terms that are a function of only SA
		#
		# This function calculates entropy, except that it does not evaluate any
		# terms that are functions of Absolute Salinity alone.  By not calculating
		# these terms, which are a function only of Absolute Salinity, several
		# unnecessary computations are avoided (including saving the computation
		# of a natural logarithm).  These terms are a necessary part of entropy,
		# but are not needed when calculating potential temperature from in-situ
		# temperature.
		#
		# VERSION NUMBER: 3.05 (27th January 2015)
		#
		sfac = 0.0248826675584615;                # sfac = 1/(40*(35.16504/35));

		x2 = sfac.*SA;
		x = sqrt(x2);
		y = t.*0.025;
		z = p.*1e-4;

		g03 = z.*(-270.983805184062 +
		    z.*(776.153611613101 + z.*(-196.51255088122 + (28.9796526294175 - 2.13290083518327.*z).*z))) +
		    y.*(-24715.571866078 + z.*(2910.0729080936 +
		    z.*(-1513.116771538718 + z.*(546.959324647056 + z.*(-111.1208127634436 + 8.68841343834394.*z)))) +
		    y.*(2210.2236124548363 + z.*(-2017.52334943521 +
		    z.*(1498.081172457456 + z.*(-718.6359919632359 + (146.4037555781616 - 4.9892131862671505.*z).*z))) +
		    y.*(-592.743745734632 + z.*(1591.873781627888 +
		    z.*(-1207.261522487504 + (608.785486935364 - 105.4993508931208.*z).*z)) +
		    y.*(290.12956292128547 + z.*(-973.091553087975 +
		    z.*(602.603274510125 + z.*(-276.361526170076 + 32.40953340386105.*z))) +
		    y.*(-113.90630790850321 + y.*(21.35571525415769 - 67.41756835751434.*z) +
		    z.*(381.06836198507096 + z.*(-133.7383902842754 + 49.023632509086724.*z)))))));

		g08 = x2.*(z.*(729.116529735046 +
		    z.*(-343.956902961561 + z.*(124.687671116248 + z.*(-31.656964386073 + 7.04658803315449.*z)))) +
		    x.*( x.*(y.*(-137.1145018408982 + y.*(148.10030845687618 + y.*(-68.5590309679152 + 12.4848504784754.*y))) -
		    22.6683558512829.*z) + z.*(-175.292041186547 + (83.1923927801819 - 29.483064349429.*z).*z) +
		    y.*(-86.1329351956084 + z.*(766.116132004952 + z.*(-108.3834525034224 + 51.2796974779828.*z)) +
		    y.*(-30.0682112585625 - 1380.9597954037708.*z + y.*(3.50240264723578 + 938.26075044542.*z)))) +
		    y.*(1760.062705994408 + y.*(-675.802947790203 +
		    y.*(365.7041791005036 + y.*(-108.30162043765552 + 12.78101825083098.*y) +
		    z.*(-1190.914967948748 + (298.904564555024 - 145.9491676006352.*z).*z)) +
		    z.*(2082.7344423998043 + z.*(-614.668925894709 + (340.685093521782 - 33.3848202979239.*z).*z))) +
		    z.*(-1721.528607567954 + z.*(674.819060538734 +
		    z.*(-356.629112415276 + (88.4080716616 - 15.84003094423364.*z).*z)))));

		return -(g03 + g08).*0.0250;

	end

	function gsw_gibbs_pt0_pt0(SA,pt0)
		# This function calculates the second derivative of the specific Gibbs
		# function with respect to temperature at zero sea pressure.  The inputs
		# are Absolute Salinity and potential temperature with reference sea
		# pressure of zero dbar.  This library function is called by both
		# "gsw_pt_from_CT(SA,CT)" ,"gsw_pt0_from_t(SA,t,p)" and
		# "gsw_pt_from_entropy(SA,entropy)".
		#
		# VERSION NUMBER: 3.05 (27th January 2015)
		#

		sfac = 0.0248826675584615;                # sfac = 1/(40*(35.16504/35));

		x2 = sfac.*SA;
		x = sqrt(x2);
		y = pt0.*0.025;

		g03 = -24715.571866078 + y.*(4420.4472249096725 +
		    y.*(-1778.231237203896 + y.*(1160.5182516851419 +
		    y.*(-569.531539542516 + y.*128.13429152494615))));

		g08 = x2.*(1760.062705994408 + x.*(-86.1329351956084 +
		    x.*(-137.1145018408982 + y.*(296.20061691375236 +
		    y.*(-205.67709290374563 + 49.9394019139016.*y))) +
		    y.*(-60.136422517125 + y.*10.50720794170734)) +
		    y.*(-1351.605895580406 + y.*(1097.1125373015109 +
		    y.*(-433.20648175062206 + 63.905091254154904.*y))));

		return (g03 + g08).*0.000625;
	end

	function gsw_pt0_from_t(SA,t,p)
		# gsw_pt0_from_t                               potential temperature with a
		#                                       reference sea pressure of zero dbar
		#
		# USAGE:
		#  pt0 = gsw_pt0_from_t(SA,t,p)
		#
		# DESCRIPTION:
		#  Calculates potential temperature with reference pressure, p_ref = 0 dbar.
		#  The present routine is computationally faster than the more general
		#  function "gsw_pt_from_t(SA,t,p,p_ref)" which can be used for any
		#  reference pressure value.
		#  This subroutine calls "gsw_entropy_part(SA,t,p)",
		#  "gsw_entropy_part_zerop(SA,pt0)" and "gsw_gibbs_pt0_pt0(SA,pt0)".
		#
		# INPUT:
		#  SA  =  Absolute Salinity                                        [ g/kg ]
		#  t   =  in-situ temperature (ITS-90)                            [ deg C ]
		#  p   =  sea pressure                                             [ dbar ]
		#         ( i.e. absolute pressure - 10.1325 dbar )
		#
		#  SA & t need to have the same dimensions.
		#  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & t are MxN.
		#
		# OUTPUT:
		#  pt0  =  potential temperature                                  [ deg C ]
		#          with reference sea pressure (p_ref) = 0 dbar.
		#  Note. The reference sea pressure of the output, pt0, is zero dbar.
		#
		# AUTHOR:
		#  Trevor McDougall, David Jackett, Claire Roberts-Thomson and Paul Barker.
		#                                                      [ help@teos-10.org ]
		#
		# VERSION NUMBER: 3.05 (27th January 2015)
		#
		# REFERENCES:
		#  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
		#   seawater - 2010: Calculation and use of thermodynamic properties.
		#   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
		#   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
		#    See section 3.1 of this TEOS-10 Manual.
		#
		#  McDougall T. J. and S. J. Wotherspoon, 2013: A simple modification of
		#   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied
		#   Mathematics Letters, 29, 20-25.
		#
		#  The software is available from http://www.TEOS-10.org
		#
		gsw_cp0=3991.86795711963;
		gsw_SSO=35.16504;
		gsw_T0=273.15;

		SSO = gsw_SSO;                    # from section 2.4 of IOC et al. (2010).

		s1 = SA*(35.0 ./ SSO);

		pt0 = t + p.*( 8.65483913395442e-6  -
		         s1.*  1.41636299744881e-6  -
		          p.*  7.38286467135737e-9  +
		          t.*(-8.38241357039698e-6  +
		         s1.*  2.83933368585534e-8  +
		          t.*  1.77803965218656e-8  +
		          p.*  1.71155619208233e-10));

		dentropy_dt = gsw_cp0./((gsw_T0 + pt0).*(1 - 0.05.*(1 - SA./SSO)));

		true_entropy_part = gsw_entropy_part(SA,t,p);

		for Number_of_iterations = 1:2
		    pt0_old = pt0;
		    dentropy = gsw_entropy_part_zerop(SA,pt0_old) - true_entropy_part;
		    pt0 = pt0_old - dentropy./dentropy_dt ; # this is half way through the modified method (McDougall and Wotherspoon, 2012)
		    pt0m = 0.5*(pt0 + pt0_old);
		    dentropy_dt = -gsw_gibbs_pt0_pt0(SA,pt0m);
		    pt0 = pt0_old - dentropy./dentropy_dt;
		end

		return pt0



	end


	function gsw_CT_from_pt(SA,pt)
		#gsw_CT_from_pt        Conservative Temperature from potential temperature
		# USAGE:
		#  CT = gsw_CT_from_pt(SA,pt)
		#
		# DESCRIPTION:
		#  Calculates Conservative Temperature of seawater from potential
		#  temperature (whose reference sea pressure is zero dbar).
		#
		# INPUT:
		#  SA  =  Absolute Salinity                                        [ g/kg ]
		#  pt  =  potential temperature (ITS-90)                          [ deg C ]
		#
		#  SA & pt need to have the same dimensions.
		#
		# OUTPUT:
		#  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
		#
		# AUTHOR:
		#  David Jackett, Trevor McDougall and Paul Barker     [ help@teos-10.org ]
		#
		# VERSION NUMBER: 3.05 (27th January 2015)
		#
		# REFERENCES:
		#  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
		#   seawater - 2010: Calculation and use of thermodynamic properties.
		#   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
		#   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
		#    See section 3.3 of this TEOS-10 Manual.
		#
		#  The software is available from http://www.TEOS-10.org
		#
		sfac = 0.0248826675584615;                  # sfac = 1/(40.*(35.16504/35)).

		x2 = sfac.*SA;
		x = sqrt(x2);
		y = pt.*0.025;                              # normalize for F03 and F08.
		gsw_cp0=3991.86795711963;
		pot_enthalpy =  61.01362420681071 + y.*(168776.46138048015 +
		    y.*(-2735.2785605119625 + y.*(2574.2164453821433 +
		    y.*(-1536.6644434977543 + y.*(545.7340497931629 +
		    (-50.91091728474331 - 18.30489878927802.*y).*y))))) +
		    x2.*(268.5520265845071 + y.*(-12019.028203559312 +
		    y.*(3734.858026725145 + y.*(-2046.7671145057618 +
		    y.*(465.28655623826234 + (-0.6370820302376359 -
		    10.650848542359153.*y).*y)))) +
		    x.*(937.2099110620707 + y.*(588.1802812170108 +
		    y.*(248.39476522971285 + (-3.871557904936333 -
		    2.6268019854268356.*y).*y)) +
		    x.*(-1687.914374187449 + x.*(246.9598888781377 +
		    x.*(123.59576582457964 - 48.5891069025409.*x)) +
		    y.*(936.3206544460336 +
		    y.*(-942.7827304544439 + y.*(369.4389437509002 +
		    (-33.83664947895248 - 9.987880382780322.*y).*y))))));

			CT = pot_enthalpy./gsw_cp0;

			return CT
	end


	function gsw_CT_from_t(SA,t,p)

		# gsw_CT_from_t           Conservative Temperature from in-situ temperature
		#
		# USAGE:
		#  CT = gsw_CT_from_t(SA,t,p)
		#
		# DESCRIPTION:
		#  Calculates Conservative Temperature of seawater from in-situ
		#  temperature.
		#
		# INPUT:
		#  SA  =  Absolute Salinity                                        [ g/kg ]
		#  t   =  in-situ temperature (ITS-90)                            [ deg C ]
		#  p   =  sea pressure                                             [ dbar ]
		#         ( i.e. absolute pressure - 10.1325 dbar )
		#
		#  SA & t need to have the same dimensions.
		#  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & t are MxN.
		#
		# OUTPUT:
		#  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
		#
		# AUTHOR:
		#  David Jackett, Trevor McDougall and Paul Barker     [ help@teos-10.org ]
		#
		# VERSION NUMBER: 3.05 (27th January 2015)
		#
		# REFERENCES:
		#  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
		#   seawater - 2010: Calculation and use of thermodynamic properties.
		#   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
		#   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
		#    See section 3.3 of this TEOS-10 Manual.
		#
		#  The software is available from http://www.TEOS-10.org
		#
		#

		pt0 = gsw_pt0_from_t(SA,t,p);
		CT = gsw_CT_from_pt(SA,pt0);
	end



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
	    #sfac = 0.0248826675584615;

        sfac = 1/(40*(35.16504/35)) ;                  # sfac = 1/(40*(35.16504/35)).
	    # offset = 5.971840214030754e-1;
		deltaS = 24;
		offset = deltaS*sfac                      # offset = deltaS*sfac.

	    x2 = sfac.*sa;
	    xs = sqrt.(x2 .+ offset);
	    ys = ct.*0.025;
	    z = p.*1e-4;

	    v000 =  1.0769995862e-3;
	    v001 = -6.0799143809e-5;
	    v002 =  9.9856169219e-6;
	    v003 = -1.1309361437e-6;
	    v004 =  1.0531153080e-7;
	    v005 = -1.2647261286e-8;
	    v006 =  1.9613503930e-9;
	    v010 = -1.5649734675e-5;
	    v011 =  1.8505765429e-5;
	    v012 = -1.1736386731e-6;
	    v013 = -3.6527006553e-7;
	    v014 =  3.1454099902e-7;
	    v020 =  2.7762106484e-5;
	    v021 = -1.1716606853e-5;
	    v022 =  2.1305028740e-6;
	    v023 =  2.8695905159e-7;
	    v030 = -1.6521159259e-5;
	    v031 =  7.9279656173e-6;
	    v032 = -4.6132540037e-7;
	    v040 =  6.9111322702e-6;
	    v041 = -3.4102187482e-6;
	    v042 = -6.3352916514e-8;
	    v050 = -8.0539615540e-7;
	    v051 =  5.0736766814e-7;
	    v060 =  2.0543094268e-7;
	    v100 = -3.1038981976e-4;
	    v101 =  2.4262468747e-5;
	    v102 = -5.8484432984e-7;
	    v103 =  3.6310188515e-7;
	    v104 = -1.1147125423e-7;
	    v110 =  3.5009599764e-5;
	    v111 = -9.5677088156e-6;
	    v112 = -5.5699154557e-6;
	    v113 = -2.7295696237e-7;
	    v120 = -3.7435842344e-5;
	    v121 = -2.3678308361e-7;
	    v122 =  3.9137387080e-7;
	    v130 =  2.4141479483e-5;
	    v131 = -3.4558773655e-6;
	    v132 =  7.7618888092e-9;
	    v140 = -8.7595873154e-6;
	    v141 =  1.2956717783e-6;
	    v150 = -3.3052758900e-7;
	    v200 =  6.6928067038e-4;
	    v201 = -3.4792460974e-5;
	    v202 = -4.8122251597e-6;
	    v203 =  1.6746303780e-8;
	    v210 = -4.3592678561e-5;
	    v211 =  1.1100834765e-5;
	    v212 =  5.4620748834e-6;
	    v220 =  3.5907822760e-5;
	    v221 =  2.9283346295e-6;
	    v222 = -6.5731104067e-7;
	    v230 = -1.4353633048e-5;
	    v231 =  3.1655306078e-7;
	    v240 =  4.3703680598e-6;
	    v300 = -8.5047933937e-4;
	    v301 =  3.7470777305e-5;
	    v302 =  4.9263106998e-6;
	    v310 =  3.4532461828e-5;
	    v311 = -9.8447117844e-6;
	    v312 = -1.3544185627e-6;
	    v320 = -1.8698584187e-5;
	    v321 = -4.8826139200e-7;
	    v330 =  2.2863324556e-6;
	    v400 =  5.8086069943e-4;
	    v401 = -1.7322218612e-5;
	    v402 = -1.7811974727e-6;
	    v410 = -1.1959409788e-5;
	    v411 =  2.5909225260e-6;
	    v420 =  3.8595339244e-6;
	    v500 = -2.1092370507e-4;
	    v501 =  3.0927427253e-6;
	    v510 =  1.3864594581e-6;
	    v600 =  3.1932457305e-5;

	    specvol = (v000 + xs.*(v100 + xs.*(v200 + xs.*(v300 + xs.*(v400 + xs.*(v500 + xs.*v600)))))
	       + ys.*(v010 + xs.*(v110 + xs.*(v210 + xs.*(v310 + xs.*(v410 + xs.*v510))))
	       + ys.*(v020 + xs.*(v120 + xs.*(v220 + xs.*(v320 + xs.*v420)))
	       + ys.*(v030 + xs.*(v130 + xs.*(v230 + xs.*v330))
	       + ys.*(v040 + xs.*(v140 + xs.*v240)
	       + ys.*(v050 + xs.*v150
	       + ys.*v060)))))
	    + z.*(    v001 + xs.*(v101 + xs.*(v201 + xs.*(v301 + xs.*(v401 + xs.*v501))))
	       + ys.*(v011 + xs.*(v111 + xs.*(v211 + xs.*(v311 + xs.*v411)))
	       + ys.*(v021 + xs.*(v121 + xs.*(v221 + xs.*v321))
	       + ys.*(v031 + xs.*(v131 + xs.*v231)
	       + ys.*(v041 + xs.*v141
	       + ys.*v051))))
	    + z.*(    v002 + xs.*(v102 + xs.*(v202 + xs.*(v302 + xs.*v402)))
	       + ys.*(v012 + xs.*(v112 + xs.*(v212 + xs.*v312))
	       + ys.*(v022 + xs.*(v122 + xs.*v222)
	       + ys.*(v032 + xs.*v132
	       + ys.*v042)))
	    + z.*(    v003 + xs.*(v103 + xs.*v203)
	       + ys.*(v013 + xs.*v113
	       + ys.*v023)
	    + z.*(    v004 + xs.*v104
	       + ys.*v014
	    + z.*(    v005
	    + z.*v006)))))
		)



	    # value = v000
	    # + xs*(v010 + xs*(v020 + xs*(v030 + xs*(v040 + xs*(v050
	    # + v060*xs))))) + ys*(v100 + xs*(v110 + xs*(v120 + xs*(v130 + xs*(v140
	    # + v150*xs)))) + ys*(v200 + xs*(v210 + xs*(v220 + xs*(v230 + v240*xs)))
	    # + ys*(v300 + xs*(v310 + xs*(v320 + v330*xs)) + ys*(v400 + xs*(v410
	    # + v420*xs) + ys*(v500 + v510*xs + v600*ys))))) + z*(v001 + xs*(v011
	    # + xs*(v021 + xs*(v031 + xs*(v041 + v051*xs)))) + ys*(v101 + xs*(v111
	    # + xs*(v121 + xs*(v131 + v141*xs))) + ys*(v201 + xs*(v211 + xs*(v221
	    # + v231*xs)) + ys*(v301 + xs*(v311 + v321*xs) + ys*(v401 + v411*xs
	    # + v501*ys)))) + z*(v002 + xs*(v012 + xs*(v022 + xs*(v032 + v042*xs)))
	    # + ys*(v102 + xs*(v112 + xs*(v122 + v132*xs)) + ys*(v202 + xs*(v212
	    # + v222*xs) + ys*(v302 + v312*xs + v402*ys))) + z*(v003 + xs*(v013
	    # + v023*xs) + ys*(v103 + v113*xs + v203*ys) + z*(v004 + v014*xs + v104*ys
	    # + z*(v005 + v006*z)))));

	    return (specvol);
	end
	function gsw_rho(sa, ct, p)

	    return (1.0./gsw_specvol.(sa,ct,p));
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
