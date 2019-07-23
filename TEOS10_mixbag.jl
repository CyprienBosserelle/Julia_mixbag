


function gsw_SP_from_SA(SA,p,long,lat)
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
    
