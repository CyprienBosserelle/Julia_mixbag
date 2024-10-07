

using Dates

function moon_age_location(jd::Float64)
    # pinched from Astro.jl
    earth_equ_radius = 6378.137
    v  = (jd - 2451550.1) / 29.530588853
    ip = v - floor(Integer, v)
    ag = ip * 29.530588853   # Moon's age from new moon in days
    ip = ip * 2pi            # Converts to radian

    # Calculate distance from anomalistic phase
    v= (jd - 2451562.2) / 27.55454988
    dp = v- floor(Integer, v)
    dp = dp * 2pi
    di = 60.4 - 3.3 * cos(dp) - .6 * cos(2 * ip - dp) - .5 * cos(2 * ip)

    # Calculate ecliptic latitude from nodal (draconic) phase
    v = (jd - 2451565.2) / 27.212220817
    np = v - floor(Integer, v)
    np = np * 2pi
    la = 5.1 * sin(np)

    # Calculate ecliptic longitude from sidereal motion
    v = (jd - 2451555.8) / 27.321582241
    rp = v - floor(Integer, v)
    lo = 360 * rp + 6.3 * sin(dp) + 1.3 * sin(2 * ip - dp) + .7 * sin(2 * ip)

    return (ag, di * earth_equ_radius, deg2rad(la), deg2rad(lo))
end


phase_date=Dates.now()



fixangle(a) = a - 360.0 * floor(a/360.0)

torad(d)=d * pi / 180.0;
todeg(r)= r * 180.0 / pi;
dsin(d)=sin(torad(d));
dcos(d)=cos(torad(d));


function kepler(m, ecc)
    #Solve the equation of Kepler.

    epsilon = 1e-6

    m = torad(m)
    e = m
    while true
        delta = e - ecc * sin(e) - m
        e = e - delta / (1.0 - ecc * cos(e))

        if abs(delta) <= epsilon
            break
        end
    end

    return e
end

function phasemoon(phase_date)
    #Calculate phase of moon as a fraction:
    #The argument is the time for which the phase is requested,
    #
    #expressed in either a DateTime
    # Translated from:
    # moon.py, based on code by John Walker (http://www.fourmilab.ch/)
    # ported to Python by Kevin Turner <acapnotic@twistedmatrix.com>
    # on June 6, 2001 (JDN 2452066.52491), under a full moon.
    #
    # This program is in the public domain: "Do what thou wilt shall be
    # the whole of the law".
    # Ported to Julia by Cyprien Bosserelle

    # Calculation of the Sun's position
    c_epoch = 2444238.5;
    # Ecliptic longitude of the Sun at epoch 1980.0
    ecliptic_longitude_epoch = 278.833540

    # Ecliptic longitude of the Sun at perigee
    ecliptic_longitude_perigee = 282.596403

    # Eccentricity of Earth's orbit
    eccentricity = 0.016718

    # Semi-major axis of Earth's orbit, in kilometers
    sun_smaxis = 1.49585e8

    # Sun's angular size, in degrees, at semi-major axis distance
    sun_angular_size_smaxis = 0.533128

    ## Elements of the Moon's orbit, epoch 1980.0

    # Moon's mean longitude at the epoch
    moon_mean_longitude_epoch = 64.975464
    # Mean longitude of the perigee at the epoch
    moon_mean_perigee_epoch = 349.383063

    # Mean longitude of the node at the epoch
    node_mean_longitude_epoch = 151.950429

    # Inclination of the Moon's orbit
    moon_inclination = 5.145396

    # Eccentricity of the Moon's orbit
    moon_eccentricity = 0.054900

    # Moon's angular size at distance a from Earth
    moon_angular_size = 0.5181

    # Semi-mojor axis of the Moon's orbit, in kilometers
    moon_smaxis = 384401.0
    # Parallax at a distance a from Earth
    moon_parallax = 0.9507

    # Synodic month (new Moon to new Moon), in days
    synodic_month = 29.53058868

    # Base date for E. W. Brown's numbered series of lunations (1923 January 16)
    lunations_base = 2423436.0

    ## Properties of the Earth

    earth_radius = 6378.16

    # date within the epoch
    day = Dates.datetime2julian(phase_date) - c_epoch

    # Mean anomaly of the Sun
    N = fixangle((360/365.2422) * day)
    # Convert from perigee coordinates to epoch 1980
    M = fixangle(N + ecliptic_longitude_epoch - ecliptic_longitude_perigee)

    # Solve Kepler's equation
    Ec = kepler(M, eccentricity)
    Ec = sqrt((1 + eccentricity) / (1 - eccentricity)) * tan(Ec/2.0)
    # True anomaly
    Ec = 2 * todeg(atan(Ec))
    # Suns's geometric ecliptic longuitude
    lambda_sun = fixangle(Ec + ecliptic_longitude_perigee)

    # Orbital distance factor
    F = ((1 + eccentricity * cos(torad(Ec))) / (1 - eccentricity^2))

    # Distance to Sun in km
    sun_dist = sun_smaxis / F
    sun_angular_diameter = F * sun_angular_size_smaxis

    ########
    #
    # Calculation of the Moon's position

    # Moon's mean longitude
    moon_longitude = fixangle(13.1763966 * day + moon_mean_longitude_epoch)

    # Moon's mean anomaly
    MM = fixangle(moon_longitude - 0.1114041 * day - moon_mean_perigee_epoch)

    # Moon's ascending node mean longitude
    # MN = fixangle(c.node_mean_longitude_epoch - 0.0529539 * day)

    evection = 1.2739 * sin(torad(2*(moon_longitude - lambda_sun) - MM))

    # Annual equation
    annual_eq = 0.1858 * sin(torad(M))

    # Correction term
    A3 = 0.37 * sin(torad(M))

    MmP = MM + evection - annual_eq - A3

    # Correction for the equation of the centre
    mEc = 6.2886 * sin(torad(MmP))

    # Another correction term
    A4 = 0.214 * sin(torad(2 * MmP))

    # Corrected longitude
    lP = moon_longitude + evection + mEc - annual_eq + A4

    # Variation
    variation = 0.6583 * sin(torad(2*(lP - lambda_sun)))

    # True longitude
    lPP = lP + variation

    #
    # Calculation of the Moon's inclination
    # unused for phase calculation.

    # Corrected longitude of the node
    # NP = MN - 0.16 * sin(torad(M))

    # Y inclination coordinate
    # y = sin(torad(lPP - NP)) * cos(torad(c.moon_inclination))

    # X inclination coordinate
    # x = cos(torad(lPP - NP))

    # Ecliptic longitude (unused?)
    # lambda_moon = todeg(atan2(y,x)) + NP

    # Ecliptic latitude (unused?)
    # BetaM = todeg(asin(sin(torad(lPP - NP)) * sin(torad(c.moon_inclination))))

    #######
    #
    # Calculation of the phase of the Moon

    # Age of the Moon, in degrees
    moon_age = lPP - lambda_sun

    # Phase of the Moon
    moon_phase = (1 - cos(torad(moon_age))) / 2.0





    # Calculate distance of Moon from the centre of the Earth
    moon_dist = (moon_smaxis * (1 - moon_eccentricity^2)) / (1 + moon_eccentricity * cos(torad(MmP + mEc)))

    # Calculate Moon's angular diameter
    moon_diam_frac = moon_dist / moon_smaxis
    moon_angular_diameter = moon_angular_size / moon_diam_frac

    return fixangle(moon_age) / 360.0, synodic_month * fixangle(moon_age) / 360.
end



#Demo
# phasemoon(DateTime(2019,8,30,11))
#
# moon_age_location(Dates.datetime2julian(Dates.now()))
