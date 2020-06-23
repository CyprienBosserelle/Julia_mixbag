
"""
 Tsunami module
 Collection of function useful for tsunami simulations.

     Cyprien Bosserelle 2020

"""
module Tsunami

    import Okada
    export faultparam,faultkm2m!,rotatexy,unrotatexy,unrotatexyCompass,rotatexyCompass,sphericDist,sphericOffset,mvBLref2centroid!,emptygrid,cartsphdist2eq,cartdistance2eq,CalcMw,Mw2slip,Calcslip!,Okadavert


    mutable struct faultparam
        lon::Float64
        lat::Float64
        length::Float64
        width::Float64
        depth::Float64
        strike::Float64
        dip::Float64
        rake::Float64
        slip::Float64
        tinit::Float64
        trise::Float64
    end

    function Okadavert(ef,nf,fault::faultparam)
        return Okada.Okada85vert(ef,nf,fault.depth,fault.strike,fault.dip,fault.length,fault.width,fault.rake,fault.slip,0);
    end

    function CalcMw(fault::faultparam)
        # Calculate Mw based on fault dimention and slip
        Mo=4.0e11*fault.slip*100*fault.length*100*fault.width*100;
        Mw=2/3*log10(Mo)-10.7;
        return Mw
    end

    function Mw2slip(Mw,length,width)
        #Calculate slip based on Mw and fault dimension
        Mo=10^((Mw+10.7).*(3.0/2.0))
        slip=Mo/(4.0e11*100*length*100*width*100);


        return slip
    end

    function Calcslip!(fault::faultparam,Mw)
        slip=Mw2slip(Mw,fault.length,fault.width)
        fault.slip=slip;
        return fault
    end


    function faultkm2m!(fault::faultparam)
        #modify fault length and width from km to m
        fault.length*=1000.0;
        fault.width*=1000.0;
        fault.depth*=1000.0;
        return fault
    end

    function rotatexy(x,y,xo,yo,theta)
        # Rotate x y coordinate by theta.
        # theta is in radian following the maths convention

        newx=cos(theta)*(x-xo)-sin(theta)*(y-yo);
        newy=sin(theta)*(x-xo)+cos(theta)*(y-yo);
        return newx,newy
    end

    function unrotatexy(x,y,xo,yo,theta)
        ## Rotate x y coordinate by theta.
        # theta is in radian following the maths convention
        newx=cos(-1*theta)*(x)-sin(-1*theta)*(y)+xo;
        newy=sin(-1*theta)*(x)+cos(-1*theta)*(y)+yo;

        return newx,newy
    end

    function unrotatexyCompass(x,y,xo,yo,alpha)
        ## Rotate x y coordinate by theta.
        # alpha is in degrees following the compass convention

        theta=-1.0*deg2rad(alpha);
        return unrotatexy(x,y,xo,yo,theta);
    end

    function rotatexyCompass(x,y,xo,yo,alpha)
        # Rotate x y coordinate by alpha.
        # alpha is in degrees following the compass convention
        theta=-1.0*deg2rad(alpha);
        return rotatexy(x,y,xo,yo,theta);

    end





    function sphericDist(lon1,lat1,lon2,lat2)
        # Calculate spherical distance between 2 pts
        # distance is in m

        earthradius=6356750.52;
        #  Determine proper longitudinal shift.
        delta=lon2-lon1;
        l=abs(delta);
        l=l>=180 ? 360-l : l;
        # Convert Decimal degrees to radians.
        beta1 = deg2rad(lat1);
        beta2 = deg2rad(lat2);
        l = deg2rad(l);
        # Calculate S/Bo subformulas.
        st = sqrt(((sin(l).*cos(beta2)).^2)+(((sin(beta2).*cos(beta1))-(sin(beta1).*cos(beta2).*cos(l))).^2));

        return asin(st)*earthradius
    end

    function sphericOffset(lon,lat,de,dn)
        # Offset lat and lon by given eats and north offsets meters
        # de: delta east;  dn: delta north  offsets are in meters



        earthradius=6356750.52;

        bearing=atan(dn,de);

        dist=hypot(dn,de);

        # this is a bit rubbish because it assumes a flat earth!
        # # Coordinate offsets in radians
        # dLat = dn/R;
        # dLon = de/(R*cos(deg2rad(lat)));
        #
        # # OffsetPosition, decimal degrees
        # latO = lat + rad2deg(dLat);
        # lonO = lon + rad2deg(dLon);
        beta1=deg2rad(lat);
        alpha1=deg2rad(lon);

        latO=asin(sin(beta1)*cos(dist/earthradius)+cos(beta1)*sin(dist/earthradius)*cos(bearing));
        lonO=lon+rad2deg(atan(sin(bearing)*sin(dist/earthradius)*cos(beta1),cos(dist/earthradius)-sin(beta1)*sin(latO)))

        return lonO,rad2deg(latO)
    end



    function mvBLref2centroid!(fault::faultparam)
        # move reference coordinates and depth
        # from bottom left corner (relative to the strike: for stroke of zero it is the south west corner, for stike of 180 it would be the north east corner)
        # to centroid of the fault plane
        fault.depth=sind(fault.dip)*fault.width;

        # Moving to the centroid is a simple rotation problem

        newX=0.5*fault.width
        newY=0.5*fault.length

        rotX,rotY=rotatexyCompass(newX,newY,0.0,0.0,fault.strike);

        # bearing=rad2deg(atan(rotY,rotX));

        # Warning this need to be in meters!!!
        newlon,newlat=sphericOffset(fault.lon,fault.lat,rotX,rotY);



        fault.lon=newlon;
        fault.lat=newlat;
        return fault


    end

    function emptygrid(xo,xmax,yo,ymax,inc)
        # porduce arrays of lat,lon/easting,northing representing an empty grid

        x=xo:inc:xmax
        y=yo:inc:ymax

        xx=collect(x);
        yy=collect(y);

        E=xx*ones(size(yy))'
        N=collect(transpose(yy*ones(size(xx))'))

        return E,N
    end

    function cartdistance2eq(e,n,eqx,eqy)
        de=e-eqx;
        dn=n-eqy;
    end

    function cartsphdist2eq(e,n,eqx,eqy)
        facE=ones(size(e));
        facE[e .< eqx] .= -1.0;
        facN=ones(size(n));
        facN[n .< eqy] .= -1.0;

        ef=sphericDist.(e,eqy,eqx,eqy).*facE;
        nf=sphericDist.(eqx,n,eqx,eqy).*facN;

        return ef,nf
    end
end
