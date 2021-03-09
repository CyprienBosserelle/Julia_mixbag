"""
    DEMproc module
    Collection of function to create DEM for flood/inundation simulation

    Available functions:
    SetEdges, patchgrids, Nestedges, Calcnxny, Checkregion, regrid

    #Example:
    using BG
    SetEdges(x,y,topo,LBRT,buffer)
"""
module DEMproc

    import Printf, NetCDF, junkinthetrunk
    export SetEdges, patchgrids, Nestedges, Calcnxny, Checkregion, regrid, getxy

    # Function to write md files (next gen?)
    """
        Function to set grid edges

        Setting edges means making the n last rows/columns the same values as the n-1 row/column

        Usage:
        SetEdges(x,y,topo,LBRT,buffer)
        Where x and y are vector of floats; topo is an array of floats
        LBRT is a Bool vect whether to set a given edges for Left, Bottom, Right and Top
        LBRT=[true,true,true,true]
        buffer is a vaue between 0 and 1 of the percentge of grid width/height to set [0.01 for 1%]


        Example:
        using BG
        using NetCDF

        infile="2019_DEM_1m_fixed_wt_drains.nc"

        topo=ncread("2019_DEM_1m_fixed_wt_drains.nc","z");
        x=ncread("2019_DEM_1m_fixed_wt_drains.nc","x");
        y=ncread("2019_DEM_1m_fixed_wt_drains.nc","y");

        newtopo=SetEdges(x,y,topo,[true,true,true,true],0.01);

    """
    function SetEdges(x,y,topo,LBRT,buffer)

        toponew =copy(topo);
        if(LBRT[1])

            leftedge=x[1]+buffer*(x[end]-x[1]);
            #1908758.0;
            indexrx=x.<fill(leftedge,size(x));
            indexz=argmin(abs.(x.-fill(leftedge,size(x))));
            indexz=max(indexz,2);
            rgedgetopo=topo[indexz,:];

            for n=1:(indexz-1)
                toponew[n,:]=rgedgetopo;
            end
        end
        if(LBRT[2])

            botedge=y[1]+buffer*(y[end]-y[1]);
            #1908758.0;
            #indexrx=y.<fill(botedge,size(y));
            indexz=argmin(abs.(y.-fill(botedge,size(y))));
            indexz=max(indexz,2);
            rgedgetopo=topo[:,indexz];

            for n=1:(indexz-1)
                toponew[:,n]=rgedgetopo;
            end
        end
        if(LBRT[3])

            rightedge=x[end]-buffer*(x[end]-x[1]);
            #1908758.0;
            #indexrx=x.>fill(rightedge,size(x));
            indexz=argmin(abs.(x.-fill(rightedge,size(x))));
            indexz=min(indexz,length(x)-1);
            rgedgetopo=topo[indexz,:];

            for n=(indexz+1):length(x)
                toponew[n,:]=rgedgetopo;
            end
        end
        if(LBRT[4])
            topedge=y[end]-buffer*(y[end]-y[1]);
            #5863100
            # indexry=y.>fill(topedge,size(y));
            indexz=argmin(abs.(y.-fill(topedge,size(y))));
            indexz=min(indexz,length(y)-1);
            toedgetopo=topo[:,indexz];
            for n=(indexz+1):length(y)
                toponew[:,n]=toedgetopo;
            end
        end

        return toponew
    end
    # """
	# XXX
	# """
    # function rotate(xin,ygrid,zgrid,xo,yo,theta)
    #     #
    #
    #
    #     #print(cos(theta)*($1-xo)-sin(theta)*($2-yo),sin(theta)*($1-xo)+cos(theta)*($2-yo),$3)
    #
    #     return xrot,yrot,zrot
    # end

    """
        #replace Nans in a grid with values from another grid even if the second grid has a different domain and resolution

        Usage:
        patchgrids(xgridprimary,ygridprimary,zgridprimary,xgridsecondary,ygridsecondary,zgridsecondary)

    """
    function patchgrids(xgridprimary,ygridprimary,zgridprimary,xgridsecondary,ygridsecondary,zgridsecondary)
        #replace Nans in a grid with values from another grid even if the second grid has a different domain and resolution
        for i=1:length(xgridprimary)
            for j=1:length(ygridprimary)
                if isnan(zgridprimary[i,j])
                    if (xgridprimary[i]>=xgridsecondary[1]) & (xgridprimary[i]<=xgridsecondary[end]) & (ygridprimary[j]>=ygridsecondary[1]) & (xgridprimary[i]>=xgridsecondary[1]) & (ygridprimary[j]<=ygridsecondary[end])
                        #Interpolate to the nearest values
                        zgridprimary[i,j]=junkinthetrunk.bilinearinterpUG(xgridsecondary,ygridsecondary,zgridsecondary,xgridprimary[i],ygridprimary[j]);


                    end
                end
            end
        end
        return zgridprimary
    end

    """
        Set the edges but for nesting conditions
        Take the secondary bathy data and put it on the edge of the primary DEM and across the buffer
    """
    function Nestedges(xgridprimary,ygridprimary,zgridprimary,xgridsecondary,ygridsecondary,zgridsecondary,LBRT,buffer)
        toponew =copy(zgridprimary);



        if(LBRT[1])

            leftedge=xgridprimary[1]+buffer*(xgridprimary[end]-xgridprimary[1]);
            indexz=argmin(abs.(xgridprimary.-fill(leftedge,size(xgridprimary))));
            indexz=max(indexz,2);


            rgedgetopo=zgridprimary[indexz,:]; #Fill first with original data
            i=1;
            for j=1:length(ygridprimary)

                    if (xgridprimary[i]>=xgridsecondary[1]) & (xgridprimary[i]<=xgridsecondary[end]) & (ygridprimary[j]>=ygridsecondary[1]) & (xgridprimary[i]>=xgridsecondary[1]) & (ygridprimary[j]<=ygridsecondary[end])
                        #Interpolate to the nearest values
                        rgedgetopo[j]=junkinthetrunk.bilinearinterpUG(xgridsecondary,ygridsecondary,zgridsecondary,xgridprimary[i],ygridprimary[j]);


                    end

            end

            for n=1:(indexz-1)
                toponew[n,:]=rgedgetopo;
            end


            #
        end
        if(LBRT[2])

            botedge=ygridprimary[1]+buffer*(ygridprimary[end]-ygridprimary[1]);

            indexz=argmin(abs.(ygridprimary.-fill(botedge,size(ygridprimary))));
            indexz=max(indexz,2);
            rgedgetopo=zgridprimary[:,indexz];

            j=1;
            for i=1:length(xgridprimary)

                    if (xgridprimary[i]>=xgridsecondary[1]) & (xgridprimary[i]<=xgridsecondary[end]) & (ygridprimary[j]>=ygridsecondary[1]) & (xgridprimary[i]>=xgridsecondary[1]) & (ygridprimary[j]<=ygridsecondary[end])
                        #Interpolate to the nearest values
                        rgedgetopo[i]=junkinthetrunk.bilinearinterpUG(xgridsecondary,ygridsecondary,zgridsecondary,xgridprimary[i],ygridprimary[j]);


                    end

            end

            for n=1:(indexz-1)
                toponew[:,n]=rgedgetopo;
            end
        end
        if(LBRT[3])

            rightedge=xgridprimary[end]-buffer*(xgridprimary[end]-xgridprimary[1]);

            indexz=argmin(abs.(xgridprimary.-fill(rightedge,size(xgridprimary))));
            indexz=min(indexz,length(xgridprimary)-1);
            rgedgetopo=zgridprimary[indexz,:];

            i=length(xgridprimary);
            for j=1:length(ygridprimary)

                    if (xgridprimary[i]>=xgridsecondary[1]) & (xgridprimary[i]<=xgridsecondary[end]) & (ygridprimary[j]>=ygridsecondary[1]) & (xgridprimary[i]>=xgridsecondary[1]) & (ygridprimary[j]<=ygridsecondary[end])
                        #Interpolate to the nearest values
                        rgedgetopo[j]=junkinthetrunk.bilinearinterpUG(xgridsecondary,ygridsecondary,zgridsecondary,xgridprimary[i],ygridprimary[j]);


                    end

            end

            for n=(indexz+1):length(xgridprimary)
                toponew[n,:]=rgedgetopo;
            end
        end
        if(LBRT[4])
            topedge=ygridprimary[end]-buffer*(ygridprimary[end]-ygridprimary[1]);

            indexz=argmin(abs.(ygridprimary.-fill(topedge,size(ygridprimary))));
            indexz=min(indexz,length(ygridprimary)-1);
            toedgetopo=zgridprimary[:,indexz];

            j=length(ygridprimary);
            for i=1:length(xgridprimary)

                    if (xgridprimary[i]>=xgridsecondary[1]) & (xgridprimary[i]<=xgridsecondary[end]) & (ygridprimary[j]>=ygridsecondary[1]) & (xgridprimary[i]>=xgridsecondary[1]) & (ygridprimary[j]<=ygridsecondary[end])
                        #Interpolate to the nearest values
                        rgedgetopo[i]=junkinthetrunk.bilinearinterpUG(xgridsecondary,ygridsecondary,zgridsecondary,xgridprimary[i],ygridprimary[j]);


                    end

            end


            for n=(indexz+1):length(ygridprimary)
                toponew[:,n]=toedgetopo;
            end
        end

        return toponew



    end
    """
    Regrid a cartesian grid to another location/resolution
    """
    function regrid(xin, yin, zin, regionout::NTuple{4,AbstractFloat}, dxout::AbstractFloat)
        #
        region=Checkregion(regionout,dxout);

        dxin=(xin[2]-xin[1]);

        if dxin<dxout
            regrid_coarsen(xin, yin, zin, regionout, dxout)
        else
            regrid_refine(xin, yin, zin, regionout, dxout)
        end
    end

    """
    Regrid to a finer resolution
    """
    function regrid_refine(xin, yin, zin, regionout, dxout)
        ## Perform bilinear interpolation

        region=Checkregion(regionout,dxout);
        nx,ny=Calcnxny(region,dxout);

        z=fill(NaN,(nx,ny));

        for i=1:nx
            for j=1:ny
                y=region[3]+(j-1)*dxout;
                x=region[1]+(i-1)*dxout;
                z[i,j]=bilinearinterpUG(xin,yin,zin,x,y);
            end
        end

        return z


    end

    function regrid_coarsen(xin, yin, zin, regionout, dxout)
        ## Block average the input grid
        region=Checkregion(regionout,dxout);
        nx,ny=Calcnxny(region,dxout);

        z=fill(NaN,(nx,ny));

        for i=1:nx
            for j=1:ny
                y=region[3]+(j-1)*dxout;
                x=region[1]+(i-1)*dxout;

                ii=(xin .> (x - 0.5*dxout)) .& (xin .<= ( x + 0.5*dxout ));
                jj=(yin .> (y - 0.5*dxout)) .& (yin .<= ( y + 0.5*dxout ));

                z[i,j]= junkinthetrunk.nanmean(zin[ii,jj]);
            end
        end

        return z;

    end

    function getxy(region::NTuple{4,AbstractFloat},res::AbstractFloat)
        nx,ny=Calcnxny(region,res);

        x=fill(NaN,(nx));
        y=fill(NaN,(ny));

        for j = 1:ny
            y[j]=region[3]+(j-1)*res;
        end
        for i=1:nx
            x[i]=region[1]+(i-1)*res;
        end

        return x,y
    end

    function Checkregion(region::NTuple{4,AbstractFloat},res::AbstractFloat)
        xmin=region[1];
        ymin=region[3];
        xmax=xmin + ceil((region[2] - xmin) / res) * res;
        ymax=ymin + ceil((region[4] - ymin) / res) * res;

        regfixed=(xmin,xmax,ymin,ymax)
        return regfixed
    end


    function Calcnxny(region::NTuple{4,AbstractFloat},res::AbstractFloat)

        xmin=region[1];
        ymin=region[3];
        xmax=region[2];
        ymax=region[4];

        nx=Int((xmax - xmin) / res) + 1
        ny=Int((ymax - ymin) / res) + 1

        return nx,ny
    end

    function bilinearinterpBase(q11,q12,q21,q22,x1,x2,y1,y2,x,y)
        x2x1 = x2 - x1;
		y2y1 = y2 - y1;
		x2x = x2 - x;
		y2y = y2 - y;
		yy1 = y - y1;
		xx1 = x - x1;
		return 1.0 / (x2x1 * y2y1) * (q11 * x2x * y2y +	q21 * xx1 * y2y + q12 * x2x * yy1 +	q22 * xx1 * yy1	);
	end

	"""
		Bilinearinterpolation on a regular (non-uniform) grid

		Usage: z=bilinearinterpUG(xgrid,ygrid,zb,x,y)
		xgrid and ygrid are expected to be vectors
	"""
	function bilinearinterpUG(xgrid,ygrid,zb,x,y)
		#Bilinearinterpolation on a regular (non-uniform) grid
		#xgrid and ygrid are expected to be vectors
		xmax=xgrid[end];
		xo=xgrid[1];
		ymax=ygrid[end];
		yo=ygrid[1];

		nx=length(xgrid);
		ny=length(ygrid);

		# Extrapolation safeguard
        # This way any extrapolation is clamped to that nearest value
		x = max(min(x, xmax), xo);
		y = max(min(y, ymax), yo);

		xdiff=x.-xgrid;
		cfi=argmin(abs.(xdiff));

		if(xdiff[cfi]<0.)
			cfi=cfi-1;
		end

		cfi=min(max(cfi, 1), nx - 1);

		cfip = cfi + 1;

		x1 = xgrid[cfi];
		x2 = xgrid[cfip];

		ydiff=y.-ygrid;
		cfj=argmin(abs.(ydiff));
		if(ydiff[cfj]<0.)
			cfj=cfj-1;
		end

		cfj=min(max(cfj, 1), ny - 1);

		cfjp = cfj + 1;

		y1 = ygrid[cfj];
		y2 = ygrid[cfjp];

		q11 = zb[cfi,cfj];
		q21 = zb[cfip,cfj];
		q12 = zb[cfi,cfjp];
		q22 = zb[cfip,cfjp];

		return bilinearinterpBase(q11,q12,q21,q22,x1,x2,y1,y2,x,y);
	end


    ## Add paintnan function

    # function SetEdges(x,y,topo)
    #     LBRT=[true,true,true,true];
    #     buffer=0.01;
    #     return SetEdges(x,y,topo,LBRT,buffer)
    #
    # end
    # function SetEdges(x,y,intopo,LBRT)
    #     buffer=0.01;
    #     return SetEdges(x,y,topo,LBRT,buffer)
    #
    # end
    # function SetEdges(x,y,intopo,buffer::bool)
    #     LBRT=[true,true,true,true];
    #     return SetEdges(x,y,topo,LBRT,buffer)
    #
    # end

    #function readBGbathy(filename)



end
