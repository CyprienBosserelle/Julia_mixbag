"""
    BG module
    Collection of function to create input and process output of the BG model

    Available functions:
    SetEdges, patchgrids

    #Example:
    using BG

"""
module BG

    import Printf, NetCDF, junkinthetrunk, GMT
    export CalcRunup, plotvar


    function scanminmax(file,var,stp)


        zmin=1/eps(eps(1.0))
        zmax=-1.0.*zmin
        levmax=NetCDF.ncgetatt(file, "Global", "maxlevel");
        levmin=NetCDF.ncgetatt(file, "Global", "minlevel");
        for level=levmin:levmax

            if level<0
                stlv="N";
            else
                stlv="P";
            end
            fvar=var*"_"*stlv*string(abs(level))

            scalefac=NetCDF.ncgetatt(file, fvar, "scale_factor");
            addoffset=NetCDF.ncgetatt(file, fvar, "add_offset");
            missval=NetCDF.ncgetatt(file, fvar, "_FillValue");



            zz=NetCDF.ncread(file,fvar,[1,1,max(stp,1)], [-1,-1,1]);
            if any([!isnothing(scalefac) !isnothing(addoffset)])
                zz=zz.*scalefac.+addoffset;
            end
            zz[zz.==missval].=NaN;

            zz=zz[.!isnan.(zz)];


            if !isempty(zz)
                zmin=min(zmin,minimum(zz[.!isnan.(zz)]));
                zmax=max(zmax,maximum(zz[.!isnan.(zz)]));
            end



        end
        if zmin==zmax
            zz=zmin
            zmin=zz-1.0;
            zmax=zz+1.0;
        end


        return zmin,zmax
    end

    function plotblocks(file,var,step; showid=false)

        # Draw the blocks
        blkxo = NetCDF.ncread(file, "blockxo")
        blkyo = NetCDF.ncread(file, "blockyo")
        blkwidth= NetCDF.ncread(file, "blockwidth") .*16.0;
        dx=NetCDF.ncread(file, "blockwidth");


        nblk=length(blkxo);

        for ib=1:nblk
            rect = [blkxo[ib].-dx[ib]*0.5 blkyo[ib].-dx[ib]*0.5; blkxo[ib]+blkwidth[ib].-dx[ib]*0.5 blkyo[ib].-dx[ib]*0.5; blkxo[ib]+blkwidth[ib].-dx[ib]*0.5 blkyo[ib]+blkwidth[ib].-dx[ib]*0.5; blkxo[ib].-dx[ib]*0.5 blkyo[ib]+blkwidth[ib].-dx[ib]*0.5; blkxo[ib].-dx[ib]*0.5 blkyo[ib].-dx[ib]*0.5];
            GMT.plot!(rect, lw=1,t=80);
            if showid
                GMT.text!(GMT.text_record([blkxo[ib] blkyo[ib]], [string(ib)]), attrib=(font=(8,"Helvetica",:black),angle=0,justify=:LM));
            end
        end
    end

    function plotvar(file,var,step,figname;region=(NaN,NaN,NaN,NaN),zrange=(NaN,NaN), plotblock=true, plotid=false)
        #

        levmax=NetCDF.ncgetatt(file, "Global", "maxlevel");
        levmin=NetCDF.ncgetatt(file, "Global", "minlevel");

        xmin=isnan(region[1]) ? NetCDF.ncgetatt(file, "Global", "xmin") : region[1];
        xmax=isnan(region[2]) ? NetCDF.ncgetatt(file, "Global", "xmax") : region[2];
        ymin=isnan(region[3]) ? NetCDF.ncgetatt(file, "Global", "ymin") : region[3];
        ymax=isnan(region[4]) ? NetCDF.ncgetatt(file, "Global", "ymax") : region[4];


        global zmin=9999.0
        global zmax=-9999.0
        #zrange[2]

        if any([isnan(zrange[1]) isnan(zrange[2])])
            zmin,zmax=scanminmax(file,var,step)
        else
            zmin=zrange[1];
            zmax=zrange[2];
        end

        #regplot=(xmin,xmax,ymin,ymax)

        topo = GMT.makecpt(color=:jet, range=(zmin,zmax), continuous=true)
        GMT.basemap(proj=:linear, region=(xmin,xmax,ymin,ymax), frame=(annot=:auto, ticks=:auto))
        #basemap(proj=:linear, figsize=(24.5,1.0), region=(xmin,xmax,ymin,ymax))
        for level=levmin:levmax

            if level<0
                stlv="N";
            else
                stlv="P";
            end

            GMT.grdimage!(file*"?"*var*"_"*stlv*string(abs(level))*"["*string(step)*"]", frame=:a, color=topo, Q=true)
        end
        if plotblock
            plotblocks(file,var,step,showid=plotid);
        end


        GMT.colorbar!(pos=(anchor=:TC,length=(15.5,0.6), horizontal=true, offset=(0,1.0)),
                  color=topo, frame=(ylabel=:m,),savefig=figname*".png" , fmt=:png, show=true)





    end


    """
    ttt
    """
    function CalcRunup(file)
        ncfile=NetCDF.open(file);


        fvar="zsmax"




        scalefac=NetCDF.ncgetatt(file, fvar, "scale_factor");
        addoffset=NetCDF.ncgetatt(file, fvar, "add_offset");
        missval=NetCDF.ncgetatt(file, fvar, "_FillValue");



        zz=NetCDF.ncread(file,fvar,[1,1,max(stp,1)], [-1,-1,1]);
        if any([!isnothing(scalefac) !isnothing(addoffset)])
            zz=zz.*scalefac.+addoffset;
        end
        zz[zz.==missval].=NaN;

        zz=zz[.!isnan.(zz)];
        #
    end


end
