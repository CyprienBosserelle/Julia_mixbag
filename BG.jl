"""
    BG module
    Collection of function to create input and process output of the BG model

    Available functions:
    SetEdges, patchgrids

    #Example:
    using BG
    SetEdges(x,y,topo,LBRT,buffer)
"""
module BG

    import Printf, NetCDF, junkinthetrunk
    export SetEdges, patchgrids

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
