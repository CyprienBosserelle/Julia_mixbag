module BG

    import Printf, NetCDF, junkinthetrunk
    export SetEdges, patchgrids

    # Function to write md files (next gen?)
    function SetEdges(x,y,topo,LBRT,buffer)

        toponew =copy(topo);
        if(LBRT[1])

            leftedge=x[1]+buffer*(x[end]-x[1]);
            #1908758.0;
            indexrx=x.<fill(leftedge,size(x));
            indexz=findmin(abs.(x.-fill(leftedge,size(x))));

            rgedgetopo=topo[indexz[2],:];

            for n=1:(length(findall(indexrx)))
                toponew[indexz[2]+n,:]=rgedgetopo;
            end
        end
        if(LBRT[2])

            botedge=y[1]+buffer*(y[end]-y[1]);
            #1908758.0;
            indexrx=y.<fill(botedge,size(y));
            indexz=findmin(abs.(y.-fill(botedge,size(y))));

            rgedgetopo=topo[:,indexz[2]];

            for n=1:(length(findall(indexrx)))
                toponew[:,indexz[2]+n]=rgedgetopo;
            end
        end
        if(LBRT[3])

            rightedge=x[end]-buffer*(x[end]-x[1]);
            #1908758.0;
            indexrx=x.>fill(rightedge,size(x));
            indexz=findmin(abs.(x.-fill(rightedge,size(x))));

            rgedgetopo=topo[indexz[2],:];

            for n=1:(length(findall(indexrx)))
                toponew[indexz[2]+n,:]=rgedgetopo;
            end
        end
        if(LBRT[4])
            topedge=y[end]-buffer*(y[end]-y[1]);
            #5863100
            indexry=y.>fill(topedge,size(y));
            indexz=findmin(abs.(y.-fill(topedge,size(y))));

            toedgetopo=topo[:,indexz[2]];
            for n=1:(length(findall(indexry)))
                toponew[:,indexz[2]+n]=toedgetopo;
            end
        end

        return toponew
    end

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
