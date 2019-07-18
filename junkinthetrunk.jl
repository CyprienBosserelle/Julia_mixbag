module junkinthetrunk

    import StatsBase, NetCDF, Interpolations, Dates, ProgressMeter
    export nanmean, write2nc, pol2cart, interp1


    #nanmean does a mean while ignoring nans
    nanmean(x) = StatsBase.mean(filter(!isnan,x))
    nanmean(x,y) = mapslices(nanmean,x,y)

    #cart2pol
    function pol2cart(theta,speed)
         ee=speed*cosd(90-theta);
         nn=speed.*sind(90-theta);
         return ee,nn;
     end

    interptime(next,prev,timenext,time)=prev + (time) / (timenext)*(next - prev);

    function interp1(xA,yA,newx)
        # Interpolate to time nx even if xA is not monotonic (although it has to be increasing)
        # This is a naive implementation and could do with improvements.
        # This function deals with extrapolation by padding the first/last known value
        #
        y=zeros(length(newx))
        for n=1:length(newx)
            index=findfirst(xA.>newx[n]);
            prev=yA[max(index[1]-1,1)];
            next=yA[index[1]];

            time=newx[n]-xA[max(index[1]-1,1)];
            timenext=xA[max(index[1],1)]-xA[max(index[1]-1,1)];


            y[n]=prev + (time) / (timenext)*(next - prev);
        end
        return y;
    end


    #write2nc writes a matrix to netcdf file
    function write2nc(x,y,z,ncfile)

        xatts = Dict("longname" => "eastings",
          "units"    => "m")
        yatts = Dict("longname" => "northings",
                  "units"    => "m")
        varatts = Dict("longname" => "Topo / Bathy above MVD53",
                  "units"    => "m")
        NetCDF.nccreate(ncfile,"z","x",x,xatts,"y",y,yatts,atts=varatts)
        NetCDF.ncwrite(x,ncfile,"x");
        NetCDF.ncwrite(y,ncfile,"y");
        NetCDF.ncwrite(z,ncfile,"z");
        NetCDF.ncclose(ncfile);
    end


# open(outname,"w") do io
#     for i=1:(nx-1)
#         for j=1:(ny-1)
#             if (mask[i,j]==1)
#                 Printf.@printf(io,">> -Z%f\n",Z[i,j]);
#                 Printf.@printf(io,"%f\t%f\n",X[i,j],Y[i,j]);
#                 Printf.@printf(io,"%f\t%f\n",X[i+1,j],Y[i+1,j]);
#                 Printf.@printf(io,"%f\t%f\n",X[i+1,j+1],Y[i+1,j+1]);
#                 Printf.@printf(io,"%f\t%f\n",X[i,j+1],Y[i,j+1]);
#             end
#        end
#     end

#end




end
