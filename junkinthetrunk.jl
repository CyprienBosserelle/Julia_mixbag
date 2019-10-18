"""
    junkinthetrunk module
    Collection of random usefull functions

    Available functions:
    nanmean, write2nc, pol2cart, interp1, nearneighb1, bilinearinterpUG, skills


"""
module junkinthetrunk

    import StatsBase, NetCDF, Dates
    export nanmean, write2nc, pol2cart, interp1, nearneighb1, bilinearinterpUG, skills


    #nanmean does a mean while ignoring nans
	"""
	    average of a vector ignoring NaNs
	    usage nanmean(x)
	"""
    nanmean(x) = StatsBase.mean(filter(!isnan,x))
    nanmean(x,y) = mapslices(nanmean,x,y)


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
		x = max(min(x, xmax), xo);
		y = max(min(y, ymax), yo);

		xdiff=x.-xgrid;
		cfi=argmin(abs.(xdiff));

		if(xdiff[cfi]>0.)
			cfi=cfi-1;
		end

		cfi=min(max(cfi, 1), nx - 1);

		cfip = cfi + 1;

		x1 = xgrid[cfi];
		x2 = xgrid[cfip];

		ydiff=y.-ygrid;
		cfj=argmin(abs.(ydiff));
		if(ydiff[cfj]>0.)
			cfj=cfj-1;
		end

		cfj=min(max(cfj, 1), ny - 1);

		cfjp = cfj + 1;

		y1 = ygrid[cfj];
		y2 = ygrid[cfjp];

		q11 = zb[cfi,cfj];
		q12 = zb[cfi,cfjp];
		q21 = zb[cfip,cfj];
		q22 = zb[cfip,cfjp];

		return bilinearinterpBase(q11,q12,q21,q22,x1,x2,y1,y2,x,y);
	end






    #cart2pol
    function pol2cart(theta,speed)
         ee=speed*cosd(90-theta);
         nn=speed.*sind(90-theta);
         return ee,nn;
     end

    interptime(next,prev,timenext,time)=prev + (time) / (timenext)*(next - prev);

    function interp1old(xA,yA,newx)
        # Interpolate to time nx even if xA is not monotonic (although it has to be increasing)
        # This is a naive implementation and could do with improvements.
        # This function deals with extrapolation by padding the first/last known value
        #
        y=zeros(length(newx))
        for n=1:length(newx)
            index=findfirst(xA.>newx[n]);
            if isnothing(index)
                if newx[n]<=xA[1]
                    index=1;
                elseif newx[n]>=xA[end]
                    index=length(xA)
                end
            end
            prev=yA[max(index[1]-1,1)];
            next=yA[index[1]];

            time=newx[n]-xA[max(index[1]-1,1)];
            timenext=xA[max(index[1],1)]-xA[max(index[1]-1,1)];


            y[n]=prev + (time) / (timenext)*(next - prev);
        end
        return y;
    end

    function interp1(xA, yA, newx)
        # Interpolate to time nx even if xA is not monotonic (although it has to be increasing)
        # This is a naive implementation and could do with improvements.
        # This function deals with extrapolation by padding the first/last known value
        #

        #selcasttimedep[indnonan], error[indnonan], ttt

        #xA=selcasttimedep[indnonan]
        #yA=error[indnonan]
        #newx=copy(ttt)


        y = zeros(length(newx))
        for n = 1:length(newx)
            index = findfirst(xA .> newx[n]);
            if index==nothing
                if newx[n] <= xA[1]
                    index = 1;
                elseif newx[n] >= xA[end]
                    index = length(xA)+1;
                end
            end
            prev = yA[max(index[1] - 1, 1)];
            next = yA[min(index[1],length(xA))];

            time = newx[n] - xA[max(index[1] - 1, 1)];
            timenext = xA[max(min(index[1],length(xA)), 1)] - xA[max(index[1] - 1, 1)];

            if max(min(index[1],length(xA)), 1) == max(index[1] - 1, 1)

                y[n] = yA[max(min(index[1],length(xA)), 1)];
            else
                y[n] = prev + (time) / (timenext) * (next - prev);

            end
        end
        return y;
    end

    function nearneighb1(xA,yA,newx,dxmax)
        # Find the nearest neighbour within and assign value
        #
        y = zeros(length(newx));
        for n = 1:length(newx)
            minim = minimum(abs.(xA .- newx[n]));
            minind = argmin(abs.(xA .- newx[n]));
            if(minim<=dxmax)
                y[n]=yA[minind];
            else
                y[n]=NaN;
            end



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


    """
        Evaluate a model skills
		usage:
        RMS,Bias,Wcorr,Bss=skills(tmeas,xmeas,tmodel,xmodel)
    """
    function skills(tmeas,xmeas,tmodel,xmodel)
        # Calculate RMS. Bias, Index of agreement, and skill score
        #

        #crop non overlapping
        index=(tmodel.>=tmeas[1]) .& (tmodel.<=tmeas[end]);

        tmod=tmodel[index];
        xmod=xmodel[index];

        indxmea=(tmeas.>=tmod[1]) .& (tmeas.<=tmod[end]);
        tmea=tmeas[indxmea];
        xmea=xmeas[indxmea];


        xint=interp1(tmod,xmod,tmea);


        RMS=sqrt(StatsBase.mean((xint.-xmea).^2));
        Bias=StatsBase.mean(xint)-StatsBase.mean(xmea);
        Wcorr=1-(sum((xmea-xint).^2)/(sum((abs.(xmea.-StatsBase.mean(xint)).+abs.(xint.-StatsBase.mean(xint))).^2)));
        Bss=1-(StatsBase.var((xmea-xint))/(StatsBase.var((xint.-StatsBase.mean(xint)))));

        return RMS,Bias,Wcorr,Bss
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
