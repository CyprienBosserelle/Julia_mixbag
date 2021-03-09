"""
    junkinthetrunk module
    Collection of random usefull functions

    Available functions:
    nanmean, nanmedian, write2nc, pol2cartCompass, cart2polCompass, interp1, nearneighb1, bilinearinterpUG, skills, sinwave, Gaussian

"""
module junkinthetrunk

    import StatsBase, NetCDF, Dates
    export nanmean, nanmedian, write2nc, pol2cartCompass, cart2polCompass, interp1, nearneighb1, bilinearinterpUG, skills, sinwave, Gaussian


    #nanmean does a mean while ignoring nans
	"""
	    average of a vector ignoring NaNs
	    usage nanmean(x)
	"""
    nanmean(x) = StatsBase.mean(filter(!isnan,x))
    nanmean(x,y) = mapslices(nanmean,x,y)

	function nanmedian(x)
		y=filter(!isnan,x);
		if !isempty(y)
			return StatsBase.median(filter(!isnan,x))
		else
			return NaN
		end
	end
	"""
	    produce a sinusoidal wave for
	    usage:
		Waterlevel=sinwave(timevector, amplitude, period, phasedalayinsec, datumshift)
		period and phase delay are in sec
		timevector is expected to be a vector of doubles
	"""
	function sinwave(timevector, amplitude, period, phasedalayinsec, datumshift)

		ttt=timevector
		A=amplitude;

		T=period;#12.8*3600.0

		B=2.0*pi/(T);
		C=phasedalayinsec; # time delay
		D=datumshift; # Datum shift


		return A.*sin.(B.*(ttt.+C)).+D;
	end
	function sinwave(timevector, period)
		return sinwave(timevector, 1.0, period, 0.0, 0.0);
	end
	function sinwave(timevector, amplitude, period)
		return sinwave(timevector, amplitude, period, 0.0, 0.0);
	end
	"""
	    produce a gaussian wave for
	    usage:
		Waterlevel=Gaussian(timevector,amplitude,position,width)
		timevector is expected to be a vector of doubles
	"""
	function Gaussian(x,a,b,c)
	    u=x-b;
	    up=-1.0*u*u/(2*c*c);
	    return a*exp.(up);
	end

    


    #cart2pol
    function pol2cartCompass(theta,speed)
         ee=speed*cosd(90.0 .- theta);
         nn=speed.*sind(90.0 .- theta);
         return ee,nn;
     end

	 function cart2polCompass(ee,nn)
		 speed=hypot.(ee,nn);
		 theta=90.0 .- rad2deg.(atan.(nn,ee));;
          return speed,theta;
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

		xA=vec(xA);
	    yA=vec(yA);


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
    function write2nc(x,y,z,ncfile,varnames)

        xatts = Dict("longname" => "eastings",
          "units"    => "m")
        yatts = Dict("longname" => "northings",
                  "units"    => "m")
        varatts = Dict("longname" => "Topo / Bathy above MVD53",
                  "units"    => "m")
        NetCDF.nccreate(ncfile,varnames[3],varnames[1],x,xatts,varnames[2],y,yatts,atts=varatts)
        NetCDF.ncwrite(x,ncfile,varnames[1]);
        NetCDF.ncwrite(y,ncfile,varnames[2]);
        NetCDF.ncwrite(z,ncfile,varnames[3]);
        NetCDF.ncclose(ncfile);
    end
	function write2nc(x,y,z,ncfile)
		varnames=["x","y","z"];
		write2nc(x,y,z,ncfile,varnames)
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

		PmO=sum(abs.(xint.-xmea));
		OmOb=sum(abs.(xmea.-StatsBase.mean(xmea)))
		Cwcorr=2.0;
		Wcorr=0.0;

		if(PmO<=(Cwcorr*OmOb))
			Wcorr=1.0-(PmO/(Cwcorr*OmOb));
		else
			Wcorr=((Cwcorr*OmOb)/PmO)-1.0;
		end

        RMS=sqrt(StatsBase.mean((xint.-xmea).^2));
        Bias=StatsBase.mean(xint)-StatsBase.mean(xmea);
        #Wcorr=1-(sum((xmea-xint).^2)/(sum((abs.(xmea.-StatsBase.mean(xint)).+abs.(xint.-StatsBase.mean(xint))).^2)));

		Bss=1-(StatsBase.var((xmea-xint))/(StatsBase.var((xint.-StatsBase.mean(xint)))));

        return RMS,Bias,Wcorr,Bss
    end

	"""
        Detrend a serie
		(Safe-ish with NaN)
		usage:
		y=detrend(x) or y=detrend(x,bp) where bp is idex of breakpoints

    """
	function detrend(s,bp)

	    x=copy(s)

	    bp=1

	    # make sure you ignore the NANs
	    mmean=nanmean(x);
	    if any(isnan.(x))
	        x[isnan.(x)].=mmean;
	    end

	    N = size(x,1);
	    bp = unique([1; bp; N]);   # Include both endpoints
	    bp = bp[bp .>= 1 .& bp .<=N];   # Should error in the future
	    lbp = length(bp);
	    #Build regressor with linear pieces + DC
	    a = zeros(N,lbp);
	    a[1:N,1] = collect(1:N)./N;
	    for k = 2:(lbp-1)
	        M = N - bp[k];
	        global a[(bp[k]+1):end,k] = (1:M)./M;
	    end
	    a[1:N,end] .= 1;
	    y = s - a*(a\x);

	    return y
	end
	function detrend(x)
	    return detrend(x,1);
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
