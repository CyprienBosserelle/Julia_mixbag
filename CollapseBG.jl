# Collapse resolution from adaptive run into a single level


using NetCDF


function ftoi(x)
	if x>=0.0
		intcast=trunc(Int, x+0.5)
	else
		intcast=trunc(Int, x-0.5)
	end
	return intcast
end

function Injection(i,j, Arr)
	
	ic=ftoi(ceil(i*0.5))
	jc=ftoi(ceil(j*0.5))

	return Arr[ic,jc]

end

function CellAvg(i,j,Arr)
	iclb=i*2-1
	jclb=j*2-1

	icrb=i*2
	jcrb=jclb

	iclt=iclb
	jclt=j*2

	icrt=icrb
	jcrt=jclt

	Vc=0.25*(Arr[iclb,jclb]+Arr[iclt,jclt]+Arr[icrb,jcrb]+Arr[icrt,jcrt])
	
	return Vc

end

function VolumeAvg(i,j,Arr,harr)
	iclb=i*2-1
	jclb=j*2-1

	icrb=i*2
	jcrb=jclb

	iclt=iclb
	jclt=j*2

	icrt=icrb
	jcrt=jclt

	havg=0.25*(harr[iclb,jclb]+harr[iclt,jclt]+harr[icrb,jcrb]+harr[icrt,jcrt])

	Vavg=(0.25/havg)*(Arr[iclb,jclb]*harr[iclb,jclb]+Arr[iclt,jclt]*harr[iclt,jclt]+Arr[icrb,jcrb]*harr[icrb,jcrb]+Arr[icrt,jcrt]*harr[icrt,jcrt])
	
	return Vavg
end

function Getlevelstr(var::String,level::Signed)
	if level<0
        stlv="N";
    else
        stlv="P";
    end
    fvar=var*"_"*stlv*string(abs(level))

	return fvar
end

function readBGFinfo(file::String)

	

	
	levmax=NetCDF.ncgetatt(file, "Global", "maxlevel");
    levmin=NetCDF.ncgetatt(file, "Global", "minlevel");

    xmin=NetCDF.ncgetatt(file, "Global", "xmin") 
    xmax=NetCDF.ncgetatt(file, "Global", "xmax") 
    ymin=NetCDF.ncgetatt(file, "Global", "ymin") 
    ymax=NetCDF.ncgetatt(file, "Global", "ymax")

	region=(xmin,xmax,ymin,ymax)

	levels=levmin:levmax

	xxstr=Getlevelstr("xx",levmin)

	xx_lmin=NetCDF.ncread(file,xxstr)

	dxmin=xx_lmin[2]-xx_lmin[1];
		
	
	return region,levels,dxmin
end



function readBGvar(file::String,var::String,level::Signed,step::Signed;istart=1,jstart=1,ni=-1,nj=-1)
	
	fvar=Getlevelstr(var::String,level::Signed)
	
	scalefac=NetCDF.ncgetatt(file, fvar, "scale_factor");
    addoffset=NetCDF.ncgetatt(file, fvar, "add_offset");
    missval=NetCDF.ncgetatt(file, fvar, "_FillValue");
	
	zz=NetCDF.ncread(file,fvar,[istart,jstart,max(step,1)], [ni,nj,1]);

	
    if any([!isnothing(scalefac) !isnothing(addoffset)])
        zz=zz.*scalefac.+addoffset;
    end
    zz[zz.==missval].=NaN;

	return dropdims(zz,dims=3)
end

function Cellavg(zref,zfine)
	nx,ny=size(zref);

	for i=1:nx
		for j=1:ny
			if(isnan(zref[i,j]))
				zref[i,j]=CellAvg(i,j,zfine)
			end
		end
	end
	return zref
end

function InjectVar(zref,zcoarse)
	nx,ny=size(zref);

	for i=1:nx
		for j=1:ny
			if(isnan(zref[i,j]))
				zref[i,j]=Injection(i,j, zcoarse)
			end
		end
	end

	return zref
end

function denanArr(zref,znonan)
	nx,ny=size(zref);
	nxnn,nynn=size(znonan);
	if (nxnn==nx) & (nynn==ny)
		for i=1:nx
			for j=1:ny
				if(isnan(zref[i,j]))
					zref[i,j]=znonan[i,j]
				end
			end
		end
	end
	return zref
end

function InjectVar(file::String,var::String,step::Signed,levref::Signed)
	zref=readBGvar(file,var,levref,step);

	zcoarse=readBGvar(file,var,levref-1,step);

	zrefined=InjectVar(zref,zcoarse)

	return zrefined
end


function Getindexregion(file::String,reflevel::Signed,region::NTuple{4, Float64})
	
	
	xxstr=Getlevelstr("xx",reflevel)
	yystr=Getlevelstr("yy",reflevel)

	xx_rlev=NetCDF.ncread(file,xxstr)
	yy_rlev=NetCDF.ncread(file,yystr)


	nx=length(xx_rlev)
	ny=length(yy_rlev)

	println(reflevel)

	

	#dxrlev==dyrlev
	dxrlev=xx_rlev[2]-xx_rlev[1]

	println(dxrlev)
	
	irefst=ftoi(min(max(floor((region[1]-xx_rlev[1])/dxrlev),1),nx))
	irefend=ftoi(min(max(ceil((region[2]-xx_rlev[1])/dxrlev),1),nx))
	jrefst=ftoi(min(max(floor((region[3]-yy_rlev[1])/dxrlev),1),ny))
	jrefend=ftoi(min(max(ceil((region[4]-yy_rlev[1])/dxrlev),1),ny))

	regij = (irefst,irefend,jrefst,jrefend)

	newBB=(xx_rlev[irefst] - dxrlev*0.5, xx_rlev[irefend] + dxrlev*0.5, yy_rlev[jrefst] - dxrlev*0.5, yy_rlev[jrefend] + dxrlev *0.5) 

	return regij,newBB
	
	
	
end

function Collapseindex(refij,reflevel,level)
	ist=refij[1]
	jst=refij[3]
	ind=refij[2]
	jnd=refij[4]

	
	if level>=reflevel
		for levi=(reflevel+1):level

			ist=ist*2-1
			ind=ind*2
	
			jst=jst*2-1
			jnd=jnd*2
		end
	else
		for levi=reflevel:-1:level
			ist=(ist+1)*0.5
			ind=ind*0.5;
			jst=(jst+1)*0.5
			jnd=jnd*0.5
		end
	end

	return (ist,ind,jst,jnd)
end

function Getxy(file::String,level::Signed,region::NTuple{4, Float64})
	xxstr=Getlevelstr("xx",level)
	yystr=Getlevelstr("yy",level)

	xx_lev=NetCDF.ncread(file,xxstr)
	yy_lev=NetCDF.ncread(file,yystr)

	

	nx=length(xx_lev)
	ny=length(yy_lev)

	indx=(xx_lev .> region[1]) .& (xx_lev .< region[2])
	indy=(yy_lev .> region[3]) .& (yy_lev .< region[4])
	x=xx_lev[indx]
	y=yy_lev[indy]


	return x,y
end

function CollapseBG(file::String,var::String;step=-1,level=999999,region=(NaN,NaN,NaN,NaN))

	regionwhole,levels,dxmin=readBGFinfo(file)


	xmin= isnan(region[1]) ? regionwhole[1] : region[1];
	xmax= isnan(region[2]) ? regionwhole[2] : region[2];
	ymin= isnan(region[3]) ? regionwhole[3] : region[3];
	ymax= isnan(region[4]) ? regionwhole[4] : region[4];

	region=(xmin,xmax,ymin,ymax)

	println(region)

	outputlevel=max(min(levels[end],level),levels[1]);

	
	refij,refBB=Getindexregion(file,levels[1],region);

	ist=refij[1]
	jst=refij[3]
	ind=refij[2]
	jnd=refij[4]

	#println(refij)

	############################
	## create the x/y axis
	
	xref,yref=Getxy(file,outputlevel,refBB)

	time=NetCDF.ncread(file,"time")
	
	if step<0
		step=length(time)
	else
		step=max(min(step,length(time)),1);
	end

	
	##################
	## DO the Injection

	
	#zref=readBGvar(file,var,outputlevel,step,);
	zcoarse=readBGvar(file,var,levels[1],step,istart=ist,jstart=jst,ni=(ind-ist+1),nj=(jnd-jst+1));
	for levi=levels[2]:outputlevel

		ist=ist*2-1
		ind=ind*2

		jst=jst*2-1
		jnd=jnd*2
		
		zfine=readBGvar(file,var,levi,step,istart=ist,jstart=jst,ni=(ind-ist+1),nj=(jnd-jst+1));
		
		zcoarse=InjectVar(zfine,zcoarse)
	end

	zref=zcoarse;

	if outputlevel<levels[end]



	
	levbb=Collapseindex(refij,levels[1],levels[end]);

	println(levbb)
	zcoarseref=readBGvar(file,var,levels[end],step,istart=levbb[1],jstart=levbb[3],ni=(levbb[2]-levbb[1]+1),nj=(levbb[4]-levbb[3]+1));

	##########################
	## Do the cell avg
	for levi=levels[end-1]:-1:(outputlevel)


		levbb=Collapseindex(refij,levels[1],levi);

		println(levbb)
		zcoarse=readBGvar(file,var,levi,step,istart=levbb[1],jstart=levbb[3],ni=(levbb[2]-levbb[1]+1),nj=(levbb[4]-levbb[3]+1));

		
		

		
		levbb=Collapseindex(refij,levels[1],levi+1);

		zfine=readBGvar(file,var,levi+1,step,istart=levbb[1],jstart=levbb[3],ni=(levbb[2]-levbb[1]+1),nj=(levbb[4]-levbb[3]+1));

		

		zcoarseref=Cellavg(zcoarse,denanArr(zcoarseref,zfine));

		

		
	end

	zref=denanArr(zref,zcoarseref)
	

	end
		

	

	return xref,yref,zref
end


function GetFlooding(file::String;step=-1,level=999999, eps=0.03, region=(NaN,NaN,NaN,NaN))
	x,y,zsmax=CollapseBG(file,"zsmax",region=region,step=step,level=level);
	x,y,hmax=CollapseBG(file,"hmax",region=region,step=step,level=level);
	
	dry = hmax .< eps;

	zsmax[dry] .= NaN;
	hmax[dry] .= NaN;

	return x,y,zsmax,hmax
end

function WriteFlooding(file::String,outroot::String;region=(NaN,NaN,NaN,NaN),step=-1,level=999999, eps=0.03)
	x,y,zsmax,hmax=GetFlooding(file,region=region,step=step,level=level,eps=eps);
	write2nc(x,y,zsmax,"zsmax",outroot*"_zsmax.nc")
	write2nc(x,y,hmax,"hmax",outroot*"_hmax.nc");
end


function write2nc(x,y,z,varname::String,ncfile::String)
    #ny,nx=size(z)
    xdimid=NcDim("x",length(x),values=x[1:end])
	ydimid=NcDim("y",length(y),values=y[1:end])
	
	
	varid=NcVar(varname,[xdimid,ydimid],t=Float32,compress=3)
	
	ncid=NetCDF.create(ncfile,[varid], mode=NC_NETCDF4)
    NetCDF.putvar(ncid,varname,Float32.(z),start=[1,1],count=[-1,-1])
		
		
	NetCDF.sync(ncid)
end


# Usage:
# WriteFlooding(BGfloodfile,outputrootname;region=(NaN,NaN,NaN,NaN),step=-1,level=999999, eps=0.03)