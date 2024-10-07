

using Dates, DelimitedFiles, Printf, NetCDF, Printf
#using Statistics

include("Moon.jl")

function weeksunday(t)
    dow=mod(Dates.dayofweek(t)+1,7);
    if dow==0
        dow=7;
    end
    return dow;
end

function plottitle(year,mon,site)
	gmt( "gmtset FONT_LABEL 10p FONT_ANNOT_PRIMARY 10p FONT_ANNOT_SECONDARY 10p MAP_FRAME_TYPE plain PS_MEDIA A3 PS_PAGE_ORIENTATION Portrait")
	
	gmtbegin(site*".ps")

	text!(J="X24.5c/1c",R="0/10/0/1",text_record([9.9 1 ], [site*" tide predictions"]), attrib=(font=(12,"Helvetica",:black),angle=0,justify=:RM),xshift=0,yshift=0,noclip=true,P=false,par=(:PS_MEDIA,:A3));

		
	basemap!(J="X24.5c/1c",R="2019-08-11T00:00:00/2019-08-18T00:00:00/0/1", frame=(axes=:N, annot="1K", ticks="1K"), xshift="3c", yshift="37c")
    
	text!(J="X24.5c/1c",R="0/10/0/1",text_record([9.9 3], [Dates.monthname(DateTime(year,mon,01))*" "*string(year)]), attrib=(font=(20,"Helvetica-Bold",:black),angle=0,justify=:RM),noclip=true);
    
	text!(J="X24.5c/1c",R="0/10/0/1",text_record([0.1 3 ], [site*" tide predictions"]), attrib=(font=(20,"Helvetica-Bold",:black),angle=0,justify=:LM),noclip=true);
end

function weekstart(year,mon)
	mmo=mon;
	t=[DateTime(year,mmo,01,0,0,0) (DateTime(year,mmo,01,0,0,0)+Dates.Month(1))];
    td=t[1]:Dates.Day(1):t[end];
	
	allsundays=findall(weeksunday.(td).==1);
    allsaturday=findall(weeksunday.(td).==7);
	ParamWeekSt=[t[1] t[1] t[1] t[1] t[1] t[1]]
    ParamWeekend=[td[allsaturday[1]] t[1] t[1] t[1] t[1] t[1]]
    nweek=1;

	for n=allsaturday[1]:(length(td)-1)
		if(weeksunday.(td[n]).==1)

			nweek=nweek+1;
			ParamWeekSt[nweek]=td[n];
		end
		if(weeksunday.(td[n]).==7)
			ParamWeekend[nweek]=td[n];
		end
	end
	if (ParamWeekend[nweek]==t[1])
		ParamWeekend[nweek]=td[end]-Dates.Day(1);
	end

	return nweek,ParamWeekSt,ParamWeekend
end

function paramshiflen(year,mon)

	mmo=mon;

	
	t=[DateTime(year,mmo,01,0,0,0) (DateTime(year,mmo,01,0,0,0)+Dates.Month(1))];
    td=t[1]:Dates.Day(1):t[end];

	
	
	Xshit_w1=(weeksunday(t[1])-1)*3.5; #in cm

	Xlen_w1=7*3.5-Xshit_w1; #in cm

	Xlen_w5=0;
	Xlen_w6=0;

	allsundays=findall(weeksunday.(td).==1);
	allsaturday=findall(weeksunday.(td).==7);

	if (Dates.isleapyear(year) && mmo==2 && weeksunday(t[1])==1)
		#No week 5
		Xlen_w5=0;
	else
		if((length(allsundays)==5) .& (length(allsaturday)==5))
			Xlen_w6=weeksunday(t[end]-Dates.Day(1))*3.5;
			Xlen_w5=7*3.5;
		else
			Xlen_w6=0;
			Xlen_w5=weeksunday(t[end]-Dates.Day(1))*3.5;
		end
	end

	ParamXShift=[Xshit_w1 -1.0*Xshit_w1 0.0 0.0 0.0 0.0] ;
	ParamYShift=[-1.0 -3 -3 -3 -3 32]
	ParamLen=[ Xlen_w1 7*3.5 7*3.5 7*3.5 Xlen_w5 Xlen_w6] ;

	return ParamXShift,ParamYShift,ParamLen
end

function tidepeaks(HT,Hlev,LowT,LowLev)
	tidepeaks=[Hlev; LowLev]
	tidetimes=[HT;LowT]

	p=sortperm([HT;LowT])
	tidepeaks=tidepeaks[p];
	tidetimes=tidetimes[p];

	# Xtide and Ytide are the location where the tide are going to be plotted
	Xtides=floor.(tidetimes,Dates.Day(1));
	Ytides=zeros(size(tidepeaks));

	m=0.8;
	Ytides[1]=m;
	dd=Xtides[1]
	for n=2:length(Ytides)
		if dd!==Xtides[n]
			dd=Xtides[n];
			m=0.8;
		else
			m=m-0.2;
		end
		Ytides[n]=m;
	end

	Xtidesdays=Dates.days.(Xtides.-floor.(tidetimes,Dates.Month(1))).+1

	tidepeaksR=round.(tidepeaks; digits=2);
	txtstring=Dates.format.(tidetimes,"HH:MM").*" ".*string.(tidepeaksR);

	return Xtidesdays,Ytides,txtstring
end
	

function GetMoon(year,month, UTCoffset)
	t=[DateTime(year,month,01,0,0,0) (DateTime(year,month,01,0,0,0)+Dates.Month(1))];
    #Warning constituents are relative to UTC time


    tt=t[1]:Dates.Minute(1):t[end];
    td=t[1]:Dates.Day(1):t[end];
	
	MPstday=first.(phasemoon.(td.-Dates.Hour(UTCoffset)));

	MPndday=first.(phasemoon.(td.-Dates.Hour(UTCoffset).+Dates.Day(1)));

	NewMoonDay=Dates.days.(td[(MPstday.>0.8) .& (MPndday.<0.1)].-td[1]).+1;
	FQmoonDay=Dates.days.(td[(MPstday.<=0.25) .& (MPndday.>0.25)].-td[1]).+1;
	FullmoonDay=Dates.days.(td[(MPstday.<=0.5) .& (MPndday.>0.5)].-td[1]).+1;
	LQmoonDay=Dates.days.(td[(MPstday.<=0.75) .& (MPndday.>0.75)].-td[1]).+1;

	

	return NewMoonDay,FQmoonDay,FullmoonDay,LQmoonDay
end

function Makesolarshade(lon,lat,year,mon,UTCoffset)

	

	t=[DateTime(year,mon,01,0,0,0) (DateTime(year,mon,01,0,0,0)+Dates.Month(1))];

	
    td=t[1]:Dates.Day(1):t[end];

	shadey=Vector{Int64}(undef,0);
	shadex=Vector{DateTime}(undef,0);

	push!(shadex,t[1],t[1]);
	push!(shadey,0,1);

	for i in td
		datestr=Dates.format.(i,"yyyy-mm-ddT01:00:00")
		
		sol=solar(sun=(pos=(lon,lat), date=datestr,  TZ=UTCoffset))
		sunrisestr=split(last(split(sol.text[6],"=")),":")
		sunsetstr=split(last(split(sol.text[7],"=")),":")

		sunrisehour=parse(Int64,first(sunrisestr));
		sunrisemin=parse(Int64,last(sunrisestr));

		sunsethour=parse(Int64,first(sunsetstr));
		sunsetmin=parse(Int64,last(sunsetstr));
		
		#println(sunrisestr," - ",sunsetstr)
		sunrise=i+Dates.Hour(sunrisehour)+Dates.Minute(sunrisemin)
		sunset=i+Dates.Hour(sunsethour)+Dates.Minute(sunsetmin)
		#println(sunrise,sunset)
		push!(shadey,1,0);
		push!(shadex,sunrise,sunrise)

		push!(shadey,0,1);
		push!(shadex,sunset,sunset)
	end

	push!(shadex,t[end],t[end]);
	push!(shadey,1,0);

	return shadex,shadey
end

function plotweek(weekst,weekend,Xshift,Yshift,width,Xtidesdays,Ytides,txtstring,ts,sl,solarx,solary,NewMoonDay,FQmoonDay,FullmoonDay,LQmoonDay;msl=0.0,tidemin=-2.0,tidemax=2.0,datumname="TGZ")

	pproj="X"*@sprintf("%f",width)*"c/2c"
	pprojsl="X"*@sprintf("%f",width)*"c/4c"

	
	if abs(msl-0.0)<=0.0001
		datumname="MSL"
	end
	ysllabel="m above "*datumname

	ndays=Dates.days.(weekend-weekst).+1;


	regionpl=Dates.format.(weekst,"yyyy-mm-dd")*"T00:00:00/"*Dates.format.(weekend,"yyyy-mm-dd")*"T23:59:59/0/1"

	td=weekst:Dates.Day(1):weekend
	
	regionpltxt=(Dates.dayofmonth(weekst),Dates.dayofmonth(weekend)+1,0,1)

	regionplsl=Dates.format.(weekst,"yyyy-mm-dd")*"T00:00:00/"*Dates.format.(weekend,"yyyy-mm-dd")*"T23:59:59/"*string(tidemin)*"/"*string(tidemax)

	
	
	## Plot the box with day of month and tide times
	basemap!(J=pproj,R=regionpl, frame=(axes=:wsen, annot="1d", grid="1d"),xshift=Xshift, yshift=Yshift,par=(:MAP_FRAME_TYPE,:inside))

	

	#text!(J=pproj,R=regionpltxt,["Hello"],x=1,y=0.5)
	# Plot the tide time and level
	text!(J=pproj,R=regionpltxt,text_record([Xtidesdays Ytides], txtstring), attrib=(font=(10,"Helvetica",:black),angle=0,justify=:LM),D="1.5c/0c");

	# Plot the day of Month
	text!(J=pproj,R=regionpltxt,text_record([Dates.dayofmonth.(td) ones(size(td))*0.75], string.(Dates.dayofmonth.(td))), attrib=(font=(20,"Helvetica-Bold",:black),angle=0,justify=:LM),D="0.2c/0c")
	
	# Plot the moon phase
	
	#Full
	scatter!(J=pproj, R=regionpltxt, FullmoonDay, ones(size(FullmoonDay)).*0.25, markersize = 0.5, marker=:circle, fill=:white,ml=(1,:black),D="0.5c/0c")
	
	#New
	scatter!(J=pproj, R=regionpltxt, NewMoonDay, ones(size(NewMoonDay)).*0.25, markersize = 0.5, marker=:circle, fill=:black,lw="1p",D="0.5c/0c")
	
	#First Quarter
	scatter!(J=pproj, R=regionpltxt, FQmoonDay, ones(size(FQmoonDay)).*0.25, markersize = 0.5, marker=:circle, fill=:white,ml=(1,:black),D="0.5c/0c")
	
	plot!(J=pproj, R=regionpltxt, [FQmoonDay.*1.0 ones(size(FQmoonDay)).*0.25  ones(size(FQmoonDay)).*-90.0 ones(size(FQmoonDay)).*90.0], markersize = 0.5, marker=:wedge, fill=:black,ml=(1,:black),D="0.5c/0c")
	
	
	#Last Quarter
	scatter!(J=pproj, R=regionpltxt, LQmoonDay, ones(size(LQmoonDay)).*0.25, markersize = 0.5, marker=:circle, fill=:white,ml=(1,:black),D="0.5c/0c")
	
	plot!(J=pproj, R=regionpltxt, [LQmoonDay.*1.0 ones(size(LQmoonDay)).*0.25  ones(size(LQmoonDay)).*90.0 ones(size(LQmoonDay)).*-90.0], markersize = 0.5, marker=:wedge, fill=:black,ml=(1,:black),D="0.5c/0c")
	
	
	## Plot the tide signal+solar
	basemap!(J=pprojsl,R=regionplsl,xaxis=(annot="6h",grid="1d"), yaxis=(annot=:auto, ticks=:auto, label=ysllabel),yshift="-4c")

	plot!(J=pprojsl,R=regionplsl,ts,sl.*0.0.+msl,pen=(0.5, :dashed))

	plot!(J=pprojsl,R=regionpl,solarx,solary,fill=:100,alpha=25)
	
	plot!(J=pprojsl,R=regionplsl,[ts[1]; ts; ts[end]],[tidemin; sl; tidemin],fill=:200,alpha=25)
	plot!(J=pprojsl,R=regionplsl,ts,sl,lw="1p",lc=:50)
	
	#println(size(sl))

end

function MakeCalendarData(year::Int64,month::Int64,coef::coef_s,opt::opt_s; UTCoffset=0,datumoffset=0.0)
    #opt=opt_s(datafolder,false,2.0,0.0,false,false,false,false,false);
    t=[DateTime(year,month,01,0,0,0) (DateTime(year,month,01,0,0,0)+Dates.Month(1))];
    #Warning constituents are relative to UTC time


    tt=t[1]:Dates.Minute(1):t[end];
    td=t[1]:Dates.Day(1):t[end];
    # remove utc offset for predicting tide time for each given month 
    HT,Hlev,LowT,LowLev=Predicttidetime(t.-Dates.Hour(UTCoffset),coef,opt);
    sl=ut_reconstr1(tt.-Dates.Hour(UTCoffset),coef,opt).+datumoffset;
    # Add UTC offset again 
    HT=HT.+Dates.Hour(UTCoffset);
    LowT=LowT.+Dates.Hour(UTCoffset);

	Hlev=Hlev.+datumoffset;
	LowLev=LowLev.+datumoffset;
	
	return tt,sl,HT,Hlev,LowT,LowLev
end


function MakeCalendar(CST,opt,year,sitename,lon,lat;UTCoffset=0,datumoffset=0.0)
	
	
	site=sitename;
	

	tidemax=sum(CST.A).+datumoffset;

	tidemin=-1.0*sum(CST.A).+datumoffset;


	#HTexcurve,LTexcurve=exceedencecurve(CST,opt;nyear=1);

	for month=1:12
	
		plottitle(year,month,site);
	
		ts,sl,HT,Hlev,LowT,LowLev=MakeCalendarData(year,month,CST,opt; UTCoffset=UTCoffset,datumoffset=datumoffset)
	
		nweek,WeekSt,Weekend=weekstart(year,month);
	
		ParamXShift,ParamYShift,ParamLen=paramshiflen(year,month)
	
		Xtidesdays,Ytides,txtstring=tidepeaks(HT,Hlev,LowT,LowLev)

		solarx,solary=Makesolarshade(lon,lat,year,month,UTCoffset)

		NewMoonDay,FQmoonDay,FullmoonDay,LQmoonDay=GetMoon(year,month, UTCoffset)
	
		#basemap!(J="X"*@sprintf("%f",Xlen_w1)*"c/1c",R="2019-01-01T00:00:00/2019-01-5T23:59:59/0/1", frame=(axes=:wsen, annot="1d", grid="1d"),xshift=Xshit_w1, yshift="0c",par=(:MAP_FRAME_TYPE,:inside),show=true)
	
		#text_record([Xtidesdays Ytides], [txtstring])
		for iw=1:nweek
			plotweek(WeekSt[iw],Weekend[iw],ParamXShift[iw],ParamYShift[iw],ParamLen[iw],Xtidesdays,Ytides,txtstring,ts,sl,solarx,solary,NewMoonDay,FQmoonDay,FullmoonDay,LQmoonDay,msl=datumoffset,tidemin=tidemin,tidemax=tidemax)
		end
		#plotweek(WeekSt[2],Weekend[2],ParamXShift[2],ParamYShift[2],ParamLen[2],Xtidesdays,Ytides,txtstring)
	
		
		gmtend(show=false)
		psconvert(site*".ps", format=:f, P=[], out_name=site*"_"*string(year)*"_"*string(month))

	
	
	end
	
	open("listcombine.tmp","w") do io
		for i=1:12
			Printf.@printf(io,"%s\n",site*"_"*string(year)*"_"*string(i)*".pdf")
		end
	end
	psconvert(list_file="listcombine.tmp",  format=:F , out_name=site*"_"*string(year)*"_Tide_Calendar")
end
