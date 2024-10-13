

# phasemoon(DateTime(2019,8,30,11))
# phasemoon(Dates.now())

# moon_age_location(Dates.datetime2julian(Dates.now()))
include("MakeTideCalendar.jl")


function get_constituentsFES(lon,lat,FESdatafolder,opt)
	
    println("Check Constituents")
	#constituents=ncread(tidemodel,"parameters",start=[1,nodeidx],count=[-1,1]);
    # read available constituents
	cstfiles=readdir(FESdatafolder);

    cstnames=first.(split.(cstfiles,".nc"));
    cstnamesU=uppercase.(cstnames);

    cstfreq=zeros(length(cstnames));
    cstlind=zeros(length(cstnames));
    # Build database of freq and lind for those constituents
    Namedict=["M2" "N2" "S2" "K1" "O1" "NU2" "MU2" "P1" "2N2" "K2" "L2" "Q1" "M3" "MKS2" "EPS2" "LDA2" "J1" "SK3" "MS4" "NO1" "M4" "OQ2" "2MS6" "OO1" "M6" "THE1" "SO1" "PHI1" "UPS1" "SIG1" "TAU1" "SO3" "BET1" "MN4" "2SK5" "2Q1" "MK3" "RHO1" "2MN6" "MK4" "ETA2" "MSN2" "2MK6" "S4" "2SM6" "MO3" "SSA" "ALP1" "SK4" "SN4" "2MK5" "MSF" "MM" "MSM" "CHI1" "MSK6" "3MK7" "MF" "M8"]
    freqdict=vec([0.0805114006717706 0.0789992486986775 0.0833333333333333 0.0417807462216577 0.0387306544501129 0.0792016199833009 0.0776894680102079 0.0415525871116757 0.0774870967255845 0.0835614924433154 0.0820235526448637 0.0372185024770199 0.120767101007656 0.0807395597817526 0.0761773160371148 0.0818211813602403 0.0432928981947507 0.125114079554991 0.163844734005104 0.0402685942485646 0.161022801343541 0.0759749447524915 0.244356134676875 0.0448308379932025 0.241534202015312 0.0430905269101274 0.0446026788832204 0.0420089053316397 0.0463429899662955 0.0359087217885502 0.0389588135600949 0.122063987783446 0.0400404351385826 0.159510649370448 0.208447412888324 0.0357063505039268 0.122292146893428 0.0374208737616432 0.240022050042219 0.164072893115086 0.0850736444164084 0.0848454853064264 0.244584293786857 0.166666666666667 0.247178067338437 0.119242055121884 0.000228159109982028 0.0343965698154571 0.166894825776649 0.162332582032011 0.202803547565199 0.00282193266156274 0.00151215197309305 0.00130978068846969 0.0404709655331880 0.247406226448419 0.283314948236969 0.00305009177154477 0.322045602687082]);
    #
    linddict=vec([48 42 57 21 13 43 40 19 39 59 54 11 69 50 35 53 25 74 86 16 82 34 110 28 106 24 27 23 29 10 14 71 15 79 99 9 72 12 103 87 61 60 111 89 113 68 3 8 90 84 96 6 5 4 17 114 120 7 125]);

    Namefreqdict=Dict{String,Float64}(Namedict .=> freqdict');
    Namelinddict=Dict{String,Int64}(Namedict .=> linddict');

    for i=1:length(cstnames)
        cstfreq[i]=get(Namefreqdict,cstnamesU[i],NaN);
        cstlind[i]=get(NNamelinddict,cstnamesU[i],0);
    end

    emptyval=zeros(length(cstfreq))

    cstA=zeros(length(cstfreq))
    cstg=zeros(length(cstfreq))


    for i=1:length(cstnames)
        ncfile=FESdatafolder*cstnames[i]*".nc";

        latmodel=ncread(ncfile,"lat");
        lonmodel=ncread(ncfile,"lon");

        nodelon=argmin(abs.(lonmodel.-lon))
        nodelat=argmin(abs.(latmodel.-lat))

        cstA[i]=ncread(ncfile,"amplitude",start=[nodelon,nodelat],count=[1,1])*0.01;#cm to m
        cstg[i]=ncread(ncfile,"phase",start=[nodelon,nodelat],count=[1,1]);
        
    end
    





	
	tempcoef=coef_s(cstnamesU,cstA,emptyval,cstg,emptyval,cstfreq,emptyval,0.0,0.0,lat,DateTime(2000,1,1))

	linds=makelind(tempcoef,opt);

	nzmodelcoef=coef_s(cstnamesU,cstA,emptyval,cstg,emptyval,cstfreq,linds,0.0,0.0,lat,DateTime(2000,1,1))

	return nzmodelcoef
end

"""
get_constituentsNZ(lon,lat,tidemodel,opt)

read original nz tide constituents
"""
function get_constituentsNZ(lon,lat,tidemodel,opt)
	println("Load lat lon")
	latlonmodel=ncread(tidemodel,"parameters",start=[2,1],count=[2,-1]);
	nodeidx=argmin(hypot.(latlonmodel[1,:].-lat,latlonmodel[2,:].-lon))'
	println("Load constituents")
	constituents=ncread(tidemodel,"parameters",start=[1,nodeidx],count=[-1,1]);

	cstnames=["M2", "S2", "N2", "K1", "O1", "Q1", "L2", "P1", "MU2", "T2", "K2", "NU2", "2N2"];
	cstfreq=[0.08051140067176844, 0.08333333333333333, 0.07899924869868231, 0.04178074622165745, 0.038730654450111, 0.03721850247702487, 0.08202355264485457, 0.04155258711167587, 0.07768946801020357, 0.08321925922953312, 0.0835614924433149, 0.0792016199832897, 0.0774870967255962];
	cstA=constituents[5:17];
	cstg=constituents[18:(end-1)];
	emptyval=zeros(length(cstfreq))
	tempcoef=coef_s(cstnames,cstA,emptyval,cstg,emptyval,cstfreq,emptyval,0.0,0.0,constituents[2],DateTime(2000,1,1))

	linds=makelind(tempcoef,opt);

	nzmodelcoef=coef_s(cstnames,cstA,emptyval,cstg,emptyval,cstfreq,linds,0.0,0.0,constituents[2],DateTime(2000,1,1))

	return nzmodelcoef
end

function GetFESCalendar(lon,lat,year,site; optdatafolder = "D:\\Projects\\Tonga\\tonga-ocean-forecasting-tools\\Tide-Calendars\\", FESdatafolder = "D:\\Data\\Wlevel\\FES2014\\fes2014b_elevations_extrapolated\\ocean_tide_extrapolated\\",UTCoffset=0,datumoffset=0.0)
	opt=opt_s(optdatafolder,false,2.0,0.0,false,false,false,false,false);

	

	CST=get_constituentsFES(lon,lat,FESdatafolder,opt)

	MakeCalendar(CST,opt,year,site,lon,lat,UTCoffset=UTCoffset,datumoffset=datumoffset)
	
end




