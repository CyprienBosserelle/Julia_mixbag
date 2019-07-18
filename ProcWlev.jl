module ProcWlev
# This modelue is a collection of functiion used to process Water Level data

    import Dates, SQLite, DataFrames
    export pressure2depth, getRBRDepthdata



    function Calcdepth(seapressure, latitude)
        # Calculate depth from sea pressure.
        # Calculated using
        # the Saunders & Fofonoff method.
        # originaly from matlab rsktools
        # Author: RBR Ltd. Ottawa ON, Canada
        #Adapted to Julia by Cyprien Bosserelle NIWA
        x = (sin(latitude/57.29578))^2;
        gr = 9.780318*(1.0 + (5.2788e-3 + 2.36e-5*x).*x) + 1.092e-6.*seapressure;
        depth = (((-1.82e-15*seapressure + 2.279e-10).*seapressure - 2.2512e-5).*seapressure + 9.72659).*seapressure;
        depth = depth./gr;
        return depth;
    end
    function pressure2depth(seapressure)=Calcdepth(seapressure, latitude)
    function pressure2depth(seapressure,lat)=Calcdepth(seapressure, lat)

    function rsktime2datetime(time)
        #Convert from ruskin (unix) time to datetime
        #
        return Dates.unix2datetime(time./1000);
    end


    function getRBRDepthdata(file, latitiude)
        #file=Folder*fff[6]
        db = SQLite.DB(file);
        SQLite.columns(db,"parameterKeys")
        #keyvalue=DataFrame(SQLite.Query(db, "SELECT value FROM parameterKeys"));
        #key=DataFrame(SQLite.Query(db, "SELECT key FROM parameterKeys"));

        meanatmpress=10.132501;


        ttime= DataFrame(SQLite.Query(db, "SELECT tstamp FROM data"));

        Pressure= DataFrame(SQLite.Query(db, "SELECT channel02 FROM data"));
        #Temp= DataFrame(SQLite.Query(db, "SELECT channel05 FROM data"));

        truetime=rsktime2datetime.(collect(skipmissing(ttime[:,1])));

        Depth=Calcdepth.(collect(skipmissing(Pressure[:,1])).-meanatmpress,latitiude);




        return [truetime, Depth]
    end

    function writespec2nc(time,freq,Pspec,outfilename)
        xatts = Dict("longname" => "time",
                  "units"    => "ms")
        yatts = Dict("longname" => "freq",
                  "units"    => "hz")
        varatts = Dict("longname" => "Power Spectrum",
                  "units"    => "m2/hz")
        nccreate(outfilename,"PSpec","time",collect(1:length(time)),xatts,"freq",freq,yatts,atts=varatts)
        ncwrite(Dates.value.(time),outfilename,"time");
        ncwrite(freq,outfilename,"freq");
        ncwrite(Pspec,outfilename,"PSpec");
        ncclose(outfilename);
        
    end

end
