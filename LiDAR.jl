

"""
Module to read and process LiDAR data
"""
module LiDAR

    import FileIO, LasIO, junkinthetrunk, LazIO

    export Readlas, Rasterise

    """
    Readlas(file::String)
    Function to read .las and .laz files
    """
    function Readlas(file::String)
        fileext=rsplit(file,".";limit=2);
        if (lowercase(fileext[end])=="laz")
            header, points = LazIO.load(file)
        else
            header, points = FileIO.load(file)
        end

        return header, points
    end

    function Rasterise(file::String, res::AbstractFloat; region=(NaN,NaN,NaN,NaN)::NTuple{4,AbstractFloat})

        header, points = Readlas(file);

        xmin = isnan(region[1]) ? floor(header.x_min / res) * res : region[1];
        ymin = isnan(region[3]) ? floor(header.y_min / res) * res : region[3];

        xmax = isnan(region[2]) ? ceil(header.x_max / res) * res : region[2];
        ymax = isnan(region[4]) ? ceil(header.y_max / res) * res : region[4];

        region=(xmin,xmax,ymin,ymax);
        regfilnal=Checkregion(region,res);

        xc,yc,zc=blockmean(points,header,regfilnal,res)

        #write2nc(xc,yc,zcheck,file*".nc");
        return xc,yc,zc
    end

    function Rasterise(filelst::Vector{String}, res::AbstractFloat; region=(NaN,NaN,NaN,NaN)::NTuple{4,AbstractFloat})

        header, points = Readlas(filelst[1]);

        xmin = isnan(region[1]) ? floor(header.x_min / res) *res : region[1];
        ymin = isnan(region[3]) ? floor(header.y_min / res) *res : region[3];

        xmax = isnan(region[2]) ? ceil(header.x_max / res) * res : region[2];
        ymax = isnan(region[4]) ? ceil(header.y_max / res) * res : region[4];

        regf=(xmin,xmax,ymin,ymax);
        region=Checkregion(region,res);

        nx,ny=Calcnxny(region,res);


        xx=fill(NaN,(nx));
        yy=fill(NaN,(ny));

        for j = 1:ny
            yy[j]=region[3]+(j-1)*res;
        end
        for i=1:nx
            xx[i]=region[1]+(i-1)*res;
        end

        z=fill(0.0,(nx,ny));
        n=fill(0,(nx,ny));

        for il=1:length(filelst)
            if il>1
                header, points = Readlas(filelst[il]);
            end

            nn,zz=blockcount(points,header,region,res);

            n +=  nn

            z += zz

        end

        indxNaN=n.==0;
        z[.!indxNaN] = z[.!indxNaN] ./ n[.!indxNaN];
        z[indxNaN] .= NaN

        return xx,yy,z
    end

    function Checkregion(region::NTuple{4,AbstractFloat},res::AbstractFloat)
        xmin=region[1];
        ymin=region[3];
        xmax=xmin + ceil((region[2] - xmin) / res) * res;
        ymax=ymin + ceil((region[4] - ymin) / res) * res;

        regfixed=(xmin,xmax,ymin,ymax)
        return regfixed
    end


    function Calcnxny(region::NTuple{4,AbstractFloat},res::AbstractFloat)

        xmin=region[1];
        ymin=region[3];
        xmax=region[2];
        ymax=region[4];

        nx=Int((xmax - xmin) / res) + 1
        ny=Int((ymax - ymin) / res) + 1

        return nx,ny
    end



    function blockij(p::LasIO.LasPoint,h::LasIO.LasHeader,region::NTuple{4,AbstractFloat},res::AbstractFloat)

        i=Int(0)
        j=Int(0)


        xp=LasIO.xcoord(p,h);
        yp=LasIO.ycoord(p,h);
        #zp=zcoord(p,h);
        if (xp>=(region[1] - 0.5 * res)) &&  (xp<=(region[2] + 0.5 * res)) && (yp>=(region[3] - 0.5 * res)) &&  (yp<=(region[4] + 0.5 * res))
            i=Int(ceil((xp - (region[1] - 0.5 * res)) / res))
            j=Int(ceil((yp - (region[3] - 0.5 * res)) / res))
        else
            i=NaN
            j=NaN
        end

        return i,j
    end

    function blockcount(p::Vector{<:LasIO.LasPoint},h::LasIO.LasHeader,region::NTuple{4,AbstractFloat},res::AbstractFloat)

        nx,ny=Calcnxny(region,res);

        zz=fill(0.0,(nx,ny));
        nn=fill(0,(nx,ny));


        for ip=1:length(p)
            i,j=blockij(p[ip],h,region,res);
            if !isnan(i) || !isnan(j)
                zz[Int(i),Int(j)]+=LasIO.zcoord(p[ip],h);
                nn[Int(i),Int(j)]+=1;
            end
        end

        return nn,zz

    end

    function blockmean(p::Vector{<:LasIO.LasPoint},h::LasIO.LasHeader,region::NTuple{4,AbstractFloat},res::AbstractFloat)

        nx,ny=Calcnxny(region,res);


        xx=fill(NaN,(nx));
        yy=fill(NaN,(ny));

        for j = 1:ny
            yy[j]=region[3]+(j-1)*res;
        end
        for i=1:nx
            xx[i]=region[1]+(i-1)*res;
        end
        nn,zz=blockcount(p,h,region,res);

        indxNaN=nn.==0;
        zz[.!indxNaN] = zz[.!indxNaN] ./ nn[.!indxNaN];
        zz[indxNaN] .= NaN

        return xx,yy,zz

    end




end
