module Delft3D

    import Printf, DelimitedFiles, NetCDF
    export readdep,writedep,Get_curv_grid_XZ,Save_curv_grid_XZ,writecurv2gmt

    # Function to write md files (next gen?)
    function writedep(z,filename)
        open(filename,"w") do io
            nx,ny=size(z);
            for jj=1:ny
                #Printf.@printf(io,"%d\n",jj);
                # for ii=1:floor(nx/12)
                #     Printf.@printf(io,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",z[Int64((ii-1)*12+1),jj],z[Int64((ii-1)*12+2),jj],z[Int64((ii-1)*12+3),jj],z[Int64((ii-1)*12+4),jj],z[Int64((ii-1)*11+5),jj],z[Int64((ii-1)*12+6),jj],z[Int64((ii-1)*12+7),jj],z[Int64((ii-1)*12+8),jj],z[Int64((ii-1)*12+9),jj],z[Int64((ii-1)*12+10),jj],z[Int64((ii-1)*12+11),jj],z[Int64((ii-1)*12+12),jj]);
                # end
                #
                #
                # if mod(nx,12)!=0
                #     for ii=(floor(nx/12)*12+1):(floor(nx/12)*12+(nx/12-floor(nx/12))*12)
                #         Printf.@printf(io,"%f\t",z[Int64(ii),jj])
                #     end
                #     Printf.@printf(io,"\n");
                # end
                for ii=1:nx
                    Printf.@printf(io,"%g",z[ii,jj])
                    if mod(ii,12)==0 || ii==nx
                        Printf.@printf(io,"\n")
                    else
                        Printf.@printf(io,"\t")
                    end
                end

            end

        end
    end


    function readdep(depfile,nx,ny)
        data=DelimitedFiles.readdlm(depfile);


        data[data.==""].=-999.0;

        A=convert(Array{Float64,2},data)

        #heatmap(A)

        B=zeros(nx,ny);

        cm=ceil(nx/12)*12;

        for i=1:nx
            for j=1:ny
                n=Int((j-1)*cm/12+ceil(i/12))
                m=Int(i-(ceil(i/12)-1)*12)
                B[i,j]=A[n,m];
            end
        end
        return B;
    end
    function writecurv2gmt(X,Y,Z,mask,outname)
        # Write curvilinear grid to gmt for plotting
        si=size(X);
        nx=si[1];
        ny=si[2];

        open(outname,"w") do io
            for i=1:(nx-1)
                for j=1:(ny-1)
                    if (mask[i,j]==1 && mask[i+1,j]==1 && mask[i,j+1]==1 && mask[i+1,j+1]==1)
                        Printf.@printf(io,">> -Z%f\n",Z[i,j]);
                        Printf.@printf(io,"%f\t%f\n",X[i,j],Y[i,j]);
                        Printf.@printf(io,"%f\t%f\n",X[i+1,j],Y[i+1,j]);
                        Printf.@printf(io,"%f\t%f\n",X[i+1,j+1],Y[i+1,j+1]);
                        Printf.@printf(io,"%f\t%f\n",X[i,j+1],Y[i,j+1]);
                    end
                end
            end

        end
    end

    function Get_curv_grid_XZ(ncfile,var,step,lev)
        #println("This code assumes you have a bathymetry file called trim-2008_TideCorr.nc in your working directory")
        #ncfile="trim-2008_TideCorr.nc" #"D:\\Projects\\Port_Otago\\2008\\D3D\\trim-2008_TideCorr.nc"
        #println("The output file will be called D3dmeshtoplot.gmt in your working directory")
        #outname="D3dmeshtoplot.gmt" #"D:\\Projects\\Port_Otago\\Bathy\\D3dmeshtoplot.gmt"

        #var="DPS0"
        nc=NetCDF.open(ncfile);

        sz=size(nc[var]);

        if length(sz)==2
            #Get bathy data
            Z=NetCDF.ncread(ncfile, var, start=[1,1], count = [-1,-1]);
        elseif length(sz)==3
            Z=NetCDF.ncread(ncfile, var, start=[1,1,step], count = [-1,-1,1]);
        elseif length(sz)==4
            Z=NetCDF.ncread(ncfile, var, start=[1,1,lev,step], count = [-1,-1,1,1]);
        end

        #Get X and Y data fo that Variable
        X=NetCDF.ncread(ncfile,"XZ", start=[1,1], count = [-1,-1]);
        Y=NetCDF.ncread(ncfile,"YZ", start=[1,1], count = [-1,-1]);

        #Get mask info
        mask=NetCDF.ncread(ncfile,"KCS", start=[1,1], count = [-1,-1]);

        return X,Y,Z,mask
    end

    function Save_curv_grid_XZ(ncfile,var,outname,step,lev)
        #println("This code assumes you have a bathymetry file called trim-2008_TideCorr.nc in your working directory")
        #ncfile="trim-2008_TideCorr.nc" #"D:\\Projects\\Port_Otago\\2008\\D3D\\trim-2008_TideCorr.nc"
        #println("The output file will be called D3dmeshtoplot.gmt in your working directory")
        #outname="D3dmeshtoplot.gmt" #"D:\\Projects\\Port_Otago\\Bathy\\D3dmeshtoplot.gmt"

        #var="DPS0"
        X,Y,Z,mask=Get_curv_grid_XZ(ncfile,var,step,lev)

        writecurv2gmt(X,Y,Z,mask,outname)
        # si=size(X);
        # nx=si[1];
        # ny=si[2];
        #
        # open(outname,"w") do io
        #     for i=1:(nx-1)
        #         for j=1:(ny-1)
        #             if (mask[i,j]==1 && mask[i+1,j]==1 && mask[i,j+1]==1 && mask[i+1,j+1]==1)
        #                 Printf.@printf(io,">> -Z%f\n",-1.0.*Z[i,j]);
        #                 Printf.@printf(io,"%f\t%f\n",X[i,j],Y[i,j]);
        #                 Printf.@printf(io,"%f\t%f\n",X[i+1,j],Y[i+1,j]);
        #                 Printf.@printf(io,"%f\t%f\n",X[i+1,j+1],Y[i+1,j+1]);
        #                 Printf.@printf(io,"%f\t%f\n",X[i,j+1],Y[i,j+1]);
        #             end
        #         end
        #     end
        #
        # end
    end

    function Save_curv_grid_XZ(ncfile,var,outname)
        Save_curv_grid_XZ(ncfile,var,outname,1,1);
    end
    function Save_curv_grid_XZ(ncfile,var,outname,step)
        Save_curv_grid_XZ(ncfile,var,outname,step,1);
    end
    function Get_curv_grid_XZ(ncfile,var)
        Get_curv_grid_XZ(ncfile,var,1,1);
    end
    function Get_curv_grid_XZ(ncfile,var,step)
        Get_curv_grid_XZ(ncfile,var,step,1);
    end









end
