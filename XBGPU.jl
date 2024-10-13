module XBGPU

    import Printf
    export writemd,writedep,curegrid

    """
    Function to cure the bathy boundary to make it suitable for XBGPU reqirements:
            1 . uniform offshore boundary (in for 5 cells wide) with smooth transition to real bathymetry (buffersize cell wide)
            2 . no bathymetry gradient on the top and bottom boundary (10 cell wide)
            3 . Wall on the land boundary (2 cells wide)
    """
    function curegrid(z;buffersize=20,offshorezs=NaN)

        if isnan(offshorezs)
            offshorezb=maximum(z[1,:]);
        end

        inshorezb=minimum(z[end,:]);
        nx,ny=size(z);
        znew = copy(z);

        # Set edges for the side of bnd

        #First the lower side
        for n=1:9
            znew[:,n]=z[:,10];
        end
        # Now the upper side
        for n=ny:-1:(ny-9)
            znew[:,n]=z[:,ny-10];
        end

        z=copy(znew);

        # Uniform offshore boundary
        znew[1:5,:] .= offshorezb;

        # smooth transition to real bathy

        for n=1:buffersize
            znew[5+n,:] = ((offshorezb*(buffersize-n)) .+ (z[5+buffersize,:].*n))./(buffersize*1.0);
        end



        # Also add a trumpian wall on the 2 cell wide side of the land boundary
        znew[end,:] .= inshorezb;
        znew[end-1,:] .= inshorezb;

        return znew
    end

    """
    """
    function readmd(filename::String)
        #read header
        
        #read md bathymety file
        #[depth nx ny dx dy sst]=readmd(filename)

        open("/usr/share/dict/words") do io
            readline(io)
            






            readlines(io)
        end
            
        #fid=fopen(filename);
            s=fgets(fid);
            A=textscan(s,'%f %f %f %f %s');
            nx=A{1};
            ny=A{2};
            dx=A{3};
            dy=A{4};
            sst=A{5};
            depth=zeros(nx,ny);
            for yy=1:ny
                s=fgets(fid);
                j=str2num(s);
                B=fscanf(fid,'%f',nx);
                depth(:,j)=B;
                s=fgets(fid);
            end
            fclose(fid);
            
            pcolor(depth'); shading flat; axis image
            
            % depth(1:27,216:240)=-9;

        return x,y,z
    end


    """
    Function to write md files ((suitable to use in XBeach_GPU))
    """
    function writemd(z,dx,nx,ny,grdalpha,filename)
        open(filename,"w") do io
            Printf.@printf(io,"%d\t%d\t%f\t%f\t%f\n",nx,ny,dx,dx,grdalpha);
            for jj=ny:-1:1
                Printf.@printf(io,"%d\n",jj);
                for ii=1:floor(nx/11)
                    Printf.@printf(io,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",z[Int64((ii-1)*11+1),jj],z[Int64((ii-1)*11+2),jj],z[Int64((ii-1)*11+3),jj],z[Int64((ii-1)*11+4),jj],z[Int64((ii-1)*11+5),jj],z[Int64((ii-1)*11+6),jj],z[Int64((ii-1)*11+7),jj],z[Int64((ii-1)*11+8),jj],z[Int64((ii-1)*11+9),jj],z[Int64((ii-1)*11+10),jj],z[Int64((ii-1)*11+11),jj]);
                end


                if mod(nx,11)!=0
                    for ii=convert(Int64, round((floor(nx/11)*11+1), digits=0)):convert(Int64, round((floor(nx/11)*11+(nx/11-floor(nx/11))*11),digits=0))
                        Printf.@printf(io,"%f\t",z[Int64(ii),jj])
                    end
                    Printf.@printf(io,"\n");
                end

            end

        end
    end
    

    """
    Function to write .dep files (suitable to use in XBeach)
    """
    function writedep(z,dx,filename)
        nx,ny=size(z)
        open(filename,"w") do io
            for jj=1:ny
                for ii=1:ceil(nx/12)
                    if ii*12<=nx
                        for n=1:12
                            ind=Int64((ii-1)*12+n);
                            Printf.@printf(io,"%f\t",z[ind,jj]);
                        end
                    else
                        for n=1:(nx-(ii-1)*12)
                            ind=Int64((ii-1)*12+n);
                            Printf.@printf(io,"%f\t",z[ind,jj]);
                        end

                    end
                    Printf.@printf(io,"\n");
                end
            end
        end
    end



end
