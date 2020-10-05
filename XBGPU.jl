module XBGPU

    import Printf
    export writemd,writedep

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
                    for ii=(floor(nx/11)*11+1):(floor(nx/11)*11+(nx/11-floor(nx/11))*11)
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
