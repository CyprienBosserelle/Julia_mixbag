module Delft3D

    import Printf, DelimitedFiles
    export readdep,writedep

    # Function to write md files (next gen?)
    function writedep(z,filename)
        open(filename,"w") do io
            nx,ny=size(z);
            for jj=1:1:ny
                #Printf.@printf(io,"%d\n",jj);
                for ii=1:floor(nx/12)
                    Printf.@printf(io,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",z[Int64((ii-1)*12+1),jj],z[Int64((ii-1)*12+2),jj],z[Int64((ii-1)*12+3),jj],z[Int64((ii-1)*12+4),jj],z[Int64((ii-1)*11+5),jj],z[Int64((ii-1)*12+6),jj],z[Int64((ii-1)*12+7),jj],z[Int64((ii-1)*12+8),jj],z[Int64((ii-1)*12+9),jj],z[Int64((ii-1)*12+10),jj],z[Int64((ii-1)*12+11),jj],z[Int64((ii-1)*12+12),jj]);
                end


                if mod(nx,12)!=0
                    for ii=(floor(nx/12)*12+1):(floor(nx/12)*12+(nx/12-floor(nx/12))*12)
                        Printf.@printf(io,"%f\t",z[Int64(ii),jj])
                    end
                    Printf.@printf(io,"\n");
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






end
