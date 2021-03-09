"""
    BG module
    Collection of function to create input and process output of the BG model

    Available functions:
    SetEdges, patchgrids

    #Example:
    using BG

"""
module BG

    import Printf, NetCDF, junkinthetrunk
    export CalcRunup
    """
    ttt
    """
    function CalcRunup(file)
        ncfile=NetCDF.open(file);


        fvar="zsmax"




        scalefac=NetCDF.ncgetatt(file, fvar, "scale_factor");
        addoffset=NetCDF.ncgetatt(file, fvar, "add_offset");
        missval=NetCDF.ncgetatt(file, fvar, "_FillValue");



        zz=ncread(file,fvar,[1,1,max(stp,1)], [-1,-1,1]);
        if any([!isnothing(scalefac) !isnothing(addoffset)])
            zz=zz.*scalefac.+addoffset;
        end
        zz[zz.==missval].=NaN;

        zz=zz[.!isnan.(zz)];
        #
    end


end
