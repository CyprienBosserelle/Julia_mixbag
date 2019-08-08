
# first build a constituent object
struct cnstit
    #constituent parameter
    Name::String
    A::Float64
    A_ci::Float64
    g::Float64
    g_ci::Float64
    freq::Float64
    lind::Int64
end

struct option
    twodim::Bool
    minsnr::Float64
    minpe::Float64
    nodsatlint::Bool# 0
    nodsatnone::Bool# 0
    gwchlint::Bool# 0
    gwchnone::Bool# 0
end


#Input should be an array of cnstit
coef=Vector{cnstit}(undef,5)

opt=option(false,2.0,0.0,false,false)


coef[1]=cnstit("M2",0.493138384657320,0.000588412284822346,206.472841113229,0.0651156148667687,0.0805114006717706,48)
coef[2]=cnstit("N2",0.122220407530187,0.000449550463417388,189.592580251570,0.248441538866563,0.0789992486986775,42)
coef[3]=cnstit("S2",0.0625455828349263,0.000578787340178285,228.577462497562,0.485431942020961,0.0833333333333333,57)
coef[4]=cnstit("K1",0.0624069393945612,0.00102801034135234,226.373348002011,0.871333797973204,0.0417807462216577,21)
coef[5]=cnstit("O1",0.0319278720784309,0.000901721376152978,234.835289824733,1.62594515196797,0.0387306544501129,13)




# single tide reconstruction based on UT
# UT_RECONSTR1()
# Reconstruction for a single record. See comments for UT_RECONSTR().
# UTide v1p0 9/2011 d.codiga@gso.uri.edu

#fprintf('ut_reconstr: ');

# parse inputs and options
#[t,opt] = ut_rcninit(tin,varargin);

#Preallocate a few vars




# determine constituents to include
if !isempty([])
    [~,ind] = ismember(cellstr(opt.cnstit),coef.name);
    if ~isequal(length(ind),length(cellstr(opt.cnstit)))
        error(['ut_reconstr: one or more of input constituents Cnstit '...
            'not found in coef.name']);
    end
else
    SNR=zeros(length(coef));
    PE=zeros(length(coef));
    #ind = collect(1:length(allcon));
    if opt.twodim
        SNR = (coef.Lsmaj.^2 +coef.Lsmin.^2)./((coef.Lsmaj_ci/1.96).^2 + (coef.Lsmin_ci/1.96).^2);
        PE = sum(coef.Lsmaj.^2 + coef.Lsmin.^2);
        PE = 100*(coef.Lsmaj.^2 + coef.Lsmin.^2)/PE;
    else
        for el=1:length(coef)
            SNR[el] = (coef[el].A.^2)./((coef[el].A_ci/1.96).^2);
            PE[el] = 100*coef[el].A.^2/sum(coef[el].A.^2);
        end

    end
    ind = (SNR.>=opt.minsnr) .& (PE.>=opt.minpe);
end



#############################
# complex coefficients
rpd = pi/180;
if opt.twodim
    ap = 0.5*(coef.Lsmaj(ind) + coef.Lsmin(ind)) .* exp(1i*(coef.theta(ind) - coef.g(ind))*rpd);
    am = 0.5*(coef.Lsmaj(ind) - coef.Lsmin(ind)) .* exp(1i*(coef.theta(ind) + coef.g(ind))*rpd);
else
    ap=zeros(ComplexF64,length(coef))
    for el=1:length(coef)
        ap[el] = 0.5*coef[el].A.*exp(-1im * coef[el].g.*rpd);
    end
        am = conj.(ap);

end




# exponentials
ngflgs = [opt.nodsatlint opt.nodsatnone opt.gwchlint opt.gwchnone];
#fprintf('prep/calcs ... ');
E = ut_E(t,coef.aux.reftime,coef.aux.frq(ind),coef.aux.lind(ind),coef.aux.lat,ngflgs,coef.aux.opt.prefilt);

# fit
fit = E*ap + conj(E)*am;

# mean (& trend)
u = nan*ones(size(tin));
whr = ~isnan(tin);
if coef.aux.opt.twodim
    v = u;
    if coef.aux.opt.notrend
        u(whr) = real(fit) + coef.umean;
        v(whr) = imag(fit) + coef.vmean;
    else
        u(whr) = real(fit) + coef.umean + ...
            coef.uslope*(t-coef.aux.reftime);
        v(whr) = imag(fit) + coef.vmean + ...
            coef.vslope*(t-coef.aux.reftime);
    end
else
    if coef.aux.opt.notrend
        u(whr) = real(fit) + coef.mean;
    else
        u(whr) = real(fit) + coef.mean + ...
            coef.slope*(t-coef.aux.reftime);
    end
    v = [];
end
fprintf('done.\n');
