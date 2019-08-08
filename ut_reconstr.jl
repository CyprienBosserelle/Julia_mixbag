
using Dates

function  ut_astron(jd)
    # % UT_ASTRON()
    # % calculate astronomical constants
    # % input
    # %   jd = time [datenum UTC] (1 x nt)
    # % outputs
    # %   astro = matrix [tau s h p np pp]T, units are [cycles] (6 x nt)
    # %   ader = matrix of derivatives of astro [cycles/day] (6 x nt)
    # % UTide v1p0 9/2011 d.codiga@gso.uri.edu
    # % (copy of t_astron.m from t_tide, Pawlowicz et al 2002)
    nt=length(jd);
    jd=collect(DateTime(2009,11,30,12,0,0):Dates.Hour(1):DateTime(2009,12,31,12,0,0));

    d=datetime2julian.(jd).-datetime2julian(DateTime(1899,12,31,12,0,0));
    D=d./10000;
    args=[ones(nt),
          d,
          D.*D,
          D.^3]
    sc= [ 270.434164,13.1763965268,-0.0000850, 0.000000039];
    hc= [ 279.696678, 0.9856473354, 0.00002267,0.000000000];
    pc= [ 334.329556, 0.1114040803,-0.0007739,-0.00000026];
    npc=[-259.183275, 0.0529539222,-0.0001557,-0.000000050];
    ppc=[ 281.220844, 0.0000470684, 0.0000339, 0.000000070];

    astro=zeros(6,nt)

    #astro=rem( [sc,hc,pc,npc,ppc].*args./360.0 ,1);
    astro[2,:]=rem.( sc*args'./360.0 ,1);
    astro[3,:]=rem( hc.*args./360.0 ,1);
    astro[4,:]=rem( pc.*args./360.0 ,1);
    astro[5,:]=rem( npc.*args./360.0 ,1);
    astro[6,:]=rem( ppc.*args./360.0 ,1);



    tau=rem(jd(:)',1)+astro[3,:]-astro[2,:];
    astro[1,:]=tau;
    dargs=[zeros(size(jd));
           ones(size(jd));
           2.0e-4.*D;
           3.0e-4.*D.*D];
    ader=[sc;hc;pc;npc;ppc]*dargs./360.0;
    dtau=1.0+ader(2,:)-ader(1,:);
    ader=[dtau;ader];

    return astro,ader
end

function [F,U,V] = ut_FUV(t,tref,lind,lat,ngflgs)
    # % UT_FUV()
    # % compute nodal/satellite correction factors and astronomical argument
    # % inputs
    # %   t = times [datenum UTC] (nt x 1)
    # %   tref = reference time [datenum UTC] (1 x 1)
    # %   lind = list indices of constituents in ut_constants.mat (nc x 1)
    # %   lat = latitude [deg N] (1 x 1)
    # %   ngflgs = [NodsatLint NodsatNone GwchLint GwchNone] each 0/1
    # % output
    # %   F = real nodsat correction to amplitude [unitless] (nt x nc)
    # %   U = nodsat correction to phase [cycles] (nt x nc)
    # %   V = astronomical argument [cycles] (nt x nc)
    # % UTide v1p0 9/2011 d.codiga@gso.uri.edu
    # % (uses parts of t_vuf.m from t_tide, Pawlowicz et al 2002)

    nt = length(t);
    nc = length(lind);
    # nodsat
    if ngflgs(2) # none
        F = ones(nt,nc);
        U = zeros(nt,nc);
    else
        if ngflgs(1) # linearized times
            #tt = tref;
        else         # exact times
            tt = t;
        end
        ntt = length(tt);
        load('ut_constants.mat');
        [astro,~]=ut_astron(tt');
        if abs(lat)<5
            lat=sign(lat).*5;
        end
        slat=sin(pi*lat/180);
        rr=sat.amprat;
        j=find(sat.ilatfac==1);
        rr(j)=rr(j).*0.36309.*(1.0-5.0.*slat.*slat)./slat;
        j=find(sat.ilatfac==2);
        rr(j)=rr(j).*2.59808.*slat;
        uu=rem( sat.deldood*astro(4:6,:)+sat.phcorr(:,ones(1,ntt)), 1);
        nfreq=length(const.isat); %#ok
        mat = rr(:,ones(1,ntt)).*exp(1i*2*pi*uu);
        F = ones(nfreq,ntt);
        ind = unique(sat.iconst);
        for i = 1:length(ind)
            F(ind(i),:) = 1+sum(mat(sat.iconst==ind(i),:),1);
        end
        U = imag(log(F))/(2*pi); % faster than angle(F)
        F=abs(F);
        for k=find(isfinite(const.ishallow))'
            ik=const.ishallow(k)+(0:const.nshallow(k)-1);
            j = shallow.iname(ik);
            exp1 = shallow.coef(ik);
            exp2 = abs(exp1);
            F(k,:)=prod(F(j,:).^exp2(:,ones(ntt,1)),1);
            U(k,:)=sum(U(j,:).*exp1(:,ones(ntt,1)),1);
        end
        F=F(lind,:)';
        U=U(lind,:)';
        if ngflgs(1) % nodal/satellite with linearized times
            F = F(ones(nt,1),:);
            U = U(ones(nt,1),:);
        end
    end
    # gwch (astron arg)
    if ngflgs(4) % none (raw phase lags not greenwich phase lags)
        if ~exist('const','var')
            load('ut_constants.mat','const');
        end
        [~,ader] = ut_astron(tref);
        ii=isfinite(const.ishallow);
        const.freq(~ii) = (const.doodson(~ii,:)*ader)/(24);
        for k=find(ii)'
            ik=const.ishallow(k)+(0:const.nshallow(k)-1);
            const.freq(k)=sum(const.freq(shallow.iname(ik)).*shallow.coef(ik));
        end
        V = 24*(t-tref)*const.freq(lind)';
    else
        if ngflgs(3)  # linearized times
            tt = tref;
        else
            tt = t;   # exact times
        end
        ntt = length(tt);
        if exist('astro','var')
            if ~isequal(size(astro,2),ntt)
                [astro,~]=ut_astron(tt');
            end
        else
            [astro,~]=ut_astron(tt');
        end
        if ~exist('const','var')
            load('ut_constants.mat');
        end
        V=rem( const.doodson*astro+const.semi(:,ones(1,ntt)), 1);
        for k=find(isfinite(const.ishallow))'
            ik=const.ishallow(k)+(0:const.nshallow(k)-1);
            j = shallow.iname(ik);
            exp1 = shallow.coef(ik);
            V(k,:) = sum(V(j,:).*exp1(:,ones(ntt,1)),1);
        end
        V=V(lind,:)';
        if ngflgs(3)    % linearized times
            [~,ader] = ut_astron(tref);
            ii=isfinite(const.ishallow);
            const.freq(~ii) = (const.doodson(~ii,:)*ader)/(24);
            for k=find(ii)'
                ik=const.ishallow(k)+(0:const.nshallow(k)-1);
                const.freq(k)=sum( const.freq(shallow.iname(ik)).* ...
                    shallow.coef(ik) );
            end
            V = V(ones(1,nt),:) + 24*(t-tref)*const.freq(lind)';
        end
    end
end

function ut_E(t,tref,frq,lind,lat,ngflgs,prefilt)
    # UT_E()
    # compute complex exponential basis function
    # inputs
    #   t = times [datenum UTC] (nt x 1)
    #   tref = reference time [datenum UTC] (1 x 1)
    #   frq = frequencies [cph] (nc x 1)
    #   lind = list indices of constituents in ut_constants.mat (nc x 1)
    #   lat = latitude [deg N] (1 x 1)
    #   ngflgs = [NodsatLint NodsatNone GwchLint GwchNone] each 0/1
    #       ([0 1 0 1] case not allowed, and not needed, in ut_E)
    #   prefilt = 'prefilt' input to ut_solv
    # output
    #   E = complex exponential basis function [unitless] (nt x nc)
    # UTide v1p0 9/2011 d.codiga@gso.uri.edu

    nt = length(t);
    nc = length(lind);
    if ngflgs(2) && ngflgs(4)
        F = ones(nt,nc);
        U = zeros(nt,nc);
        V = 24*(t-tref)*frq';
    else
        [F,U,V] = ut_FUV(t,tref,lind,lat,ngflgs);
    end
    E = F.*exp(1im*(U+V)*2*pi);
    if ~isempty(prefilt)
        # P=interp1(prefilt.frq,prefilt.P,frq)';
        # P( P>max(prefilt.rng) | P<min(prefilt.rng) | isnan(P) )=1;
        # E = E.*P(ones(nt,1),:);
    end
    return E
end


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

opt=option(false,2.0,0.0,false,false,false,false)


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
    # [~,ind] = ismember(cellstr(opt.cnstit),coef.name);
    # if ~isequal(length(ind),length(cellstr(opt.cnstit)))
    #     error(['ut_reconstr: one or more of input constituents Cnstit '...
    #         'not found in coef.name']);
    # end
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
    # ap = 0.5*(coef.Lsmaj(ind) + coef.Lsmin(ind)) .* exp(1i*(coef.theta(ind) - coef.g(ind))*rpd);
    # am = 0.5*(coef.Lsmaj(ind) - coef.Lsmin(ind)) .* exp(1i*(coef.theta(ind) + coef.g(ind))*rpd);
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
