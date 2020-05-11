###########################
## Okada 1985

using NetCDF
using Plots
File="Bathy_model.nc";


EE=ncread(File,"x");
NN=ncread(File,"y");


E=EE*ones(size(NN))'
N=collect(transpose(NN*ones(size(EE))'))

xf=389529;#[389529 346996 346996];
yf=8284999#[8284999 8229431 8229431];
ddepth=18.0#[18.0 18.0 18.0];



slipTI=10.0#[10.00 4.62 4.71];
ddip=65.0#[65 29 29];
# rrake=[265 90 90];
# sstrike=[324 175 180];
rrake=-91#[-91 90 90];
sstrike=144#[144 175 180];
xinc=130#[130 50 50];
yinc=30#[30 75 75];

ef=(E.-xf)./1000;#spheriq_dist(E,lat,lon,lat).*facef/1000;
nf=(N.-yf)./1000;#spheriq_dist(lon,N,lon,lat).*facnf/1000;%-5000:10:5000;


#[uE,uN,uZ] = okada85(ef,nf,Depth,strike,dip,lengthg,width,rake,slip,open)
uE,uN,uZ=Okada85(ef,nf,ddepth,sstrike,ddip,xinc,yinc,rrake,slipTI,0)

heatmap(uZ)



#	Translated to Julia By Cyprien Bosserelle 2020
#	References:
#	   Aki K., and P. G. Richards, Quantitative seismology, Freemann & Co,
#	      New York, 1980.
#	   Okada Y., Surface deformation due to shear and tensile faults in a
#	      half-space, Bull. Seismol. Soc. Am., 75:4, 1135-1154, 1985.
#
#	Acknowledgments: Dmitry Nicolsky, Qian Yao, Halldor Geirsson
#	Development history:
#	   [2014-05-24]: fixes a bug for tilt calculation (K1) when DIP=90.
#	      Detected by Halldor Geirsson.
#	   [2012-11-08]: solves partially mathematical singularities in
#	      specific cases like DIP=90, STRIKE=0, and fault reaching surface.
#	      Detected by Qian Yao.
#	   [2012-08-29]: allows vectorization of RAKE, SLIP and OPEN.
#	   [2011-03-08]: help review.
#	   [2011-03-06]: new optional argument to plot fault geometry with
#	      output arguments, and bug correction for the fault centroid position
#	      (in calculation and plot).
#	   [2010-11-29]: change coordinates and depth to fault centroid
#	      (instead of middle top edge).
#	   [2010-09-24]: bugs correction in the syntax of I1, K2 and uyy_tf
#	      functions, affecting some components. Detected by Dmitry Nicolsky.
#
#	Copyright (c) 1997-2012, François Beauducel, covered by BSD License.
#	All rights reserved.
#
#	Redistribution and use in source and binary forms, with or without
#	modification, are permitted provided that the following conditions are
#	met:
#
#	   * Redistributions of source code must retain the above copyright
#	     notice, this list of conditions and the following disclaimer.
#	   * Redistributions in binary form must reproduce the above copyright
#	     notice, this list of conditions and the following disclaimer in
#	     the documentation and/or other materials provided with the distribution
#
#	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#	AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#	IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
#	ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
#	LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
#	CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
#	SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
#	INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
#	CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
#	ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#	POSSIBILITY OF SUCH DAMAGE.

"""
OKADA85 Surface deformation due to a finite rectangular source.
	[uE,uN,uZ,uZE,uZN,uNN,uNE,uEN,uEE] = OKADA85(...
	   E,N,DEPTH,STRIKE,DIP,LENGTH,WIDTH,RAKE,SLIP,OPEN)
	computes displacements, tilts and strains at the surface of an elastic
	half-space, due to a dislocation defined by RAKE, SLIP, and OPEN on a
	rectangular fault defined by orientation STRIKE and DIP, and size LENGTH and
	WIDTH. The fault centroid is located (0,0,-DEPTH).

	   E,N    : coordinates of observation points in a geographic referential
	            (East,North,Up) relative to fault centroid (units are described below)
	   DEPTH  : depth of the fault centroid (DEPTH > 0)
	   STRIKE : fault trace direction (0 to 360° relative to North), defined so
	            that the fault dips to the right side of the trace
	   DIP    : angle between the fault and a horizontal plane (0 to 90°)
	   LENGTH : fault length in the STRIKE direction (LENGTH > 0)
	   WIDTH  : fault width in the DIP direction (WIDTH > 0)
	   RAKE   : direction the hanging wall moves during rupture, measured relative
	            to the fault STRIKE (-180 to 180°).
	   SLIP   : dislocation in RAKE direction (length unit)
	   OPEN   : dislocation in tensile component (same unit as SLIP)

	returns the following variables (same matrix size as E and N):
	   uN,uE,uZ        : displacements (unit of SLIP and OPEN)
	   uZE,uZN         : tilts (in rad * FACTOR)
	   uNN,uNE,uEN,uEE : horizontal strains POSITIVE = COMPRESSION (unit of FACTOR)

	Length unit consistency: E, N, DEPTH, LENGTH, and WIDTH must have the same
	unit (e.g. km) which can be different from that of SLIP and OPEN (e.g. m) but
	with a possible FACTOR on tilt and strain results (in this case, an
	amplification of km/m = 1000). To have FACTOR = 1 (tilt in radians and
	correct strain unit), use the same length unit for all aforesaid variables.

	[...] = OKADA85(...,NU) specifies Poisson's ratio NU (default is 0.25 for
	an isotropic medium).

	Formulas and notations from Okada [1985] solution excepted for strain
	convention (here positive strain means compression), and for the fault
	parameters after Aki & Richards [1980], e.g.:
	      DIP=90, RAKE=0   : left lateral (senestral) strike slip
	      DIP=90, RAKE=180 : right lateral (dextral) strike slip
	      DIP=70, RAKE=90  : reverse fault
	      DIP=70, RAKE=-90 : normal fault

	It is also possible to produce partial outputs, with following syntax:
	   [uE,uN,uZ] = OKADA85(...) for displacements only;
	   [uE,uN,uZ,uZE,uZN] = OKADA85(...) for displacements and tilts;
	   [uE,uN,uZ,uNN,uNE,uEN,uEE] = OKADA85(...) for displacements and strains;
	   [uZE,uZN] = OKADA85(...) for tilts only;
	   [uZE,uZN,uNN,uNE,uEN,uEE] = OKADA85(...) for tilts and strains;
	   [uNN,uNE,uEN,uEE] = OKADA85(...) for strains only.

	Note that vertical strain components can be obtained with following equations:
	   uNZ = -uZN;
	   uEZ = -uZE;
	   uZZ = -(uEE + uNN)*NU/(1-NU);

	[...] = OKADA85(...,'plot') or OKADA85(...) without output argument
	produces a 3-D figure with fault geometry and dislocation at scale (if
	all of the fault parameters are scalar).

	Equations are all vectorized excepted for argument DIP which must be
	a scalar (beacause of a singularity in Okada's equations); all other
	arguments can be scalar or matrix of the same size.

	Example:

	   [E,N] = meshgrid(linspace(-10,10,50));
	   [uE,uN,uZ] = okada85(E,N,2,30,70,5,3,-45,1,1,'plot');
	   figure, surf(E,N,uN)

	considers a 5x3 fault at depth 2, N30°-strike, 70°-dip, and unit dislocation
	in all directions (reverse, senestral and open). Displacements are computed
	on a regular grid from -10 to 10, and North displacements are plotted as a
	surface.


	Author: François Beauducel <beauducel@ipgp.fr>
	   Institut de Physique du Globe de Paris
	Created: 1997
	Updated: 2014-05-24
	"""
function Okada85(e,n,depth,strike,dip,L,W,rake,slip,U3)

	strike = strike*pi/180;	# converting STRIKE in radian
	dip = dip*pi/180;	# converting DIP in radian ('delta' in Okada's equations)
	rake = rake*pi/180;	#converting RAKE in radian
	nu = 0.25;

	# Defines dislocation in the fault plane system
	U1 = cos(rake).*slip;
	U2 = sin(rake).*slip;
	# Converts fault coordinates (E,N,DEPTH) relative to centroid
	# into Okada's reference system (X,Y,D)
	d = depth + sin(dip).*W/2;	# d is fault's top edge
	ec = e .+ cos(strike).*cos(dip).*W/2;
	nc = n .- sin(strike).*cos(dip).*W/2;
	x = cos(strike).*nc .+ sin(strike).*ec .+ L/2;
	y = sin(strike).*nc .- cos(strike).*ec .+ cos(dip).*W;
	# Variable substitution (independent from xi and eta)
	p = y.*cos(dip) .+ d.*sin(dip);
	q = y.*sin(dip) .- d.*cos(dip);

	# Displacements
	ux = ( -U1/(2*pi) .* chinnery_ux_ss.(x,p,L,W,q,dip,nu) # strike-slip
		.- U2/(2*pi) .* chinnery_ux_ds.(x,p,L,W,q,dip,nu) # dip-slip
		.+ U3/(2*pi) .* chinnery_ux_tf.(x,p,L,W,q,dip,nu)); # tensile fault

	uy = ( -U1/(2*pi) .* chinnery_uy_ss.(x,p,L,W,q,dip,nu) # strike-slip
		.- U2/(2*pi) .* chinnery_uy_ds.(x,p,L,W,q,dip,nu) # dip-slip
		.+ U3/(2*pi) .* chinnery_uy_tf.(x,p,L,W,q,dip,nu)); # tensile fault

	uz = ( -U1/(2*pi) .* chinnery_uz_ss.(x,p,L,W,q,dip,nu) # strike-slip
		.- U2/(2*pi) .* chinnery_uz_ds.(x,p,L,W,q,dip,nu) # dip-slip
		.+ U3/(2*pi) .* chinnery_uz_tf.(x,p,L,W,q,dip,nu)); # tensile fault
	# Rotation from Okada's axes to geographic
	ue = sin(strike).*ux .- cos(strike).*uy;
	un = cos(strike).*ux .+ sin(strike).*uy;


	# Tilt
	uzx =( -U1/(2*pi) .* chinnery_uzx_ss.(x,p,L,W,q,dip,nu) # strike-slip
		 .- U2/(2*pi) .* chinnery_uzx_ds.(x,p,L,W,q,dip,nu) # dip-slip
		 .+ U3/(2*pi) .* chinnery_uzx_tf.(x,p,L,W,q,dip,nu)); # tensile fault

	uzy =( -U1/(2*pi) .* chinnery_uzy_ss.(x,p,L,W,q,dip,nu) # strike-slip
		 .- U2/(2*pi) .* chinnery_uzy_ds.(x,p,L,W,q,dip,nu) # dip-slip
		 .+ U3/(2*pi) .* chinnery_uzy_tf.(x,p,L,W,q,dip,nu)); # tensile fault

	# Rotation from Okada's axes to geographic
	uze = -sin(strike).*uzx .+ cos(strike).*uzy;
	uzn = -cos(strike).*uzx .- sin(strike).*uzy;``

	# Strain
	uxx = ( -U1/(2*pi) .* chinnery_uxx_ss.(x,p,L,W,q,dip,nu) # strike-slip
			- U2/(2*pi) .* chinnery_uxx_ds.(x,p,L,W,q,dip,nu) # dip-slip
			+ U3/(2*pi) .* chinnery_uxx_tf.(x,p,L,W,q,dip,nu)); # tensile fault

	uxy = ( -U1/(2*pi) .* chinnery_uxy_ss.(x,p,L,W,q,dip,nu) # strike-slip
			.- U2/(2*pi) .* chinnery_uxy_ds.(x,p,L,W,q,dip,nu) # dip-slip
			.+ U3/(2*pi) .* chinnery_uxy_tf.(x,p,L,W,q,dip,nu)); # tensile fault

	uyx = ( -U1/(2*pi) .* chinnery_uyx_ss.(x,p,L,W,q,dip,nu) # strike-slip
		 .- U2/(2*pi) .* chinnery_uyx_ds.(x,p,L,W,q,dip,nu) # dip-slip
		.+ U3/(2*pi) .* chinnery_uyx_tf.(x,p,L,W,q,dip,nu)); # tensile fault

	uyy = ( -U1/(2*pi) .* chinnery_uyy_ss.(x,p,L,W,q,dip,nu) # strike-slip
			.- U2/(2*pi) .* chinnery_uyy_ds.(x,p,L,W,q,dip,nu) # dip-slip
			.+ U3/(2*pi) .* chinnery_uyy_tf.(x,p,L,W,q,dip,nu)); # tensile fault

	# Rotation from Okada's axes to geographic
	unn = cos(strike) .^ 2 .* uxx .+ sin(2*strike).*(uxy .+ uyx)/2 + sin(strike).^2 .*uyy;
	une = sin(2*strike) .* (uxx .- uyy)/2 + sin(strike).^2 .* uyx .- cos(strike).^2 .*uxy;
	uen = sin(2*strike) .* (uxx .- uyy)/2 - cos(strike).^2 .* uyx .+ sin(strike).^2 .*uxy;
	uee = sin(strike).^2 .* uxx .- sin(2*strike).*(uyx .+ uxy)./2 .+ cos(strike).^2 .*uyy;


	return ue,un,uz



end

"""
# Notes for I... and K... subfunctions:
#
#	1. original formulas use Lame's parameters as mu/(mu+lambda) which
#	   depends only on the Poisson's ratio = 1 - 2*nu
#	2. tests for cos(dip) == 0 are made with "cos(dip) > eps"
#	   because cos(90*pi/180) is not zero but = 6.1232e-17 (!)
#	   NOTE: don't use cosd and sind because of incompatibility
#	   with Matlab v6 and earlier...
"""


# ###################
# ##### TEST
# ####################
# x=2
# y=3
# d=4
# dip=70
# L=3
# W=2
# u3=0
# rake=0
# slip=1
#
#
# Okada85(x-L/2,y-cosd(dip)*W/2,d-sind(dip)*W/2,90,dip,L,W,rake,slip,u3)
#
# xi=x
# eta=1
# q=1
# dip=dip*pi/180
# nu=0.25
#
# ux_ss(xi,eta,q,dip*pi/180,0.25)

"""
 Chinnery's notation [equation (24) p. 1143]

 In Julia this seem like a very slow way of doing things becasue it may prevent the functions to be compiled
"""
function chinnery(f,x,p,L,W,q,dip,nu)
	u = eval(Symbol(f)).(x,p,q,dip,nu) - eval(Symbol(f)).(x,p-W,q,dip,nu) - eval(Symbol(f)).(x-L,p,q,dip,nu) + eval(Symbol(f))(x-L,p-W,q,dip,nu);
	return u
end

function chinnery_ux_ss(x,p,L,W,q,dip,nu)
	u =ux_ss.(x,p,q,dip,nu) - ux_ss.(x,p-W,q,dip,nu) - ux_ss.(x-L,p,q,dip,nu) + ux_ss.(x-L,p-W,q,dip,nu);
	return u
end

function chinnery_ux_ds(x,p,L,W,q,dip,nu)
	u =ux_ds.(x,p,q,dip,nu) - ux_ds.(x,p-W,q,dip,nu) - ux_ds.(x-L,p,q,dip,nu) + ux_ds.(x-L,p-W,q,dip,nu);
	return u
end

function chinnery_ux_tf(x,p,L,W,q,dip,nu)
	u =ux_tf.(x,p,q,dip,nu) - ux_tf.(x,p-W,q,dip,nu) - ux_tf.(x-L,p,q,dip,nu) + ux_tf.(x-L,p-W,q,dip,nu);
	return u
end

function chinnery_uy_ss(x,p,L,W,q,dip,nu)
	u =uy_ss.(x,p,q,dip,nu) - uy_ss.(x,p-W,q,dip,nu) - uy_ss.(x-L,p,q,dip,nu) + uy_ss.(x-L,p-W,q,dip,nu);
	return u
end

function chinnery_uy_ds(x,p,L,W,q,dip,nu)
	u =uy_ds.(x,p,q,dip,nu) - uy_ds.(x,p-W,q,dip,nu) - uy_ds.(x-L,p,q,dip,nu) + uy_ds.(x-L,p-W,q,dip,nu);
	return u
end

function chinnery_uy_tf(x,p,L,W,q,dip,nu)
	u =uy_tf.(x,p,q,dip,nu) - uy_tf.(x,p-W,q,dip,nu) - uy_tf.(x-L,p,q,dip,nu) + uy_tf.(x-L,p-W,q,dip,nu);
	return u
end


function chinnery_uz_ss(x,p,L,W,q,dip,nu)
	u =uz_ss.(x,p,q,dip,nu) - uz_ss.(x,p-W,q,dip,nu) - uz_ss.(x-L,p,q,dip,nu) + uz_ss.(x-L,p-W,q,dip,nu);
	return u
end

function chinnery_uz_ds(x,p,L,W,q,dip,nu)
	u =uz_ds.(x,p,q,dip,nu) - uz_ds.(x,p-W,q,dip,nu) - uz_ds.(x-L,p,q,dip,nu) + uz_ds.(x-L,p-W,q,dip,nu);
	return u
end

function chinnery_uz_tf(x,p,L,W,q,dip,nu)
	u =uz_tf.(x,p,q,dip,nu) - uz_tf.(x,p-W,q,dip,nu) - uz_tf.(x-L,p,q,dip,nu) + uz_tf.(x-L,p-W,q,dip,nu);
	return u
end

function chinnery_uzx_ss(x,p,L,W,q,dip,nu)
	u =uzx_ss.(x,p,q,dip,nu) - uzx_ss.(x,p-W,q,dip,nu) - uzx_ss.(x-L,p,q,dip,nu) + uzx_ss.(x-L,p-W,q,dip,nu);
	return u
end

function chinnery_uzx_ds(x,p,L,W,q,dip,nu)
	u =uzx_ds.(x,p,q,dip,nu) - uzx_ds.(x,p-W,q,dip,nu) - uzx_ds.(x-L,p,q,dip,nu) + uzx_ds.(x-L,p-W,q,dip,nu);
	return u
end

function chinnery_uzx_tf(x,p,L,W,q,dip,nu)
	u =uzx_tf.(x,p,q,dip,nu) - uzx_tf.(x,p-W,q,dip,nu) - uzx_tf.(x-L,p,q,dip,nu) + uzx_tf.(x-L,p-W,q,dip,nu);
	return u
end

function chinnery_uzy_ss(x,p,L,W,q,dip,nu)
	u =uzy_ss.(x,p,q,dip,nu) - uzy_ss.(x,p-W,q,dip,nu) - uzy_ss.(x-L,p,q,dip,nu) + uzy_ss.(x-L,p-W,q,dip,nu);
	return u
end

function chinnery_uzy_ds(x,p,L,W,q,dip,nu)
	u =uzy_ds.(x,p,q,dip,nu) - uzy_ds.(x,p-W,q,dip,nu) - uzy_ds.(x-L,p,q,dip,nu) + uzy_ds.(x-L,p-W,q,dip,nu);
	return u
end

function chinnery_uzy_tf(x,p,L,W,q,dip,nu)
	u =uzy_tf.(x,p,q,dip,nu) - uzy_tf.(x,p-W,q,dip,nu) - uzy_tf.(x-L,p,q,dip,nu) + uzy_tf.(x-L,p-W,q,dip,nu);
	return u
end

function chinnery_uxx_ss(x,p,L,W,q,dip,nu)
	u =uxx_ss.(x,p,q,dip,nu) - uxx_ss.(x,p-W,q,dip,nu) - uxx_ss.(x-L,p,q,dip,nu) + uxx_ss.(x-L,p-W,q,dip,nu);
	return u
end

function chinnery_uxx_ds(x,p,L,W,q,dip,nu)
	u =uxx_ds.(x,p,q,dip,nu) - uxx_ds.(x,p-W,q,dip,nu) - uxx_ds.(x-L,p,q,dip,nu) + uxx_ds.(x-L,p-W,q,dip,nu);
	return u
end

function chinnery_uxx_tf(x,p,L,W,q,dip,nu)
	u =uxx_tf.(x,p,q,dip,nu) - uxx_tf.(x,p-W,q,dip,nu) - uxx_tf.(x-L,p,q,dip,nu) + uxx_tf.(x-L,p-W,q,dip,nu);
	return u
end

function chinnery_uyy_ss(x,p,L,W,q,dip,nu)
	u =uyy_ss.(x,p,q,dip,nu) - uyy_ss.(x,p-W,q,dip,nu) - uyy_ss.(x-L,p,q,dip,nu) + uyy_ss.(x-L,p-W,q,dip,nu);
	return u
end

function chinnery_uyy_ds(x,p,L,W,q,dip,nu)
	u =uyy_ds.(x,p,q,dip,nu) - uyy_ds.(x,p-W,q,dip,nu) - uyy_ds.(x-L,p,q,dip,nu) + uyy_ds.(x-L,p-W,q,dip,nu);
	return u
end

function chinnery_uyy_tf(x,p,L,W,q,dip,nu)
	u =uyy_tf.(x,p,q,dip,nu) - uyy_tf.(x,p-W,q,dip,nu) - uyy_tf.(x-L,p,q,dip,nu) + uyy_tf.(x-L,p-W,q,dip,nu);
	return u
end

function chinnery_uxy_ss(x,p,L,W,q,dip,nu)
	u =uxy_ss.(x,p,q,dip,nu) - uxy_ss.(x,p-W,q,dip,nu) - uxy_ss.(x-L,p,q,dip,nu) + uxy_ss.(x-L,p-W,q,dip,nu);
	return u
end

function chinnery_uxy_ds(x,p,L,W,q,dip,nu)
	u =uxy_ds.(x,p,q,dip,nu) - uxy_ds.(x,p-W,q,dip,nu) - uxy_ds.(x-L,p,q,dip,nu) + uxy_ds.(x-L,p-W,q,dip,nu);
	return u
end

function chinnery_uxy_tf(x,p,L,W,q,dip,nu)
	u =uxy_tf.(x,p,q,dip,nu) - uxy_tf.(x,p-W,q,dip,nu) - uxy_tf.(x-L,p,q,dip,nu) + uxy_tf.(x-L,p-W,q,dip,nu);
	return u
end

function chinnery_uyx_ss(x,p,L,W,q,dip,nu)
	u =uyx_ss.(x,p,q,dip,nu) - uyx_ss.(x,p-W,q,dip,nu) - uyx_ss.(x-L,p,q,dip,nu) + uyx_ss.(x-L,p-W,q,dip,nu);
	return u
end

function chinnery_uyx_ds(x,p,L,W,q,dip,nu)
	u =uyx_ds.(x,p,q,dip,nu) - uyx_ds.(x,p-W,q,dip,nu) - uyx_ds.(x-L,p,q,dip,nu) + uyx_ds.(x-L,p-W,q,dip,nu);
	return u
end

function chinnery_uyx_tf(x,p,L,W,q,dip,nu)
	u =uyx_tf.(x,p,q,dip,nu) - uyx_tf.(x,p-W,q,dip,nu) - uyx_tf.(x-L,p,q,dip,nu) + uyx_tf.(x-L,p-W,q,dip,nu);
	return u
end
"""
Displacement subfunctions
strike-slip displacement subfunctions [equation (25) p. 1144]
"""
function ux_ss(xi,eta,q,dip,nu)
	R = sqrt.(xi.^2 + eta.^2 + q.^2);
	u = xi.*q./(R.*(R + eta)) + I1(xi,eta,q,dip,nu,R).*sin(dip);
	k = findall(q.!=0);
	if q!=0
		u = u + atan(xi.*eta./(q.*R));
	end
	return u
end

function uy_ss(xi,eta,q,dip,nu)
	R = sqrt.(xi.^2 + eta.^2 + q.^2);
	u = (eta.*cos.(dip) + q.*sin.(dip)).*q./(R.*(R + eta)) + q.*cos.(dip)./(R + eta)+ I2.(eta,q,dip,nu,R).*sin.(dip);
	return u
end

function uz_ss(xi,eta,q,dip,nu)
	R = sqrt.(xi.^2 + eta.^2 + q.^2);
	db = eta.*sin.(dip) - q.*cos.(dip);
	u = (eta.*sin.(dip) - q.*cos.(dip)).*q./(R.*(R + eta)) + q.*sin.(dip)./(R + eta) + I4.(db,eta,q,dip,nu,R).*sin.(dip);
	return u
end

"""
dip-slip displacement subfunctions [equation (26) p. 1144]
"""
function ux_ds(xi,eta,q,dip,nu)
	R = sqrt.(xi.^2 + eta.^2 + q.^2);
	u = q./R - I3.(eta,q,dip,nu,R).*sin.(dip).*cos.(dip);
	return u
end

function uy_ds(xi,eta,q,dip,nu)
	R = sqrt.(xi.^2 + eta.^2 + q.^2);
	u = (eta.*cos.(dip) + q.*sin.(dip)).*q./(R.*(R + xi)) - I1.(xi,eta,q,dip,nu,R).*sin.(dip).*cos.(dip);
	#k = findall(q.!=0);
	if q!=0
		u = u + cos.(dip).*atan.(xi.*eta./(q.*R));
	end
	#u[k] = u[k] + cos.(dip).*atan.(xi[k].*eta[k]./(q[k].*R[k]));
	return u
end
function uz_ds(xi,eta,q,dip,nu)
	R = sqrt.(xi.^2 + eta.^2 + q.^2);
	db = eta.*sin.(dip) - q.*cos.(dip);
	u = db.*q./(R.*(R + xi)) - I5.(xi,eta,q,dip,nu,R,db).*sin.(dip).*cos.(dip);
	#k = findall(q.!=0);
	if q!=0
		u = u + sin.(dip).*atan.(xi.*eta./(q.*R));
	end
	#u[k] = u[k] + sin.(dip).*atan.(xi[k].*eta[k]./(q[k].*R[k]));
	return u
end

"""
tensile fault displacement subfunctions [equation (27) p. 1144]
"""
function ux_tf(xi,eta,q,dip,nu)
	R = sqrt.(xi.^2 + eta.^2 + q.^2);
	u = q.^2 ./(R.*(R + eta)) - I3.(eta,q,dip,nu,R).*sin.(dip).^2;
	return u
end

function uy_tf(xi,eta,q,dip,nu)
	R = sqrt.(xi.^2 + eta.^2 + q.^2);
	u = -(eta.*sin.(dip) - q.*cos.(dip)).*q./(R.*(R + xi)) - sin.(dip).*xi.*q./(R.*(R + eta)) - I1.(xi,eta,q,dip,nu,R).*sin.(dip).^2;
	#k = findall(q.!=0);
	if q!=0
		u = u + sin.(dip).*atan.(xi.*eta./(q.*R));
	end
	#u[k] = u[k] + sin.(dip).*atan.(xi[k].*eta[k]./(q[k].*R[k]));
	return u
end

function uz_tf(xi,eta,q,dip,nu)
	R = sqrt.(xi.^2 + eta.^2 + q.^2);
	db = eta.*sin.(dip) - q.*cos.(dip);
	u = (eta.*cos.(dip) + q.*sin.(dip)).*q./(R.*(R + xi)) + cos.(dip).*xi.*q./(R.*(R + eta)) - I5.(xi,eta,q,dip,nu,R,db).*sin.(dip).^2;
	#k = findall(q.!=0);
	if q!=0
		u = u + cos.(dip).*atan.(xi.*eta./(q.*R));
	end
	#u[k] = u[k] - cos.(dip).*atan.(xi[k].*eta[k]./(q[k].*R[k]));
	return u
end

"""
I... displacement subfunctions [equations (28) (29) p. 1144-1145]
"""
function I1(xi,eta,q,dip,nu,R)
	db = eta.*sin.(dip) - q.*cos(dip);
	if cos(dip) > eps()
		I = (1 - 2*nu) * (-xi./(cos(dip).*(R + db)))- sin(dip)./cos(dip).*I5(xi,eta,q,dip,nu,R,db);
	else
		I = -(1 - 2*nu)/2 * xi.*q./(R + db).^2;
	end
	return I
end


function I2(eta,q,dip,nu,R)
	return (1 - 2*nu) * (-log(R + eta)) - I3(eta,q,dip,nu,R);
end

function I3(eta,q,dip,nu,R)
	yb = eta.*cos(dip) + q.*sin(dip);
	db = eta.*sin(dip) - q.*cos(dip);
	if cos(dip) > eps()
		I = (1 - 2*nu) * (yb./(cos(dip)*(R + db)) - log(R + eta)) + sin(dip)./cos(dip) * I4(db,eta,q,dip,nu,R);
	else
		I = (1 - 2*nu)/2 * (eta./(R + db) + yb.*q./(R + db).^2 - log(R + eta));
	end
	return I
end

function I4(db,eta,q,dip,nu,R)
	if cos(dip) > eps()
		I = (1 - 2*nu) * 1 ./ cos(dip) * (log(R + db) - sin(dip).*log(R + eta));
	else
		I = -(1 - 2*nu) * q ./ (R + db);
	end
	return I
end

function I5(xi,eta,q,dip,nu,R,db)
	X = sqrt(xi.^2 + q.^2);
	if cos(dip) > eps()
		I = (1 - 2*nu) * 2 ./ cos(dip) .* atan((eta.*(X + q.*cos(dip)) + X.*(R + X).*sin(dip)) ./(xi.*(R + X).*cos(dip)));
		if xi == 0
			I = 0;
		end
	else
		I = -(1 - 2*nu) * xi.*sin(dip)./(R + db);
	end
	return I
end

"""
 Tilt subfunctions
 strike-slip tilt subfunctions [equation (37) p. 1147]
"""
function uzx_ss(xi,eta,q,dip,nu)
	R = sqrt.(xi.^2 + eta.^2 + q.^2);
	u = -xi .* q .^ 2 .* A.(eta,R).*cos.(dip) + ((xi.*q)./R.^3 - K1.(xi,eta,q,dip,nu,R)).*sin.(dip);
	return u
end

function uzy_ss(xi,eta,q,dip,nu)
	R = sqrt.(xi.^2 + eta.^2 + q.^2);
	db = eta.*sin(dip) - q.*cos(dip);
	yb = eta.*cos(dip) + q.*sin(dip);
	u = (db.*q./R.^3).*cos(dip) + (xi .^ 2 .* q.*A(eta,R).*cos(dip) - sin(dip)./R + yb.*q./R.^3 - K2(xi,eta,q,dip,nu,R)).*sin(dip);
	return u
end

"""
 dip-slip tilt subfunctions [equation (38) p. 1147]
"""
function uzx_ds(xi,eta,q,dip,nu)
	R = sqrt(xi.^2 + eta.^2 + q.^2);
	db = eta.*sin(dip) - q.*cos(dip);
	u = db.*q./R.^3 + q.*sin(dip)./(R.*(R + eta)) + K3(xi,eta,q,dip,nu,R).*sin(dip).*cos(dip);
	return u
end

function uzy_ds(xi,eta,q,dip,nu)
	R = sqrt(xi.^2 + eta.^2 + q.^2);
	db = eta.*sin(dip) - q.*cos(dip);
	yb = eta.*cos(dip) + q.*sin(dip);
	u = yb.*db.*q.*A(xi,R) 	- (2*db./(R.*(R + xi)) + xi.*sin(dip)./(R.*(R + eta))).*sin(dip) + K1(xi,eta,q,dip,nu,R).*sin(dip).*cos(dip);
	return u
end
"""
 tensile fault tilt subfunctions [equation (39) p. 1147]
"""
function uzx_tf(xi,eta,q,dip,nu)
	R = sqrt(xi.^2 + eta.^2 + q.^2);
	u = q .^ 2 ./R .^ 3 .* sin(dip) - q .^ 3 .*A(eta,R).*cos(dip) + K3(xi,eta,q,dip,nu,R).*sin(dip).^2;
	return u
end

function uzy_tf(xi,eta,q,dip,nu)
	R = sqrt(xi.^2 + eta.^2 + q.^2);
	db = eta.*sin(dip) - q.*cos(dip);
	yb = eta.*cos(dip) + q.*sin(dip);
	u = (yb.*sin(dip) + db.*cos(dip)).*q .^ 2 .*A(xi,R) + xi.*q .^ 2 .*A(eta,R).*sin(dip).*cos(dip) - (2*q./(R.*(R + xi)) - K1(xi,eta,q,dip,nu,R)).*sin(dip).^2;
	return u
end

function A(x,R)
	a = (2*R + x)./(R .^ 3 .*(R + x).^2);
	return a
end

"""
 K... tilt subfunctions [equations (40) (41) p. 1148]
"""
function K1(xi,eta,q,dip,nu,R)
	db = eta.*sin(dip) - q.*cos(dip);
	if cos(dip) > eps()
		K = (1 - 2*nu) * xi./cos(dip) .* (1 ./ (R.*(R + db)) - sin(dip)./(R.*(R + eta)));
	else
		K = (1 - 2*nu) * xi.*q./(R.*(R + db).^2);
	end
	return K
end

function K2(xi,eta,q,dip,nu,R)
	K = (1 - 2*nu) * (-sin(dip)./R + q.*cos(dip)./(R.*(R + eta))) - K3(xi,eta,q,dip,nu,R);
	return K
end

function K3(xi,eta,q,dip,nu,R)
	db = eta.*sin(dip) - q.*cos(dip);
	yb = eta.*cos(dip) + q.*sin(dip);
	if cos(dip) > eps()
		K = (1 - 2*nu) * 1 ./cos(dip) .* (q./(R.*(R + eta)) - yb./(R.*(R + db)));
	else
		K = (1 - 2*nu) * sin(dip)./(R + db) .* (xi.^ 2 ./(R.*(R + db)) - 1);
	end
	return K
end


"""
 Strain subfunctions
 strike-slip strain subfunctions [equation (31) p. 1145]
"""
function uxx_ss(xi,eta,q,dip,nu)
	R = sqrt(xi.^2 + eta.^2 + q.^2);
	u = xi .^ 2 .* q .* A(eta,R) - J1(xi,eta,q,dip,nu,R).*sin(dip);
	return u
end

function uxy_ss(xi,eta,q,dip,nu)
	R = sqrt(xi.^2 + eta.^2 + q.^2);
	db = eta.*sin(dip) - q.*cos(dip);
	u = xi .^3 .* db./(R .^ 3 .*(eta .^ 2 + q .^ 2))- (xi .^ 3 .*A(eta,R) + J2(xi,eta,q,dip,nu,R)).*sin(dip);
	return u
end

function uyx_ss(xi,eta,q,dip,nu)
	R = sqrt(xi.^2 + eta.^2 + q.^2);
	u = xi .* q ./ R .^ 3 .* cos(dip) + (xi .* q .^ 2 .* A(eta,R) - J2(xi,eta,q,dip,nu,R)).*sin(dip);
	return u
end

function uyy_ss(xi,eta,q,dip,nu)
	R = sqrt(xi.^2 + eta.^2 + q.^2);
	yb = eta.*cos(dip) + q.*sin(dip);
	u = yb .* q ./ R .^ 3 .* cos(dip) + (q .^ 3 .* A(eta,R).*sin(dip) - 2*q.*sin(dip)./(R.*(R + eta)) - (xi.^2 + eta.^2)./ R .^ 3 .* cos(dip) - J4(xi,eta,q,dip,nu,R)).*sin(dip);
	return u
end

"""
 dip-slip strain subfunctions [equation (32) p. 1146]
"""
function uxx_ds(xi,eta,q,dip,nu)
	R = sqrt(xi.^2 + eta.^2 + q.^2);
	u = xi.*q./R.^3 + J3(xi,eta,q,dip,nu,R).*sin(dip).*cos(dip);
	return u
end

function uxy_ds(xi,eta,q,dip,nu)
	R = sqrt(xi.^2 + eta.^2 + q.^2);
	yb = eta.*cos(dip) + q.*sin(dip);
	u = yb.*q./R.^3 - sin(dip)./R + J1(xi,eta,q,dip,nu,R).*sin(dip).*cos(dip);
	return u
end

function uyx_ds(xi,eta,q,dip,nu)
	R = sqrt(xi.^2 + eta.^2 + q.^2);
	yb = eta.*cos(dip) + q.*sin(dip);
	u = yb.*q./R.^3 + q.*cos(dip)./(R.*(R + eta)) + J1(xi,eta,q,dip,nu,R).*sin(dip).*cos(dip);
	return u
end

function uyy_ds(xi,eta,q,dip,nu)
	R = sqrt(xi.^2 + eta.^2 + q.^2);
	yb = eta.*cos(dip) + q.*sin(dip);
	u = yb .^ 2 .* q .* A(xi,R) - (2*yb./(R.*(R + xi)) + xi.*cos(dip)./(R.*(R + eta))).*sin(dip) + J2(xi,eta,q,dip,nu,R).*sin(dip).*cos(dip);
	return u
end

"""
 tensile fault strain subfunctions [equation (33) p. 1146]
"""
function uxx_tf(xi,eta,q,dip,nu)
	R = sqrt(xi.^2 + eta.^2 + q.^2);
	u = xi .* q .^ 2 .* A(eta,R) + J3(xi,eta,q,dip,nu,R).*sin(dip).^2;
	return u
end

function uxy_tf(xi,eta,q,dip,nu)
	R = sqrt(xi.^2 + eta.^2 + q.^2);
	db = eta.*sin(dip) - q.*cos(dip);
	u = -db.*q./ R .^ 3 - xi.^ 2 .*q.*A(eta,R).*sin(dip) + J1(xi,eta,q,dip,nu,R).*sin(dip).^2;
	return u
end

function uyx_tf(xi,eta,q,dip,nu)
	R = sqrt(xi.^2 + eta.^2 + q.^2);
	u = q .^ 2 ./R .^ 3 .*cos(dip) + q .^ 3 .* A(eta,R).*sin(dip) + J1(xi,eta,q,dip,nu,R).*sin(dip).^2;
	return u
end

function uyy_tf(xi,eta,q,dip,nu)
	R = sqrt(xi.^2 + eta.^2 + q.^2);
	db = eta.*sin(dip) - q.*cos(dip);
	yb = eta.*cos(dip) + q.*sin(dip);
	u = (yb .* cos(dip) - db.*sin(dip)) .* q .^ 2 .* A(xi,R) - q.*sin(2*dip)./(R.*(R + xi)) - (xi .* q .^ 2 .*A(eta,R) - J2(xi,eta,q,dip,nu,R)).*sin(dip).^2;
	return u
end

"
 J... tensile fault subfunctions [equations (34) (35) p. 1146-1147]
"
function J1(xi,eta,q,dip,nu,R)
	db = eta.*sin(dip) - q.*cos(dip);
	if cos(dip) > eps()
		J = (1 - 2*nu) * 1 ./cos(dip) * (xi .^ 2 ./(R.*(R + db).^2) - 1 ./(R + db)) - sin(dip)./cos(dip)*K3(xi,eta,q,dip,nu,R);
	else
		J = (1 - 2*nu)/2 * q./(R + db).^2 .* (2 * xi .^ 2 ./(R.*(R + db)) - 1);
	end
	return J
end

function J2(xi,eta,q,dip,nu,R)
	db = eta.*sin(dip) - q.*cos(dip);
	yb = eta.*cos(dip) + q.*sin(dip);
	if cos(dip) > eps()
		J = (1 - 2*nu) * 1 ./cos(dip) * xi.*yb./(R.*(R + db).^2) - sin(dip)./cos(dip)*K1(xi,eta,q,dip,nu,R);
	else
		J = (1 - 2*nu)/2 * xi.*sin(dip)./(R + db).^ 2 .* (2*q .^ 2 ./(R.*(R + db)) - 1);
	end
	return J
end

function J3(xi,eta,q,dip,nu,R)
	J = (1 - 2*nu) * -xi./(R.*(R + eta)) - J2(xi,eta,q,dip,nu,R);
	return J
end

function J4(xi,eta,q,dip,nu,R)
	J = (1 - 2*nu) * (-cos(dip)./R - q.*sin(dip)./(R.*(R + eta))) - J1(xi,eta,q,dip,nu,R);
	return J
end
