function [x, y] = dolorentz_general(x0,y0,Gamma1, Gamma2, E1, E2, E3, E4, Npoints)
%dolorentz_general : do lorentzian convolution with variable HWHM.
% Formula:
%==========
% the distribution function is:
%                    1          Gamma
% f(x; x0, Gamma) = ----  -------------------------
%                    pi    (x-x0)^2 + Gamma^2
%(formulae comes from http://en.wikipedia.org/wiki/Lorentzian_function)

%   Scheme:
%==========
%   E1......E2......E2......E4  
%     Gamma1          Gamma2
%
%   [E1 E2]; Gamma1
%   (E2 E3); linear increasing from Gamma1 to Gamma2
%   [E3 E4]; Gamma2

%  INPUT VARIABLES:
%==========
%	x0, y0                -- bar values,should be COLUMN VECTOR
% 	Gamma1, Gamma2        -- half-width at half-maximum (HWHM)  
%                            Gamma = FWHM/2
% 	E1, E2, E3, E4        -- the area you want do calculate the Lorentz y
% 	Npoints               -- Npoints you want to calulate in this area(prefer an ODD number)

%
%  OUTPUT VARIABLES:
%===================
% 	x       -- sample points , ROW VECTOR
% 	y       -- value at the sample point,ROW VECTOR

%   Update:
%==========
%2011-03-03 extend to use different Gamma values in different range.

%------------------------------------------------------------------
% Copyright: Weijie HUA 
%            Department of Theoretical Chemistry and Biology
%            Royal Institute of Technology, SWEDEN
%            huaweijie@gmail.com 
%-----------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the sample points x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
siz = size(x0,1); 
dE = (E4 - E1)/(Npoints - 1);
x = [E1:dE:E4];           % x    :     1*Npoints

X = repmat(x, siz, 1); % X:  siz*Npoints

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate BIG_Gamma 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gamma = ones (size(x0));
ind_a = find (x0 <= E2 ); Gamma(ind_a) = Gamma(ind_a) * Gamma1;
ind_c = find (x0 >= E3 ); Gamma(ind_c) = Gamma(ind_c) * Gamma2;


ind_b = find ( E2 < x0  & x0 < E3 );
dGamma = Gamma2 - Gamma1;
Gamma(ind_b) = Gamma(ind_b) .* ( Gamma1 + (Gamma(ind_b) - Gamma1) .* dGamma );
BIG_Gamma = repmat(Gamma, 1, Npoints);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate f
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lorentz y at the sample points
X0 = repmat(x0, 1, Npoints);
Y0 = repmat(y0, 1, Npoints);
f = BIG_Gamma./pi./((X - X0).^2.0 + BIG_Gamma.^2.0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = Y0 .* f; 
y = sum(Y, 1);


%------------------------------------
%[driver to test function dolorentz_general]
%clear all,clc,format long ;
%x0=0; y0=1;E1=-10;E2=-5; E3=5; E2=10;Npoints=21;Gamma1=0.5; Gamma2=1.5;
%[x, y] = dolorentz_general(x0,y0, Gamma1, Gamma2, E1, E2,E3, E4, Npoints);
%ans=[x' y']
