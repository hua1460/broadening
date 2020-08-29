function [x, y] = dogaussian_general(x0,y0,Gamma1, Gamma2, E1, E2, E3, E4, Npoints)
% dogaussian_general : do gaussian convolution with variable hwhm

% Formula:
%==========
% the distribution function is:
%                                        (x - x0)^2
% f_g(x; x0, sigma)  =  a *  exp{ -  -------------------------}
%                                        2*sigma^2 
%  Here: 
%  constant =  2 * sqrt ( 2 * ln (2) ) = 2.35482004503095 
%  sigma is standard deviation: sigma = FWHM/constant       
%  a is normalization factor:  a = 1 / (sigma *   sqrt(2*pi) )
%  (formula from http://en.wikipedia.org/wiki/Gaussian_function
% See also http://hyperphysics.phy-astr.gsu.edu/hbase/math/gaufcn.html)
%
%   Scheme:
%==========
%   E1......E2......E2......E4  
%     Gamma1          Gamma2
%
%   [E1 E2]; Gamma1
%   (E2 E3); linear increasing from Gamma1 to Gamma2
%   [E3 E4]; Gamma2

%  INPUT VARIABLES:
%===================
%	x0, y0       -- bar values,should be COLUMN VECTOR
% 	Gamma        -- 1/2 * FWHM 
% 	E1, E2 -- the area you want do calculate the Lorentz y
% 	Npoints      -- Npoints you want to calulate in this area(prefer an ODD number)
%
%  OUTPUT VARIABLES:
%===================
% 	x       -- sample points , ROW VECTOR
% 	y       -- value at the sample point,ROW VECTOR
%
%2011-10-17 extend to use different Gamma values in different range.

% ###########################
% ## CopyRight Weijie Hua ###
% ## huaweijie@gmail.com  ###
% ###########################

constant = 2.35482004503095;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the sample points x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
siz = size(x0,1); 
dE = (E4 - E1)/(Npoints - 1);
x = [E1:dE:E4];           % x:  1*Npoints

X = repmat(x, siz, 1);    % X:  siz*Npoints

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
% const
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FWHM =  2 * BIG_Gamma; sigma = FWHM/constant;
csq = sigma.*sigma; a = 1./(sigma *sqrt(2*pi));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate f
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lorentz y at the sample points
X0 = repmat(x0, 1, Npoints);
Y0 = repmat(y0, 1, Npoints);

f = a.* exp(- (X - X0).^2.0 ./2.0./csq);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = Y0 .* f; 
y = sum(Y, 1);


%------------------------------------
%[driver to test function dogaussian]
%clear all,clc,format long ;
%x0=0; y0=1;E1=-10;E2=10;Npoints=21;FWHW=1;
% Gamma= 1/2 * FWHW;
%[x, y] = dogaussian(x0,y0,Gamma, E1, E2, Npoints);
%ans=[x' y']

