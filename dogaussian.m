function [x, y] = dogaussian(x0, y0, Gamma, E1, E2, Npoints)
%dogaussian : do gaussian convolution
% the distribution function is:
%                                        (x - x0)^2
% f_g(x; x0, sigma)  =  a *  exp{ -  -------------------------}
%                                        2*sigma^2 
%  Here: 
%  constant1 =  2 * sqrt ( 2 * ln (2) ) = 2.35482004503095 
%  sigma is standard deviation: sigma = FWHM/constant1       
%  a is normalization factor:  a = 1 / (sigma *   sqrt(2*pi) )
%  (formula from http://en.wikipedia.org/wiki/Gaussian_function
% See also http://hyperphysics.phy-astr.gsu.edu/hbase/math/gaufcn.html)
%
%  INPUT VARIABLES:
%	x0, y0       -- bar values,should be COLUMN VECTOR
% 	Gamma        -- 1/2 * FWHM 
% 	E1, E2 -- the area you want do calculate the gaussian y
% 	Npoints      -- Npoints you want to calulate in this area(prefer an ODD number)
%
%  OUTPUT VARIABLES:
% 	x       -- sample points , ROW VECTOR
% 	y       -- value at the sample point,ROW VECTOR
%

%------------------------------------------------------------------
% Copyright: Weijie HUA 
% 			 Department of Theoretical Chemistry and Biology
%            Royal Institute of Technology, SWEDEN
%            huaweijie@gmail.com 
%-----------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% const
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FWHM =  2 * Gamma;
constant1 = 2.35482004503095;
sigma = FWHM/constant1;
csq = sigma*sigma;
a = 1/(sigma*sqrt(2*pi));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the sample points x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
siz = size(x0,1); 
dE = (E2 - E1)/(Npoints - 1);
x = [E1:dE:E2];           % x    :     1*Npoints

X = repmat(x, siz, 1); % X:  siz*Npoints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate f
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lorentz y at the sample points
X0 = repmat(x0, 1, Npoints);
Y0 = repmat(y0, 1, Npoints);

f = a.* exp(- (X - X0).^2.0 /2.0/csq);
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


