function [lamda, y] = dogaussian_ev2nm(x0, y0,  Gamma, lamda1, lamda2, Npoints)
%function [lamda, y] = dogaussian_ev2nm(x0, y0,  Gamma, lamda1, lamda2, Npoints)
%           |    |                      |   |       |    |       | 
%           nm  arb.unit                eV  arb.unit eV   nm      nm
% DO GAUSSIAN CONVOLUTION for bar data in eV and change into data in nm
% (FOR, e.g., UV ABSORPTION SPECTRUM) 
% 
% THE INPUT SETS [x0 y0] are the raw "energy (eV)"-"oscillator strnegth" data
% read from quantum chemistry software. 
% THE OUTPUT SETS [lamda y] are the wave length (nm)-"oscillator strength"
% data that can be compared directly with experiment
% NOTE: The difference between dogaussian and dogaussian_ev2nm lies in that
% output variable "lamda" and input variables  "lamda1" "lamda2" are in nm
% units
% 2018-02-27 Nanjing 

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
% 	lamda1, lamda2   -- the area you want do calculate the gaussian y
% 	Npoints      -- Npoints you want to calulate in this area(prefer an ODD number)
%
%  OUTPUT VARIABLES:
% 	lamda   -- sample points , ROW VECTOR
% 	y       -- value at the sample point,ROW VECTOR
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% const
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
constant1 = 2.35482004503095;

% E = hc/lamda   1nm ---> 1239.84193 eV
% http://halas.rice.edu/conversions
constant2 = 1239.84193;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FWHM =  2 * Gamma;
sigma = FWHM/constant1;
csq = sigma*sigma;
a = 1/(sigma*sqrt(2*pi));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the sample points lamda
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
siz = size(x0,1); 
dlamda = (lamda2 - lamda1)/(Npoints - 1);
lamda = [lamda1:dlamda:lamda2];    % lamda :     1*Npoints      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lamda in unit of nm,  x is in unit of eV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
x =  constant2./lamda;  
X = repmat(x, siz, 1); % X:  siz*Npoints

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate f
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gaussian y at the sample points
X0 = repmat(x0, 1, Npoints);
Y0 = repmat(y0, 1, Npoints);

f = a.* exp(- (X - X0).^2.0 /2.0/csq);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = Y0 .* f; 
y = sum(Y, 1);




