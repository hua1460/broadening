function [x, y] = dolorentz(x0,y0,Gamma, E1, E2, Npoints)
%updated 2012.05.14
%dolorentz : do lorentzian convolution
% the distribution function is:
%                    1          Gamma
% f(x; x0, Gamma) = ----  -------------------------
%                    pi    (x-x0)^2 + Gamma^2
%(formulae comes from http://en.wikipedia.org/wiki/Lorentzian_function)

%  INPUT VARIABLES:
%	x0, y0       -- bar values,should be COLUMN VECTOR
% 	Gamma        -- half-width at half-maximum (HWHM)  
%                    Gamma = FWHM/2
% 	E1, E2 -- the area you want do calculate the Lorentz y
% 	Npoints      -- Npoints you want to calulate in this area(prefer an ODD number)
%
%  OUTPUT VARIABLES:
% 	x       -- sample points , ROW VECTOR
% 	y       -- value at the sample point,ROW VECTOR

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
dE = (E2 - E1)/(Npoints - 1);
x = [E1:dE:E2];           % x    :     1*Npoints
y = zeros(size(x));
factor = Gamma/pi;

for j = 1:1:siz 
	y = y + y0(j) * factor./( ( x - x0 (j) ).^2.0 + Gamma^2.0);
end

%X = repmat(x, siz, 1); % X:  siz*Npoints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate f
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lorentz y at the sample points
%X0 = repmat(x0, 1, Npoints);
%factor = Gamma/pi;
%f = factor./((X - X0).^2.0 + Gamma^2.0); clear X0 X;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Y0 = repmat(y0, 1, Npoints);
%Y = Y0 .* f; clear Y0 f;
%y = sum(Y, 1); clear Y;
%
%------------------------------------
%[driver to test function dolorentz]
%clear all,clc,format long ;
%x0=0; y0=1;E1=-10;E2=10;Npoints=21;Gamma=0.5;
%[x, y] = dolorentz(x0,y0,Gamma, E1, E2, Npoints);
%ans=[x' y']
