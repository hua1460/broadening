function convoluted = broadening (barspectra, spin, convtype, Gamma, ...
                                    varargin)

% do broadening of barspectra

% Variables:
% ==================
%  barspectra can have 10, 7, 5, or 2 columns (see below)
%  Gamma : Gamma is HWHM, i.e., FWHW/2.
%          [Gamma]
%          [Gamma1 Gamma2]
%  spin:  (will only influence the filename of  the  output  file)
%        =1
%        =2
%  convtype:  g, l, G
%
% Optional variables for varargin
% ===============================
% (string)    postfix:  postfix of the saved data files 
% (numerical) Elimit: [0]
%                     [E1 E2]
%                     [E1 E2 E3 E4]
%
% They can both appear (no limitation of order), or only appear one of them.
% The program can determine automatically according to string or numerical types.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Updated 2019.06.24 NJUST
%updated 2018.10.06
%Modified 2015.09.09
%Updated 2013.11.27
% = Copyright: Weijie HUA, KTH, Stockholm  =
% =     Email: huaweijie@gmail.com         =
% =        All rights reserved.            =
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%  E1, E2: left and right limits
%  dE    : stepsize

subsection('Convoluting the barspectra');

%. Define parameters
default_dE         = 0.01;
default_range      = 40;
default_percentL    = 0.20; 
default_percentR    = 0.20; 

dE              = default_dE       ; 
range           = default_range    ; 
percentL         = default_percentL; 
percentR         = default_percentR; 

%.
nargopt = nargin - 4;

if (nargopt == 0);

	postfix = '';
	Elimit = 0;

elseif (nargopt == 1 || nargopt == 2)

	Elimit = 0;

	for k = 1:1:nargopt

		if(isstr( varargin{k} ))
			postfix = varargin{k};
		elseif(isnumeric ( varargin{k}))
			Elimit = varargin{k};
		end

	end

else
	error ('Wrong number of optional argumnemts, can only be 0,1,2!');
end


%. Define output file names
[File1, File2] = define_output_filenames(spin, convtype, postfix);
File1
File2

if (strcmp(File1,File2))
    error('The same File1 and File2 specified');
end
%. Output Bars
ind_end = print_bars(barspectra, File1);

%elapsed_time = etime (clock, globalt0);
%fprintf('bar barspectra printed out, elapsed_time = %12.1f (mins)\n',  elapsed_time/60);


%. Plot the convoluted barspectra
x0         = barspectra (:,1);
minx0      = min(x0);
maxx0      = max(x0);
heightL    = (maxx0 - minx0 ) * percentL;
heightR    = (maxx0 - minx0 ) * percentR;

% i.e. the input is Elimit = any_scalar
if (size(Elimit,2) == 1) 

	if(size(Gamma) ~= [1 1])	
		error('Variable Gamma here is expected to be a scalar') ;
	end

	E1                 = floor( (minx0 - heightL)/10 ) * 10;
	E2_scheme1         =  ceil( (maxx0 + heightR)/10 ) * 10;
	%E2_scheme2         = E1 + range
	%E2 = min(E2_scheme1, E2_scheme2)
    E2 = E2_scheme1;

	Npoints = (E2 - E1) / dE + 1;
    fprintf('The x-limit of data is [%20.10f %20.10f] \n', minx0, maxx0);
    fprintf('Broadening is done  in [%20.10f %20.10f] by default range setting scheme\n', E1, E2);
    
% i.e. the input is Elimit = [E1 E2]
elseif (size(Elimit,2) == 2)

    if(size(Gamma) ~= [1 1])
        error('Variable Gamma here is expected to be a scalar') ;
    end

	E1 = Elimit(1);
	E2 = Elimit(2);
	Npoints = (E2 - E1) / dE + 1;
    if (minx0 < E1)
        fprintf('The minumum of bardata %10.4f is lower than the specified E1=%10.4f\n', minx0, E1);
        fprintf('Please reconsider to set a smaller E1 to cover all data');
    else
        fprintf('The x-limit of data is [%20.10f %20.10f] \n', minx0, maxx0);
        fprintf('Input broadening range:[%20.10f %20.10f] \n', E1, E2);
        
    end

% i.e. the input is Elimit = [E1 E2  E3  E4];
elseif (size(Elimit,2) == 4)
	
	if(size(Gamma) ~= [1 2])
        error('Variable Gamma here is expected to be a row vector [Gamma1 Gamma2]') ;
    end

	E1 =  Elimit(1);	E2 =  Elimit(2);
	E3 =  Elimit(3);	E4 =  Elimit(4);
	Npoints = (E4 - E1) / dE + 1;

	Gamma1=Gamma(1); 
	Gamma2=Gamma(2);
else
	error(' Wrong dimension of variable ''Elimit''.')
end


convoluted = zeros(Npoints, ind_end);

for k = 2 : 1 : ind_end

	%. Convolution
	y0 = barspectra(:,k);

	if (convtype == 'g') 
		[x, y] = dogaussian(x0, y0, Gamma, E1, E2, Npoints);
	elseif (convtype == 'l')
		[x, y] = dolorentz (x0, y0, Gamma, E1, E2, Npoints);
	elseif (convtype == 'G')
		[x, y] = dogaussian_general (x0, y0, Gamma1, Gamma2, E1, E2, E3, E4, Npoints);
	else
		error ('Only convtype=g,l,G are possible (updated 2019-06-24)!')
	end

	if(k == 2)
		convoluted (:, 1) = x';
	end

	convoluted(:, k)= y';

end

%elapsed_time = etime (clock, globalt0);
%fprintf('Convolution done, elapsed_time = %12.1f (mins)\n',  elapsed_time/60);

print_convoluted(convoluted, File2);




%-----%. Plot
%-----figure
%-----plot(x, y,'-','LineWidth',2); 
%-----set(gca,'Fontsize',14);
%-----hold on;
%-----bar(x0, y0);
%-----set(gca, 'xlim', [E1, E2]);
%-----xlabel('Energy(eV)'); ylabel('Intensity(arb. unit)');
%-----titleline=['Raw data of ',pfxname,', spin=',spin,', HWHW = ', num2str(Gamma),' eV'];
%-----title(titleline);
%-----name_tif=['Raw',pfxname,'_nexafs_',spin,'.tif'];
%-----name_eps=['Raw',pfxname,'_nexafs_',spin,'.eps'];
%-----set(gcf,'PaperPositionMode','auto');
%-----saveas(gcf,name_tif,'tif');
%-----saveas(gcf,name_eps,'psc2');

