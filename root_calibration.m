% FAST Calibration v1.0 - GPLv3
% root calibration
% Moritz O. Ziegler, mziegler@gfz-potsdam.de
% Manual:  	http://doi.org/10.2312/wsm.2018.003
% Download:	http://github.com/MorZieg/FAST_Calibration
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script helps to perform a calibration of a regional model with 
% stress data records.
% 
% folder: Directions to the folder in which the scripts are located. 
% shmin and shmax: The coordinates (X,Y,Z), magnitude, and confidence
% (0 to 1) of stress data used for the calibration.
% name: Name for the macro and datafile.
% x and y: These are the push and pull boundary conditions for each of the
% modelled stress states. (Three sets of boundary conditions are required.)
%
% NOTE: The function write_macro.m is required.

clear all
close all
addpath('calib_functions')

% User input
folder = 'C:\Documents\Geom-Num-Model\Calibration';

%           X       Y      Z    Value   Confidence
shmax = [ 3000  3000    -800    17.16   1;
          5000  6500    -1800   36.0 0.6];

%           X       Y      Z    Value   Confidence
shmin = [ 3000  3000    -800    11.6   1;
          5000  6500    -3000   53.0 0.9];

name = 'root';

x = [ -5 -5 -3 ];
y = [ 2 4 2 ];

mod = length(x);
disp(write_macro(shmax,shmin,name,mod,folder))
disp(['Use macro ' name '.mcr in Tecplot'])

%% After running the macro continue here:

cd data/
shmax_calib = csvread(strcat(name,'_shmax.csv'));
shmax_calib = shmax_calib(:,1);
shmin_calib = csvread(strcat(name,'_shmin.csv'));
shmin_calib = shmin_calib(:,1);
cd ../

for i = 1:size(shmax,1)
    dshmax(:,i) = shmax_calib((((i-1)*mod)+1):(mod*i)) - shmax(i,4);
end

for i = 1:size(shmin,1)
    dshmin(:,i) = shmin_calib((((i-1)*mod)+1):(mod*i)) - shmin(i,4);
end

dshmin = bsxfun(@times,dshmin,shmin(:,5)')/sum(shmin(:,5)',1);
dshmax = bsxfun(@times,dshmax,shmax(:,5)')/sum(shmax(:,5)',1);

for i = 1:2
    if i == 1
        stress = dshmax;
    else
        stress = dshmin;
    end

    v = [x; y; stress' ]';
    r1 = v(2,:) - v(1,:);
    r2 = v(3,:) - v(1,:);

    test = r1./r2;
    if (test(1) == test(2)) && (test(2) == test(3))
        fprintf('ERROR! Planes are linearly dependent, i=%i\n',i);
    end

    n = cross(r1,r2);
    d = (n*v(1,:)');
    if i == 1
        d_x = d;
        n1_x = n(1);
        n2_x = n(2);
    else
        d_i = d;
        n1_i = n(1);
        n2_i = n(2);
    end
end

if n2_x == 0
    n2_x = 0.0001;
end

bcx = ( n2_i * d_x - n2_x * d_i ) / ( n1_x * n2_i - n2_x * n1_i );
bcy = (d_x - n1_x * bcx) / n2_x;

disp(['Root model boundary condition X: ',num2str(bcx)])
disp(['Root model boundary condition Y: ',num2str(bcy)])
