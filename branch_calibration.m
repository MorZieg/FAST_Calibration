% FAST Calibration v1.0 - GPLv3
% branch calibration
% Moritz O. Ziegler, mziegler@gfz-potsdam.de
% Manual:  	http://doi.org/10.2312/wsm.2018.003
% Download:	http://github.com/MorZieg/FAST_Calibration
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script helps to perform a calibration of a local model with a
% regional model. (Nested modelling)
% 
% branch_corner: The X and Y coordinates of the local model.
% root_nodes: The path and filename of the nodes of the local model. 
% type: type of calibration point distribution as well as the name for
% the macro and datafile. (random, border, or corner)
% num: number of calibration points.
% distrib: This is a factor for the distribution of calibration points. The
% larger the number the more evenly distributed are the calibration points.
% This variable only applies to type corner. It is ignored (but
% has to be present) for type random and border.
% minelem: The distance between the model border and the closest element
% (only for type border).
% x and y: These are the push and pull boundary conditions for each of the
% modelled stress states.
% Zmax: The lowest topographic elevation in the branch model.
% Zmin: The bottom of the branch model.
%
% NOTE: The function nodes2calibrationpoints.m is required.
% NOTE: The function write_macro.m is required.

clear all
close all
addpath('calib_functions')

% User input
folder = 'C:\Documents\Geom-Num-Model\Calibration';

branch_corner = [  2000 2000;
                   2000 6000;
                   6000 6000;
                   6000 2000];
                        
root_nodes = 'root_nodes.csv';
type = 'border';
num = 100;
distrib = 50;
minelem = 100;
Zmax = 0;
Zmin = -4000;

x = [ -3 -3 -1 ];
y = [ 1 3 1 ];

% Comment the following line if the calibration points are already
% available. For large models nespoints.m takes some time to run.
calib_points = nodes2calibrationpoints(branch_corner,root_nodes,num,distrib,type,minelem,Zmax,Zmin);

rootname = strcat('root_branch_calib_',type);
branchname = strcat('branch_calib_',type);

mod = length(x);
disp(write_macro(calib_points,calib_points,rootname,1,folder))
disp(write_macro(calib_points,calib_points,branchname,mod,folder))

disp(['Use macro ' rootname '.mcr on the model that provides the reference stress state.'])
disp(['Use macro ' branchname '.mcr on the smaller model with test scenarios that shall be calibrated.'])

%% After running the macro continue here:

cd data/
shmax = csvread(strcat(rootname,'_shmax.csv'));
shmax = shmax(:,1);
shmin = csvread(strcat(rootname,'_shmin.csv'));
shmin = shmin(:,1);

shmax_calib = csvread(strcat(branchname,'_shmax.csv'));
shmax_calib = shmax_calib(:,1);
shmin_calib = csvread(strcat(branchname,'_shmin.csv'));
shmin_calib = shmin_calib(:,1);
cd ../

local = [ shmax_calib(:,1) shmin_calib(:,1) ];
regional = [ shmax(:,1) shmin(:,1) ];

dshmin = zeros(num,mod);
dshmax = dshmin;
count = 1;

for i = 1:length(shmax)
    for k = 1:mod
        dshmin(i,k) = local(count,2) - regional(i,2);
        dshmax(i,k) = local(count,1) - regional(i,1);
        count = count + 1;
    end
end

dshmin = mean(dshmin,1)';
dshmax = mean(dshmax,1)';

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

disp(['Branch model boundary condition X: ',num2str(bcx)])
disp(['Branch model boundary condition Y: ',num2str(bcy)])


