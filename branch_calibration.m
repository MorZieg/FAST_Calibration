%% Branch model calibration
% Part of FAST Calibration v2.0 - GPLv3
% Moritz O. Ziegler, mziegler@gfz-potsdam.de
% DOI:      http://doi.org/10.5880/wsm.2021.002
% Manual:  	http://doi.org/10.48440/wsm.2021.002
% Download:	http://github.com/MorZieg/FAST_Calibration
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script helps to perform a calibration of a local model with a
% regional model. (Nested modelling)
% 
% branch_corner: The X and Y coordinates of the local model.
% root_nodes: The path and filename of the nodes of the local model. 
% type: type of calibration point distribution as well as the name for
% the macro and datafile. (random, border, corner, or user)
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

clear all
close all
addpath('calib_functions')

% User input
folder = 'F:\Users\mziegler\Documents\Software\FAST\Scripting';

branch_corner = [  2000 2000;
                   2000 6000;
                   6000 6000;
                   6000 2000];
                        
root_nodes = 'root_nodes.csv';
type = 'user';
num = 10;
distrib = 50;
minelem = 100;
Zmax = 0;
Zmin = -4000;

x = [ -3 -3 -1 ];
y = [ 1 3 1 ];

% Comment the following line if the calibration points are already
% available. For large models nespoints.m takes some time to run.
calib_points = nodes2calibrationpoints(branch_corner,root_nodes,num,distrib,type,minelem,Zmax,Zmin);
calib_points = { 1 calib_points; 1 calib_points };

rootname = strcat('root_branch_calib_',type);
branchname = strcat('branch_calib_',type);

mod = length(x);

disp(write_macro(calib_points,rootname,1,folder))
disp(write_macro(calib_points,branchname,mod,folder))

disp(['Use macro ' rootname '.mcr on the root model.'])
disp(['Use macro ' branchname '.mcr on the branch model that shall be calibrated.'])

%% After running the macro continue here:
% Load the modelled stress state from the calibration scenarios.
cd data/
shmax_calib = csvread(strcat(branchname,'_shmax.csv'));
shmax_calib = shmax_calib(:,1);
shmin_calib = csvread(strcat(branchname,'_shmin.csv'));
shmin_calib = shmin_calib(:,1);

shmax_model = csvread(strcat(rootname,'_shmax.csv'));
shmax_model = shmax_model(:,1);
shmin_model = csvread(strcat(rootname,'_shmin.csv'));
shmin_model = shmin_model(:,1);
cd ../

% Assign the modelled stress states to the calibration locations.
t = 1;
for i = 1:length(calib_points)
    for j = 1:length(calib_points{i}(:,1))
        calib_data{i}{j}(:,1) = shmax_calib(t:(t+(mod-1)));
        calib_data{i}{j}(:,2) = shmin_calib(t:(t+(mod-1)));
        t = t + mod;
    end
end

t = 1;
calib_model_max{1} = "shmax";
calib_model_min{1} = "shmin";
for j = 1:length(calib_points{1}(:,1))
    calib_model_max{2}(j,4) = shmax_model(t);
    calib_model_max{2}(j,5) = 1;
    calib_model_min{2}(j,4) = shmin_model(t);
    calib_model_min{2}(j,5) = 1;
    t = t + 1;
end

% Compute the differences between root model and branch model stress state.
dshmin = model_deviation(calib_model_min,calib_data,"shmin");
dshmax = model_deviation(calib_model_max,calib_data,"shmax");

% Compute the calibrated boundary conditions.
[bcx, bcy] = calibrate(dshmin,dshmax,x,y);    

disp(['Branch model boundary condition X: ',num2str(bcx)])
disp(['Branch model boundary condition Y: ',num2str(bcy)])
disp([' '])

