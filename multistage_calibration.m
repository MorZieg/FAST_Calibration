%% FAST Calibration v2.4 - Multistage Modelling
% Part of FAST Calibration v2.4 - GPLv3
% Moritz O. Ziegler, mziegler@gfz-potsdam.de
% DOI:      http://doi.org/10.5880/wsm.2023.002
% Manual:  	http://doi.org/10.48440/wsm.2023.002
% Download:	http://github.com/MorZieg/FAST_Calibration
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script helps to perform a calibration of a geomechanical-numerical
% model on the stress state provided by a larger scale model in a so-called
% multistage approach.
%
% For more information consult the manual.

clear all
close all
addpath('calib_functions')

% User input
folder = 'C:\Users\Documents\FAST';

% Corner coordinates of branch model in according coordinate system
branch_corner = [   4450000   5520000;
                    4450000   5540000;
                    4475000   5540000;
                    4475000   5520000;];

root_nodes = 'root_nodes.csv';
%root_nodes = 'calib_points.mat';
distrib = 20;
type = 'border'; 
num = 40;
minelem = 200;
Zmax = 0;
Zmin = -8000;

% Displacement boundary conditions branch model
x = [ -10 -10 -5 ];
y = [ 10 5 10 ];

% Relation between material IDs in root and branch model.
root_litho =    [1 3 2];
branch_litho =  [1 5 2];

% If branch model and root model are in a different coordinate system:
p_root = projcrs(23032);
p_branch = projcrs(5684);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of user input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP I
% Check for the independence of test boundary conditions
disp('FAST Calibration v2.4 - Multistage approach')
disp('STEP I')
disp(' ')
check(x,y)

% Coordinate systems are adjusted.
if exist("p_branch")
    disp('Coordinate systems of root and branch model are adjusted.')

    [lat,lon] = projinv(p_branch,branch_corner(:,1),branch_corner(:,2));
    [east,north] = projfwd(p_root,lat,lon);
    branch_corner = [east north];
end

% Create calibration points or read calibration points from file.
if isequal(root_nodes(end-2:end),'csv')
    disp('Generating new calibration points')
    disp('This may take a while for large models an/or many calibration points.')
    calib_points_root = nodes2calibrationpoints(branch_corner,root_nodes,num,distrib,type,minelem,Zmax,Zmin);
    time = string(datetime('now','Format','uuuuMMdd''T''HHmmss'));
    save(strcat("calib_points_initial_",time,".mat"),"calib_points_root")

elseif isequal(root_nodes(end-2:end),'mat')
    disp(strcat('Reading calibration points from file ',root_nodes))
    load(root_nodes)

else
    disp('ERROR! No root nodes file or calibration points specified.')
end

calib_points_root = { 1 calib_points_root; 1 calib_points_root };

% Calibration points are transferred from root model CRS to branch model
% CRS
if exist("p_branch")
    [lat,lon] = projinv(p_root,calib_points_root{2,2}(:,1),calib_points_root{2,2}(:,2));
    [east,north] = projfwd(p_branch,lat,lon);
    calib_points_branch = { 1 [east,north,calib_points_root{2,2}(:,3)]; 1 [east,north,calib_points_root{2,2}(:,3)] };
else
    calib_points_branch = calib_points_root;
end

rootname = strcat('root_matID_',type);
branchname = strcat('branch_matID_',type);

disp(write_macro_mID(calib_points_root,rootname,1,folder))
disp(write_macro_mID(calib_points_branch,branchname,3,folder))

disp(['Use macro ' rootname '.mcr on the root model.'])
disp(['Use macro ' branchname '.mcr on the branch model.'])
disp(' ')

%% STEP II
% Load the material IDs at potential calibrations for root and branch model.
disp('STEP II')
cd data/
check(branchname,'_matID.csv',rootname,'_matID.csv')

branch_calib = csvread(strcat(branchname,'_matID.csv'));
branch_calib = branch_calib(:,1);

root_model = csvread(strcat(rootname,'_matID.csv'));
root_model = root_model(:,1);
cd ../

check(length(branch_calib),length(calib_points_root{1,2}(:,1)),1)

% Assign the material IDs to the calibration locations.
for i = 1:length(calib_points_root{2,2}(:,1))
    mat_root(i) = root_model(i);
    mat_branch(i) = branch_calib(i);
end

% Compare lithologies:
k = 1;
for i = 1:length(mat_root)
    if ismember(mat_root(i),root_litho)
        r = find(mat_root(i) == root_litho);
        b = branch_litho(r);

        if ismember(mat_branch(i),b)
             final_calib_points(k,:) = calib_points_root{1,2}(i,:);
             k = k + 1;
        end
    end
end

scatter3(final_calib_points(:,1),final_calib_points(:,2),final_calib_points(:,3));
hold on
scatter3(branch_corner(:,1),branch_corner(:,2),[0,0,0,0],50,'black','LineWidth',1.5);
plot3([branch_corner(:,1);branch_corner(1,1)],[branch_corner(:,2);branch_corner(1,2)],[0,0,0,0,0],'k','LineWidth',1.5)
hold off

time = string(datetime('now','Format','uuuuMMdd''T''HHmmss'));
save(strcat("calib_points_final_",time,"_litho.mat"),"final_calib_points")

disp('Calibration points tested for same lithology.')
disp([num2str(size(final_calib_points,1)) ' calibration points remaining.'])
disp(' ')


%% STEP III
disp('STEP III')

if isequal(root_nodes(end-8:end),'litho.mat')
    load(root_nodes)
end

calib_points_root = { 1 final_calib_points; 1 final_calib_points };

rootname = strcat('root_',type);
branchname = strcat('branch_',type);

% If required: Coordinate transformation.
if exist("p_branch")
    [lat,lon] = projinv(p_root,calib_points_root{2,2}(:,1),calib_points_root{2,2}(:,2));
    [east,north] = projfwd(p_branch,lat,lon);
    calib_points_branch = { 1 [east,north,calib_points_root{2,2}(:,3)]; 1 [east,north,calib_points_root{2,2}(:,3)] };
else
    calib_points_branch = calib_points_root;
end

disp(write_macro(calib_points_root,rootname,1,folder))
disp(write_macro(calib_points_branch,branchname,3,folder))

disp(['Use macro ' rootname '.mcr on the model that provides the reference stress state.'])
disp(['Use macro ' branchname '.mcr on the smaller model with test scenarios that shall be calibrated.'])
disp(' ')


%% STEP IV
% Load the modelled stress state from the calibration scenarios.
disp('STEP IV')
cd data/

check(branchname,'_shmax.csv',rootname,'_shmax.csv')

shmax_calib = csvread(strcat(branchname,'_shmax.csv'));
shmax_calib = shmax_calib(:,1);
shmin_calib = csvread(strcat(branchname,'_shmin.csv'));
shmin_calib = shmin_calib(:,1);

shmax_model = csvread(strcat(rootname,'_shmax.csv'));
shmax_model = shmax_model(:,1);
shmin_model = csvread(strcat(rootname,'_shmin.csv'));
shmin_model = shmin_model(:,1);
cd ../

check(length(shmax_calib),length(calib_points_root{1,2}(:,1))*2,3)

% Assign the modelled stress states to the calibration locations.
for i = 1:2
    t = 1;
    for j = 1:length(calib_points_root{1,2}(:,1))
        calib_data{i}{j}(:,1) = shmax_calib(t:(t+(3-1)));
        calib_data{i}{j}(:,2) = shmin_calib(t:(t+(3-1)));
        t = t + 3;
    end
end

t = 1;
calib_model_max{1} = "shmax";
calib_model_min{1} = "shmin";
for j = 1:length(calib_points_root{1,2}(:,1))
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
disp(' ')
