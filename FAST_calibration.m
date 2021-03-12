%% FAST Calibration v2.0
% Distributed under licence GPLv3
% Moritz O. Ziegler, mziegler@gfz-potsdam.de
% DOI:      http://doi.org/10.5880/wsm.2021.002
% Manual:  	http://doi.org/10.48440/wsm.2021.002
% Download:	http://github.com/MorZieg/FAST_Calibration
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script helps to perform a calibration of a geomechanical-numerical
% model on stress data records.
% 
% folder: Directions to the folder in which the scripts are located. 
% name: Name for the macro and datafile.
% x and y: Displacement boundary conditions for each of the test
% model stress states. (Three sets of boundary conditions are required.)
% stress_data: A nx2 cell-variable in which the data records for
% calibration are stored sorted according to types of stress information.

clear all
close all
addpath('calib_functions')

% User input
folder = 'F:\Users\mziegler\Documents\Software\FAST\Scripting';
name = 'root';
x = [ -5 -5 -3 ];
y = [ 2 4 2 ];

stress_data  = {
%%%% Data used for calibration
% SHmax
%    X     Y      Z    Value   Confidence
"shmax"...
[ 3000  3000    -800    17.16   1;
  5000  6500    -1800   36.0 0.6];

% Shmin
%    X     Y      Z    Value   Confidence
"shmin"...
[ 3000  3000    -800    11.6   1;
  5000  6500    -3000   53.0 0.9];

% k-ratio used for SHmax calibration
%    X     Y      Z    Value   Confidence
"k_shmax"...
[ 5000  5000    -4000   1.0 0.3];
  
% Values defining faults that should be used for critical shmax
%  X    Y     Z  strike  dip  cohesion  friction
%"critical_shmax_detail"...
%[ 5000 4000 -3000 310 30 10 0.6 10 10 5 0.2];

%%%% Indicators and limits
% FITs
%    X     Y      Z    Value   Confidence
"fit_detail"...
[ 3000  3000    -400    5   1;
  3000  3000    -700    14   1;
  5000  6500    -1000   12 1;
  5000  6500    -3000   45.0 0.1];
  
% Borehole Breakouts
%    X     Y      Z    Value   Confidence
"bbo_detail"...
[ 3000  3000    -750    40   1;
  3000  3000    -430    20   1];
  
% Intact borehole section (no breakouts)
%    X     Y      Z    Value   Confidence
"nbo_detail"...
[ 3000  3000    -750    40   1;
  3000  3000    -430    20   1];
  
% Drilling induced fractures
%    X     Y      Z    Value   Confidence
"dif_detail"...
[ 5000  6500    -800   24 1;
  5000  6500    -1500   30 1];

% Intact borehole section (no drilling induced fractures)
%    X     Y      Z    Value   Confidence
"nif_detail"...
[ 5000  6500    -800   24 1;
  5000  6500    -1500   30 1];
  
% Regime Stress Ratio
%    X     Y      Z    Value   Confidence
"k-ratio_detail"...
[ 5000  5000    -4000   0.9 1;
  5000  5000    -3000   1.0 1];

% Regime Stress Ratio
%    X     Y      Z    Value   Confidence
"rsr_detail"...
[ 8500  9200    -4500   0.9 0.8;
  8400  8900    -4800   1.5 0.8];
};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check user input
if sum(contains([stress_data{:,1}],'shmin')) ==  0
    disp(['ATTENTION! Calibration point(s) for Shmin missing.'])
    disp([' '])
    return
end

if (sum(contains([stress_data{:,1}],'shmax')) + sum(contains([stress_data{:,1}],'critical_shmax'))) == 0
    disp(['ATTENTION! Calibration point(s) for SHmax missing.'])
    disp([' '])
    return
end

% Write Tecplot macro
mod = length(x);
disp(write_macro(stress_data,name,mod,folder))
disp(['Use macro ' name '.mcr in Tecplot on the calibration scenarios'])
disp([' '])

%% After running the macro continue here:

% Load the modelled stress state from the calibration scenarios.
cd data/
shmax_calib = csvread(strcat(name,'_shmax.csv'));
shmax_calib = shmax_calib(:,1);
shmin_calib = csvread(strcat(name,'_shmin.csv'));
shmin_calib = shmin_calib(:,1);
sv_calib = csvread(strcat(name,'_sv.csv'));
sv_calib = sv_calib(:,1);
shazi_calib = csvread(strcat(name,'_shazi.csv'));
shazi_calib = shazi_calib(:,1);
cd ../

% Assign the modelled stress states to the calibration locations.
t = 1;
for i = 1:size(stress_data,1)
    for j = 1:length(stress_data{i,2}(:,1))
        calib_data{i}{j}(:,1) = shmax_calib(t:(t+(mod-1)));
        calib_data{i}{j}(:,2) = shmin_calib(t:(t+(mod-1)));
        calib_data{i}{j}(:,3) = sv_calib(t:(t+(mod-1)));
        calib_data{i}{j}(:,4) = shazi_calib(t:(t+(mod-1)));
        t = t + mod;
    end
end

% Compute the differences between modelled and observed stress state.
dshmin = model_deviation(stress_data,calib_data,"shmin");

% Check if direct SHmax measurements are available
if sum(contains([stress_data{:,1}],'shmax')) == 0
    dshmax = zeros(length(x),1);
else
    dshmax = model_deviation(stress_data,calib_data,"shmax");
end

disp(['RESULTS:'])
% Check for other SHmax indicators or if critical shmax is assumed
if sum(contains([stress_data{:,1}],'critical_shmax')) > 0
    [bcx, bcy] = calibrate(dshmin,dshmax,x,y);
    dshmax = critical_shmax(stress_data,calib_data,x,y,bcx,bcy);
elseif sum(contains([stress_data{:,1}],'k_shmax')) > 0
    [bcx, bcy] = calibrate(dshmin,dshmax,x,y);
    dshmax = kratio(stress_data,calib_data,x,y,bcx,bcy);
end

% Provide the final best-fit boundary conditions
[bcx, bcy] = calibrate(dshmin,dshmax,x,y);

% Evaluate the indicators and limits according to the best-fit boundary
% conditions
if sum(contains([stress_data{:,1}],'fit')) > 0 || sum(contains([stress_data{:,1}],'rsr')) > 0 || sum(contains([stress_data{:,1}],'k')) > 0 ||sum(contains([stress_data{:,1}],'bbo')) > 0 ||sum(contains([stress_data{:,1}],'nbo')) > 0 ||sum(contains([stress_data{:,1}],'dif')) > 0 ||sum(contains([stress_data{:,1}],'nif')) > 0
     accuracy(stress_data,calib_data,x,y,bcx,bcy)
end
    
disp(['BOUNDARY CONDITIONS:'])
disp(['Root model boundary condition X: ',num2str(bcx)])
disp(['Root model boundary condition Y: ',num2str(bcy)])
disp([' '])