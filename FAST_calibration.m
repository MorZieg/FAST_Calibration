%% FAST Calibration v2.4
% Distributed under licence GPLv3
% Moritz O. Ziegler, mziegler@gfz-potsdam.de
% DOI:      http://doi.org/10.5880/wsm.2023.002
% Manual:  	http://doi.org/10.48440/wsm.2023.002
% Download:	http://github.com/MorZieg/FAST_Calibration
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script helps to perform a calibration of a geomechanical-numerical
% model on provided stress data records and preliminarily assesses the
% plausibility of the modelled stress state using indirect data.
%
% For more information consult the manual.

clear all
close all
addpath('calib_functions')

% User input
folder = 'C:\Users\Documents\Software\FAST';
name = 'root';
x = [ -50 -50 -30 ];
y = [ 20 40 20 ];

stress_data  = {
%%%% Data used for calibration
% SHmax
%    X     Y      Z    Value   Confidence
"shmax"...
[ 700000  5530000    -800    17.16   1;
  680000  5540000    -1800   36.0  0.6];

% Shmin
%    X     Y      Z    Value   Confidence
"shmin"...
[ 700000  5530000    -800    11.6   1;
  680000  5540000    -3000   53.0 0.9];

% k-ratio used for SHmax calibration
%    X     Y      Z    Value   Confidence
"k_shmax"...
[ 700000  5530000    -4000   1.0 0.3];
  
% Values defining faults that should be used for critical shmax
%  X    Y     Z  strike  dip  cohesion  friction
%"critical_shmax_detail"...
%[ 700000 5530000 -6000 20 60 5 0.6 10 10 5 0.2];

%%%% Indicators and limits
% FITs
%    X     Y      Z    FIT-pressure   Confidence
"fit_detail"...
[ 700000  5530000    -400    5    1;
  700000  5530000    -700    14   1;
  680000  5540000    -1000   12   1;
  680000  5540000    -3000   45   0.1];
  
% Borehole Breakouts
%    X     Y      Z    Compressive_strength   Confidence
"bbo_detail"...
[ 700000  5530000    -750    40   1;
  700000  5530000    -430    50   1];
  
% Intact borehole section (no breakouts)
%    X     Y      Z    Compressive_strength   Confidence
"nbo_detail"...
[ 680000  5540000    -750    40   1;
  680000  5540000    -430    50   1];
  
% Drilling induced fractures
%    X     Y      Z    Tensile_strength   Confidence
"dif_detail"...
[ 700000 5530000    -800   24 1;
  700000 5530000    -1500   30 1];

% Intact borehole section (no drilling induced fractures)
%    X     Y      Z    Tensile_strength   Confidence
"nif_detail"...
[ 680000  5540000    -800   24 1;
  680000  5540000    -1500   30 1];
  
% K-ratio
%    X     Y      Z    k_ratio   Confidence
"k-ratio_detail"...
[ 700000 5530000    -4000   0.9 1;
  700000 5530000    -3000   1.0 1];

% Regime Stress Ratio
%    X     Y      Z    RSR   Confidence
"rsr_detail"...
[ 700000 5530000    -4500   0.9 0.8;
  700000 5530000    -4800   1.5 0.8];
};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check user input
check(stress_data)
check(x,y)

% Write Tecplot macro
mod = length(x);
disp(write_macro(stress_data,name,mod,folder))
disp(['Use macro ' name '.mcr in Tecplot on the calibration scenarios'])
disp(' ')

%% After running the macro continue here:

% Load the modelled stress state from the calibration scenarios.
cd data/

check(name,'_shmax.csv',name,'_sv.csv')

shmax_calib = csvread(strcat(name,'_shmax.csv'));
shmax_calib = shmax_calib(:,1);
shmin_calib = csvread(strcat(name,'_shmin.csv'));
shmin_calib = shmin_calib(:,1);
sv_calib = csvread(strcat(name,'_sv.csv'));
sv_calib = sv_calib(:,1);
shazi_calib = csvread(strcat(name,'_shazi.csv'));
shazi_calib = shazi_calib(:,1);
cd ../

% Sanity check
lsd = 0;
for i = 1:size(stress_data,1)
    lsd = lsd + length(stress_data{i,2}(:,1));
end
check(length(shmax_calib),lsd,3)

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

disp('RESULTS:')
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

disp(' ')
disp('BOUNDARY CONDITIONS:')
disp(['Root model boundary condition X: ',num2str(bcx)])
disp(['Root model boundary condition Y: ',num2str(bcy)])
disp(' ')