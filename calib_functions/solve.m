function [ stress ] = solve(stress_scenario,x,y,bcx,bcy)
% Part of FAST Calibration v2.0 - GPLv3
% Moritz O. Ziegler, mziegler@gfz-potsdam.de
% DOI:      http://doi.org/10.5880/wsm.2021.002
% Manual:  	http://doi.org/10.48440/wsm.2021.002
% Download:	http://github.com/MorZieg/FAST_Calibration
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function computes a horizontal stress component from the difference
% between the modelled and observed stress state and the best-fit boundary
% conditions.
%
% stress_scenario: calib_data for the stress component that is to be
% solved.
% x: Displacement in x’ direction prescribed at different test scenarios.
% y: Displacement in y’ direction prescribed at different test scenarios.
% bcx: Best-fit boundary conditions in x' direction.
% bcy: Best-fit boundary conditions in y' direction.
%
v = [x; y; stress_scenario' ]';

p = v(1,:);
u = v(2,:) - v(1,:);
v = v(3,:) - v(1,:);

test = u./v;
if test(1) == test(2) == test(3)
   fprintf('ERROR! Planes are linearly dependent.\n')
end

a = u(2) * v(3) - u(3) * v(2);
b = u(3) * v(1) - u(1) * v(3);
c = u(1) * v(2) - u(2) * v(1);
d = p(1) * a + p(2) * b + p(3) * c;

stress = ( d - a * bcx - b * bcy) / c;

end

