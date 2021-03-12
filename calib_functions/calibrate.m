function [ bcx, bcy ] = calibrate(dshmin,dshmax,x,y)
% Part of FAST Calibration v2.0 - GPLv3
% Moritz O. Ziegler, mziegler@gfz-potsdam.de
% DOI:      http://doi.org/10.5880/wsm.2021.002
% Manual:  	http://doi.org/10.48440/wsm.2021.002
% Download:	http://github.com/MorZieg/FAST_Calibration
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is the heart of the FAST calibration tool and computes the
% boundary condition of the best-fit model.
%
% dshmin: Difference between the modelled calibration stress scenarios and
% the actual stress data records for Shmin.
% dshmax: Difference between the modelled calibration stress scenarios and
% the actual stress data records for SHmax.
% x: Displacement in x’ direction prescribed at different test scenarios.
% y: Displacement in y’ direction prescribed at different test scenarios.
%

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
    if test(1) == test(2) == test(3)
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

bcx = ( n2_i * d_x - n2_x * d_i ) / ( n1_x * n2_i - n2_x * n1_i );
bcy = (d_x - n1_x * bcx) / n2_x;

end