function [ dshmax ] = shmax(stress_data,calib_data,x,y,bcx,bcy)
% Part of FAST Calibration v2.0 - GPLv3
% Moritz O. Ziegler, mziegler@gfz-potsdam.de
% DOI:      http://doi.org/10.5880/wsm.2021.002
% Manual:  	http://doi.org/10.48440/wsm.2021.002
% Download:	http://github.com/MorZieg/FAST_Calibration
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function provides SHmax calibration from either SHmax magnitude data or
% derives SHmax magnitudes from k-ratios.
%
% stress_data: A nx2 cell-variable in which the data records for
% calibration are stored sorted according to types of stress indicators.
% calib_data: A 1xn cell-variable in which the calibration scenarios stress
% states are stored for each location of a calibration point.
% x: Displacement in x’ direction prescribed at different test scenarios.
% y: Displacement in y’ direction prescribed at different test scenarios.
% bcx: Best-fit boundary conditions in x' direction.
% bcy: Best-fit boundary conditions in y' direction.
%

m = 1;

% k-value Shmax computation
for i = 1:size(stress_data,1)
    if contains(stress_data{i,1},"k_shmax")
        for j = 1:length(calib_data{i})
            shmin = solve(calib_data{i}{j}(:,2),x,y,bcx,bcy);
            shmax_temp = 2 * stress_data{i,2}(j,4) * calib_data{i}{j}(1,3) - shmin;
            if shmin > shmax_temp
                disp(['k-value #',num2str(j),' not reasonable'])
            else
                calib{1}{m} = calib_data{i}{j};
                shmaxv(m) = shmax_temp;
                weight(m) = stress_data{i,2}(j,5);
                m = m + 1;
            end
        end
    end
end

% Shmax values
for i = 1:size(stress_data,1)
    if contains(stress_data{i,1},"k_shmax")
        continue
    end
    if contains(stress_data{i,1},"shmax")
        stress_data{i,1};
        for j = 1:length(calib_data{i})
            shmaxv(m) = stress_data{i,2}(j,4);
            weight(m) = stress_data{i,2}(j,5);
            calib{1}{m} =  calib_data{i}{j};
            m = m + 1;
        end
    end
end

% Compute dshmax
stress{1} = "shmax";
stress{2}(:,4) = shmaxv;
stress{2}(:,5) = weight;

dshmax = model_deviation(stress,calib,"shmax");

end