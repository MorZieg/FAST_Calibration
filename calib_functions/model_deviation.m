function ds_out = model_deviation(stress_data,calib_data,type)
% Part of FAST Calibration v2.0 - GPLv3
% Moritz O. Ziegler, mziegler@gfz-potsdam.de
% DOI:      http://doi.org/10.5880/wsm.2021.002
% Manual:  	http://doi.org/10.48440/wsm.2021.002
% Download:	http://github.com/MorZieg/FAST_Calibration
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model_deviation is used to compare modelled calibration stress scenarios
% to the actual stress data records.
% 
% stress_data: A nx2 cell-variable in which the data records for
% calibration are stored sorted according to types of stress indicators.
% calib_data: A 1xn cell-variable in which the calibration scenarios stress
% states are stored for each location of a calibration point.
% type: The type of indicator considered. Either "shmax" or "shmin".
%

if contains(type,"shmax")
    t = 1;
elseif contains(type,"shmin")
    t = 2;
end

for j = 1:size(stress_data,1)
    if contains(stress_data{j,1},type)
        for i = 1:size(stress_data{j,2},1)
            ds_out(:,i) = calib_data{j}{i}(:,t) - stress_data{j,2}(i,4);
        end
        ds_out = bsxfun(@times,ds_out,stress_data{j,2}(:,5)')/sum(stress_data{j,2}(:,5)',1);
    end
end

end