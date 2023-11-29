function dshmax = critical_shmax(stress_data,calib_data,x,y,bcx,bcy)
% Part of FAST Calibration v2.4 - GPLv3
% Moritz O. Ziegler, mziegler@gfz-potsdam.de
% DOI:      http://doi.org/10.5880/wsm.2023.002
% Manual:  	http://doi.org/10.48440/wsm.2023.002
% Download:	http://github.com/MorZieg/FAST_Calibration
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SHmax magnitude is calibrated on the assumption of a fault that
% is about to fail. Therefore, for a given fault orientation and properties
% the amgnitude of SHmax is computed that is sufficient to lead to a
% failure.
%
% stress_data: A nx2 cell-variable in which the data records for
% calibration are stored sorted according to types of stress indicators.
% calib_data: A 1xn cell-variable in which the calibration scenarios stress
% states are stored for each location of a calibration point.
% x: Displacement in x' direction prescribed at different test scenarios.
% y: Displacement in y' direction prescribed at different test scenarios.
% bcx: Best-fit boundary conditions in x' direction.
% bcy: Best-fit boundary conditions in y' direction.
%

% Determine position of critical shmax in stress_data
for i = 1:size(stress_data,1)
    if contains(stress_data{i,1},"critical_shmax")
        pos = i;
    end
end

tensor = [];

% Perform for each critical shmax location
for i = 1:length(stress_data{pos,2}(:,1))
    strike = round(stress_data{pos,2}(i,4) - stress_data{pos,2}(i,8)):1:round(stress_data{pos,2}(i,4) + stress_data{pos,2}(i,8));
    strike_unaltered = strike;
    dip = round(stress_data{pos,2}(i,5) - stress_data{pos,2}(i,9)):1:round(stress_data{pos,2}(i,5) + stress_data{pos,2}(i,9));
    C = round(stress_data{pos,2}(i,6) - stress_data{pos,2}(i,10)):1:round(stress_data{pos,2}(i,6) + stress_data{pos,2}(i,10));
    mu = round(stress_data{pos,2}(i,7) - stress_data{pos,2}(i,11),1):0.1:round(stress_data{pos,2}(i,7) + stress_data{pos,2}(i,11),1);
    
    % Compute Shmin and Sv
    Shmin = solve(calib_data{pos}{i}(:,2),x,y,bcx,bcy);
    Sv = calib_data{pos}{i}(1,3);
    
    prel_tensor = [];
    
    % Loop over the strike variability
    for s = 1:length(strike)
        % Adjust strike angle
        strike(s) = strike(s) - mean(calib_data{pos}{i}(:,4)) + 90;
        if strike(s) < 0
            strike(s) = 360 + strike(s);
        elseif strike(s) > 360
            strike(s) = strike(s) - 360;
        end
    
        %Loop over the dip
        for d = 1:length(dip)
            % Compute normal of the fault plane
            n(1) = sind(dip(d)).*cosd(strike(s)); % Y is north!
            n(2) = -sind(dip(d)).*sind(strike(s)); % X is east!
            n(3) = -cosd(dip(d)); % vertical component

            % Compute Shmax for a slip tendency of 1
            A_temp = n(2)^2 * Shmin + n(3)^2 * Sv;
            B_temp = n(2)^2 * Shmin^2 + n(3)^2 * Sv^2;
    
            % Loop over the friction coefficient
            for m = 1:length(mu)
                % Loop over the cohesion
                for co = 1:length(C)
                    a = n(1)^4 + mu(m)^2 * n(1)^4 - n(1)^2;
                    b = 2 * n(1)^2 * A_temp + 2 * mu(m)^2 * n(1)^2 * A_temp + 2 * mu(m) * C(co) * n(1)^2;
                    c = A_temp^2 + mu(m)^2 * A_temp^2 + 2 * mu(m) * C(co) * A_temp + C(co)^2 - B_temp;
    
                    Shmax_temp(1) = ( -b + sqrt( b^2 - 4 * a * c)) / ( 2 * a );
                    Shmax_temp(2) = ( -b - sqrt( b^2 - 4 * a * c)) / ( 2 * a );
    
                    % Find S3
                    prel_tensor(end+1,1) = max(Shmax_temp);
                    prel_tensor(end,2) = strike_unaltered(s);
                    prel_tensor(end,3) = dip(d);
                    prel_tensor(end,4) = C(co);
                    prel_tensor(end,5) = mu(m);
                end
            end
        end
    end
    
    % Set S3 depending on tectonic stress regime it is Sv or Shmin.
    if Shmin < Sv
        s3 = Shmin;
    else
        s3 = Sv;
    end
    
    % Remove results in which Shmax is smaller than S3.
    t = 0;
    for j = 1:length(prel_tensor(:,1))
        if prel_tensor(j,1) < s3
            t = [ t j];
        end
    end
    t(1) = [];
    prel_tensor(t,:) = [];
    
    % Find the smallest SHmax which leads to failure.
    if length(prel_tensor(:,1)) > 0
        [~,minshmax] = min(prel_tensor(:,1));
        tensor(end+1,:) = prel_tensor(minshmax,:);
    end

end

% Find the fault with the smallest SHmax required for failure.
[~,i] = min(tensor(:,1));
stress{1} = "shmax";
stress{2}(:,4) = tensor(i,1);
stress{2}(:,5) = 1;
calib{1}{1} = calib_data{pos}{i};

disp(['Orientation #',num2str(i),' (',num2str(tensor(i,2)),'|',num2str(tensor(i,3)),') with a cohesion of ',num2str(tensor(i,4)),' MPa and a friction coefficient of ',num2str(tensor(i,5)),' fails for SHmax = ',num2str(tensor(i,1)),' MPa.'])
disp([' '])
dshmax = model_deviation(stress,calib,"shmax");

end
