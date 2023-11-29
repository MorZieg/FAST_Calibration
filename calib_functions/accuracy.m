function [ ] = accuracy(stress_data,calib_data,x,y,bcx,bcy)
% Part of FAST Calibration v2.4 - GPLv3
% Moritz O. Ziegler, mziegler@gfz-potsdam.de
% DOI:      http://doi.org/10.5880/wsm.2023.002
% Manual:  	http://doi.org/10.48440/wsm.2023.002
% Download:	http://github.com/MorZieg/FAST_Calibration
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function provides estimates on the validity of the calibrated stress
% field according to expectations of formation integrity tests, k-ratios
% and regime stress ratios (RSR) at different locations.
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

disp(' ')
disp('Assessment of indirect data')
disp(' ')
disp(' ')

for i = 1:size(stress_data,1)
    % FITs
    if contains(stress_data{i,1},"fit")
       failed = 0;
       failed_sv = 0;
       failed_shmin = 0;
       for j = 1:size(stress_data{i,2},1)
           dshlb = calib_data{i}{j}(:,2);
           dsvlb = calib_data{i}{j}(:,3);
           
           mod_shmin = solve(dshlb,x,y,bcx,bcy);
           shmin_lb_mod(j) = mod_shmin;
           
           if mod_shmin < stress_data{i,2}(j,4)
                failed = failed + stress_data{i,2}(j,5);
                failed_shmin = failed_shmin + 1;
                if contains(stress_data{i,1},'detail')
                    diff = stress_data{i,2}(j,4) - mod_shmin;
                    disp(['FIT #',num2str(j),' is ',num2str(diff),' MPa larger than Shmin.'])
                end
           end
           if dsvlb(1) < stress_data{i,2}(j,4)
                failed = failed + stress_data{i,2}(j,5);
                failed_sv = failed_sv + 1;
                if contains(stress_data{i,1},'detail')
                    diff = stress_data{i,2}(j,4) - dsvlb;
                    disp(['FIT #',num2str(j),' is ',num2str(diff(1)),' MPa larger than Sv.'])
                end
           end
       end
       failed_percent = failed / sum(stress_data{i,2}(:,5));
       
       disp(['Summary:'])
       disp([num2str(round(failed_percent * 100)),'% of the weighted lower boundary constraints (FITs) failed.'])
       disp([num2str(failed_shmin),' of ',num2str(size(stress_data{i,2},1)),' lower boundary constraints failed for Shmin.'])
       disp([num2str(failed_sv),' of ',num2str(size(stress_data{i,2},1)),' lower boundary constraints failed for Sv.'])
       disp([' '])
    end
    
    % RSR
    if contains(stress_data{i,1},"rsr")
        for j = 1:size(stress_data{i,2},1)
            dshminlb = calib_data{i}{j}(:,2);
            dshmaxlb = calib_data{i}{j}(:,1);
            
            shmin = solve(dshminlb,x,y,bcx,bcy);
            shmax = solve(dshmaxlb,x,y,bcx,bcy);
            sv = calib_data{i}{j}(1,3);
            
            if shmin < sv < shmax
                n = 1;
            elseif sv < shmin < shmax
                n = 2;
            else
                n = 0;
            end
            
            s = sort([shmin; shmax; sv],'descend');
            sr = (s(2) - s(3))/(s(1) - s(3));
            rsr = (n + 0.5) + (-1)^n * (sr - 0.5);
            rsr_diff(j) = rsr - stress_data{i,2}(j,4);
            
            if contains(stress_data{i,1},'detail')
                disp(['Position #',num2str(j),': Difference in RSR: ',num2str(rsr_diff(j)),'.'])
            end  
            
        end
        rsr_mean = mean(rsr_diff);
        rsr_std = std(rsr_diff);
        
        disp(['Summary: Mean deviation of RSR is ',num2str(rsr_mean),' Standard Deviation ',num2str(rsr_std),'.']);
        disp([' '])
    end
    
    % k-value
    if contains(stress_data{i,1},"k-ratio")
        for j = 1:size(stress_data{i,2},1)
            dshminlb = calib_data{i}{j}(:,2);
            dshmaxlb = calib_data{i}{j}(:,1);
            
            shmin = solve(dshminlb,x,y,bcx,bcy);
            shmax = solve(dshmaxlb,x,y,bcx,bcy);
            sv = calib_data{i}{j}(1,3);
            
            k = (shmax + shmin)/(2*sv);
            k_diff(j) = k - stress_data{i,2}(j,4);
            
            if contains(stress_data{i,1},'detail')
                disp(['Position #',num2str(j),': modelled k-ratio is ',num2str(k),'. Difference in k-ratio: ',num2str(k_diff(j)),'.'])
            end            
        end
        k_mean = mean(k_diff);
        k_std = std(k_diff);
        
        disp(['Summary: Mean deviation of k is ',num2str(k_mean),' Standard Deviation ',num2str(k_std),'.']);
        disp([' '])
    end
    
    % Observed borehole breakouts
    if contains(stress_data{i,1},'bbo')
        bo = 0;
        for j = 1:size(stress_data{i,2},1)
            dshminlb = calib_data{i}{j}(:,2);
            dshmaxlb = calib_data{i}{j}(:,1);
            
            shmin = solve(dshminlb,x,y,bcx,bcy);
            shmax = solve(dshmaxlb,x,y,bcx,bcy);
            
            sigma_circ_max = 3 * shmax - shmin;
            if sigma_circ_max > stress_data{i,2}(j,4)
                if contains(stress_data{i,1},'detail')
                    disp(['Breakout at observed borehole position #',num2str(j),'.'])
                end
                bo = bo + 1;
            else
                if contains(stress_data{i,1},'detail')
                    disp(['No breakout at observed borehole position #',num2str(j),'.'])
                end
            end
        end
        disp(['Summary: ',num2str(bo),' breakouts out of ',num2str(j),' observed locations.'])
        disp([' '])
    end
    
    % Intact borehole section (no breakouts)
    if contains(stress_data{i,1},'nbo')
        nobo = 0;
        for j = 1:size(stress_data{i,2},1)
            dshminlb = calib_data{i}{j}(:,2);
            dshmaxlb = calib_data{i}{j}(:,1);
            
            shmin = solve(dshminlb,x,y,bcx,bcy);
            shmax = solve(dshmaxlb,x,y,bcx,bcy);
            
            sigma_circ_max = 3 * shmax - shmin;
            if sigma_circ_max > stress_data{i,2}(j,4)
                if contains(stress_data{i,1},'detail')
                    disp(['Breakout at intact borehole section #',num2str(j),'.'])
                end
            else
                if contains(stress_data{i,1},'detail')
                    disp(['No breakout at intact borehole section #',num2str(j),'.'])
                end
                nobo = nobo + 1;
            end
        end
        disp(['Summary: ',num2str(nobo),' intact borehole walls out of ',num2str(j),' observed intact locations.'])
        disp([' '])
    end
    
    % Drilling induced fracture
    if contains(stress_data{i,1},'dif')
        dif = 0;
        for j = 1:size(stress_data{i,2},1)
            dshminlb = calib_data{i}{j}(:,2);
            dshmaxlb = calib_data{i}{j}(:,1);
            
            shmin = solve(dshminlb,x,y,bcx,bcy);
            shmax = solve(dshmaxlb,x,y,bcx,bcy);
            
            sigma_circ_min = 3 * shmin - shmax;
            if sigma_circ_min < stress_data{i,2}(j,4)
                if contains(stress_data{i,1},'detail')
                    disp(['Drilling induced fracture at observed position #',num2str(j),'.'])
                end
                dif = dif + 1;
            else
                if contains(stress_data{i,1},'detail')
                    disp(['No drilling induced fracture at observed position #',num2str(j),'.'])
                end
            end
        end
        disp(['Summary: ',num2str(dif),' drilling induced fractures out of ',num2str(j),' possible locations.'])
        disp([' '])
    end
    
    % Intact borehole wall (no Drilling induced fracture)
    if contains(stress_data{i,1},'nif')
        nodif = 0;
        for j = 1:size(stress_data{i,2},1)
            dshminlb = calib_data{i}{j}(:,2);
            dshmaxlb = calib_data{i}{j}(:,1);
            
            shmin = solve(dshminlb,x,y,bcx,bcy);
            shmax = solve(dshmaxlb,x,y,bcx,bcy);
            
            sigma_circ_min = 3 * shmin - shmax;
            if sigma_circ_min < stress_data{i,2}(j,4)
                if contains(stress_data{i,1},'detail')
                    disp(['Drilling induced fracture at observed intact section #',num2str(j),'.'])
                end
            else
                if contains(stress_data{i,1},'detail')
                    disp(['No drilling induced fracture at observed intact section #',num2str(j),'.'])
                end
                nodif = nodif + 1;
            end
        end
        disp(['Summary: ',num2str(nodif),' intact borehole walls out of ',num2str(j),' observed intact locations.'])
        disp([' '])
    end
end
disp(' ')

end