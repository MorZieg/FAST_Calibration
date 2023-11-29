function check(varargin)
% Part of FAST Calibration v2.4 - GPLv3
% Moritz O. Ziegler, mziegler@gfz-potsdam.de
% DOI:      http://doi.org/10.5880/wsm.2023.002
% Manual:  	http://doi.org/10.48440/wsm.2023.002
% Download:	http://github.com/MorZieg/FAST_Calibration
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function provides some sanity checks on the user input variables and
% the data read from files for the Calibration and Multistage approaches.

    if nargin == 1
        % Check whether data are available for SHmax and Shmin
        stress_data = varargin{1};
        if sum(contains([stress_data{:,1}],'shmin')) ==  0
            disp('ERROR! Not enough calibration data specified for Shmin.')
            pause;
        end
        
        if (sum(contains([stress_data{:,1}],'shmax')) + sum(contains([stress_data{:,1}],'critical_shmax'))) == 0
            disp('ERROR! Not enough calibration data specified for SHmax and/or Shmin.')
            pause;
        end
    
    elseif nargin == 2
        % Check whether the initial boundary consitions are linearly dependent.
        x = varargin{1};
        y = varargin{2};

        m = (y(2)-y(1))/(x(2)-x(1));
        n = y(1) - m * x(1);
        
        test = m * x(3) + n;
        
        if round(test) == y(3)
            disp('ERROR! The test boundary conditions are linearly dependent. Errors or wrong results are expected.')
            pause;
        end

    elseif nargin == 3
        % Check whether the number of data/synthetic calibration points is
        % in agreement with the provided data by the macro.
        model = varargin{1};
        data = varargin{2};

        test = model / varargin{3};

        if test ~= data
            disp('ERROR! Mismatch between loaded data and number of calibration points. Rerun the macro in Tecplot.')
            pause;
        end

    elseif nargin == 4
        % Check the age of the files in data and prompt if they are older
        % than 1 hour.
        file1 = strcat(varargin{1},varargin{2});
        file2 = strcat(varargin{3},varargin{4});

        FileInfo = dir(file1);
        t1 = FileInfo.date;

        FileInfo = dir(file2);
        t2 = FileInfo.date;

        out1 = diff(datenum([t1;datetime('now')]))*24;
        out2 = diff(datenum([t2;datetime('now')]))*24;

        if (out1 > 1) || (out2 > 1)
            disp('ATTENTION! Old files in /data. Tecplot macro has been run longer than 1 hour ago. Press any key to continue.')
            pause;
        end        

    end

end