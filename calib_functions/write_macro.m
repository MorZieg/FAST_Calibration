function [ text ] = write_macro(stress,name,mod,folder)
% Part of FAST Calibration v2.4 - GPLv3
% Moritz O. Ziegler, mziegler@gfz-potsdam.de
% DOI:      http://doi.org/10.5880/wsm.2023.002
% Manual:  	http://doi.org/10.48440/wsm.2023.002
% Download:	http://github.com/MorZieg/FAST_Calibration
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write_macro creates a Tecplot 360 EX Macro to extract variables at
% certain locations (nodes) from a calibrated regional model.
% 
% stress: A nx1 cell-variable in which the data records for calibration and
% their coordinates are stored sorted according to types of stress
% indicators.
% name: is the name for the macro (e.g.: regional_name.mcr)
% mod: The number of test scenarios
% folder: The full path to the current folder
%

coords = [];
for i = 1:size(stress,1)
    coords = [coords; stress{i,2}(:,1:3)];
end
num = length(coords(:,1));

% Macro Text
header = '#!MC 1410\n';
line1 = '$!VarSet |MFBD| = ''E:\\Program Files\\Tecplot\\Tecplot 360 EX 2015 R2''\n\n';

zonecheck = '$!IF |NUMZONES| != %i\n$!PAUSE "The expected number of %i zones is not encountered. The results can be wrong or an error may occur."\n$!ENDIF\n\n';

zone1 = '$!CREATERECTANGULARZONE\nIMAX = 1\nJMAX = 1\nKMAX = 1\n';
zone2 = 'X1 = %i\nY1 = %i\nZ1 = %i\n';
zone3 = 'X2 = %i\nY2 = %i\nZ2 = %i\n\n';

varnum = '$!GETVARNUMBYNAME |SHMAX|\nNAME = "SHmax"\n$!GETVARNUMBYNAME |SHMIN|\nNAME = "Shmin"\n$!GETVARNUMBYNAME |SV|\nNAME = "Sv"\n$!GETVARNUMBYNAME |SHAZI|\nNAME = "SHazi"\n\n';

interp1 = '$!LINEARINTERPOLATE\nSOURCEZONES =  [%i]\n';
interp2 = 'DESTINATIONZONE = %i\nVARLIST =  [|SHMAX|,|SHMIN|,|SV|,|SHAZI|]\nLINEARINTERPCONST = 0\nLINEARINTERPMODE = DONTCHANGE\n\n';

export1 = '$!EXTENDEDCOMMAND\nCOMMANDPROCESSORID = ''excsv''\n';
export2 = 'COMMAND = ''FrOp=1:ZnCount=%i:ZnList=[%i-%i]:VarCount=1:VarList=[|%s|]:ValSep=",":FNAME="%s\\data\\%s.csv"''\n\n';

del = '$!DELETEZONES [%i-%i]\n\n';

footer = '$!RemoveVar |MFBD|';

% Write Macro File
filename = strcat(name,'.mcr');
fid = fopen(filename,'w');

fprintf(fid,header);
fprintf(fid,line1);

fprintf(fid,zonecheck,mod,mod);

fprintf(fid,varnum);

% Create Zones
for i = 1:num
    for j = 1:mod
        a = coords(i,1);
        b = coords(i,2);
        c = coords(i,3);

        fprintf(fid,zone1);
        fprintf(fid,zone2,a,b,c);
        fprintf(fid,zone3,a,b,c);
    end
end

% Interpolate
for i = 1:num
    for j = 1:mod
        fprintf(fid,interp1,j);
        fprintf(fid,interp2,(mod+((i-1)*mod)+j));
    end
end

% Export
for i = 1:4
    if i == 1
        inst = 'SHMAX';
        id = '_shmax';
    elseif i == 2
        inst = 'SHMIN';
        id = '_shmin';
    elseif i == 3
        inst = 'SV';
        id = '_sv';
    elseif i == 4
        inst = 'SHAZI';
        id = '_shazi';
    end
    
    fprintf(fid,export1);
    fprintf(fid,export2,(num*mod),(mod+1),(mod+(num*mod)),inst,folder,strcat(name,id));
end

fprintf(fid,del,mod+1,((num*mod)+mod));
fprintf(fid,footer);

fclose(fid);
text = 'Macro File created';
end