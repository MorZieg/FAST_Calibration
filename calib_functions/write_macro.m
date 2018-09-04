function [ text ] = write_macro(shmax,shmin,name,mod,folder)
% Part of FAST Calibration v1.0 - GPLv3
% Moritz O. Ziegler, mziegler@gfz-potsdam.de
% Manual:  	http://doi.org/10.2312/wsm.2018.003
% Download:	http://github.com/MorZieg/FAST_Calibration
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write_macro creates a Tecplot Macro to extract variables at certain
% locations (nodes) from a calibrated regional model.
% 
% shmax: is an at least n x 3 Matrix with the coordinates of the SHmax
% datapoints.
% shmin: is the same for Shmin datapoints.
% name: is the name for the macro (e.g.: regional_name.mcr)
% folder: The full path to the current folder
%

% Macro Text
header = '#!MC 1410\n';
line1 = '$!VarSet |MFBD| = ''E:\\Program Files\\Tecplot\\Tecplot 360 EX 2015 R2''\n\n';

zone1 = '$!CREATERECTANGULARZONE\nIMAX = 1\nJMAX = 1\nKMAX = 1\n';
zone2 = 'X1 = %i\nY1 = %i\nZ1 = %i\n';
zone3 = 'X2 = %i\nY2 = %i\nZ2 = %i\n\n';

varnum = '$!GETVARNUMBYNAME |SHMAX|\nNAME = "SHmax"\n$!GETVARNUMBYNAME |SHMIN|\nNAME = "Shmin"\n\n';

interp1 = '$!LINEARINTERPOLATE\nSOURCEZONES =  [%i]\n';
interp2 = 'DESTINATIONZONE = %i\nVARLIST =  [|%s|]\nLINEARINTERPCONST = 0\nLINEARINTERPMODE = DONTCHANGE\n\n';

export1 = '$!EXTENDEDCOMMAND\nCOMMANDPROCESSORID = ''excsv''\n';
export2 = 'COMMAND = ''FrOp=1:ZnCount=%i:ZnList=[%i-%i]:VarCount=1:VarList=[|%s|]:ValSep=",":FNAME="%s\\data\\%s.csv"''\n\n';

del = '$!DELETEZONES [%i-%i]\n\n';

footer = '$!RemoveVar |MFBD|';

% Write Macro File

filename = strcat(name,'.mcr');
fid = fopen(filename,'w');

fprintf(fid,header);
fprintf(fid,line1);

fprintf(fid,varnum);

% SHmax
num = length(shmax(:,1)); % Number of datapoints
% Create Zones
for i = 1:num
    for j = 1:mod
        a = shmax(i,1);
        b = shmax(i,2);
        c = shmax(i,3);
   
        fprintf(fid,zone1);
        fprintf(fid,zone2,a,b,c);
        fprintf(fid,zone3,a,b,c);
    end
end

% Interpolate
for i = 1:num
    for j = 1:mod
        fprintf(fid,interp1,j);
        fprintf(fid,interp2,(mod+((i-1)*mod)+j),'SHMAX');
    end
end
    
fprintf(fid,export1);
fprintf(fid,export2,(num*mod),(mod+1),(mod+(num*mod)),'SHMAX',folder,strcat(name,'_shmax'));
    
fprintf(fid,del,mod+1,((num*mod)+mod));

% Shmin
num = length(shmin(:,1)); % Number of datapoints
% Create Zones
for i = 1:num
    for j = 1:mod
        a = shmin(i,1);
        b = shmin(i,2);
        c = shmin(i,3);
   
        fprintf(fid,zone1);
        fprintf(fid,zone2,a,b,c);
        fprintf(fid,zone3,a,b,c);
    end
end

% Interpolate
for i = 1:num
    for j = 1:mod
        fprintf(fid,interp1,j);
        fprintf(fid,interp2,(mod+((i-1)*mod)+j),'SHMIN');
    end
end
    
fprintf(fid,export1);
fprintf(fid,export2,(num*mod),(mod+1),(mod+(num*mod)),'SHMIN',folder,strcat(name,'_shmin'));

fprintf(fid,del,mod+1,((num*mod)+mod));


fprintf(fid,footer);

fclose(fid);


text = 'Macro File created';

end

