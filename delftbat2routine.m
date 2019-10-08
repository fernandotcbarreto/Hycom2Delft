function [bat]=delftbat2routine(depth,mdf)
%convert .dep file to hycomtoriemann bat routine

dep=fopen(depth);
fid=fopen(mdf, 'r');
%dep=fopen('L3_maxV.dep', 'r');
%fid=fopen('L3_maxV.mdf', 'r')
counter = 0;                                % line counter (sloppy but works)
while 1                                     % infinite loop
    tline = fgetl(fid);                     % read a line
    counter = counter + 1;                  % we are one line further
    if ischar(tline)                        % if the line is string 
        U = strfind(tline, 'MNKmax'); % where the string start (if at all)
        if isfinite(U) == 1;                % if it is a number actually
            break                           % we found it, lets go home
        end
     end
end
sep=textscan(tline, '%s');
sep1=sep{1};
nlon=sep1{2};
nlon=str2double(nlon);
nlat=sep1{3};
nlat=str2double(nlat);
sumc=textscan(dep, '%f', nlon);
for i=2:nlat
    c{i}=textscan(dep, '%f', nlon);
%    c{i}=c{i}';
    sumc=[sumc,c{i}];
end
d=cell2mat(sumc);
bat=d';
bat=bat(1:end-1,1:end-1);
fclose(dep)
%frewind(fileID);
clc

    