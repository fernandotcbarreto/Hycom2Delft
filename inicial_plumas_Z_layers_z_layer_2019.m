clear all
close all
clc

% este programa faz as condições iniciais para o delft3d
dirs='C:\Users\fernando\Desktop\back_up\delft3d_FM_method_routines\d-data_backup';

cd (dirs)


fid = fopen('hycom.ini','wt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HYCOM
fileh='2016-01-01.nc';

lath=ncread(fileh,'Latitude');
lonh=ncread(fileh,'Longitude')-360;
lat=permute(lath, [2,1]);    
lon=permute(lonh, [2,1]);
[nn,mm2] = size(lat);

ssh=zeros(size(lat));  %remo nao tem SSH

[la lb]=size(lat)
nnn=cell(la,lb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Delft3d
ncfileb='trim-meso.nc';
bat=ncread(ncfileb,'depth');
bator=ncread(ncfileb,'depth');
X = ncread(ncfileb,'longitude');
Y = ncread(ncfileb,'latitude');
t=isnan(X);
jedi=griddata(X(~t),Y(~t),bat(~t),lon(:),lat(:));
jedi=reshape(jedi,size(lat));
t=isnan(jedi);
jedi(t)=griddata(lat(~t),lon(~t),jedi(~t),lat(t),lon(t), 'nearest'); 
bat=jedi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bat(bat<0)=0;
bat=bat+ssh;


zmax=max(max(bat)); %%%% this variable must be put in delft mdf



prop=[            0.001111 
                  0.001111
                  0.002222
                  0.002222
                  0.004444
                  0.005556
                  0.005556
                  0.005556
                  0.005556
                  0.011111
                  0.011111
                  0.011111
                  0.022222
                  0.022222
                  0.022222
                  0.022222
                  0.022222
                  0.022222
                  0.022222
                  0.022222
                  0.022222
                  0.022222
                  0.022222
                  0.022222
                  0.055556
                  0.055556
                  0.111111
                  0.111111
                  0.111111
                  0.111112
                  0.111112];  
              
for i=1:la
     for j=1:lb
       propbat=zmax*prop;
       [a b]=size(propbat);
       propbats=zeros(a,b);
            for t=1:length(prop)
              if t==1         
               propbats(t)=propbat(t);
              else    
               propbats(t)=propbat(t)+propbats(t-1);
              end
            end
       nnn{i,j}=propbats;
     end 
end              



  uh=ncread(fileh, 'u');
  uh=uh(:,:,:,1);
  vh=ncread(fileh, 'v');
  vh=vh(:,:,:,1);
  temph=ncread(fileh,'temperature');
  temph = temph(:,:,:,1);
  salth=ncread(fileh,'salinity');
  salth = salth(:,:,:,1);
  
  z = ncread(fileh, 'Depth');

  u=permute(uh, [3,2,1]); u(u>1000)=NaN;
  v=permute(vh, [3,2,1]); v(v>1000)=NaN;  
  temp=permute(temph, [3,2,1]); temp(temp>1000)=NaN;
  salt=permute(salth, [3,2,1]); salt(salt>1000)=NaN;

  [nn,mm2] = size(lat);
  
   siz=size(u);
  
  for i=1:siz(1)
    kt=squeeze(temp(i,:,:));
    ks=squeeze(salt(i,:,:));
    ku=squeeze(u(i,:,:));
    kv=squeeze(v(i,:,:));
    
    tn=~isnan(ku);
    
    snans=sum(sum(tn));
    
    if snans~=0
        
      u(i,~tn)=griddata(lat(tn),lon(tn), ku(tn), lat(~tn), lon(~tn), 'nearest');
    
      v(i,~tn)=griddata(lat(tn),lon(tn), kv(tn), lat(~tn), lon(~tn), 'nearest');
      
      temp(i,~tn)=griddata(lat(tn),lon(tn), kt(tn), lat(~tn), lon(~tn), 'nearest');
    
      salt(i,~tn)=griddata(lat(tn),lon(tn), ks(tn), lat(~tn), lon(~tn), 'nearest');      
    
    else 
      u(i,:,:)=u(i-1,:,:);
    
      v(i,:,:)=v(i-1,:,:);
      
      temp(i,:,:)=temp(i-1,:,:);
    
      salt(i,:,:)=salt(i-1,:,:);
    end 

  end
  
  ui=zeros(length(prop), la, lb);
  vi=zeros(length(prop), la, lb);
  ti=zeros(length(prop), la, lb);
  si=zeros(length(prop), la, lb);

    for i=1:la
    for j=1:lb      

        ui(:,i,j)=interp1(z,u(:,i,j), nnn{i,j});
        
        tn=isnan(ui(:,i,j));
        
        snans=sum(sum(tn));
        
        if snans~=0
         ui(tn,i,j) = interp1(z,u(:,i,j), nnn{i,j}(tn), 'pchip');  %Extrapolation
        end
        
        %vi
        vi(:,i,j)=interp1(z,v(:,i,j), nnn{i,j});
        
        tn=isnan(vi(:,i,j));
        
        snans=sum(sum(tn));
        
        if snans~=0
         vi(tn,i,j) = interp1(z,v(:,i,j), nnn{i,j}(tn), 'pchip');  %Extrapolation
        end        
        
        
        ti(:,i,j)=interp1(z,temp(:,i,j), nnn{i,j});
        
        tn=isnan(ti(:,i,j));
        
        snans=sum(sum(tn));
        
        if snans~=0
         ti(tn,i,j) = interp1(z,temp(:,i,j), nnn{i,j}(tn), 'pchip');   %Extrapolation
        end        
        
        %si
        si(:,i,j)=interp1(z,salt(:,i,j), nnn{i,j});
        
        tn=isnan(si(:,i,j));
        
        snans=sum(sum(tn));
        
        if snans~=0
         si(tn,i,j) = interp1(z,salt(:,i,j), nnn{i,j}(tn), 'pchip');  %Extrapolation
        end          
                
    end 
  end 
  
  


ssh = griddata(lon,lat,ssh,X,Y);

[niw njw] = size(ssh);
ssh(isnan(ssh)) = 0;

tstart = tic;
fprintf('%s\t','Writing initial conditions for ssh...')

%ssh=ssh*0;
for i = 1:niw
    for j = 1:njw
        fprintf(fid,'%20.13e',ssh(i,j));
        fprintf(fid,'%s',' ');
    end
    fprintf(fid,'\n');
    %fprintf(fid,'%20.13e\n',ssh(i,j));
end
% ssh_aux = ssh(end,:);
% ssh_aux(end+1) = ssh_aux(end);
% 
% for aux = 1:nj+1
%     fprintf(fid,'%20.13e',ssh_aux(aux));
%     fprintf(fid,'%s',' ');
% end
%fprintf(fid,'\n');
% fprintf(fid,'\n'); 
telapsed = toc(tstart);
fprintf('%s\n',['done! Elapsed time: ' sprintf('%6.4f',telapsed) ' s'])





 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                     UUUUU
  %
%ui=ui*0;

for k = length(prop):-1:1 %remo
%for k = [30 29 28 27 26 25 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 0]
    tstart = tic;
    fprintf('%s\t',['Writing initial conditions for velocity u, layer ' sprintf('%02g',k) '...'])
   

    uint = griddata(lon,lat,squeeze(ui(k,:,:)),X,Y); % or u = griddata(lat,lon,u,Y,X)you choose
     
    uint(isnan(uint)) = 0;
    
    
    for i = 1:niw
        for j = 1:njw
            fprintf(fid,'%20.13e',uint(i,j));
            fprintf(fid,'%s',' ');
        end
        fprintf(fid,'\n');
    end
    
    telapsed = toc(tstart);
    fprintf('%s\n',['done! Elapsed time: ' sprintf('%6.4f',telapsed) ' s'])
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%% V

%vi=vi*0;


for k =length(prop):-1:1
    tstart = tic;
    fprintf('%s\t',['Writing initial conditions for velocity v, layer ' sprintf('%02g',k) '...'])

        vint = griddata(lon,lat,squeeze(vi(k,:,:)),X,Y);
        
        vint(isnan(vint)) = 0;

    
    for i = 1:niw
        for j = 1:njw
            fprintf(fid,'%20.13e',vint(i,j));
            fprintf(fid,'%s',' ');
        end
        fprintf(fid,'\n');
    end
    
    telapsed = toc(tstart);
    fprintf('%s\n',['done! Elapsed time: ' sprintf('%6.4f',telapsed) ' s'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%% Salinity
   




for k = length(prop):-1:1
%for k = [30 29 28 27 26 25 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 0]
    tstart = tic;
    fprintf('%s\t',['Writing initial conditions for salinity, layer ' sprintf('%02g',k) '...'])
%    if k >= 0
%    S = nc_varget(hycom_file,'salinity',[0 k 0 0],[1 1 nii njj]);
    
    
        
    sint = griddata(lon,lat,squeeze(si(k,:,:)),X,Y);
    
    sint(isnan(sint)) = mode(sint(:));
    
    sint(sint<0) = mode(sint(:));


    
    for i = 1:niw
        for j = 1:njw
            fprintf(fid,'%20.13e',sint(i,j));
            fprintf(fid,'%s',' ');
        end
        fprintf(fid,'\n');
    end
    
    telapsed = toc(tstart);
    fprintf('%s\n',['done! Elapsed time: ' sprintf('%6.4f',telapsed) ' s'])
end


  %%%%%%%%%%%%%%%%%%%%%        Temperature
  
    

for k = length(prop):-1:1
%for k = [30 29 28 27 26 25 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 0]
    tstart = tic;
    fprintf('%s\t',['Writing initial conditions for temperature, layer ' sprintf('%02g',k) '...'])
 %   if k >= 0
 %   T = nc_varget(hycom_file,'temperature',[0 k 0 0],[1 1 nii njj]);
        
        
    tint = griddata(lon,lat,squeeze(ti(k,:,:)),X,Y);
    
    
    tint(isnan(tint)) = mode(tint(:));
    tint(tint<0) = mode(tint(:));
    
    
    for i = 1:niw
        for j = 1:njw
            fprintf(fid,'%20.13e',tint(i,j));
            fprintf(fid,'%s',' ');
        end
        fprintf(fid,'\n');
    end
    
    telapsed = toc(tstart);
    fprintf('%s\n',['done! Elapsed time: ' sprintf('%6.4f',telapsed) ' s'])
    
end


% %EFLUENT
% 
% for k = [18    12     8     7     6     5     4     3     2     1     0]
% %for k = [30 29 28 27 26 25 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 0]
%     tstart = tic;
%     fprintf('%s\t',['Writing initial conditions for temperature, layer ' sprintf('%02g',k) '...'])
%     if k >= 0
%     T = nc_varget(hycom_file,'temperature',[0 k 0 0],[1 1 nii njj]);
%         
%     
%     T(isnan(T)) = mode(T(:));
%     T = griddata(lon,lat,T,X,Y);
%     T(isnan(T)) = mode(T(:));
%     T(T<0) = mode(T(:));
%     T=T*0;
%     
%     else
%         T = zeros(ni,nj);
%     end
%     
%     for i = 1:ni
%         for j = 1:nj
%             fprintf(fid,'%20.13e',T(i,j));
%             fprintf(fid,'%s',' ');
%         end
%         fprintf(fid,'\n');
%     end
%     
%     telapsed = toc(tstart);
%     fprintf('%s\n',['done! Elapsed time: ' sprintf('%6.4f',telapsed) ' s'])
%     
% end
% 

% Generate layers to put in .MDF



 fclose(fid);