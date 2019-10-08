clear all; ; clc
cd C:\Users\fernando\Desktop\back_up\delft3d_FM_method_routines\d-data_backup

%/////////////////
%put .dep e o .mdf in the directory above !!!!!!

%put the name of .mdf file
variable='Riemann';

file_tide='great_2_ref_mod_L2_elev.bct';


fidbnd=fopen('meso.bnd')


fileh='2016-01-01.nc';
ncfilebgreat='trim-l1.nc';


t_aux = [datenum(2016,01,01,0,0,0),datenum(2016,01,26,0,0,0)];


thetarad=0;
alfarad=0;
%t=datestr(t_aux, 'yyyy-mm-dd')

range=3;

new_wl=true

withtide = false

curv_interp = false;

old_bat=false
% 
% name='bilo2';
% 
% bat=delftbat2routine(sprintf([name '.dep']),sprintf([name '.mdf']));
% maxx=max(max(bat));
% depth=ncread(ncfile,'Depth');
% diff=abs(depth-maxx);
% [i j]=min(diff);
% lay=j;


modo='REMO';
%modo='HYCOM';

west=false;


north_bnd = load('north_acu.dat'); %north_bnd = north_bnd(:,2);
south_bnd = load('south_acu.dat'); %south_bnd = south_bnd(:,2);
east_bnd = load('east_acu.dat'); %east_bnd = east_bnd(:,1);


if west
  west_bnd = load('westr.dat'); %west_bnd = west_bnd(:,2);
else 
  west_bnd = cat(1,north_bnd,south_bnd);
end
% coef_north = load('coef_north.dat');
% coef_south = load('coef_south.dat');
% coef_east = load('coef_east.dat');
% coef_west = load('coef_west.dat');


% 
% if strcmp(modo,'REMO')
%   disp('ENTRANDO MODO REMO')
%   [bat,lon,lat] = get_hycom_bat_remo(ncfile);
% else
%   disp('ENTRANDO MODO HYCOM')
%   [bat,lon,lat] = get_hycom_bat(ncfile);
% end 

% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Reorganize REMO data
% lonr=ncread(ncfile,'Longitude');
% latr=ncread(ncfile,'Latitude');
% [lon,lat]=meshgrid(lonr,latr);
% u=permute(squeeze(ncread(ncfile,'u')), [3,2,1]);  u(u>1000)=NaN;
% v=permute(squeeze(ncread(ncfile,'v')), [3,2,1]);  v(v>1000)=NaN;
% temp=permute(squeeze(ncread(ncfile,'temperature')), [3,2,1]);  temp(temp>1000)=NaN;
% salt=permute(squeeze(ncread(ncfile,'salinity')), [3,2,1]);     salt(salt>1000)=NaN;
% ssh=zeros(size(lat));  %remo nao tem SSH
% z=squeeze(ncread(ncfile,'Depth'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NEW WATER LEVEL

if new_wl
    lonG = nc_read(ncfilebgreat,'longitude');
    latG = nc_read(ncfilebgreat,'latitude');
    WL = nc_read(ncfilebgreat, 'waterlevel');
    [nnG,mm2G] = size(latG);
    WL(isnan(WL))=0;
    lonG(isnan(lonG))=0;
    latG(isnan(latG))=0;
end    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Reorganize HYCOM data
lath=ncread(fileh,'Latitude');
lonh=ncread(fileh,'Longitude')-360;
lat=permute(lath, [2,1]);    
lon=permute(lonh, [2,1]);
[nn,mm2] = size(lat);

ssh=zeros(size(lat));  %remo nao tem SSH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INDEXES TO CROP HYCOM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA

one=ones(nn,mm2);
northl=one*max(max(north_bnd(:,1)));
northl=northl(:,1);
latn=lat(:,1);
difn=abs(latn-northl);
[i limlatN]=min(difn);      %limlatN +range  since north is at the bottom of the matrice

one=ones(nn,mm2);
northl=one*min(min(south_bnd(:,1)));
northl=northl(:,1);
latn=lat(:,1);
difn=abs(latn-northl);
[i limlatS]=min(difn);    %limlatS -range
    
one=ones(nn,mm2);
northl=one*max(max(east_bnd(:,2)));
northl=northl(1,:);
latn=lon(1,:);
difn=abs(latn-northl);
[i limlonE]=min(difn);    %limlonE + range

one=ones(nn,mm2);
northl=one*min(min(west_bnd(:,2)));
northl=northl(1,:);
latn=lon(1,:);
difn=abs(latn-northl);
[i limlonW]=min(difn);   %limlonW -range


lat=lat(limlatS-range:limlatN+range, limlonW-range:limlonE+range);

lon=lon(limlatS-range:limlatN+range, limlonW-range:limlonE+range);

ssh=ssh(limlatS-range:limlatN+range, limlonW-range:limlonE+range);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
%%%%BATHYMETRY FROM HYCOM
% 
% nc='depth_GLBa0.08_09.nc';
% bati=ncread(nc,'bathymetry');
% latet=ncread(nc,'Latitude');
% lonet=ncread(nc,'Longitude');
% 
% bat=griddata(latet,lonet, bati, lat,lon);
% bat(bat>6000)=NaN;
% bat(bat>max(z))=max(z);  as=isnan(bat);
% bat=griddata(lat(~as), lon(~as), bat(~as), lat, lon, 'nearest');

%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% DELFT3D FILE

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bat(bat<0)=0;
bat=bat+ssh;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BND sections name
sizlat=size(Y);
n=0 
tline='a'
while ischar(tline)
   n=n+1;
   tline=fgetl(fidbnd);
end
n=n-1;
frewind(fidbnd)
sections=cell(n, 5);
tline=fgetl(fidbnd);
i=0;
j=1;
while i~=n
    d=textscan(tline, '%s');
    m1=str2double(d{1}{4})
    n1=str2double(d{1}{5})
    m2=str2double(d{1}{6})
    n2=str2double(d{1}{7})
    if n1==1
       n1=n1+1
    end
    if n2==1
       n2=n2+1
    end
    if m1==1
       m1=m1+1
    end     
    if m2==1
       m2=m2+1
    end     
    if n1==sizlat(1)
       n1=n1-1
    end
    if n2==sizlat(1)
        n2=n2-1
    end           
    if m1==sizlat(2)
       m1=m1-1
    end     
    if m2==sizlat(2)
       m2=m2-1
    end              
    sect{j,1} = d{1}{1};
    sect{j,2} = m1;
    sect{j,3} = n1;
    sect{j,4} = m2;
    sect{j,5} = n2;
    tline=fgetl(fidbnd);    
    i=i+1;  
    j=j+1;
end

%/////////////////////



[la lb]=size(lat)
nnn=cell(la,lb);

bat(bat>5000)=5000;

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
                           


zmax=max(max(bat)); %%%% this variable must be put in delft mdf
fprintf('  Top Depth is  %f m, press Enter to continue', zmax);


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indext=1

  ui=zeros(length(prop), la, lb, length(t_aux(1):t_aux(2)));
  vi=zeros(length(prop), la, lb, length(t_aux(1):t_aux(2)));
  ti=zeros(length(prop), la, lb, length(t_aux(1):t_aux(2)));
  si=zeros(length(prop), la, lb, length(t_aux(1):t_aux(2)));

for v=t_aux(1):t_aux(2)
    
  name=datestr(v, 'yyyy-mm-dd');
  fileh = [name '.nc'];
  
  lath=ncread(fileh,'Latitude');
  lonh=ncread(fileh,'Longitude')-360;
  uh=ncread(fileh, 'u');
  uh=uh(:,:,:,1);
  vh=ncread(fileh, 'v');
  vh=vh(:,:,:,1);
%  temph=ncread(fileh,'temperature');
%  temph = temph(:,:,:,1);
%  salth=ncread(fileh,'salinity');
%  salth = salth(:,:,:,1);
  z = ncread(fileh, 'Depth');
  lat=permute(lath, [2,1]);    
  lon=permute(lonh, [2,1]);
  u=permute(uh, [3,2,1]); u(u>1000)=NaN;
  v=permute(vh, [3,2,1]); v(v>1000)=NaN;  
%  temp=permute(temph, [3,2,1]); temp(temp>1000)=NaN;
%  salt=permute(salth, [3,2,1]); salt(salt>1000)=NaN;
  %a=rot90(lath);
  [nn,mm2] = size(lat);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INDEXES TO CROP HYCOM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA

  one=ones(nn,mm2);
  northl=one*max(max(north_bnd(:,1)));
  northl=northl(:,1);
  latn=lat(:,1);
  difn=abs(latn-northl);
  [i limlatN]=min(difn);      %limlatN +5

  one=ones(nn,mm2);
  northl=one*min(min(south_bnd(:,1)));
  northl=northl(:,1);
  latn=lat(:,1);
  difn=abs(latn-northl);
  [i limlatS]=min(difn);    %limlatS -5
    
  one=ones(nn,mm2);
  northl=one*max(max(east_bnd(:,2)));
  northl=northl(1,:);
  latn=lon(1,:);
  difn=abs(latn-northl);
  [i limlonE]=min(difn);    %limlonE + 5

  one=ones(nn,mm2);
  northl=one*min(min(west_bnd(:,2)));
  northl=northl(1,:);
  latn=lon(1,:);
  difn=abs(latn-northl);
  [i limlonW]=min(difn);   %limlonW -5

  lat=lat(limlatS-range:limlatN+range, limlonW-range:limlonE+range);
  
  lon=lon(limlatS-range:limlatN+range, limlonW-range:limlonE+range);

  u=u(:,limlatS-range:limlatN+range, limlonW-range:limlonE+range);

  v=v(:,limlatS-range:limlatN+range, limlonW-range:limlonE+range);

%  temp=temp(:,limlatS-range:limlatN+range, limlonW-range:limlonE+range);

%  salt=salt(:,limlatS-range:limlatN+range, limlonW-range:limlonE+range);
  
  siz=size(u);
  
  for i=1:siz(1)
%    kt=squeeze(temp(i,:,:));
%    ks=squeeze(salt(i,:,:));
    ku=squeeze(u(i,:,:));
    kv=squeeze(v(i,:,:));
    
    tn=~isnan(ku);
    
    snans=sum(sum(tn));
    
    if snans~=0
      u(i,~tn)=griddata(lat(tn),lon(tn), ku(tn), lat(~tn), lon(~tn), 'nearest');
    
      v(i,~tn)=griddata(lat(tn),lon(tn), kv(tn), lat(~tn), lon(~tn), 'nearest');
    
    else 
      u(i,:,:)=u(i-1,:,:);
    
      v(i,:,:)=v(i-1,:,:);
    
    end 

end

%%%%%%%%%%%%%%%%%%%%%%
% nc='depth_GLBa0.08_09.nc';
% bati=ncread(nc,'bathymetry');
% latet=ncread(nc,'Latitude');
% lonet=ncread(nc,'Longitude');
% 
% bat=griddata(latet,lonet, bati, lat,lon);
% bat(bat>6000)=NaN; as=isnan(bat);
% bat=griddata(lat(~as), lon(~as), bat(~as), lat, lon, 'nearest');

%y=zeros(size(u));

%for i=1:la
%   for j=1:lb
%        y(:,i,j)=y(:,i,j)+ubar(i,j);
%    end 
%end 

%u=u+y;

%for i=1:la
%    for j=1:lb
%        y(:,i,j)=y(:,i,j)+vbar(i,j);
%    end 
%end 

%v=v+y;

  for i=1:la
    for j=1:lb      
%         bef=interp1(z,u(:,i,j), nnn{i,j});
%         a1=isnan(bef);
%         b1=nnn{i,j};
%         bef(isnan(bef))=interp1(z,u(:,i,j),b1(a1),'pchip');
%         ui(:,i,j)=bef;
%         %vi
%         bef=interp1(z,v(:,i,j), nnn{i,j});
%         a1=isnan(bef);
%         b1=nnn{i,j};
%         bef(isnan(bef))=interp1(z,v(:,i,j),b1(a1),'pchip');        
%         vi(:,i,j)=bef;
        %ti
        ui(:,i,j,indext)=interp1(z,u(:,i,j), nnn{i,j});
        
        tn=isnan(ui(:,i,j,indext));
        
        snans=sum(sum(tn));
        
        if snans~=0
         ui(tn,i,j,indext) = interp1(z,u(:,i,j), nnn{i,j}(tn), 'pchip');  %Extrapolation
        end
        
        %vi
        vi(:,i,j,indext)=interp1(z,v(:,i,j), nnn{i,j});
        
        tn=isnan(vi(:,i,j,indext));
        
        snans=sum(sum(tn));
        
        if snans~=0
         vi(tn,i,j,indext) = interp1(z,v(:,i,j), nnn{i,j}(tn), 'pchip');  %Extrapolation
        end        
                
    end 
  end 
  
  indext=indext+1

end

latfm=lat;
lonfm=lon;
numb=length(prop);
vii=vi;

if new_wl
    
    save('transpose.mat', 'latfm', 'lonfm', 'numb', 'bat', 'ssh', 'lonG', 'latG');  
    
else 
    
    save('transpose.mat', 'latfm', 'lonfm', 'numb', 'bat', 'ssh');  
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END SOUTELINO'S FM



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END SOUTELINO'S FM

%bat=delftbat2routine('l1.dep','l1.mdf');
%bat(:,:)=2000
ni=1;
nj=1;


% bat_north = bat(end,:);
% bat_east = bat(:,end);
% bat_south = bat(1,:);
% bat_west = bat(:,1);

% lon = lon(1,:);
% lat = lat(:,1); lat = lat';

timevec = datevec(t_aux(1));

dd = timevec(3);
mm = timevec(2);
yyyy = timevec(1);

fid = fopen('water_level_old_1.bct','wt');


[nn,mm2] = size(lat);


c = 1;

bat(bat>5000)=5000;

lay=numb;  

bdsec=1;

depth_i=zeros(2,1);

numtimes=diff(t_aux)*24 + 1;

for bd = 1:length(north_bnd)/2
    
    lati = north_bnd(c:c+1,1);
    loni = north_bnd(c:c+1,2);
    coef=[1,1];
    tref = 60;    
    tin = 1;
    
    if curv_interp
      deltax=abs(diff(loni));
      deltay=abs(diff(lati));
      
      thetarad = atan(deltax/deltay);
      thetadeg = rad2deg(atan(deltax/deltay));
      alfadeg=abs(thetadeg + 90 - 180);
      alfarad = deg2rad(alfadeg);
        
     
    end
    
    %    coef = coef_north(bd,:);
    
    one=ones(nn,mm2);
    northl=one*lati(1,1);
    difn=abs(lat-northl);
    
    one=ones(nn,mm2);
    eastl=one*loni(1,1);
    dife=abs(lon-eastl);
    
    sumdif=difn+dife;    
    
    [limlat, limlon]=find(sumdif==min(sumdif(:)));    
    
    if new_wl
        
      one=ones(nnG,mm2G);
      northl=one*lati(1,1);
      difn=abs(latG-northl);
    
      one=ones(nnG,mm2G);
      eastl=one*loni(1,1);
      dife=abs(lonG-eastl);
    
      sumdif=difn+dife;    
    
      [limlatG, limlonG]=find(sumdif==min(sumdif(:)));   
              
      tide=0;
     
    end
    
    
    write_delft_header_mod(fid,sect{bdsec},bd,dd,mm,yyyy,length(prop),numtimes,variable);
    
    fprintf(fid,'%15.2f %s',0,' ');fprintf_bound(fid,prop*0); fprintf_bound(fid,prop*0); fprintf(fid,'\n');
 
    
    disp(sect{bdsec})
    
    if ~new_wl && withtide       
      eta_tide = get_bnd_tide(file_tide,bd,'North');
      
    end
    
   
    for t = 1:t_aux(2) - t_aux(1)
        
        hycom_file_0=t;
        hycom_file_1=hycom_file_0+1;
        
%        fprintf('%s\n',['Escrevendo dados - Boundary North' sprintf('%03g',bd)...
%            ' - hydro ---> tempo ' sprintf('%02g',t)])
        
        if ~new_wl && withtide        
          tide = eta_tide(tin:tin+24,1:2); tide = tide';
        else
          tide=0;
        end

        
      depth_i(1) = bator(sect{bdsec,3},sect{bdsec,2});
      depth_i(2) = bator(sect{bdsec,5},sect{bdsec,4});
     

     
      if new_wl
       [interp_property] = interp_delft_boundary_hydro_REF_2019(hycom_file_0,hycom_file_1,...
            'North',variable,ni,nj,lat,lon,lati,loni,tide,coef(1),coef(2), depth_i, old_bat, ui,vii, bat,lay,modo,...
            limlat,limlon, new_wl,limlatG, limlonG, tin, WL, curv_interp, thetarad, alfarad);          

       else
       [interp_property] = interp_delft_boundary_hydro_REF_2019(hycom_file_0,hycom_file_1,...
            'North',variable,ni,nj,lat,lon,lati,loni,tide,coef(1),coef(2), depth_i, old_bat,ui,vii, bat,lay,modo,...
            limlat,limlon, new_wl, curv_interp, thetarad, alfarad);
          
      end
         
      
        for tinp = 2:25
            
             fprintf(fid,'%15.2f %s',tref,' ');
             
             if strcmp(variable, 'WaterLevel')
               fprintf_bound(fid,interp_property(1,tinp));
               fprintf_bound(fid,interp_property(2,tinp));                                 
             else
               fprintf_bound(fid,interp_property(:,1,tinp));
               fprintf_bound(fid,interp_property(:,2,tinp));     
             end 
             tref = tref + 60;
             fprintf(fid,'\n');
        end
        tin = tin + 24;
    end
    c = c + 2;
    
    bdsec = bdsec + 1;
    
end




c = 1;

for bd = 1:length(south_bnd)/2
    
    lati = south_bnd(c:c+1,1);
    loni = south_bnd(c:c+1,2);
    coef=[1,1];
    tref = 60;    
    tin = 1;    
%    coef = coef_south(bd,:);

    if curv_interp
      deltax=abs(diff(loni));
      deltay=abs(diff(lati));
      
      thetarad = atan(deltax/deltay); 
      thetadeg = rad2deg(atan(deltax/deltay)); 
      alfadeg=abs(thetadeg + 90 - 180);
      alfarad = deg2rad(alfadeg);
     
    end
    
    
    one=ones(nn,mm2);
    northl=one*lati(1,1);
    difn=abs(lat-northl);
    
    one=ones(nn,mm2);
    eastl=one*loni(1,1);
    dife=abs(lon-eastl);
    
    sumdif=difn+dife;    
    
    [limlat, limlon]=find(sumdif==min(sumdif(:))); 
    
    if new_wl
        
      one=ones(nnG,mm2G);
      northl=one*lati(1,1);
      difn=abs(latG-northl);
    
      one=ones(nnG,mm2G);
      eastl=one*loni(1,1);
      dife=abs(lonG-eastl);
    
      sumdif=difn+dife;    
    
      [limlatG, limlonG]=find(sumdif==min(sumdif(:)));   
              
      tide=0
     
    end
       
    
    write_delft_header_mod(fid,sect{bdsec},bd,dd,mm,yyyy,length(prop),numtimes,variable);
    fprintf(fid,'%15.2f %s',0,' ');fprintf_bound(fid,prop*0); fprintf_bound(fid,prop*0); fprintf(fid,'\n');

        disp(sect{bdsec})

    
  if ~new_wl  && withtide      
  
    eta_tide = get_bnd_tide(file_tide,bd,'South');

  end 
  
    for t = 1:t_aux(2) - t_aux(1)
        
        hycom_file_0=t;
        hycom_file_1=hycom_file_0+1;
        
%        fprintf('%s\n',['Escrevendo dados - Boundary South' sprintf('%03g',bd)...
%            ' - hydro ---> tempo ' sprintf('%02g',t)])   
        
        if ~new_wl && withtide
          tide = eta_tide(tin:tin+24,1:2); tide = tide';
        else
          tide=0;
        end
        
%bdsec
      depth_i(1) = bator(sect{bdsec,3},sect{bdsec,2}); 
      depth_i(2) = bator(sect{bdsec,5},sect{bdsec,4});
      

      if new_wl
       [interp_property] = interp_delft_boundary_hydro_REF_2019(hycom_file_0,hycom_file_1,...
            'South',variable,ni,nj,lat,lon,lati,loni,tide,coef(1),coef(2),depth_i, old_bat,ui,vii, bat,lay,modo,...
            limlat,limlon, new_wl,limlatG, limlonG, tin, WL, curv_interp, thetarad, alfarad);          

       else
       [interp_property] = interp_delft_boundary_hydro_REF_2019(hycom_file_0,hycom_file_1,...
            'South',variable,ni,nj,lat,lon,lati,loni,tide,coef(1),coef(2),depth_i, old_bat,ui,vii, bat,lay,modo,...
            limlat,limlon, new_wl, curv_interp, thetarad, alfarad);
          
      end
        
        for tinp = 2:25
            
            fprintf(fid,'%15.2f %s',tref,' ');
            
             if strcmp(variable, 'WaterLevel')
               fprintf_bound(fid,interp_property(1,tinp));
               fprintf_bound(fid,interp_property(2,tinp));                                 
             else
               fprintf_bound(fid,interp_property(:,1,tinp));
               fprintf_bound(fid,interp_property(:,2,tinp));     
             end            
            
            tref = tref + 60;
            fprintf(fid,'\n');
        end
        tin = tin + 24;
    end
    c = c + 2;

    bdsec = bdsec + 1;

end


c = 1;

for bd = 1:length(east_bnd)/2
    lati = east_bnd(c:c+1,1);
    loni = east_bnd(c:c+1,2);
    coef=[1,1];
    tref = 60;    
    tin = 1;  
 %   coef = coef_east(bd,:);
    
    if curv_interp
      deltax=abs(diff(lati));
      deltay=abs(diff(loni));
      
      thetarad = atan(deltax/deltay); 
      thetadeg = rad2deg(atan(deltax/deltay)); 
      alfadeg=abs(thetadeg + 90 - 180);
      alfarad = deg2rad(alfadeg);
     
    end
    
    one=ones(nn,mm2);
    northl=one*lati(1,1);
    difn=abs(lat-northl);
    
    one=ones(nn,mm2);
    eastl=one*loni(1,1);
    dife=abs(lon-eastl);
    
    sumdif=difn+dife;    
    
    [limlat, limlon]=find(sumdif==min(sumdif(:))); 
    
    if new_wl
        
      one=ones(nnG,mm2G);
      northl=one*lati(1,1);
      difn=abs(latG-northl);
    
      one=ones(nnG,mm2G);
      eastl=one*loni(1,1);
      dife=abs(lonG-eastl);
    
      sumdif=difn+dife;    
    
      [limlatG, limlonG]=find(sumdif==min(sumdif(:)));   
              
      tide=0
     
    end   
       
    
    write_delft_header_mod(fid,sect{bdsec},bd,dd,mm,yyyy,length(prop),numtimes,variable);
    fprintf(fid,'%15.2f %s',0,' ');fprintf_bound(fid,prop*0); fprintf_bound(fid,prop*0); fprintf(fid,'\n');

        disp(sect{bdsec})

  if ~new_wl && withtide
        
    eta_tide = get_bnd_tide(file_tide,bd,'East');

  end
  
    for t = 1:t_aux(2) - t_aux(1)
        
        hycom_file_0=t;
        hycom_file_1=hycom_file_0+1;
        
%        fprintf('%s\n',['Escrevendo dados - Boundary East' sprintf('%03g',bd)...
%            ' - hydro ---> tempo ' sprintf('%02g',t)])   
        
        if ~new_wl && withtide            
        
          tide = eta_tide(tin:tin+24,1:2); tide = tide';
        else
          tide=0;
        
        end
        
 
%  bdsec
      depth_i(1) = bator(sect{bdsec,3},sect{bdsec,2});
      depth_i(2) = bator(sect{bdsec,5},sect{bdsec,4});
      

      
      if new_wl
       [interp_property] = interp_delft_boundary_hydro_REF_2019(hycom_file_0,hycom_file_1,...
            'East',variable,ni,nj,lat,lon,lati,loni,tide,coef(1),coef(2),depth_i, old_bat,ui,vii, bat,lay,modo ...
            , limlat,limlon, new_wl,limlatG, limlonG, tin, WL, curv_interp, thetarad, alfarad);          

       else
       [interp_property] = interp_delft_boundary_hydro_REF_2019(hycom_file_0,hycom_file_1,...
            'East',variable,ni,nj,lat,lon,lati,loni,tide,coef(1),coef(2),depth_i, old_bat,ui,vii, bat,lay,modo...
            , limlat,limlon, new_wl, curv_interp, thetarad, alfarad);
          
      end
        
        for tinp = 2:25
            
            fprintf(fid,'%15.2f %s',tref,' ');
            
             if strcmp(variable, 'WaterLevel')
               fprintf_bound(fid,interp_property(1,tinp));
               fprintf_bound(fid,interp_property(2,tinp));                                 
             else
               fprintf_bound(fid,interp_property(:,1,tinp));
               fprintf_bound(fid,interp_property(:,2,tinp));     
             end            
            
            tref = tref + 60;
            fprintf(fid,'\n');
        end
        tin = tin + 24;
    end
    c = c + 2;
    
    bdsec = bdsec + 1;

end

if ~west
  return
end

c = 1;

for bd = 1:length(west_bnd)/2
    lati = west_bnd(c:c+1,1);
    loni = west_bnd(c:c+1,2);
    coef=[1,1]
    tref = 60;    
    tin = 1;      
%    coef = coef_west(bd,:);

    if curv_interp
      deltax=abs(diff(lati));
      deltay=abs(diff(loni));
      
      thetarad = atan(deltax/deltay); 
      thetadeg = rad2deg(atan(deltax/deltay)); 
      alfadeg=abs(thetadeg + 90 - 180);
      alfarad = deg2rad(alfadeg);
     
    end
    
    one=ones(nn,mm2);
    northl=one*lati(1,1);
    difn=abs(lat-northl);
    
    one=ones(nn,mm2);
    eastl=one*loni(1,1);
    dife=abs(lon-eastl);
    
    sumdif=difn+dife;    
    
    [limlat, limlon]=find(sumdif==min(sumdif(:))); 
 
    if new_wl
        
      one=ones(nnG,mm2G);
      northl=one*lati(1,1);
      difn=abs(latG-northl);
    
      one=ones(nnG,mm2G);
      eastl=one*loni(1,1);
      dife=abs(lonG-eastl);
    
      sumdif=difn+dife;    
    
      [limlatG, limlonG]=find(sumdif==min(sumdif(:)));   
              
      tide=0
     
    end   
       
    
    write_delft_header_mod(fid,sect{bdsec},bd,dd,mm,yyyy,length(prop),numtimes,variable);
    fprintf(fid,'%15.2f %s',0,' ');fprintf_bound(fid,prop*0); fprintf_bound(fid,prop*0); fprintf(fid,'\n');    
    
    if ~new_wl && withtide
        
      eta_tide = get_bnd_tide(file_tide,bd,'West');
      
    end 
    
    for t = 1:t_aux(2) - t_aux(1)
        
        hycom_file_0=t
        hycom_file_1=hycom_file_0+1
        
        fprintf('%s\n',['Escrevendo dados - Boundary West' sprintf('%03g',bd)...
            ' - hydro ---> tempo ' sprintf('%02g',t)])  
        
      if ~new_wl && withtide        
         
        tide = eta_tide(tin:tin+24,1:2); tide = tide';
      else
        tide=0
        
      end 
        
      depth_i(1) = bator(sect{bdsec,3},sect{bdsec,2});
      depth_i(2) = bator(sect{bdsec,5},sect{bdsec,4});       
        
      if new_wl
       [interp_property] = interp_delft_boundary_hydro_REF_2019(hycom_file_0,hycom_file_1,...
            'West',variable,ni,nj,lat,lon,lati,loni,tide,coef(1),coef(2),depth_i, old_bat,ui,vii, bat,lay,modo...
            , limlat,limlon, new_wl,limlatG, limlonG, tin, WL, curv_interp, thetarad, alfarad);          

       else
       [interp_property] = interp_delft_boundary_hydro_REF_2019(hycom_file_0,hycom_file_1,...
            'West',variable,ni,nj,lat,lon,lati,loni,tide,coef(1),coef(2),depth_i, old_bat,ui,vii, bat,lay,modo...
            , limlat,limlon, new_wl, curv_interp, thetarad, alfarad);
          
      end
        
        for tinp = 2:25
            
            fprintf(fid,'%15.2f %s',tref,' ');
            
             if strcmp(variable, 'WaterLevel')
               fprintf_bound(fid,interp_property(1,tinp));
               fprintf_bound(fid,interp_property(2,tinp));                                 
             else
               fprintf_bound(fid,interp_property(:,1,tinp));
               fprintf_bound(fid,interp_property(:,2,tinp));     
             end                
            
            tref = tref + 60;
            fprintf(fid,'\n');
        end
        tin = tin + 24;
    end
    c = c + 2;
    
    bdsec = bdsec + 1

end

fclose(fid)