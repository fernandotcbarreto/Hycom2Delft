clear all; ; clc
cd C:\Users\fernando\Desktop\back_up\delft3d_FM_method_routines\d-data_backup

%ncfile = 'archv.2004_001_00_3ztsuv_nc4.nc';

fileh='2016-01-01.nc';

fidbnd=fopen('meso.bnd')

t_aux = [datenum(2016,01,01,0,0,0),datenum(2016,01,26,0,0,0)];

%t=datestr(t_aux, 'yyyy-mm-dd')

range=3;

%put .dep e o .mdf in the directory above !!!!!!

%put the name of .mdf file
% 
% name='island';
% 
% bat=delftbat2routine(sprintf([name '.dep']),sprintf([name '.mdf']));
% maxx=max(max(bat));
% depth=ncread(ncfile,'Depth');
% diff=abs(depth-maxx);
% [i j]=min(diff);
% lay=j;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BND sections name
n=0 
tline='a'
while ischar(tline)
   n=n+1;
   tline=fgetl(fidbnd);
end
n=n-1
frewind(fidbnd)
sections=cell(n);
tline=fgetl(fidbnd);
i=0;
j=1;
while i~=n
    d=textscan(tline, '%s');
    sections{j} = d{1}{1};
    tline=fgetl(fidbnd);    
    i=i+1;  
    j=j+1;
end


%/////////////////////
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
% if strcmp(modo,'REMO')
%   disp('ENTRANDO MODO REMO')
%   [bat,lon,lat] = get_hycom_bat_remo(ncfile);
% else
%   disp('ENTRANDO MODO HYCOM')
%   [bat,lon,lat] = get_hycom_bat(ncfile);
% end 
% 
% 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Reorganize REMO data
% lonr=ncread(ncfile,'Longitude');
% latr=ncread(ncfile,'Latitude');
% [lon,lat]=meshgrid(lonr,latr);
% u=permute(squeeze(ncread(ncfile,'u')), [3,2,1]);  u(u>1000)=NaN;
% v=permute(squeeze(ncread(ncfile,'v')), [3,2,1]);  v(v>1000)=NaN;
% temp=permute(squeeze(ncread(ncfile,'temperature')), [3,2,1]);  temp(temp>1000)=NaN;
% salt=permute(squeeze(ncread(ncfile,'salinity')), [3,2,1]);     salt(salt>1000)=NaN;
% ssh=zeros(size(lat));  %remo nao tem SSH
% z=squeeze(ncread(ncfile,'Depth'));


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


%%%%%%%%%%%%%%%%%%%%%% DELFT3D FILE

ncfileb='trim-meso.nc';
bat=ncread(ncfileb,'depth');
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
%  uh=ncread(fileh, 'u');
%  uh=uh(:,:,:,1);
%  vh=ncread(fileh, 'v');
%  vh=vh(:,:,:,1);
  temph=ncread(fileh,'temperature');
  temph = temph(:,:,:,1);
  salth=ncread(fileh,'salinity');
  salth = salth(:,:,:,1);
  z = ncread(fileh, 'Depth');
 % u=permute(uh, [3,2,1]); u(u>1000)=NaN;
 % v=permute(vh, [3,2,1]); v(v>1000)=NaN;
  lat=permute(lath, [2,1]);    
  lon=permute(lonh, [2,1]);
  temp=permute(temph, [3,2,1]); temp(temp>1000)=NaN;
  salt=permute(salth, [3,2,1]); salt(salt>1000)=NaN;
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

 % u=u(:,limlatS-range:limlatN+range, limlonW-range:limlonE+range);

  %v=v(:,limlatS-range:limlatN+range, limlonW-range:limlonE+range);

  temp=temp(:,limlatS-range:limlatN+range, limlonW-range:limlonE+range);

  salt=salt(:,limlatS-range:limlatN+range, limlonW-range:limlonE+range);
  
  siz=size(temp);
  
  for i=1:siz(1)
    kt=squeeze(temp(i,:,:));
    ks=squeeze(salt(i,:,:));
%    ku=squeeze(u(i,:,:));
%    kv=squeeze(v(i,:,:));
    
    tn=~isnan(kt);
    
    snans=sum(sum(tn));
    
    if snans~=0
      temp(i,~tn)=griddata(lat(tn),lon(tn), kt(tn), lat(~tn), lon(~tn), 'nearest');
    
      salt(i,~tn)=griddata(lat(tn),lon(tn), ks(tn), lat(~tn), lon(~tn), 'nearest');
    
    else 
      temp(i,:,:)=temp(i-1,:,:);
    
      salt(i,:,:)=salt(i-1,:,:);
    
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
        ti(:,i,j,indext)=interp1(z,temp(:,i,j), nnn{i,j});
        
        tn=isnan(ti(:,i,j,indext));
        
        snans=sum(sum(tn));
        
        if snans~=0
         ti(tn,i,j,indext) = interp1(z,temp(:,i,j), nnn{i,j}(tn), 'pchip');   %Extrapolation
        end        
        
        %si
        si(:,i,j,indext)=interp1(z,salt(:,i,j), nnn{i,j});
        
        tn=isnan(si(:,i,j,indext));
        
        snans=sum(sum(tn));
        
        if snans~=0
         si(tn,i,j,indext) = interp1(z,salt(:,i,j), nnn{i,j}(tn), 'pchip');  %Extrapolation
        end  
        
                
    end 
  end 
  
  indext=indext+1

end

latfm=lat;
lonfm=lon;
numb=length(prop);
vii=vi;
save('transpose_ts.mat','latfm', 'lonfm', 'numb');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END SOUTELINO'S FM


%bat=delftbat2routine('l1.dep','l1.mdf');
%[ni,nj] = size(lat);

ni=1

nj=1

%bat_west = bat(:,1);


% lon = lon(1,:);
% lat = lat(:,1); lat = lat';

% variaveis auxiliares de tempo

timevec = datevec(t_aux(1));

dd = timevec(3);
mm = timevec(2);
yyyy = timevec(1);

fid = fopen('AAA_3.bcc','wt');


[nn,mm2] = size(lat);


c_north = 1;

bat(bat>5000)=5000;

lay=numb   

bdsec=1


for bd = 1:length(north_bnd)/2

    lati = north_bnd(c_north:c_north+1,1);
    loni = north_bnd(c_north:c_north+1,2);
    tref = 0;
    tin = 1;
    coef0 = 0;
    coef1 = 0;
    
    
    one=ones(nn,mm2);
    northl=one*lati(1,1);
    difn=abs(lat-northl);
    
    one=ones(nn,mm2);
    eastl=one*loni(1,1);
    dife=abs(lon-eastl);
    
    sumdif=difn+dife;    
    
    [limlat, limlon]=find(sumdif==min(sumdif(:)));
    
    write_delft_header_mod(fid,sections{bdsec},bd,dd,mm,yyyy,length(prop),2,'Salinity');


    for t = 1:t_aux(2) - t_aux(1)
        hycom_file_0=t;
        hycom_file_1=hycom_file_0+1;
        
%        fprintf('%s\n',['Escrevendo dados - Boundary North' sprintf('%03g',bd)...
%            ' - Salinidade ---> tempo ' sprintf('%02g',t)])          
        
        %fim dos nomes para remo
        
%         hycom_file_0 = ['Hycom-noda_' num2str(tstep(1)) sprintf('%02g',tstep(2))...
%             sprintf('%02g',tstep(3)) '.nc'];
%         hycom_file_1 = ['Hycom-noda_' num2str(tstepf(1)) sprintf('%02g',tstepf(2))...
%             sprintf('%02g',tstepf(3)) '.nc'];
        
        coef1 = coef1 + 0.01;
        
         [interp_property] = interp_delft_boundary_TS_REF_2019(hycom_file_0,hycom_file_1,...
             'North','salinity',ni,nj,lat,lon,lati,loni,0,coef0,coef1,ti, si,bat,lay, modo, limlat,limlon);
         
%          aux = zeros(31,2,25);
%          aux(5:end,:,:) = interp_property;
        
       if tref==0       
          for tinp = 1:25
            
              fprintf(fid,'%15.2f %s',tref,' ');
            
              fprintf_bound(fid,interp_property(:,1,tinp));
              fprintf_bound(fid,interp_property(:,2,tinp));
              tref = tref + 60;
              fprintf(fid,'\n');
          end           
           

       else 
          for tinp = 2:25
            
              fprintf(fid,'%15.2f %s',tref,' ');
            
              fprintf_bound(fid,interp_property(:,1,tinp));
              fprintf_bound(fid,interp_property(:,2,tinp));
              tref = tref + 60;
              fprintf(fid,'\n');
          end
        
       end 
                    
        tin = tin + 24;
    end
    

    lati = north_bnd(c_north:c_north+1,1);
    loni = north_bnd(c_north:c_north+1,2);
    tref = 0;
    tin = 1;
    coef0 = 0;
    coef1 = 0;
    write_delft_header_mod(fid,sections{bdsec},bd,dd,mm,yyyy,length(prop),2,'Temperature');

    for t = 1:t_aux(2) - t_aux(1)
        hycom_file_0=t;
        hycom_file_1=hycom_file_0+1;        
        
        
        
%        fprintf('%s\n',['Escrevendo dados - Boundary North' sprintf('%03g',bd)...
%            ' - Temperatura ---> tempo ' sprintf('%02g',t)])
        
       
        %fim dos nomes para remo
        
%         hycom_file_0 = ['Hycom-noda_' num2str(tstep(1)) sprintf('%02g',tstep(2))...
%             sprintf('%02g',tstep(3)) '.nc'];
%         hycom_file_1 = ['Hycom-noda_' num2str(tstepf(1)) sprintf('%02g',tstepf(2))...
%             sprintf('%02g',tstepf(3)) '.nc'];
        
        coef1 = coef1 + 0.01;
        
         [interp_property] = interp_delft_boundary_TS_REF_2019(hycom_file_0,hycom_file_1,...
             'North','temperature',ni,nj,lat,lon,lati,loni,0,coef0,coef1,ti, si,bat,lay, modo, limlat,limlon);
         
%          aux = zeros(33,2,25);
%          aux(5:end,:,:) = interp_property;
%          interp_property = aux;       

        
       if tref==0       
          for tinp = 1:25
            
              fprintf(fid,'%15.2f %s',tref,' ');
            
              fprintf_bound(fid,interp_property(:,1,tinp));
              fprintf_bound(fid,interp_property(:,2,tinp));
              tref = tref + 60;
              fprintf(fid,'\n');
          end           
           

       else 
          for tinp = 2:25
            
              fprintf(fid,'%15.2f %s',tref,' ');
            
              fprintf_bound(fid,interp_property(:,1,tinp));
              fprintf_bound(fid,interp_property(:,2,tinp));
              tref = tref + 60;
              fprintf(fid,'\n');
          end
        
       end 
       
        
        tin = tin + 24;
    end
    
    
    c_north = c_north + 2;
    
    bdsec = bdsec + 1

end

c_south = 1;
 

for bd = 1:length(south_bnd)/2

    lati = south_bnd(c_south:c_south+1,1);
    loni = south_bnd(c_south:c_south+1,2);
    tref = 0;
    tin = 1;
    coef0 = 0;
    coef1 = 0;

    one=ones(nn,mm2);
    northl=one*lati(1,1);
    difn=abs(lat-northl);
    
    one=ones(nn,mm2);
    eastl=one*loni(1,1);
    dife=abs(lon-eastl);
    
    sumdif=difn+dife;    
    
    [limlat, limlon]=find(sumdif==min(sumdif(:)));
    
    
    write_delft_header_mod(fid,sections{bdsec},bd,dd,mm,yyyy,length(prop),2,'Salinity');


    for t = 1:t_aux(2) - t_aux(1)
        
        hycom_file_0=t;
        hycom_file_1=hycom_file_0+1;
        
 %       fprintf('%s\n',['Escrevendo dados - Boundary South' sprintf('%03g',bd)...
 %           ' - Salinidade ---> tempo ' sprintf('%02g',t)])
               
        
        
        %fim dos nomes para remo
        
%         hycom_file_0 = ['Hycom-noda_' num2str(tstep(1)) sprintf('%02g',tstep(2))...
%             sprintf('%02g',tstep(3)) '.nc'];
%         hycom_file_1 = ['Hycom-noda_' num2str(tstepf(1)) sprintf('%02g',tstepf(2))...
%             sprintf('%02g',tstepf(3)) '.nc'];
        
        coef1 = coef1 + 0.01;
        
         [interp_property] = interp_delft_boundary_TS_REF_2019(hycom_file_0,hycom_file_1,...
             'South','salinity',ni,nj,lat,lon,lati,loni,0,coef0,coef1,ti, si, bat,lay,modo,limlat,limlon);
         
%          aux = zeros(33,2,25);
%          aux(5:end,:,:) = interp_property;
%          interp_property = aux;
        
       if tref==0       
          for tinp = 1:25
            
              fprintf(fid,'%15.2f %s',tref,' ');
            
              fprintf_bound(fid,interp_property(:,1,tinp));
              fprintf_bound(fid,interp_property(:,2,tinp));
              tref = tref + 60;
              fprintf(fid,'\n');
          end           
           

       else 
          for tinp = 2:25
            
              fprintf(fid,'%15.2f %s',tref,' ');
            
              fprintf_bound(fid,interp_property(:,1,tinp));
              fprintf_bound(fid,interp_property(:,2,tinp));
              tref = tref + 60;
              fprintf(fid,'\n');
          end
        
       end 
     
        tin = tin + 24;
    end
    

    lati = south_bnd(c_south:c_south+1,1);
    loni = south_bnd(c_south:c_south+1,2);
    tref = 0;
    tin = 1;
    coef0 = 0;
    coef1 = 0;
    write_delft_header_mod(fid,sections{bdsec},bd,dd,mm,yyyy,length(prop),2,'Temperature');

    for t = 1:t_aux(2) - t_aux(1)
        hycom_file_0=t;
        hycom_file_1=hycom_file_0+1;
        
%        fprintf('%s\n',['Escrevendo dados - Boundary South' sprintf('%03g',bd)...
%            ' - Temperatura ---> tempo ' sprintf('%02g',t)])
        
      
        %fim dos nomes para remo
        
%         hycom_file_0 = ['Hycom-noda_' num2str(tstep(1)) sprintf('%02g',tstep(2))...
%             sprintf('%02g',tstep(3)) '.nc'];
%         hycom_file_1 = ['Hycom-noda_' num2str(tstepf(1)) sprintf('%02g',tstepf(2))...
%             sprintf('%02g',tstepf(3)) '.nc'];
        
        coef1 = coef1 + 0.01;
        
         [interp_property] = interp_delft_boundary_TS_REF_2019(hycom_file_0,hycom_file_1,...
             'South','temperature',ni,nj,lat,lon,lati,loni,0,coef0,coef1,ti, si, bat,lay,modo,limlat,limlon);
         
%          aux = zeros(33,2,25);
%          aux(5:end,:,:) = interp_property;
%          interp_property = aux;       
%         
%         
       if tref==0       
          for tinp = 1:25
            
              fprintf(fid,'%15.2f %s',tref,' ');
            
              fprintf_bound(fid,interp_property(:,1,tinp));
              fprintf_bound(fid,interp_property(:,2,tinp));
              tref = tref + 60;
              fprintf(fid,'\n');
          end           
           

       else 
          for tinp = 2:25
            
              fprintf(fid,'%15.2f %s',tref,' ');
            
              fprintf_bound(fid,interp_property(:,1,tinp));
              fprintf_bound(fid,interp_property(:,2,tinp));
              tref = tref + 60;
              fprintf(fid,'\n');
          end
        
       end 

        
        tin = tin + 24;
    end
       
    c_south = c_south + 2;
    
    bdsec = bdsec + 1
    
end


c_east = 1;

for bd = 1:length(east_bnd)/2
    
    lati = east_bnd(c_east:c_east+1,1);
    loni = east_bnd(c_east:c_east+1,2);
    tref = 0;
    tin = 1;
    coef0 = 0;
    coef1 = 0;
    
     
    one=ones(nn,mm2);
    northl=one*lati(1,1);
    difn=abs(lat-northl);
    
    one=ones(nn,mm2);
    eastl=one*loni(1,1);
    dife=abs(lon-eastl);
    
    sumdif=difn+dife;    
    
    [limlat, limlon]=find(sumdif==min(sumdif(:)));
   
    write_delft_header_mod(fid,sections{bdsec},bd,dd,mm,yyyy,length(prop),2,'Salinity');
    
    for t = 1:t_aux(2) - t_aux(1)
        
        hycom_file_0=t;
        hycom_file_1=hycom_file_0+1;
        
%        fprintf('%s\n',['Escrevendo dados - Boundary East' sprintf('%03g',bd)...
%            ' - Salinidade ---> tempo ' sprintf('%02g',t)])
         
        %fim dos nomes para remo
        
%         hycom_file_0 = ['Hycom-noda_' num2str(tstep(1)) sprintf('%02g',tstep(2))...
%             sprintf('%02g',tstep(3)) '.nc'];
%         hycom_file_1 = ['Hycom-noda_' num2str(tstepf(1)) sprintf('%02g',tstepf(2))...
%             sprintf('%02g',tstepf(3)) '.nc'];
        
        coef1 = coef1 + 0.01;
        
         [interp_property] = interp_delft_boundary_TS_REF_2019(hycom_file_0,hycom_file_1,...
            'East','salinity',ni,nj,lat,lon,lati,loni,0,coef0,coef1,ti, si,bat, lay, modo,limlat,limlon);
        
%         aux = zeros(33,2,25);
%         aux(5:end,:,:) = interp_property;
%         interp_property = aux;
%         
       if tref==0       
          for tinp = 1:25
            
              fprintf(fid,'%15.2f %s',tref,' ');
            
              fprintf_bound(fid,interp_property(:,1,tinp));
              fprintf_bound(fid,interp_property(:,2,tinp));
              tref = tref + 60;
              fprintf(fid,'\n');
          end           
           

       else 
          for tinp = 2:25
            
              fprintf(fid,'%15.2f %s',tref,' ');
            
              fprintf_bound(fid,interp_property(:,1,tinp));
              fprintf_bound(fid,interp_property(:,2,tinp));
              tref = tref + 60;
              fprintf(fid,'\n');
          end
        
       end 

        
        tin = tin + 24;
    end
    
    tref = 0;
    tin = 1;
    coef0 = 0;
    coef1 = 0;
    write_delft_header_mod(fid,sections{bdsec},bd,dd,mm,yyyy,length(prop),2,'Temperature');

    for t = 1:t_aux(2) - t_aux(1)
        
        hycom_file_0=t;
        hycom_file_1=hycom_file_0+1;
        
        fprintf('%s\n',['Escrevendo dados - Boundary East' sprintf('%03g',bd)...
            ' - Temperatura ---> tempo ' sprintf('%02g',t)])
         
        
        %fim dos nomes para remo
        
%         hycom_file_0 = ['Hycom-noda_' num2str(tstep(1)) sprintf('%02g',tstep(2))...
%             sprintf('%02g',tstep(3)) '.nc'];
%         hycom_file_1 = ['Hycom-noda_' num2str(tstepf(1)) sprintf('%02g',tstepf(2))...
%             sprintf('%02g',tstepf(3)) '.nc'];
        
        coef1 = coef1 + 0.01;
        
         [interp_property] = interp_delft_boundary_TS_REF_2019(hycom_file_0,hycom_file_1,...
            'East','temperature',ni,nj,lat,lon,lati,loni,0,coef0,coef1,ti, si,bat, lay, modo,limlat,limlon);
        
%         aux = zeros(33,2,25);
%         aux(5:end,:,:) = interp_property;
%         interp_property = aux;
        
       if tref==0       
          for tinp = 1:25
            
              fprintf(fid,'%15.2f %s',tref,' ');
            
              fprintf_bound(fid,interp_property(:,1,tinp));
              fprintf_bound(fid,interp_property(:,2,tinp));
              tref = tref + 60;
              fprintf(fid,'\n');
          end           
           

       else 
          for tinp = 2:25
            
              fprintf(fid,'%15.2f %s',tref,' ');
            
              fprintf_bound(fid,interp_property(:,1,tinp));
              fprintf_bound(fid,interp_property(:,2,tinp));
              tref = tref + 60;
              fprintf(fid,'\n');
          end
        
       end 

        
        tin = tin + 24;
    end
    
    c_east = c_east + 2;
    
    bdsec = bdsec + 1
end

if ~west
  return
end

c_west = 1;

for bd = 1:length(west_bnd)/2
    
    lati = west_bnd(c_west:c_west+1,1);
    loni = west_bnd(c_west:c_west+1,2);
    tref = 0;
    tin = 1;
    coef0 = 0;
    coef1 = 0;
    
    one=ones(nn,mm2);
    northl=one*lati(1,1);
    difn=abs(lat-northl);
    
    one=ones(nn,mm2);
    eastl=one*loni(1,1);
    dife=abs(lon-eastl);
    
    sumdif=difn+dife;    
    
    [limlat, limlon]=find(sumdif==min(sumdif(:)));
    
    write_delft_header_mod(fid,sections{bdsec},bd,dd,mm,yyyy,length(prop),2,'Salinity');
    
    for t = 1:t_aux(2) - t_aux(1)
        
        hycom_file_0=t;
        hycom_file_1=hycom_file_0+1;
        
        
        fprintf('%s\n',['Escrevendo dados - Boundary West' sprintf('%03g',bd)...
            ' - Salinidade ---> tempo ' sprintf('%02g',t)])
        %fim dos nomes para remo
        
%         hycom_file_0 = ['Hycom-noda_' num2str(tstep(1)) sprintf('%02g',tstep(2))...
%             sprintf('%02g',tstep(3)) '.nc'];
%         hycom_file_1 = ['Hycom-noda_' num2str(tstepf(1)) sprintf('%02g',tstepf(2))...
%             sprintf('%02g',tstepf(3)) '.nc'];
        
        coef1 = coef1 + 0.01;
        
         [interp_property] = interp_delft_boundary_TS_REF_2019(hycom_file_0,hycom_file_1,...
            'West','salinity',ni,nj,lat,lon,lati,loni,0,coef0,coef1,ti, si,bat, lay, modo,limlat,limlon);
        
%         aux = zeros(33,2,25);
%         aux(5:end,:,:) = interp_property;
%         interp_property = aux;
%         
        for tinp = 1:24
            
            fprintf(fid,'%15.2f %s',tref,' ');
            
            fprintf_bound(fid,interp_property(:,1,tinp));
            fprintf_bound(fid,interp_property(:,2,tinp));
            tref = tref + 60;
            fprintf(fid,'\n');
        end
        tin = tin + 24;
    end
    
        tref = 0;
    tin = 1;
    coef0 = 0;
    coef1 = 0;
    write_delft_header_mod(fid,sections{bdsec},bd,dd,mm,yyyy,length(prop),2,'Temperature');

    for t = 1:t_aux(2) - t_aux(1)
        
        hycom_file_0=t
        hycom_file_1=hycom_file_0+1
        
        fprintf('%s\n',['Escrevendo dados - Boundary West' sprintf('%03g',bd)...
            ' - Temperatura ---> tempo ' sprintf('%02g',t)])
        
        
        %fim dos nomes para remo
        
%         hycom_file_0 = ['Hycom-noda_' num2str(tstep(1)) sprintf('%02g',tstep(2))...
%             sprintf('%02g',tstep(3)) '.nc'];
%         hycom_file_1 = ['Hycom-noda_' num2str(tstepf(1)) sprintf('%02g',tstepf(2))...
%             sprintf('%02g',tstepf(3)) '.nc'];
        
        coef1 = coef1 + 0.01;
        
         [interp_property] = interp_delft_boundary_TS_REF_2019(hycom_file_0,hycom_file_1,...
            'West','temperature',ni,nj,lat,lon,lati,loni,0,coef0,coef1,ti, si,bat, lay, modo,limlat,limlon);
        
%         aux = zeros(33,2,25);
%         aux(5:end,:,:) = interp_property;
%         interp_property = aux;
        
        for tinp = 1:24
            
            fprintf(fid,'%15.2f %s',tref,' ');
            
            fprintf_bound(fid,interp_property(:,1,tinp));
            fprintf_bound(fid,interp_property(:,2,tinp));
            tref = tref + 60;
            fprintf(fid,'\n');
        end
        tin = tin + 24;
    end
    
    c_west = c_west + 2;
    
    bdsec = bdsec + 1

    
end 

fclose(fid);