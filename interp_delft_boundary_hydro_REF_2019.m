function [interp_property] = interp_delft_boundary_mod_climatolico_sigma(hycom_file_0,hycom_file_1,...
    bnd_id,type,ni,nj,lat,lon,lati,loni,tide,coef0,coef1,  depth_i, old_bat,ui,vii, varargin)
%askljdflkasjhdflkjshdfgjksdhkf
% This subroutine controls the interpolation for delft boundary conditions,
% it is based on hycom ncoda operational results. The interpolation is set
% to one hour time step and for three kinds of boundary conditions Reimann,
% Current and Elevation, the last ones for the Dirichlet boundary
% conditions.
%
% Usage: First the user must to define the ncoda's hycom file to perform
% the interpolation (hycom_file_0 and hycom_file_1), so the user is
% requested for the location of the boundary (bnd_id -
% 'North','South','East' or 'South'), its kind ('Riemann','Current' or
% 'Elevation'), the specifications of number of i cells (ni) and j cells
% (nj), location of the boundary coordinates (lon_lat, provided by hycom)
% and location where it must be interpolated (lon_lat_int), at the last the
% user need to provide the time series of tidal heights at the A and B 
% points of the boundary section and the coefficients, in case of
% 'Elevation' and 'Current' this works as slow start to this variables, for
% 'Riemann' this works as smoothing of neighbouring points (see at delft3d
% manual). Additionally in case of Riemann boundary condition the user must
% specify the bathymetry at the boundary where the intporlation is
% performed.


modo=varargin{3};
lt=varargin{4};
ll=varargin{5};

new_wl=varargin{6};

if new_wl
  ltG = varargin{7};
  llG = varargin{8};
  tin = varargin{9};
  WL=varargin{10};
  curv_interp=varargin{11};
  if curv_interp
      thetarad=varargin{12};
      alfarad=varargin{13};
  end 
end

if ~new_wl
  curv_interp=varargin{7};
  if curv_interp
      thetarad=varargin{8};
      alfarad=varargin{9};
  end 
end


range=2;


% dep=ncread(hycom_file_0,'Depth');
% u=ncread(hycom_file_0, 'u');
% 
% if strcmp(modo,'REMO')
%  dep(1)=0;
%  lats=ncread(hycom_file_0, 'Latitude'); 
%  lons=ncread(hycom_file_0, 'Longitude');
%  [latz lonz]=meshgrid(lats,lons);
%  latz=permute(latz, [2 1]);
%  lonz=permute(lonz, [2 1]);
%  [nis njs]=size(latz);
% else
%  lats=ncread(hycom_file_0, 'Latitude'); 
%  lons=ncread(hycom_file_0, 'Longitude');
%  latz=permute(lats, [2 1]);
%  lonz=permute(lons, [2 1]);
%  [nis njs]=size(latz);
% end 
 
 

%disp('LOADING FEATURES')
load transpose;

%disp('loaded')


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5 bat gecbo from delft
% ncfileb='trim-tst.nc';
% name='tst';
% bat=delftbat2routine(sprintf([name '.dep']),sprintf([name '.mdf']));
% X = nc_varget(ncfileb,'longitude');
% Y = nc_varget(ncfileb,'latitude');
% X = rot90(X);   % ROT90 E FLIPUD é iqual permute
% X = flipud(X);    
% Y = rot90(Y);     
% Y = flipud(Y);   
% 
% X(isnan(X))=0;
% Y(isnan(Y))=0;
% 
% bat=griddata(Y(1:425,1:252),X(1:425, 1:252),bat,latz,lonz);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%teste bat HYCOM
% nc='depth_GLBa0.08_09.nc';
% bati=ncread(nc,'bathymetry');
% latet=ncread(nc,'Latitude');
% lonet=ncread(nc,'Longitude');
% bat=griddata(latet,lonet, bati, latz,lonz);
%%%%%%%%%%%%%%%%%%%%%55hycom

%%%%%%%%%%%%%%%  BAT z layers from hycom
%   for i= 1:nis
%   for j= 1:njs
%    ub=u(i,j,:);
%    ub=ub<10000;
%    c=sum(ub);
%    if c ~= 0
%               aux_dep = dep(c);
%               bat(i,j) = aux_dep;
%    else
%               bat(i,j) = -999;
%    end
%   end
%   end

switch type
    case 'Current'
        
      if ~curv_interp
          
        switch bnd_id
            case {'North' , 'South'}
                u_normal_0 = vii(:,:,:,hycom_file_0);
                u_normal_1 = vii(:,:,:,hycom_file_1);
%                u_parallel_id='u';
            otherwise
                u_normal_0 = ui(:,:,:,hycom_file_0);
                u_normal_1 = ui(:,:,:,hycom_file_1);                
%                u_parallel_id='v';
        end
      else 
            u_normal_0 = ui(:,:,:,hycom_file_0);
            u_normal_1 = ui(:,:,:,hycom_file_1);       
          
            v_normal_0 = vii(:,:,:,hycom_file_0);
            v_normal_1 = vii(:,:,:,hycom_file_1); 
            
      end 
      
      
        lay=varargin{2};        
        interp_property = zeros(numb,2,25);    
        
        
        for k = 0:numb-1
%            u_normal_0 = nc_varget(hycom_file_0,u_normal_id,[0 k 0 0],[1 1 ni nj]);
%            u_normal_1 = nc_varget(hycom_file_1,u_normal_id,[0 k 0 0],[1 1 ni nj]);
            
            
%            u_normal_teste = nc_varget(hycom_file_0,u_normal_id,[0 0 0 0],[1 33 ni nj]);
%            assignin('base', 'u_teste',u_normal_teste);
            
           
            
%            u_parallel_0 = nc_varget(hycom_file_0,u_parallel_id,[0 k 0 0],[1 1 ni nj]);
%            u_parallel_1 = nc_varget(hycom_file_1,u_parallel_id,[0 k 0 0],[1 1 ni nj]);            
%             switch bnd_id
%                 case 'North'
%                     u_normal_0 = u_normal_0(end,:);
%                     u_normal_1 = u_normal_1(end,:);
%                 case 'East'
%                     u_normal_0 = u_normal_0(:,end);
%                     u_normal_1 = u_normal_1(:,end);
%                 case 'South'
%                     u_normal_0 = u_normal_0(1,:);
%                     u_normal_1 = u_normal_1(1,:);
%             end

         if ~ curv_interp         
            u_normal_0_int = griddata(latfm(lt-range:lt+range,ll-range:ll+range),lonfm(lt-range:lt+range,ll-range:ll+range),squeeze(u_normal_0(k+1,lt-range:lt+range,ll-range:ll+range)),lati,loni);
           
            u_normal_1_int = griddata(latfm(lt-range:lt+range,ll-range:ll+range),lonfm(lt-range:lt+range,ll-range:ll+range),squeeze(u_normal_1(k+1,lt-range:lt+range,ll-range:ll+range)),lati,loni);            
 
         else 
             
            u_normal_0_int = griddata(latfm(lt-range:lt+range,ll-range:ll+range),lonfm(lt-range:lt+range,ll-range:ll+range),squeeze(u_normal_0(k+1,lt-range:lt+range,ll-range:ll+range)),lati,loni);
           
            u_normal_1_int = griddata(latfm(lt-range:lt+range,ll-range:ll+range),lonfm(lt-range:lt+range,ll-range:ll+range),squeeze(u_normal_1(k+1,lt-range:lt+range,ll-range:ll+range)),lati,loni);
            
            v_normal_0_int = griddata(latfm(lt-range:lt+range,ll-range:ll+range),lonfm(lt-range:lt+range,ll-range:ll+range),squeeze(v_normal_0(k+1,lt-range:lt+range,ll-range:ll+range)),lati,loni);
           
            v_normal_1_int = griddata(latfm(lt-range:lt+range,ll-range:ll+range),lonfm(lt-range:lt+range,ll-range:ll+range),squeeze(v_normal_1(k+1,lt-range:lt+range,ll-range:ll+range)),lati,loni); 
            
              switch bnd_id
                case {'North' , 'South'}           
                    
                    
                  u_normal_0_int = u_normal_0_int * cos(thetarad);
                              
                  u_normal_1_int = u_normal_1_int * cos(thetarad);
            
                  v_normal_0_int = v_normal_0_int * cos(alfarad);
                              
                  v_normal_1_int = v_normal_0_int * cos(alfarad);
                  
                  
                  
                  otherwise
                      
                  u_normal_0_int = u_normal_0_int * cos(alfarad);
            
                  u_normal_1_int = u_normal_1_int * cos(alfarad);
            
                  v_normal_0_int = v_normal_0_int * cos(thetarad);
            
                  v_normal_1_int = v_normal_0_int * cos(thetarad);                      
                      
              end 
            
            u_normal_0_int = u_normal_0_int + v_normal_0_int;
            
            u_normal_1_int = u_normal_1_int + v_normal_1_int;
            
         end
         
            
            u_normal_0_int(isnan(u_normal_0_int)) = 0;
            u_normal_1_int(isnan(u_normal_1_int)) = 0;
                    
            
            interp_u_normal(1,1:25) = linspace(u_normal_0_int(1),u_normal_1_int(1),25);
            interp_u_normal(2,1:25) = linspace(u_normal_0_int(2),u_normal_1_int(2),25);
            
                        
           interp_property(k+1,1:2,1:25) = interp_u_normal;
        end                                      
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        
    case 'WaterLevel'

        interp_property = zeros(2,25);    

        ssh_0_interp = griddata(latfm(lt-range:lt+range,ll-range:ll+range),lonfm(lt-range:lt+range,ll-range:ll+range),ssh(lt-range:lt+range,ll-range:ll+range),lati,loni);
        ssh_1_interp = griddata(latfm(lt-range:lt+range,ll-range:ll+range),lonfm(lt-range:lt+range,ll-range:ll+range),ssh(lt-range:lt+range,ll-range:ll+range),lati,loni);
        
        ssh_0_interp(isnan(ssh_0_interp)) = 0;
        ssh_1_interp(isnan(ssh_1_interp)) = 0;
        
        interp_property(1,1:25) = linspace(ssh_0_interp(1),ssh_1_interp(1),25);
        interp_property(2,1:25) = linspace(ssh_0_interp(2),ssh_1_interp(2),25);
        
        if new_wl
            
           ixi=1;  
           for tix=tin:tin+24
             tide_great=griddata(latG(ltG-range:ltG+range,llG-range:llG+range),lonG(ltG-range:ltG+range,llG-range:llG+range),WL(ltG-range:ltG+range,llG-range:llG+range, tix),lati,loni);
             
              
             if isnan(tide_great(1))
               tide_great(1)=griddata(latG(ltG-range:ltG+range,llG-range:llG+range),lonG(ltG-range:ltG+range,llG-range:llG+range),WL(ltG-range:ltG+range,llG-range:llG+range,...
                 tix),lati(1),loni(1), 'nearest');   
            
             end
   
             if isnan(tide_great(2))
               tide_great(2)=griddata(latG(ltG-range:ltG+range,llG-range:llG+range),lonG(ltG-range:ltG+range,llG-range:llG+range),WL(ltG-range:ltG+range,llG-range:llG+range,...
                 tix),lati(2),loni(2), 'nearest'); 
             end
             
            interp_property(1,ixi) =  interp_property(1,ixi) + tide_great(1);
            interp_property(2,ixi) =  interp_property(2,ixi) + tide_great(2);
            
            ixi=ixi+1;           
           end          
                  
        else
           interp_property = interp_property + tide;        
        end
        
        
 %       if coef1 <= 1
  %          coef = linspace(coef0,coef1,25);
  %          interp_property(1,1:25) = interp_property(1,1:25).*coef;
  %          interp_property(2,1:25) = interp_property(2,1:25).*coef;
   %     end
        
    case 'Riemann'
        %get parameters from varargin
      if ~curv_interp
          
        switch bnd_id
            case {'North' , 'South'}
                u_normal_0 = vii(:,:,:,hycom_file_0);
                u_normal_1 = vii(:,:,:,hycom_file_1);
%                u_parallel_id='u';
            otherwise
                u_normal_0 = ui(:,:,:,hycom_file_0);
                u_normal_1 = ui(:,:,:,hycom_file_1);                
%                u_parallel_id='v';
        end
      else 
            u_normal_0 = ui(:,:,:,hycom_file_0);
            u_normal_1 = ui(:,:,:,hycom_file_1);       
          
            v_normal_0 = vii(:,:,:,hycom_file_0);
            v_normal_1 = vii(:,:,:,hycom_file_1); 
            
      end 
      
        lay=varargin{2};        
        
        interp_property = zeros(numb,2,25);    
        
        disp('ELEVATOIN..........................................')
        ssh_0_interp = griddata(latfm(lt-range:lt+range,ll-range:ll+range),lonfm(lt-range:lt+range,ll-range:ll+range),ssh(lt-range:lt+range,ll-range:ll+range),lati,loni);
        ssh_1_interp = griddata(latfm(lt-range:lt+range,ll-range:ll+range),lonfm(lt-range:lt+range,ll-range:ll+range),ssh(lt-range:lt+range,ll-range:ll+range),lati,loni); 
        
       
        ssh_0_interp(isnan(ssh_0_interp)) = 0;
        ssh_1_interp(isnan(ssh_1_interp)) = 0;
        
        eta(1,1:25) = linspace(ssh_0_interp(1),ssh_1_interp(1),25);
        eta(2,1:25) = linspace(ssh_0_interp(2),ssh_1_interp(2),25);
        
        
        if new_wl
            
           ixi=1;  
           for tix=tin:tin+24
             tide_great=griddata(latG(ltG-range:ltG+range,llG-range:llG+range),lonG(ltG-range:ltG+range,llG-range:llG+range),WL(ltG-range:ltG+range,llG-range:llG+range, tix),lati,loni);
            
             
             if isnan(tide_great(1))
               tide_great(1)=griddata(latG(ltG-range:ltG+range,llG-range:llG+range),lonG(ltG-range:ltG+range,llG-range:llG+range),WL(ltG-range:ltG+range,llG-range:llG+range,...
                 tix),lati(1),loni(1), 'nearest');   
             end
   
             if isnan(tide_great(2))
               tide_great(2)=griddata(latG(ltG-range:ltG+range,llG-range:llG+range),lonG(ltG-range:ltG+range,llG-range:llG+range),WL(ltG-range:ltG+range,llG-range:llG+range,...
                 tix),lati(2),loni(2), 'nearest'); 
             end
             
            eta(1,ixi) =  eta(1,ixi) + tide_great(1);
            eta(2,ixi) =  eta(2,ixi) + tide_great(2);
            
            ixi=ixi+1;           
           end          
                  
        else
           eta = eta + tide;        
        end
                
        if old_bat
          depth_i = griddata(latfm(lt-range:lt+range,ll-range:ll+range),lonfm(lt-range:lt+range,ll-range:ll+range),bat(lt-range:lt+range,ll-range:ll+range),lati,loni);
        end
        
  depth_i(1);
  depth_i(2);
  
        depthr(1,1:25) = depth_i(1);
        depthr(2,1:25) = depth_i(2);  
                
        for k = 0:numb-1            
            
                      
          if ~ curv_interp         
            u_normal_0_int = griddata(latfm(lt-range:lt+range,ll-range:ll+range),lonfm(lt-range:lt+range,ll-range:ll+range),squeeze(u_normal_0(k+1,lt-range:lt+range,ll-range:ll+range)),lati,loni);
           
            u_normal_1_int = griddata(latfm(lt-range:lt+range,ll-range:ll+range),lonfm(lt-range:lt+range,ll-range:ll+range),squeeze(u_normal_1(k+1,lt-range:lt+range,ll-range:ll+range)),lati,loni);            
 
         else 
             
            u_normal_0_int = griddata(latfm(lt-range:lt+range,ll-range:ll+range),lonfm(lt-range:lt+range,ll-range:ll+range),squeeze(u_normal_0(k+1,lt-range:lt+range,ll-range:ll+range)),lati,loni)
           
            u_normal_1_int = griddata(latfm(lt-range:lt+range,ll-range:ll+range),lonfm(lt-range:lt+range,ll-range:ll+range),squeeze(u_normal_1(k+1,lt-range:lt+range,ll-range:ll+range)),lati,loni);
            
            v_normal_0_int = griddata(latfm(lt-range:lt+range,ll-range:ll+range),lonfm(lt-range:lt+range,ll-range:ll+range),squeeze(v_normal_0(k+1,lt-range:lt+range,ll-range:ll+range)),lati,loni)
           
            v_normal_1_int = griddata(latfm(lt-range:lt+range,ll-range:ll+range),lonfm(lt-range:lt+range,ll-range:ll+range),squeeze(v_normal_1(k+1,lt-range:lt+range,ll-range:ll+range)),lati,loni); 
            
              switch bnd_id
                case {'North' , 'South'}            
                  u_normal_0_int = u_normal_0_int * cos(thetarad);
            
                  u_normal_1_int = u_normal_1_int * cos(thetarad);
            
                  v_normal_0_int = v_normal_0_int * cos(alfarad);
            
                  v_normal_1_int = v_normal_0_int * cos(alfarad);
                  
                  otherwise
                      
                  u_normal_0_int = u_normal_0_int * cos(alfarad)
            
                  u_normal_1_int = u_normal_1_int * cos(alfarad);
            
                  v_normal_0_int = v_normal_0_int * cos(thetarad)
            
                  v_normal_1_int = v_normal_0_int * cos(thetarad);                      
                      
              end 
            
            u_normal_0_int = u_normal_0_int + v_normal_0_int
            
            u_normal_1_int = u_normal_1_int + v_normal_1_int;
            
            thhhh
                         
         end

    
                                     
           u_normal_0_int(isnan(u_normal_0_int)) = 0;
           u_normal_1_int(isnan(u_normal_1_int)) = 0;
                    
            
           interp_u_normal(1,1:25) = linspace(u_normal_0_int(1),u_normal_1_int(1),25)*coef0;
           interp_u_normal(2,1:25) = linspace(u_normal_0_int(2),u_normal_1_int(2),25)*coef1;
   
         if strcmp(bnd_id,'North') || strcmp(bnd_id,'East')
                interp_property(k+1,1:2,1:25) = interp_u_normal - eta.*sqrt(9.81./depthr);
           else
                interp_property(k+1,1:2,1:25) = interp_u_normal + eta.*sqrt(9.81./depthr);
         end
         
        end
        
    case 'salinity'
        
         lay=varargin{2};
         interp_property = zeros(numb,2,2);
        
        
         
        for k = 0:numb-1
            salinity_0 = si;
            salinity_1 = si;
            

            
%             switch bnd_id
%                 case 'North'
%                     salinity_bound_0 = salinity_0(end,:);
%                     salinity_bound_1 = salinity_1(end,:);
%                 case 'East'
%                     salinity_bound_0 = salinity_0(:,end);
%                     salinity_bound_1 = salinity_1(:,end);
%                 case 'South'
%                     salinity_bound_0 = salinity_0(1,:);
%                     salinity_bound_1 = salinity_1(1,:);
%             end
            
            salinity_0_int = griddata(latfm,lonfm,squeeze(salinity_0(k+1,:,:)),lati,loni);
            salinity_1_int = griddata(latfm,lonfm,squeeze(salinity_1(k+1,:,:)),lati,loni);
                        
            salinity_0_int(isnan(salinity_0_int)) = 0;
            salinity_1_int(isnan(salinity_1_int)) = 0;
            
            interp_salinity(1,1:2) = linspace(salinity_0_int(1),salinity_1_int(1),2);
            interp_salinity(2,1:2) = linspace(salinity_0_int(2),salinity_1_int(2),2);
            
            
            interp_property(k+1,1:2,1:2) = interp_salinity;            
        end 
        
    case 'temperature'
        
         lay=varargin{2};
         interp_property = zeros(numb,2,2);
        
         
         
        for k = 0:numb-1
            temperature_0 = ti;
            temperature_1 = ti;
            
%             switch bnd_id
%                 case 'North'
%                     temperature_bound_0 = temperature_0(end,:);
%                     temperature_bound_1 = temperature_1(end,:);
%                 case 'East'
%                     temperature_bound_0 = temperature_0(:,end);
%                     temperature_bound_1 = temperature_1(:,end);
%                 case 'South'
%                     temperature_bound_0 = temperature_0(1,:);
%                     temperature_bound_1 = temperature_1(1,:);
%             end
            
            temperature_0_int = griddata(latfm,lonfm,squeeze(temperature_0(k+1,:,:)),lati,loni);
            temperature_1_int = griddata(latfm,lonfm,squeeze(temperature_1(k+1,:,:)),lati,loni);

            temperature_0_int(isnan(temperature_0_int)) = 0;
            temperature_1_int(isnan(temperature_1_int)) = 0;
            
            interp_temperature(1,1:2) = linspace(temperature_0_int(1),temperature_1_int(1),2);
            interp_temperature(2,1:2) = linspace(temperature_0_int(2),temperature_1_int(2),2);
            
            
            interp_property(k+1,1:2,1:2) = interp_temperature;
        
        
        end
                     
        
end
end
