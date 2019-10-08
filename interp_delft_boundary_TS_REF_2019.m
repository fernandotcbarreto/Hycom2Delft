function [interp_property] = interp_delft_boundary_mod_climatolico_sigma(hycom_file_0,hycom_file_1,...
    bnd_id,type,ni,nj,lat,lon,lati,loni,tide,coef0,coef1,ti, si,varargin)
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
%  
 

% disp('LOADING FEATURES')
 load transpose_ts;


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
 

    case 'salinity'
        
         lay=varargin{2};
         
         interp_property = zeros(numb,2,25);
         salinity_0 = si(:,:,:,hycom_file_0);
         salinity_1 = si(:,:,:,hycom_file_1);        
        
         
        for k = 0:numb-1
                      
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
            
            salinity_0_int = griddata(latfm(lt-range:lt+range,ll-range:ll+range),lonfm(lt-range:lt+range,ll-range:ll+range),squeeze(salinity_0(k+1,lt-range:lt+range,ll-range:ll+range)),lati,loni);
            salinity_1_int = griddata(latfm(lt-range:lt+range,ll-range:ll+range),lonfm(lt-range:lt+range,ll-range:ll+range),squeeze(salinity_1(k+1,lt-range:lt+range,ll-range:ll+range)),lati,loni);
                                    
            salinity_0_int(isnan(salinity_0_int)) = 0;
            salinity_1_int(isnan(salinity_1_int)) = 0;
            
            interp_salinity(1,1:25) = linspace(salinity_0_int(1),salinity_1_int(1),25);
            interp_salinity(2,1:25) = linspace(salinity_0_int(2),salinity_1_int(2),25);
            
            
            interp_property(k+1,1:2,1:25) = interp_salinity;            
        end 
        
    case 'temperature'
        
         lay=varargin{2};
         interp_property = zeros(numb,2,25);
         temperature_0 = ti(:,:,:,hycom_file_0);
         temperature_1 = ti(:,:,:,hycom_file_1);        
         
         
        for k = 0:numb-1
            
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
            
            temperature_0_int = griddata(latfm(lt-range:lt+range,ll-range:ll+range),lonfm(lt-range:lt+range,ll-range:ll+range),squeeze(temperature_0(k+1,lt-range:lt+range,ll-range:ll+range)),lati,loni);
            temperature_1_int = griddata(latfm(lt-range:lt+range,ll-range:ll+range),lonfm(lt-range:lt+range,ll-range:ll+range),squeeze(temperature_1(k+1,lt-range:lt+range,ll-range:ll+range)),lati,loni);

            temperature_0_int(isnan(temperature_0_int)) = 0;
            temperature_1_int(isnan(temperature_1_int)) = 0;
            
            interp_temperature(1,1:25) = linspace(temperature_0_int(1),temperature_1_int(1),25);
            interp_temperature(2,1:25) = linspace(temperature_0_int(2),temperature_1_int(2),25);
            
            
            interp_property(k+1,1:2,1:25) = interp_temperature;
        
        
        end
                     
        
end
end
