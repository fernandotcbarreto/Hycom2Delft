cd C:\Users\fernando\Desktop\meso_scale

fid=fopen('meso.bnd')

ncfile='trim-meso.nc';

lat=ncread(ncfile, 'latitude');
lon=ncread(ncfile, 'longitude');

sizlat=size(lat);

n=0 
tline='a'
while ischar(tline)
   n=n+1;
   tline=fgetl(fid);
end
n=n-1

frewind(fid)
i=0

tline='a'
q=0;
i=0
tline=fgetl(fid);

%COLOCAR NA ORDEM Q ESTIVER NO .bnd 

while i~=n
    if tline(1:1) == 'N'
       fidd=fopen('north_acu.dat', 'wt'); 
       while tline(1:1) == 'N'
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
           fprintf(fidd,'%f %f \n', lat(n1,m1), lon(n1,m1));
           fprintf(fidd,'%f %f \n', lat(n2,m2), lon(n2,m2));       
           i=i+1;          
           tline=fgetl(fid)    
       end 
    end     
    q=0
   
    if tline(1:1) == 'S'
        fidd=fopen('south_acu.dat', 'wt');
        while tline(1:1) == 'S'
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
           fprintf(fidd,'%f %f \n', lat(n1,m1), lon(n1,m1));
           fprintf(fidd,'%f %f \n', lat(n2,m2), lon(n2,m2));       
           i=i+1;
           tline=fgetl(fid)    
       end 
    end
    
    q=0
   
    if tline(1:1) == 'E'
        fidd=fopen('east_acu.dat', 'wt');
        while tline(1:1) == 'E'
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
           fprintf(fidd,'%f %f \n', lat(n1,m1), lon(n1,m1));
           fprintf(fidd,'%f %f \n', lat(n2,m2), lon(n2,m2));       

           i=i+1;
          tline=fgetl(fid)    
       end 
    end 
    
    q=0
    if tline(1:1) == 'W'
        fidd=fopen('west_acu.dat', 'wt');
        while tline(1:1) == 'W'
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

           
           fprintf(fidd,'%f %f \n', lat(n1,m1), lon(n1,m1));
           fprintf(fidd,'%f %f \n', lat(n2,m2), lon(n2,m2));       

           i=i+1;
           tline=fgetl(fid)    
       end 
    end
       
    
end 
    