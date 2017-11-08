% This script computes both sides of the thermal wind relation using NCEP
% Reanalysis climatological data. The results are plotted as a function of latitude
% and are averaged over a specified longitudinal range (Mike Byrne, March 20 2012).
%------------------------------------------------------
%USER DEFINED INPUT VARIABLES

%INPUT DATA FILE NAMES
fname1='theta_500.nc';
fname2='u_400.nc'; %For the wind, put the lower pressure level first
fname3='u_600.nc';

%INPUT PRESSURE LEVEL (mb) FOR TEMPERATURE CROSS-SECTION

p = 500;
p_0 = 1000; % The reference pressure level used to calculate pot. temp.

%INPUT PRESSURE DIFFERENCE (mb) BETWEEN TWO WIND LEVELS

dp=-200;

%INPUT LONGITUDE LIMITS FOR CROSS-SECTION
%Quantities are averaged between the specified longitude limits

lon_limit_1=110;
lon_limit_2=150;

%------------------------------------------------------

%Read NetCDF files

nc=netcdf.open(fname1,'nowrite');
%ncdump(nc);

nc2=netcdf.open(fname2,'nowrite');
%ncdump(nc2);

nc3=netcdf.open(fname3,'nowrite');
%ncdump(nc3);

%------------------------
% We now read in the data from the netdf files. 
% First, the potential temperature:

% theta = nc{'pottmp'}(:); % Dimensions are (time,level,lat,lon)
theta = ncread(fname1, 'pottmp')';

%theta = squeeze(mean(theta,1)); % Avg in time
% Now the winds:

% uwind_lev_1 = nc2{'uwnd'}(:); % Dimensions are (time,level,lat,lon)
% uwind_lev_2 = nc3{'uwnd'}(:);

uwind_lev_1 = ncread(fname2, 'uwnd')';
uwind_lev_2 = ncread(fname3, 'uwnd')';

% And the axes:
% lat = nc{'lat'}(:);
% lon = nc{'lon'}(:);
lat = ncread(fname1, 'lat');
lon = ncread(fname1, 'lon');

% Find the indices corresponding to the longitude limits:
[val lon_index_1] = min(abs(lon - lon_limit_1));
[val lon_index_2] = min(abs(lon - lon_limit_2));

% Now average the data between these indices:
theta = squeeze(mean(theta(:,lon_index_1:lon_index_2),2));
uwind_lev_1 = squeeze(mean(uwind_lev_1(:,lon_index_1:lon_index_2),2));
uwind_lev_2 = squeeze(mean(uwind_lev_2(:,lon_index_1:lon_index_2),2));

% Get the winds between the latitude points:
uwind_lev_1_mid = uwind_lev_1(1:end-1) + ( uwind_lev_1(2:end)-uwind_lev_1(1:end-1) )/2;
uwind_lev_2_mid = uwind_lev_2(1:end-1) + ( uwind_lev_2(2:end)-uwind_lev_2(1:end-1) )/2;
lat_mid = lat(1:end-1) + ( lat(2:end) - lat(1:end-1) )/2;

% Define some parameters for the thermal wind relation:
f = 2*(2*pi/(24*60*60))*sin(lat_mid.*pi/180); % Coriolis parameter, units are in 1/s
R = 287.06; % Gas constant, units are J/kg/K
kappa = 2/7; % kappa = R/c_p

% Calculate the pot. temperature gradient with respect to latitude:
% 1 degree latitude = 69 miles = 111,000 meters
theta_grad = ( theta(2:end) - theta(1:end-1) ) ./ ( (lat(2) - lat(1))*111000 );

% Calculate the RHS of the thermal wind relation:
tw_rhs=R*((p/p_0)^kappa)*theta_grad./f/p;

% Now calculate the LHS:

tw_lhs=(1/dp)*(uwind_lev_1_mid - uwind_lev_2_mid);

resid=(tw_lhs-tw_rhs);

%-------------------------
%PLOTS

%Plot temperature cross-section
figure(1);
plot(lat,theta);
xlabel('Latitude (deg. N)');
ylabel('Pot. Temperature (deg. C)');
title([num2str(p),'mb surface for January, ',num2str(lon_limit_1),'deg E to ' ,num2str(lon_limit_2),'deg E']);

%Plot temperature gradient at each latitude
figure(2);
plot(lat_mid,theta_grad);
xlabel('Latitude (deg. N)');
ylabel('Pot. Temperature gradient (deg. C / meter)');
title([num2str(p),'mb surface for January, ',num2str(lon_limit_1),'deg E to ' ,num2str(lon_limit_2),'deg E']);

%Plot temperature gradient side of thermal wind
figure(3);
plot(lat_mid,tw_rhs);
xlabel('Latitude (deg. N)');
ylabel('Pot. Temperature gradient side of TW (m/s/mb)');
title([num2str(p),'mb surface for January, ',num2str(lon_limit_1),'deg E to ' ,num2str(lon_limit_2),'deg E']);

%Plot vertical wind shear side of thermal wind
figure(4);
plot(lat_mid,tw_lhs);
xlabel('Latitude (deg. N)');
ylabel('Vertical wind shear side of TW (m/s/mb)');
title([num2str(p),'mb surface for January, ',num2str(lon_limit_1),'deg E to ' ,num2str(lon_limit_2),'deg E']);
%Plot residual
figure(5);
plot(lat_mid,resid);
xlabel('Latitude (deg. N)');
ylabel('Residual (m/s/mb)');
title([num2str(p),'mb surface for January, ',num2str(lon_limit_1),'deg E to ' ,num2str(lon_limit_2),'deg E']);
