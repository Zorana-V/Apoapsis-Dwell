clc; clear;
%-----------------------------------------------------------------------%
% Task 1: Using Kepler Solver to compute the orbit at five minute       % 
% intervals for one full period at varying values of lomega and bomega  %
%-----------------------------------------------------------------------%


% Given Values:
tau = 24; % period of the orbit [hours]
t0 = 0; %initial time [hours]
deltat = 5; % time increment [mins]
e = 0.25; % eccentricity of the orbit [unitless]
inc = 63.4; % orbital inclination [degrees]
lomega = [90 270]; % arguement of the periapsis [degrees]
bomega = [0 45 90 135 180]; % longitude of the ascending node [degrees]
mu = 398600; % gravitational parameter [km^3 * s^-2]
nu0 = 0; % initial true anomaly [degrees]
OmegaE = 2*pi/(60*60*24); % rotational rate of the Earth [rad * s^-1]
marker = ["s"; "d"; "v"; ">"]; % string of markers used in geotrack plots
az = -35; % for the 3d plot view
el = 40; %for the 3d plot view


% Unit Conversions of Given Values:
tf = tau*60*60; % final time [sec]
t0 = t0*60*60; % initial time [sec]
deltat = deltat*60; % time increment [sec]
inc = inc*(pi/180); % orbital inclination [rad]
lomega = lomega*(pi/180); % arguement of the periapsis [rad]
bomega = bomega*(pi/180); % longitude of the ascending node [rad]
nu0 = nu0*(pi/180); % initial true anomaly [rad]


% Tabulation of Important Parameters:
N = ((tf - t0)/deltat)+1; % number of data points we will have
t = linspace(t0,tf,N).'; % column vector of time increments [sec]
a = (((tf/(2*pi))^2)*mu)^(1/3); % semi-major axis [km]


% Initialization of Vector to Decrease Cumputational Expense of the Code
nu = zeros(N,1);
rECI = zeros(N,3*length(lomega)*length(bomega)); % initializing the position vector
vECI = zeros(N,3*length(lomega)*length(bomega)); % initializing the velocity vector
L = length(lomega)*length(bomega); % number of different combinations of lomega and b omega
thetae = zeros(N,L); % initializing a longitude vector
phie = zeros(N,L); %initializing a latitude vector
long = 0; %counter variable for for loops


% Setting the Initial Conditions at t = t_o
nu(1) = nu0; % true anomaly at t_o
oe = [a; e; 0; inc; 0; nu(1)]; % orbital elements at initial time t_o


% For loops which will find the [x, y, z] position of the satelite at the initial 
% time on each orbit, runs 10 iterations for all possible omega combinations
for ii = 1:length(lomega)
 for jj = 1:length(bomega)
 oe(3) = bomega(jj);
 oe(5) = lomega(ii);
 [rPCIf,~] = oe2rv_Visic_Zorana(oe,mu);
 rECI(1,(1+long):(3+long)) = rPCIf.';
 long = long + 3;
 end
end


% For loop which uses the initial values for each combination of omegas to
% determine the rest of the coordinates using kepler solver
long = 0; %reset counter variable
for ii = 1:length(lomega)
 for jj = 1:length(bomega)
 for kk = 2:N
 oe(3) = bomega(jj);
 oe(5) = lomega(ii);
 [E,nu(kk)] = KeplerSolver_Visic_Zorana(a,e,t(kk-1),t(kk),nu(kk-1),mu);
 oe(6) = nu(kk);
 [rPCIf,vPCIf] = oe2rv_Visic_Zorana(oe,mu);
 rECI(kk,(1+long):(3+long)) = rPCIf.';
 end
 long = long + 3;
 end 
end


% Plotting the position results on two seperate figures
[xx,yy,zz] = earthSphere(50);
figure(1)
plot3(rECI(:,1),rECI(:,2),rECI(:,3),rECI(:,4),rECI(:,5),rECI(:,6),rECI(:,7),rECI(:,8),rECI(:,9),rECI(:,10),rECI(:,11),rECI(:,12),rECI(:,13),rECI(:,14),rECI(:,15),xx,yy,zz)
title('\omega = 90 degrees')
legend('\Omega = 0 degrees','\Omega = 45 degrees','\Omega = 90 degrees','\Omega = 135 degrees','\Omega = 180 degrees')
view(az,el)
figure(2)
plot3(rECI(:,16),rECI(:,17),rECI(:,18),rECI(:,19),rECI(:,20),rECI(:,21),rECI(:,22),rECI(:,23),rECI(:,24),rECI(:,25),rECI(:,26),rECI(:,27),rECI(:,28),rECI(:,29),rECI(:,30),xx,yy,zz)
title('\omega = 270 degrees')
legend('\Omega = 0 degrees','\Omega = 45 degrees','\Omega = 90 degrees','\Omega = 135 degrees','\Omega = 180 degrees')
view(az,el)


%-----------------------------------------------------------------------%
% Task 2: Using the previously determined values for position to find   %
% the longitude and latitude of the orbit. Then, plotting the results   %
% on a map of the Earth for each combination of omegas.                 %
%-----------------------------------------------------------------------%


long = 0; % rest counter variable
count = 0; % counter variable


% Converting all of the position values from ECI to ECEF coordinates, then
% using the ECEF coordinated to calculate the longitude (thetae) and
% latitude (phie) of the orbits
for ii = 1:length(lomega)
 for jj = 1:length(bomega)
 rvValuesECEF = eci2ecef_Visic_Zorana(t,rECI(:,(1+long):(3+long)),OmegaE);
 xe = rvValuesECEF(:,1);
 ye = rvValuesECEF(:,2);
 ze = rvValuesECEF(:,3);
 for kk = 1:N
 thetae(kk,jj+count) = atan2(ye(kk),xe(kk));
 phie(kk,jj+count) = atan2(ze(kk),sqrt((xe(kk)^2)+(ye(kk)^2)));
 end
 long = long + 3;
 end
 count = count + length(bomega);
end


% Plotting the results for the orbit on two seperate graphs: one for each
% lomega value
figure(3)
clf
mercatorDisplay(thetae(:,1),phie(:,1))
for ii = 2:(L/2)
 plot(thetae(:,ii)*180/pi,phie(:,ii)*180/pi,marker(ii-1))
 hold on
end
title('\omega = 90 degrees')
legend('\Omega = 0 degrees','\Omega = 45 degrees','\Omega = 90 degrees','\Omega = 135 degrees','\Omega = 180 degrees')
hold off
figure(4)
clf
mercatorDisplay(thetae(:,(L/2)+1),phie(:,(L/2)+1))
for ii = ((L/2)+2):L
 plot(thetae(:,ii)*180/pi,phie(:,ii)*180/pi,marker(ii-6))
 hold on
end
title('\omega = 270 degrees')
legend('\Omega = 0 degrees','\Omega = 45 degrees','\Omega = 90 degrees','\Omega = 135 degrees','\Omega = 180 degrees')
hold off
fprintf('----------------------------------------------------------------------------------------------------------\n')
fprintf('Interesting feature of the geotrack orbits: \n')
fprintf(' \t For each angle of the argument of periapsis, the geotrack orbit shifts by a value of 90 degrees \n \t in longitude for a corresponding 90 degree shift in the longitude of the ascending node. \n')
fprintf('----------------------------------------------------------------------------------------------------------\n')
fprintf('The key difference between the geotrack orbits when the \n argument of the periapsis = 90,270 degrees: \n')
fprintf(' \t The orbit is reflected across the equatorial line when the argument of the periapsis is shifted \n \t by a value of 180 degrees \n')
fprintf('----------------------------------------------------------------------------------------------------------\n')


%-----------------------------------------------------------------------%
% Task 3: Using the longitude and latitude value to find the crossover  % 
% point of the orbit then, calculating the time the satelite spends in  % 
% between the two cross over points. This is done by finding the local  % 
% minimum or maximum (which corresponds to which hemisphere the satelite% 
% spends the most amount of time in) and stepping away from it in both  %
% directions until two points are equal to one another in longitude.    %
% This is the cross over point.                                         %
%-----------------------------------------------------------------------%
thetae = unwrap(thetae);
phie = unwrap(phie);
for ii = 1:L
 if ii <= (L/2)
 om = 1;
 OM = ii;
 else
 om = 2;
 OM = ii - 5;
 end
 arg = lomega(om)*180/(pi);
 long = bomega(OM)*180/(pi);
 [~, ~, nhemi] = find(phie(:,ii) > 0);
 if length(nhemi) >= ceil(N/2)
 local = max(phie(:,ii));
 index = find(phie(:,ii) == local);
 for ii = 1:100
 if (thetae(index+ii)-thetae(index-ii))<=0.01
 dummy = ii;
 break
 end
 end
 time = ((dummy*deltat - 2)/(60*60))*2;
 fprintf('---------------------------------------------------------------------------------------------------\n')
 fprintf('The orbit with the characteristics: \n \t argument of the periapsis = %g degrees \n \t longitude of the ascending node = %g degrees \n ',arg,long)
 fprintf("spends most of its time in the Northern Hemisphere. \n This orbit spends %g hours in between crossover points. \n",time)
 fprintf('---------------------------------------------------------------------------------------------------\n')
 else
 local = min(phie(:,ii));
 index = find(phie(:,ii) == local);
 for ii = 1:100
 if (thetae(index+ii)-thetae(index-ii))<=0.01
 dummy = ii;
 break
 end
 end
 time = ((dummy*deltat - 2)/(60*60))*2;
 fprintf('---------------------------------------------------------------------------------------------------\n')
 fprintf('The orbit with the characteristics: \n \t arguement of the periapsis = %g degrees \n \t longitude of the ascending node = %g degrees \n ',arg,long)
 fprintf("spends most of its time in the Southern Hemisphere. \n This orbit spends %g hours in between crossover points. \n",time)
 fprintf('---------------------------------------------------------------------------------------------------\n')
 end
end
