% %{
clc
clear all 
close all

mu_e = 398600;                  %km^3/s^2   (Gravitational Parameter Earth)
r_e = 6378;                     %km
g0 = 9.81*10^-3;                %km/s

mu_s = 1.32712*10^11;           %km^3/s^2   (Gravitational Parameter Sun)

mu_c = 6.26325*10^1;            %km^3/s^2   (Gravitational Parameter Ceres)
r_c = 473;                      %km
r_p = r_c + 50;                 %km         (Radius of paragee, when the sc comes into orbit of Ceres)

I_sp = 462;                     %Isp of the second stage Delta IV

% Grabbing the earth ephemeride file
fid = fopen('earth_launch_2025_2035_ephemerides.txt', 'r');
if fid == -1 % -1 if file could not be found
    disp('Error in Earth file ID');
else 
    temp = textscan(fid, '%f %s %f %f %f %f %f %f', 'Delimiter', ',');
    e_lt = temp{1};                         % Earth Launch Times
    e_ld = temp{2};                         % Earth Launch Dates (Use for lables on plot)
    e_lp = [temp{3}, temp{4}, temp{5}];     % Earth Launch Position
    e_lv = [temp{6}, temp{7}, temp{8}];     % Earth Launch Velocity
    clear temp;
    fclose(fid);
end

% Grabbing the ceres ephemeride file
fid = fopen('ceres_arrival_2025_2042_ephemerides.txt', 'r');
if fid == -1 % -1 if file could not be found
    disp('Error in Ceres file ID');
else 
    temp = textscan(fid, '%f %s %f %f %f %f %f %f', 'Delimiter', ',');
    c_at = temp{1};                         % Ceres Arrival Times
    c_ad = temp{2};                         % Ceres Arrival Dates (Use for lables on plot)
    c_ap = [temp{3}, temp{4}, temp{5}];     % Ceres Arrival Position
    c_av = [temp{6}, temp{7}, temp{8}];     % Ceres Arrival Velocity
    clear temp;
    fclose(fid);
end
clear fid;

%{  
PROCESS FOR CALCULATING DELIVERED MASS:
    1. Curve fit C3 for rocket (this case using SLS)
    2. Use lambert to calculte the departure and arrival velociy 
    3. Calculate delta V
    4. Calculate C3 from Earth
    5. Calculate launch mass for given C3
    6. Calculate insertion delta V around Ceres
    7. Calculate the delivered mass (del_m) using the rocket equation
%}
 
% Step 1: C3 for rocket
% Delta IV Heavy Performance
C3_Delta_IV = [10; 20; 30; 40; 50; 60; 70; 80; 90; 100];
mass = [8460; 6995; 5755; 4700; 3790; 3000; 2315; 1710; 1180; 705] ;

% Fit the data using a 5th-order polynomial
f=fit(C3_Delta_IV,mass,'poly4');


% % SLS Block-2 performance (From NASA/JPL NEO Deflection App: 
% % https://cneos.jpl.nasa.gov/nda/)
% C3_SLS = [15.040; 19.397; 30.569; 50.662; 76.424; 99.743; 130.376] ;
% mass = [39335.88; 36856.18; 30977.73; 22111.90; 13701.39; 8493.88; 4200.52] ;
% 
% % Fit the data using a 5th-order polynomial
% f=fit(C3_SLS,mass,'poly4');

data = [];

% Preallocate to speed up
del_m = zeros(length(e_ld),length(c_ad)); % Delivered Mass
tic
for i = 1:length(e_ld)
    for j = 1:length(c_ad)
        if c_at(j) > e_lt(i)
            
            % Step 2: Use lamber problem to solve for the velocities for the transfer orbit
            [Ve_1, Vc_1] = lambert(e_lp(i,:), c_ap(j,:), (c_at(j)-e_lt(i))*(24*60*60), "retro", mu_s);
            
            % Step 3: Compute v infinity
            V_inf_e = Ve_1 - e_lv(i, :);
            
            % Step 4: Compute C3 from Earth to Ceres at Earth
            C3 = norm(V_inf_e)^2;
            
            if (C3 > 0 && C3 < 80) % Get rid of unrealistic values
                % Step 5: Get Launch Mass
                m0 = f(C3); %limit c3 between 0-150
                
                % Step 6: Calculate insertion delta V
                V_inf_c = Vc_1-c_av(i, :);
                V_peri = sqrt(norm(V_inf_c)^2+2*mu_c/r_p);
                V_circ_c_at_peri = sqrt(mu_c/r_p);
                dV_c_arrival = V_peri - V_circ_c_at_peri;
                
                % Step 7: Calculate the del_m
                del_m(i, j) = m0*exp(-dV_c_arrival/(g0*I_sp));
                
                % data = ['launch date', 'arrival date', 'm_delivered', 'C3', 'ToF']
                data(end+1,1) = e_lt(i);
                data(end, 2) = c_at(j);
                data(end, 3) = del_m(i, j);
                data(end, 4) = C3;
                data(end, 5) = c_at(j)-e_lt(i);
            else
%                 del_m(i,j) = 0;
            end   
        end
    end
    fprintf('Iteration %g/%g\n', i, length(e_ld));
end
toc
%}
% 
% Graph the data
ax = gca;
startLD = datenum('01-01-2025');
endLD = datenum('12-31-2035');
yDates = linspace(startLD, endLD, length(e_ld));
ax.YTick = yDates;

startAD = datenum('01-01-2025');
endAD = datenum('12-31-2042');
xDates = linspace(startAD, endAD, length(c_ad));
ax.XTick = xDates;

[X,Y] = meshgrid(c_at, e_lt);
% [X,Y] = meshgrid(xDates, yDates);
contourf(X,Y,del_m,'ShowText','on');
xlabel('Arrival Dates from Earth to Ceres (dd-mmm-yyyy)') 
ylabel('Launch Dates from Earth to Ceres (dd-mmm-yyyy)') 
% 
% Find and Label the minimum on the graph
max_del_m = max(del_m(del_m>0));
[max_x_i, max_y_i] = find(del_m == max_del_m);
disp(max_del_m);
disp(X(1, max_x_i));
disp(Y(max_y_i, 1));
format long;
julianDateArrive = c_at(max_y_i)
julianDateLaunch = e_lt(max_x_i)
text(Y(max_y_i, 1), X(1, max_x_i), 0.1, "\leftarrowMax: "+max_del_m, 'Color','red')
%{
datetick('x','dd-mmm-yyyy','keepticks')
datetick('y','dd-mmm-yyyy','keepticks')
%}