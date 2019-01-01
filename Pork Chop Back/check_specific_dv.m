% %{
clc
clear all 
close all

mu_e = 398600;                  %km^3/s^2   (Gravitational Parameter Earth)
r_e = 6378;                     %km

mu_s = 1.32712*10^11;           %km^3/s^2   (Gravitational Parameter Sun)

mu_c = 6.26325*10^1;            %km^3/s^2   (Gravitational Parameter Ceres)
r_c = 473;                      %km
r_p = r_c + 50;                 %km         (Radius of paragee, when the sc comes into orbit of Ceres)


% Grabbing the ceres ephemeride file
fid = fopen('ceres_launch_2025_2050_ephemerides.txt', 'r');
if fid == -1 % -1 if file could not be found
    disp('Error in Earth file ID');
else 
    temp = textscan(fid, '%f %s %f %f %f %f %f %f', 'Delimiter', ',');
    c_Jt = temp{1};                         % Ceres Launch Times
    c_Cd = temp{2};                         % Ceres Launch Dates (Use for lables on plot)
    c_Pos = [temp{3}, temp{4}, temp{5}];     % Ceres Launch Position
    c_Vel = [temp{6}, temp{7}, temp{8}];     % Ceres Launch Velocity
    clear temp;
    fclose(fid);
end

% Grabbing the earth ephemeride file
fid = fopen('earth_arrival_2025_2050_ephemerides.txt', 'r');
if fid == -1 % -1 if file could not be found
    disp('Error in Ceres file ID');
else 
    temp = textscan(fid, '%f %s %f %f %f %f %f %f', 'Delimiter', ',');
    e_Jt = temp{1};                         % Earth Arrival Times
    e_Cd = temp{2};                         % Earth Arrival Dates (Use for lables on plot)
    e_Pos = [temp{3}, temp{4}, temp{5}];     % Earth Arrival Position
    e_Vel = [temp{6}, temp{7}, temp{8}];     % Earth Arrival Velocity
    clear temp;
    fclose(fid);
end
clear fid;

data = [];

% Initialize dates (Julian)
depart_date     = 2463806.5;
arrival_date    = 2464306.5;

% Grab the index of each of the signifigant dates
[irow_ceres, icol_ceres] = find(c_Jt == depart_date);
[irow_earth, icol_earth] = find(e_Jt == arrival_date);

    
% Step 1: Use lamber problem to solve for the velocities for the transfer orbit
[Vc_1, Ve_1] = lambert(c_Pos(irow_ceres,:), e_Pos(irow_earth,:), (e_Jt(irow_earth) - c_Jt(irow_ceres))*(24*60*60), "Pro", mu_s);

% Step 2: Compute v infinity
V_inf_c = Vc_1 - c_Vel(irow_ceres, :);

% Compute v infinity at earth
V_inf_e = norm(Ve_1 - e_Vel(irow_earth, :));

% Step 3: Check V_inf_e
if V_inf_e < 10
    % Calculate the delta_V
    delta_V = abs(calc_DV_from_planet_and_velocity(c_Vel, mu_c, V_inf_c, r_p));

   % data = ['launch date', 'arrival date', 'delta_V', 'v_inf']
   data(1, 1) = delta_V
   data(1, 2) = V_inf_e;
else
    disp("V_inf at Earth is too high: " + V_inf_e);
end
