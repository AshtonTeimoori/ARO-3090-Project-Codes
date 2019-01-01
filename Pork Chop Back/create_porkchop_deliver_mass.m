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

%{  
PROCESS FOR CALCULATING MIN DELTA V:
    1. Use lambert to calculte the departure and arrival velociy 
    2. Calculate both V_inf at Ceres and Earth
    3. Check if the V_inf at Earth is < 10km/s
    4. Calculate delta_V based on step 3
        4.1 If V_inf at Earth is less than 10km/s, store delta_V_1
        (Deprecated) 4.2 If V_inf is greater than 10km/s, calculate the delta_V_2
        needed. Store the sum of the delta_V's
    6. Find min delta_V
%}


% Preallocate to speed up
delta_V = zeros(length(c_Cd),length(e_Cd)); % Delivered Mass

data = [];

for i = 1:length(c_Cd)
    for j = 1:length(e_Cd)
        if e_Jt(j) > c_Jt(i)
    
            % Step 1: Use lamber problem to solve for the velocities for the transfer orbit
            [Vc_1, Ve_1] = lambert(c_Pos(i,:), e_Pos(j,:), (e_Jt(j) - c_Jt(i))*(24*60*60), "Pro", mu_s);

            % Step 2: Compute v infinity
            V_inf_c = Vc_1 - c_Vel(i, :);
            
            % Compute v infinity at earth
            V_inf_e = norm(Ve_1 - e_Vel(j, :));

            % Step 3: Check V_inf_e
            if V_inf_e < 10
                % Calculate the delta_V
                delta_V(i,j) = abs(calc_DV_from_planet_and_velocity(c_Vel, mu_c, V_inf_c, r_p));
                
               % data = ['launch date', 'arrival date', 'delta_V', 'v_inf']
               data(end+1,1) = e_Jt(i);
               data(end, 2) = c_Jt(j);
               data(end, 3) = delta_V(i,j);
               data(end, 4) = V_inf_e;
            end                
        end
    end
    fprintf('Iteration %g/%g\n', i, length(c_Cd));
end

%}

% %{
[X,Y] = meshgrid(c_Jt, e_Jt);
contourf(X,Y,delta_V,'ShowText','on');
xlabel('Arrival Dates from Earth to Ceres (dd-mmm-yyyy)') 
ylabel('Launch Dates from Earth to Ceres (dd-mmm-yyyy)')

% Find and Label the minimum on the graph
min_dV_m = min(delta_V(delta_V>0));
[max_x_i, max_y_i] = find(delta_V == min_dV_m);

% [data_launch, data_flyby, data_arrive, data_mass, data_deltaV, data_C3] = find(m_delivered == max_del_m);
disp(min_dV_m);
disp(X(1, max_x_i));
disp(Y(max_y_i, 1));
format long;
julianDateArrive = c_Jt(max_y_i)
julianDateLaunch = e_Jt(max_x_i)
text(Y(max_y_i, 1), X(1, max_x_i), 0.1, "\leftarrowMin: "+min_dV_m, 'Color','red')
%}

% Graph with Calendar Dates
%{
ax = gca;
startLD = datenum('01-01-2025');
endLD = datenum('12-31-2035');
yDates = linspace(startLD, endLD, length(e_Jt));
ax.YTick = yDates;

startAD = datenum('01-01-2025');
endAD = datenum('12-31-2042');
xDates = linspace(startAD, endAD, length(c_Jt));
ax.XTick = xDates;

[X,Y] = meshgrid(xDates, yDates);
contourf(X,Y,delta_V,'ShowText','on');
xlabel('Arrival Dates from Earth to Ceres (dd-mmm-yyyy)') 
ylabel('Launch Dates from Earth to Ceres (dd-mmm-yyyy)') 

datetick('x','dd-mmm-yyyy','keepticks')
datetick('y','dd-mmm-yyyy','keepticks')
%}