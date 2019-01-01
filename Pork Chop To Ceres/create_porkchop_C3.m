%%{
clc
clear all 
close all

mu_s = 1.32712*10^11;           %km^3/s^2   (Gravitational Parameter Sun)

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
PROCESS FOR CALCULATING PORKCHOP PLOT:
    1. Use lambert to calculte the departure and arrival velociy 
    2. Calculate delta V
    3. Calculate C3 and store in matrix
%}

% Preallocate to speed up
C3 = zeros(length(e_ld),length(c_ad));

for i = 1:length(e_ld)
    for j = 1:length(c_ad)
        if c_at(j) > e_lt(i)
            
            % Step 1: Use lamber problem to solve for the velocities for the transfer orbit
            [Ve_1, Vc_1] = lambert(e_lp(i,:), c_ap(j,:), (c_at(j)-e_lt(i))*(24*60*60), "Pro", mu_s);
            
            % Step 2: Compute v infinity
            V_inf_e = Ve_1 - e_lv(i, :);
            
            % Step 3: Compute C3 from Earth to Ceres at Earth
            tempC3 = norm(V_inf_e)^2;
            if (tempC3 < 100)
                C3(i,j) = tempC3;
            else
                C3(i,j) = 0;
            end   
        end
    end
    fprintf('Iteration %g/%g\n', i, length(e_ld));
end

%}

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

[X,Y] = meshgrid(xDates, yDates);
contourf(X,Y,C3,'ShowText','on');
xlabel('Arrival Dates from Earth to Ceres (dd-mmm-yyyy)') 
ylabel('Launch Dates from Earth to Ceres (dd-mmm-yyyy)') 

% Find and Label the minimum on the graph
min_C3 = min(C3(C3>0));
[min_x_i, min_y_i] = find(C3 == min_C3);
disp(X(1, min_x_i));
disp(Y(min_y_i, 1));
text(Y(min_y_i, 1), X(1, min_x_i), 0.1, "\leftarrowMin: "+min_C3, 'Color','red')

% datetick('x','dd-mmm-yyyy','keepticks')
% datetick('y','dd-mmm-yyyy','keepticks')