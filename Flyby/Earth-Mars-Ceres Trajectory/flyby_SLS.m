% %{
clc
clear all 
close all

% Sun
mu_s = 1.32712*10^11;           %km^3/s^2   (Gravitational Parameter Sun)

% Earth
mu_e = 398600;                  %km^3/s^2   (Gravitational Parameter Earth)
r_e = 6378;                     %km
g0 = 9.81*10^-3;                %km/s

% Mars
mu_m = 42828;                   %km^3/s^2   (Gravitational Parameter Mars)
r_m = 3389;                     %km

% Ceres
mu_c = 6.26325*10^1;            %km^3/s^2   (Gravitational Parameter Ceres)
r_c = 473;                      %km
r_p = r_c + 50;                 %km         (Radius of paragee, when the sc comes into orbit of Ceres)

I_sp = 465;                     %Isp of the second stage Delta IV

% Grabbing the earth ephemeride file
fid = fopen('earth_launch_2025_2035_ephemerides.txt', 'r');
if fid == -1 % -1 if file could not be found
    disp('Error in Earth file ID');
else 
    temp = textscan(fid, '%f %s %f %f %f %f %f %f', 'Delimiter', ',');
    e_Jt = temp{1};                         % Earth Launch Times
    e_CD = temp{2};                         % Earth Launch Dates (Use for lables on plot)
    e_Pos = [temp{3}, temp{4}, temp{5}];    % Earth Launch Position
    e_Vel = [temp{6}, temp{7}, temp{8}];    % Earth Launch Velocity
    clear temp;
    fclose(fid);
end

% Grabbing the mars ephemeride file
fid = fopen('mars_flyby_2025_2035_ephemerides.txt', 'r');
if fid == -1 % -1 if file could not be found
    disp('Error in mars file ID');
else 
    temp = textscan(fid, '%f %s %f %f %f %f %f %f', 'Delimiter', ',');
    m_Jt = temp{1};                         % Mars Arrival Times
    m_Cd = temp{2};                         % Mars Arrival Dates (Use for lables on plot)
    m_Pos = [temp{3}, temp{4}, temp{5}];    % Mars Arrival Position
    m_Vel = [temp{6}, temp{7}, temp{8}];    % Mars Arrival Velocity
    clear temp;
    fclose(fid);
end

% Grabbing the ceres ephemeride file
fid = fopen('ceres_arrival_2025_2042_ephemerides.txt', 'r');
if fid == -1 % -1 if file could not be found
    disp('Error in Ceres file ID');
else 
    temp = textscan(fid, '%f %s %f %f %f %f %f %f', 'Delimiter', ',');
    c_Jt = temp{1};                         % Ceres Arrival Times
    c_Cd = temp{2};                         % Ceres Arrival Dates (Use for lables on plot)
    c_Pos = [temp{3}, temp{4}, temp{5}];    % Ceres Arrival Position
    c_Vel = [temp{6}, temp{7}, temp{8}];    % Ceres Arrival Velocity
    clear temp;
    fclose(fid);
end
clear fid;

%{
PROCESS FOR CALCULATING C3:
    1. Curve fit C3 for rocket (this case using Falcon Heavy)
    2. Use lambert problem to calculate the departure velocity at Earth 
    (V_1) and arrival velocity at Mars (V_arr).
    3. Calculate V_inf at Earth
    4. Calculate C3
    5. Check if C3 is <= 80 (this will give us a min delevered mass of
    1710kg)
    6. Use the lambert problem at the arrival date at Mars (from Earths) as
    the launch date (to Ceres). [V_dep, V_2]
    7. Calculate V_inf_pos and V_inf_neg
    8. Compare V_inf_pos and V_inf_neg and make sure they are less than 
    1km/s different.
    9. Calculate del_max for the hyperbolic trajectory
    10. Calculate del for the hyperbolic trajectory
    11. Check if del > del_max (makes sure the sc doesn't crash into the
    planet
    12. Calculate delta V that the planet gives
    (2*norm(v_inf_pos)*sin(del/s))
    13. See if trailing or leading by comparing the angle between V_ceres
    adn V_if_neg. If angle < 90 -> leading (negative), angle > 90 trailing
    (positive)
    14. Calculate delta V needed to slow down into orbit around Ceres
    15. Calculate delta V needed (dV_insertion - dV_assist)
    17. Calculate mass delivered
    18. Plot


%}

% Step 1: C3 for Falcon Heavy rocket
% Falcon Heavy Performance
C3_SLS = [15.040; 19.397; 30.569; 50.662; 76.424; 99.743; 130.376] ;
mass = [39335.88; 36856.18; 30977.73; 22111.90; 13701.39; 8493.88; 4200.52] ;

% Fit the data using a 5th-order polynomial
f=fit(C3_SLS,mass,'poly4');

% Preallocate to save speed
m_delivered = zeros(length(e_Jt), length(c_Jt));
dV = zeros(length(e_Jt), length(c_Jt));
C3 = zeros(length(e_Jt), length(c_Jt));
data = [];

% Where we save the xlsx file to
filename = 'SLS Optimal Trajectories.xlsx';
data_titles = {'i', 'j', 'k', 'Launch Date', 'Flyby Date', 'Arrival Date', 'Mass Delivered', 'Delta V', 'C3', 'Delta V gained from flyby'};
data_write = table(data_titles);
writetable(data_write, filename, 'Sheet', 1, 'Range', 'A1');
tic
% Step 2: Lambert from Earth to Mars
for i = 1:length(e_Jt)
   fprintf('Iteration %g/%g\n', i, length(e_Jt));
   for j = 1:length(m_Jt)
        if m_Jt(j) > e_Jt(i)
            [V_1, V_2_arrive] = lambert(e_Pos(i,:), m_Pos(j,:), (m_Jt(j)-e_Jt(i))*(24*60*60), "Pro", mu_s);
            
            % Step 3: Calculate V_inf at Earth
            e_V_inf = V_1-e_Vel(i,:);
            
            % Step 4: Calculate C3 from Earth to Mars
            tempC3 = norm(e_V_inf)^2;
            
            % Step 5: Make sure C3 is reasonable
            if (tempC3 > 0 && tempC3 <= 80)
                
                % Step 6: Lamber from Mars to Ceres
                for k = 1:length(c_Jt)
                    if c_Jt(k) > m_Jt(j)
                        [V_1_depart, V_2] = lambert(m_Pos(j,:), c_Pos(k,:), (c_Jt(k)-m_Jt(j))*(24*60*60), "Pro", mu_s);
                        
                        % Step 7: Calculate V_infinities
                        V_inf_pos = V_1_depart - m_Vel(j, :);
                        V_inf_neg = V_2_arrive - m_Vel(j, :);
                        
                        % Step 8: Make sure V_inf's are close
                        if abs(norm(V_inf_pos)-norm(V_inf_neg)) <= .1 %km/s
                            
                            % Step 9: Calculate del_max
                            del_max = 2*asin(1/(1+((r_m*norm(V_inf_pos)^2)/mu_m)));
                            
                            % Step 10: Calculate current del
                            del = acos(dot(V_inf_pos, V_inf_neg)/(norm(V_inf_pos)*norm(V_inf_neg)));
                            
                            % Step 11: Make sure del < del_max
                            if del < del_max
                               
                                % Step 12: Calculate delta_V from Mars
                                dV_assist = 2*norm(V_inf_pos)*sin(del/2);
                                
                                % Step 13: Check if trailing or leading
                                if norm(V_1_depart) < norm(V_2_arrive)  % Checks if leading
                                    dV_assist = -2*norm(V_inf_pos)*sin(del/2);
                                end
                                
                                % Step 14: Calculate the insertion delta V at Ceres
                                c_V_inf = V_2 -  c_Vel(k, :);
                                V_peri = sqrt(norm(c_V_inf)^2+2*mu_c/r_p);
                                V_circ_c_at_peri = sqrt(mu_c/r_p);
                                dV_insertion = V_peri - V_circ_c_at_peri;
                                
                                % Step 15: Calculate total dV needed
                                dV(i,k) = dV_insertion;
                                C3(i,k) = tempC3;
                                
                                % Step 16: Calculate the mass delivered
                                % Get launch mass
                                m0 = f(tempC3); %limit tempC3 between 0-80
                                m_delivered(i,k) = m0*exp(-dV_insertion/(g0*I_sp));
                                
                                % data = ['launch date'. 'flyby date', 'arrival date', 'm_delivered', 'deltaV Needed', 'C3', 'dV gained from flyby']
                                data(end+1,1) = i;
                                data(end, 2) = j;
                                data(end, 3) = k;
                                data(end, 4) = e_Jt(i);
                                data(end, 5) = m_Jt(j);
                                data(end, 6) = c_Jt(k);
                                data(end, 7) = m_delivered(i,k);
                                data(end, 8) = dV_insertion;
                                data(end, 9) = tempC3;
                                data(end, 10) = dV_assist;
                                
                                table_write = table(data);
                                writetable(table_write, filename, 'Sheet', 1, 'Range', 'A2');
                            end
                        end
                    end
                end                
            end            
        end
   end
   fprintf('Iteration %g/%g\n', i, length(e_Jt));
end
toc
%}

% Graph with Julian Dates and max
% %{
[X,Y] = meshgrid(c_Jt, e_Jt);
contourf(X,Y,m_delivered,'ShowText','on');
xlabel('Arrival Dates from Earth to Ceres (dd-mmm-yyyy)') 
ylabel('Launch Dates from Earth to Ceres (dd-mmm-yyyy)')

% Find and Label the minimum on the graph
max_del_m = max(m_delivered(m_delivered>0));
[max_x_i, max_y_i] = find(m_delivered == max_del_m);

% [data_launch, data_flyby, data_arrive, data_mass, data_deltaV, data_C3] = find(m_delivered == max_del_m);
disp(max_del_m);
disp(X(1, max_x_i));
disp(Y(max_y_i, 1));
format long;
julianDateArrive = c_Jt(max_y_i)
julianDateLaunch = e_Jt(max_x_i)
text(Y(max_y_i, 1), X(1, max_x_i), 0.1, "\leftarrowMax: "+max_del_m, 'Color','red')
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
contourf(X,Y,m_delivered,'ShowText','on');
xlabel('Arrival Dates from Earth to Ceres (dd-mmm-yyyy)') 
ylabel('Launch Dates from Earth to Ceres (dd-mmm-yyyy)') 

datetick('x','dd-mmm-yyyy','keepticks')
datetick('y','dd-mmm-yyyy','keepticks')
%}
