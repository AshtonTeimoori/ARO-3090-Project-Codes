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

I_sp = 470;                     %Isp of the second stage Delta IV

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

% Step 1: C3 for Falcon Heavy rocket
% Falcon Heavy Performance
C3_Falcon_Heavy = [10; 20; 30; 40; 50; 60; 70; 80; 90; 100];
mass = [12345; 10115; 8230; 6640; 5280; 4100; 3080; 2195; 1425; 755];

% Fit the data using a 5th-order polynomial
f=fit(C3_Falcon_Heavy,mass,'poly4');

% Initialize dates (Julian)
depart_date     = 2462516.5;
flyby_date     = 2463256.5;
arrival_date    = 2463796.5;

% Grab the index of each of the signifigant dates
[irow_earth, icol_earth] = find(e_Jt == depart_date);
[irow_mars, icol_mars] = find(m_Jt == flyby_date);
[irow_ceres, icol_ceres] = find(c_Jt == arrival_date);


if m_Jt(irow_mars) > e_Jt(irow_earth)
    [V_1, V_2_arrive] = lambert(e_Pos(irow_earth,:), m_Pos(irow_mars,:), (m_Jt(irow_mars)-e_Jt(irow_earth))*(24*60*60), "Pro", mu_s);

    % Step 3: Calculate V_inf at Earth
    e_V_inf = V_1-e_Vel(irow_earth,:);

    % Step 4: Calculate C3 from Earth to Mars
    tempC3 = norm(e_V_inf)^2;

    % Step 5: Make sure C3 is reasonable
    if (tempC3 > 0 && tempC3 <= 80)

        % Step 6: Lamber from Mars to Ceres
        if c_Jt(irow_ceres) > m_Jt(irow_mars)
            [V_1_depart, V_2] = lambert(m_Pos(irow_mars,:), c_Pos(irow_ceres,:), (c_Jt(irow_ceres)-m_Jt(irow_mars))*(24*60*60), "Pro", mu_s);

            % Step 7: Calculate V_infinities
            V_inf_pos = V_1_depart - m_Vel(irow_mars, :);
            V_inf_neg = V_2_arrive - m_Vel(irow_mars, :);

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
                    c_V_inf = V_2 -  c_Vel(irow_ceres, :);
                    V_peri = sqrt(norm(c_V_inf)^2+2*mu_c/r_p);
                    V_circ_c_at_peri = sqrt(mu_c/r_p);
                    dV_insertion = V_peri - V_circ_c_at_peri;

                    % Step 16: Calculate the mass delivered
                    % Get launch mass
                    m0 = f(tempC3); %limit tempC3 between 0-80
                    m_delivered = m0*exp(-dV_insertion/(g0*I_sp));
                    disp("Mass Delivered: " + m_delivered);
                else
                   disp("del is larger than del_max. Del: " + del + ", del_max: " +  del_max);
                end
            else
                disp("V_infinities are not similar");
            end
        else
            disp("Flyby Date < Arrival Date"); 
        end
    else
        disp("C3 is too high: " + tempC3);
    end
else 
    disp("Launch Date < Flyby Date");
end


