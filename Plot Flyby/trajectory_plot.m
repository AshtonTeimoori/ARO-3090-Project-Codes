clc
clear all
close all

% Sun
mu	 = 1.32712*10^11; %km^3/s^2     (Gravitational Parameter Sun)

% Grabbing the Earth's ephemeride file
fid = fopen('earth_2025_2050_ephemerides.txt', 'r');
if fid == -1 % -1 if file isn't found
    disp("Error: Cannot locate Earth's ephemeride file");
else
    temp = textscan(fid, '%f %s %f %f %f %f %f %f', 'Delimiter', ',');
    e_Jt = temp{1};                         %   Earth Launch Times
    e_CD = temp{2};                         %   Earth Launch Dates
    e_Pos = [temp{3}, temp{4}, temp{5}];    %   Earth Launch Position
    e_Vel = [temp{6}, temp{7}, temp{8}];    %   Earth Launch Velocity
    clear temp;
    fclose(fid);
end

% Grabbing the Mars's ephemeride file
fid = fopen('mars_2025_2050_ephemerides.txt', 'r');
if fid == -1 % -1 if file isn't found
    disp("Error: Cannot locate Mars's ephemeride file");
else
    temp = textscan(fid, '%f %s %f %f %f %f %f %f', 'Delimiter', ',');
    m_Jt = temp{1};                         %   Mars Launch Times
    m_CD = temp{2};                         %   Mars Launch Dates
    m_Pos = [temp{3}, temp{4}, temp{5}];    %   Mars Launch Position
    m_Vel = [temp{6}, temp{7}, temp{8}];    %   Mars Launch Velocity
    clear temp;
    fclose(fid);
end

% Grabbing the Ceres's ephemeride file
fid = fopen('ceres_2025_2050_ephemerides.txt', 'r');
if fid == -1 % -1 if file isn't found
    disp("Error: Cannot locate Ceres's ephemeride file");
else
    temp = textscan(fid, '%f %s %f %f %f %f %f %f', 'Delimiter', ',');
    c_Jt = temp{1};                         %   Ceres Launch Times
    c_CD = temp{2};                         %   Ceres Launch Dates
    c_Pos = [temp{3}, temp{4}, temp{5}];    %   Ceres Launch Position
    c_Vel = [temp{6}, temp{7}, temp{8}];    %   Ceres Launch Velocity
    clear temp;
    fclose(fid);
end

% Initialize dates (Julian)
depart_date     = 2462516.5;
flyby_date      = 2463256.5;
arrival_date    = 2463796.5;

% Grab the index of each of the signifigant dates
[irow_earth, icol_earth] = find(e_Jt == depart_date);
[irow_mars, icol_mars] = find(m_Jt == flyby_date);
[irow_ceres, icol_ceres] = find(c_Jt == arrival_date);

% Integration options
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8) ;


% Transfer from Earth to Mars
% Solve lambert problem
[VL_earth, VL_mars_arrival] = lambert(e_Pos(irow_earth, :), m_Pos(irow_mars, :), (m_Jt(irow_mars, 1)-e_Jt(irow_earth, 1))*(86400), "Pro", mu);


% Initialize initial conditions for SC from Earth to Mars
Y0_earth = [e_Pos(irow_earth, :)'; VL_earth'];

% Set integration time from Earth to Mars
ToF_earth_mars = m_Jt(irow_mars, 1)-e_Jt(irow_earth, 1);
t0_earth_mars = 0; 
tf_earth_mars = ToF_earth_mars*86400; % 1 day (or 86400 s)

% call ode45 for Earth to Mars
[T_em,Y_em] = ode45(@twobodyEOM3D, [t0_earth_mars tf_earth_mars], Y0_earth, options, mu);


% Transfer from Mars to Ceres
% Solve lambert problem
[VL_mars_depart, VL_ceres] = lambert(m_Pos(irow_mars, :), c_Pos(irow_ceres, :), (c_Jt(irow_ceres, 1)-m_Jt(irow_mars, 1))*(86400), "Pro", mu);

% Initialize initial conditions for SC from Mars to Ceres
Y0_mars = [m_Pos(irow_mars, :)'; VL_mars_depart'];

% Integration time
ToF_mars_ceres = c_Jt(irow_ceres, 1)-m_Jt(irow_mars, 1);
t0_mars_ceres = 0; 
tf_mars_ceres = ToF_mars_ceres*86400; % 1 day (or 86400 s)

% call ode45 for Mars to Ceres
[T_mc,Y_mc] = ode45(@twobodyEOM3D, [t0_mars_ceres tf_mars_ceres], Y0_mars, options, mu);

color = '.g';

% Plot radius over time
figure("Name", "Trajectory")
plot3(e_Pos(:,1), e_Pos (:,2), e_Pos(:,3), '-b');
hold on
plot3(m_Pos(:,1), m_Pos(:,2), m_Pos(:,3), '-r');
hold on 
plot3(c_Pos(:,1), c_Pos(:,2), c_Pos(:,3), '-w');
hold on 
plot3(0,0,0,'*y')
set(gca,'Color','k')
hold on;

% Set level orbit plane
t=0:pi/2:2*pi;
x=cos(t);
y=sin(t);
z=zeros(1,5);

plot3(x*1e9, z, y*1e9, '-k');
hold on;
plot3(z, x*1e9, y*1e9, '-k');
hold on;
plot3(x*1e9, y*1e9, z, '-k');
hold on;

% Plot where the sc will intersect the planets
plot3(e_Pos(irow_earth, 1), e_Pos(irow_earth, 2), e_Pos(irow_earth, 3), '.b');
hold on;

plot3(m_Pos(irow_mars, 1), m_Pos(irow_mars, 2), m_Pos(irow_mars, 3), '.r');
hold on;

plot3(c_Pos(irow_ceres, 1), c_Pos(irow_ceres, 2), c_Pos(irow_ceres, 3), '.w');
hold on;

% Plot from Earth to Mars
for i=1:length(T_em)
    inc_val = floor(T_em(i)/86400);
    if i > 150
        plot3(e_Pos(irow_earth-1+inc_val, 1), e_Pos(irow_earth-1+inc_val, 2), e_Pos(irow_earth-1+inc_val, 3), '.b');
        plot3(m_Pos(irow_earth-1+inc_val, 1), m_Pos(irow_earth-1+inc_val, 2), m_Pos(irow_earth-1+inc_val, 3), '.r');
        plot3(c_Pos(irow_earth-1+inc_val, 1), c_Pos(irow_earth-1+inc_val, 2), c_Pos(irow_earth-1+inc_val, 3), '.w');
    end
    plot3(Y_em(i, 1), Y_em(i, 2), Y_em(i, 3), color);
    hold on;
    pause(0.01);
end

plot3(m_Pos(irow_earth-1+inc_val, 1), m_Pos(irow_earth-1+inc_val, 2), m_Pos(irow_earth-1+inc_val, 3), 'Or');
hold on;

% Plot from Mars to Ceres
for i=1:length(T_mc)
    inc_val = floor(T_mc(i)/86400);
    plot3(e_Pos(irow_mars-1+inc_val, 1), e_Pos(irow_mars-1+inc_val, 2), e_Pos(irow_mars-1+inc_val, 3), '.b');
    plot3(m_Pos(irow_mars-1+inc_val, 1), m_Pos(irow_mars-1+inc_val, 2), m_Pos(irow_mars-1+inc_val, 3), '.r');
    plot3(c_Pos(irow_mars-1+inc_val, 1), c_Pos(irow_mars-1+inc_val, 2), c_Pos(irow_mars-1+inc_val, 3), '.w');
    plot3(Y_mc(i, 1), Y_mc(i, 2), Y_mc(i, 3), color);
    hold on;
    pause(0.01);
end

plot3(c_Pos(irow_mars-1+inc_val, 1), c_Pos(irow_mars-1+inc_val, 2), c_Pos(irow_mars-1+inc_val, 3), 'Ow');
hold on;
