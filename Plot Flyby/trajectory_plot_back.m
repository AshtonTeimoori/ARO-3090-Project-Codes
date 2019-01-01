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
depart_date     = 2466306.5;
arrival_date    = 2466686.5;

% Grab the index of each of the signifigant dates
[irow_ceres, icol_ceres] = find(c_Jt == depart_date);
[irow_earth, icol_earth] = find(e_Jt == arrival_date);

% Integration options
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8) ;


% Transfer from Ceres to Earth
% Solve lambert problem
disp((e_Jt(irow_earth, 1)-c_Jt(irow_ceres, 1)));
[VL_ceres, VL_earth] = lambert(c_Pos(irow_ceres, :), e_Pos(irow_earth, :), (e_Jt(irow_earth, 1)-c_Jt(irow_ceres, 1))*(86400), "Pro", mu);


% Initialize initial conditions for SC from Ceres to Earth
Y0_ceres = [c_Pos(irow_ceres, :)'; VL_ceres'];

% Set integration time from Ceres to Earth
ToF_ce = e_Jt(irow_earth, 1)-c_Jt(irow_ceres, 1);
t0_ce = 0; 
tf_ce = ToF_ce*86400; % 1 day (or 86400 s)

% call ode45 for Earth to Mars
[T_ce,Y_ce] = ode45(@twobodyEOM3D, [t0_ce tf_ce], Y0_ceres, options, mu);


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

% Less than 0.6 AU
oneAU = 149597870.700;
auBorder = 0:0.05:2*pi;
xAU = oneAU*0.6*cos(auBorder);
yAU = oneAU*0.6*sin(auBorder);
plot3(xAU, yAU, zeros(size(auBorder)), '-y');
hold on 

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

% Plot from Ceres to Earth
for i=1:length(T_ce)
    inc_val = floor(T_ce(i)/86400);
    plot3(e_Pos(irow_ceres-1+inc_val, 1), e_Pos(irow_ceres-1+inc_val, 2), e_Pos(irow_ceres-1+inc_val, 3), '.b');
    plot3(m_Pos(irow_ceres-1+inc_val, 1), m_Pos(irow_ceres-1+inc_val, 2), m_Pos(irow_ceres-1+inc_val, 3), '.r');
    plot3(c_Pos(irow_ceres-1+inc_val, 1), c_Pos(irow_ceres-1+inc_val, 2), c_Pos(irow_ceres-1+inc_val, 3), '.w');
    plot3(Y_ce(i, 1), Y_ce(i, 2), Y_ce(i, 3), color);
    hold on;
    pause(0.01);
end

plot3(e_Pos(irow_ceres-1+inc_val, 1), e_Pos(irow_ceres-1+inc_val, 2), e_Pos(irow_ceres-1+inc_val, 3), 'Ob');
hold on;
