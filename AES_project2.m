clear; clc; close all;

%% Load US06 Data
us06_data = readtable('us06col.txt', 'Delimiter', '\t', 'HeaderLines', 2);
time_US06 = us06_data{:, 1}; % Time (s)
speed_mph_US06 = us06_data{:, 2}; % Speed (mph)
v_US06 = speed_mph_US06 * 0.44704; % Convert mph to m/s

%% Load UDDS Data
udds_data = readtable('uddscol.txt', 'Delimiter', '\t', 'HeaderLines', 2);
time_UDDS = udds_data{:, 1}; % Time (s)
speed_mph_UDDS = udds_data{:, 2}; % Speed (mph)
v_UDDS = speed_mph_UDDS * 0.44704; % Convert mph to m/s

%% Compute Acceleration
a_US06 = gradient(v_US06, mean(diff(time_US06))); 
a_UDDS = gradient(v_UDDS, mean(diff(time_UDDS)));

%% Vehicle and Battery Parameters (From Table 2 and Table 1)
g = 9.81;       % Gravity (m/s^2)
theta = 0;      % Slope (deg)
cr = 0.013;     % Rolling resistance coeff
Af = 2.65;      % Frontal area (m^2)
rho_air = 1.2;  % Air density (kg/m^3)
cd = 0.23;      % Drag coefficient
eta = 0.8;      % Powertrain efficiency
ca = 0.12;      % Auxiliary consumption factor
m_vehicle = 1300; % Vehicle mass excluding battery (kg)

% Tesla (2170) Battery Pack: 96S46P configuration
x_2170 = 96; 
y_2170 = 46;
Vcell_2170 = 3.65;   % V
Qcell_2170 = 4.60;   % Ah
num_cells_2170 = x_2170 * y_2170;
battery_mass_2170 = num_cells_2170 * 68.6e-3;  % cell mass * number of cells
m_total_2170 = m_vehicle + battery_mass_2170;
Vpack_2170 = x_2170 * Vcell_2170;
Qpack_2170 = y_2170 * Qcell_2170; % Ah
Epack_2170 = x_2170 * y_2170 * Vcell_2170 * Qcell_2170; % Wh

% Tesla (4680) Battery Pack: 92S9P configuration
x_4680 = 92; 
y_4680 = 9;
Vcell_4680 = 3.6;  % V
Qcell_4680 = 22;   % Ah
num_cells_4680 = x_4680 * y_4680;
battery_mass_4680 = num_cells_4680 * 358e-3; % kg
m_total_4680 = m_vehicle + battery_mass_4680;
Vpack_4680 = x_4680 * Vcell_4680;
Qpack_4680 = y_4680 * Qcell_4680; % Ah
Epack_4680 = x_4680 * y_4680 * Vcell_4680 * Qcell_4680; % Wh

%% Function to Analyze Cycle
function [P, I, s, E_vehicle, E_con, range] = analyze_cycle(time, v, a, m_total, Vpack, Qpack, Epack, eta, ca, g, theta, cr, rho_air, cd, Af)
    dt = mean(diff(time)); 
    Qpack_As = Qpack * 3600; % Convert Ah to As

    % Compute Power P(t)
    P = (m_total*g*sind(theta) + cr*m_total*g*cosd(theta) + ...
         0.5*rho_air*cd*Af.*v.^2 + m_total.*a).*v; 

    % Compute Current I(t)
    I = P / Vpack;

    % Compute SOC s(t) using Coulomb counting (Corrected sign)
    s = zeros(size(time));
    s(1) = 0.8; % Initial SOC = 80%
    for k = 2:length(time)
        s(k) = s(k-1) - (I(k)*dt)/Qpack_As;  % Use minus sign for discharge
    end

    % Compute Energy Used and Range
    Ed = sum(P(P > 0))*dt/3600; % Discharge energy in Wh
    Ec = sum(P(P < 0))*dt/3600; % Regen energy in Wh
    E_vehicle = (eta*Ec + Ed/eta)*(1 - ca); % Net vehicle energy considering efficiency and auxiliaries

    % Distance traveled
    D = sum(v)*dt;    % meters
    D_km = D/1000;    % km

    % Energy consumption (Wh/km)
    E_con = E_vehicle/D_km;

    % Range (km)
    range = Epack / E_con;
end

%% Analyze Both Cycles for Tesla (2170)
[P_US06_2170, I_US06_2170, s_US06_2170, E_vehicle_US06_2170, E_con_US06_2170, range_US06_2170] = ...
    analyze_cycle(time_US06, v_US06, a_US06, m_total_2170, Vpack_2170, Qpack_2170, Epack_2170, eta, ca, g, theta, cr, rho_air, cd, Af);

[P_UDDS_2170, I_UDDS_2170, s_UDDS_2170, E_vehicle_UDDS_2170, E_con_UDDS_2170, range_UDDS_2170] = ...
    analyze_cycle(time_UDDS, v_UDDS, a_UDDS, m_total_2170, Vpack_2170, Qpack_2170, Epack_2170, eta, ca, g, theta, cr, rho_air, cd, Af);

%% Analyze Both Cycles for Tesla (4680)
[P_US06_4680, I_US06_4680, s_US06_4680, E_vehicle_US06_4680, E_con_US06_4680, range_US06_4680] = ...
    analyze_cycle(time_US06, v_US06, a_US06, m_total_4680, Vpack_4680, Qpack_4680, Epack_4680, eta, ca, g, theta, cr, rho_air, cd, Af);

[P_UDDS_4680, I_UDDS_4680, s_UDDS_4680, E_vehicle_UDDS_4680, E_con_UDDS_4680, range_UDDS_4680] = ...
    analyze_cycle(time_UDDS, v_UDDS, a_UDDS, m_total_4680, Vpack_4680, Qpack_4680, Epack_4680, eta, ca, g, theta, cr, rho_air, cd, Af);

%% Presenting Results According to Questions

% Q1(a): Compute and Plot P(t) for Both Packs
figure('Name','Q1(a): P(t) Plots for Both Packs');
subplot(2,2,1);
plot(time_US06, P_US06_2170/1000);
xlabel('Time (s)'); ylabel('Power (kW)');
title('Q1(a): US06 Power - Tesla (2170)');

subplot(2,2,2);
plot(time_UDDS, P_UDDS_2170/1000);
xlabel('Time (s)'); ylabel('Power (kW)');
title('Q1(a): UDDS Power - Tesla (2170)');

subplot(2,2,3);
plot(time_US06, P_US06_4680/1000);
xlabel('Time (s)'); ylabel('Power (kW)');
title('Q1(a): US06 Power - Tesla (4680)');

subplot(2,2,4);
plot(time_UDDS, P_UDDS_4680/1000);
xlabel('Time (s)'); ylabel('Power (kW)');
title('Q1(a): UDDS Power - Tesla (4680)');

% Q1(b): Compute and Plot I(t) for Both Packs
figure('Name','Q1(b): I(t) Plots for Both Packs');
subplot(2,2,1);
plot(time_US06, I_US06_2170);
xlabel('Time (s)'); ylabel('Current (A)');
title('Q1(b): US06 Current - Tesla (2170)');

subplot(2,2,2);
plot(time_UDDS, I_UDDS_2170);
xlabel('Time (s)'); ylabel('Current (A)');
title('Q1(b): UDDS Current - Tesla (2170)');

subplot(2,2,3);
plot(time_US06, I_US06_4680);
xlabel('Time (s)'); ylabel('Current (A)');
title('Q1(b): US06 Current - Tesla (4680)');

subplot(2,2,4);
plot(time_UDDS, I_UDDS_4680);
xlabel('Time (s)'); ylabel('Current (A)');
title('Q1(b): UDDS Current - Tesla (4680)');

% Q1(c): Compute and Plot s(t) for Both Packs
figure('Name','Q1(c): s(t) Plots for Both Packs');
subplot(2,2,1);
plot(time_US06, s_US06_2170);
xlabel('Time (s)'); ylabel('SOC');
title('Q1(c): US06 SOC - Tesla (2170)');

subplot(2,2,2);
plot(time_UDDS, s_UDDS_2170);
xlabel('Time (s)'); ylabel('SOC');
title('Q1(c): UDDS SOC - Tesla (2170)');

subplot(2,2,3);
plot(time_US06, s_US06_4680);
xlabel('Time (s)'); ylabel('SOC');
title('Q1(c): US06 SOC - Tesla (4680)');

subplot(2,2,4);
plot(time_UDDS, s_UDDS_4680);
xlabel('Time (s)'); ylabel('SOC');
title('Q1(c): UDDS SOC - Tesla (4680)');

% Q2: Range for Both Packs and Cycles
disp('Q2: Range Estimations');
disp('--- Tesla (2170) Results ---');
disp(['US06 Range: ', num2str(range_US06_2170), ' km']);
disp(['UDDS Range: ', num2str(range_UDDS_2170), ' km']);

disp('--- Tesla (4680) Results ---');
disp(['US06 Range: ', num2str(range_US06_4680), ' km']);
disp(['UDDS Range: ', num2str(range_UDDS_4680), ' km']);
