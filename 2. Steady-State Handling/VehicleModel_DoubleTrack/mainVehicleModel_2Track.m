% ----------------------------------------------------------------
%% Main script for a basic simulation framework with a double track vehcile model
%  authors: 
%  rev. 1.0 Mattia Piccinini & Gastone Pietro Papini Rosati
%  rev. 2.0 Edoardo Pagot
%  date:
%  rev 1.0:    13/10/2020
%  rev 2.0:    16/05/2022
%  rev 2.1:    08/07/2022 (Biral)
%       - added Fz saturation. Correceted error in Fx
%       - initial condition is now parametric in initial speed
%       - changed the braking torque parameters to adapt to a GP2 model
% ----------------------------------------------------------------

% ----------------------------
%% Initialization
% ----------------------------
initialize_environment;


% ----------------------------
%% Load vehicle data
% ----------------------------

% test_tyre_model; % some plot to visualize the curvers resulting from the
% loaded data

vehicle_data = getVehicleDataStruct();
% pacejkaParam = loadPacejkaParam();

% ----------------------------
%% Define initial conditions for the simulation
% ----------------------------
V0 = 0/3.6; % Initial speed
X0 = loadInitialConditions(V0);

% ----------------------------
%% Define the desired speed
% ----------------------------
V_des = 50/3.6; % Initial speed

% ----------------------------
%% Simulation parameters
% ----------------------------
simulationPars = getSimulationParams(); 
Ts = simulationPars.times.step_size;  % integration step for the simulation (fixed step)
T0 = simulationPars.times.t0;         % starting time of the simulation
Tf = simulationPars.times.tf;         % stop time of the simulation

% ----------------------------
%% Speed Ramp Test
% ----------------------------
fprintf('Starting Simulation\n')
test_type = 1; % speed_ramp
tic;
model_sim = sim('Vehicle_Model_2Track');
elapsed_time_simulation = toc;
fprintf('Simulation completed\n')
fprintf('The total simulation time was %.2f seconds\n',elapsed_time_simulation)

% ----------------------------
%% Post-Processing
% ----------------------------
dataAnalysis(model_sim,vehicle_data,Ts,test_type);
% vehicleAnimation(model_sim,vehicle_data,Ts);

%% Camber Angle, Toe Angle and Roll Stiffness Effects
% Camber
camber = -6:2:6; % [deg]

for i=camber
    vehicle_data.front_wheel.static_camber = i; % deg2rad(i);
    if i<0
        fieldname = strcat('camber_min', num2str(abs(i)));
    else
        fieldname = strcat('camber_', num2str(i));
    end

    model_sim_camber.(fieldname) = sim('Vehicle_Model_2Track');
    fprintf('Test with Camber = %2d [deg] completed\n', i)
end
% Camber Post-Processing
camber_effect(model_sim_camber, vehicle_data, Ts, camber);

% Toe
toe = -1:0.5:1;  % [deg]

vehicle_data = getVehicleDataStruct(); 

for i = toe
    vehicle_data.front_wheel.delta_f_0 = deg2rad(i); 
    if i<0 
        fieldname = strcat('toe_min', num2str(abs(i)*10)); 
    else 
        fieldname = strcat('toe_', num2str(i*10)); 
    end  
    model_sim_toe.(fieldname) = sim('Vehicle_Model_2Track');
    fprintf('Test with Toe = %4.1f [deg] completed\n', i)

end 
% Toe Post-Processing
toe_effect(model_sim_toe, vehicle_data, Ts, toe);

% Roll Stiffness
roll_stiff = 0.2:0.1:0.8;

vehicle_data = getVehicleDataStruct(); 
Ks_r = vehicle_data.rear_suspension.Ks_r; 
Ks_f = vehicle_data.front_suspension.Ks_f; 
Ks_tot = Ks_r + Ks_f; 

for i=roll_stiff 
    Ks_f = i*Ks_tot;
    Ks_r = (1-i)*Ks_tot;
    vehicle_data.front_suspension.Ks_f = Ks_f;
    vehicle_data.rear_suspension.Ks_r  = Ks_r;

    fieldname = strcat('roll_stiff_', num2str(i*10)); 
    model_sim_roll_stiff.(fieldname) = sim('Vehicle_Model_2Track');
    fprintf('Test with Epsilon = %3.1f [deg] completed\n', i)

end 
% Roll Stiffness Post-Processing
roll_stiff_effect(model_sim_roll_stiff, vehicle_data, Ts, roll_stiff);

% ----------------------------
%% Steer Ramp Test
% ----------------------------
vehicle_data = getVehicleDataStruct(); 

fprintf('Starting Simulation\n')
test_type = 2; % steer_ramp
tic;
model_sim = sim('Vehicle_Model_2Track');
elapsed_time_simulation = toc;
fprintf('Simulation completed\n')
fprintf('The total simulation time was %.2f seconds\n',elapsed_time_simulation)

% ----------------------------
%% Post-Processing
% ----------------------------
dataAnalysis(model_sim,vehicle_data,Ts,test_type);
% vehicleAnimation(model_sim,vehicle_data,Ts);
