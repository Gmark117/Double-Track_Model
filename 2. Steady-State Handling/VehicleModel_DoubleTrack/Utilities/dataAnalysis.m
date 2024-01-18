function dataAnalysis(model_sim,vehicle_data,Ts, test_type)

    % ----------------------------------------------------------------
    %% Post-Processing and Data Analysis
    % ----------------------------------------------------------------
    dataset_start = 50000;


    % ---------------------------------
    %% Load vehicle data
    % ---------------------------------
    Lf = vehicle_data.vehicle.Lf;  % [m] Distance between vehicle CoG and front wheels axle
    Lr = vehicle_data.vehicle.Lr;  % [m] Distance between vehicle CoG and front wheels axle
    L  = vehicle_data.vehicle.L;   % [m] Vehicle length
    Wf = vehicle_data.vehicle.Wf;  % [m] Width of front wheels axle 
    Wr = vehicle_data.vehicle.Wr;  % [m] Width of rear wheels axle                   
    m  = vehicle_data.vehicle.m;   % [kg] Vehicle Mass
    g  = vehicle_data.vehicle.g;   % [m/s^2] Gravitational acceleration
    tau_D = vehicle_data.steering_system.tau_D;  % [-] steering system ratio (pinion-rack)
    h_rc_f = vehicle_data.front_suspension.h_rc_f; % roll center height front
    h_rc_r = vehicle_data.rear_suspension.h_rc_r; % roll center height rear
    hGs = vehicle_data.vehicle.hGs;
    h_r = h_rc_r + (h_rc_f - h_rc_r)*Lr/L;
    h_s = hGs - h_r;
    Ks_r=vehicle_data.rear_suspension.Ks_r;
    Karb_r = vehicle_data.rear_suspension.Karb_r;
    Ks_f=vehicle_data.front_suspension.Ks_f;
    Karb_f = vehicle_data.front_suspension.Karb_f;
    % ---------------------------------
    %% Extract data from simulink model
    % ---------------------------------
    time_sim = model_sim.states.u.time;
    % dt = time_sim(2)-time_sim(1);

    % -----------------
    % Inputs
    % -----------------
    ped_0      = model_sim.inputs.ped_0.data;
    delta_D    = model_sim.inputs.delta_D.data;

    % -----------------
    % States
    % -----------------
    x_CoM      = model_sim.states.x.data;
    y_CoM      = model_sim.states.y.data;
    psi        = model_sim.states.psi.data;
    u          = model_sim.states.u.data;
    v          = model_sim.states.v.data;
    Omega      = model_sim.states.Omega.data;
    Fz_rr      = model_sim.states.Fz_rr.data;
    Fz_rl      = model_sim.states.Fz_rl.data;
    Fz_fr      = model_sim.states.Fz_fr.data;
    Fz_fl      = model_sim.states.Fz_fl.data;
    % delta      = model_sim.states.delta.data;
    omega_rr   = model_sim.states.omega_rr.data;
    omega_rl   = model_sim.states.omega_rl.data;
    omega_fr   = model_sim.states.omega_fr.data;
    omega_fl   = model_sim.states.omega_fl.data;
    alpha_rr   = model_sim.states.alpha_rr.data;
    alpha_rl   = model_sim.states.alpha_rl.data;
    alpha_fr   = model_sim.states.alpha_fr.data;
    alpha_fl   = model_sim.states.alpha_fl.data;
    kappa_rr   = model_sim.states.kappa_rr.data;
    kappa_rl   = model_sim.states.kappa_rl.data;
    kappa_fr   = model_sim.states.kappa_fr.data;
    kappa_fl   = model_sim.states.kappa_fl.data;

    % -----------------
    % Extra Parameters
    % -----------------
    Tw_rr      = model_sim.extra_params.Tw_rr.data;
    Tw_rl      = model_sim.extra_params.Tw_rl.data;
    Tw_fr      = model_sim.extra_params.Tw_fr.data;
    Tw_fl      = model_sim.extra_params.Tw_fl.data;
    Fx_rr      = model_sim.extra_params.Fx_rr.data;
    Fx_rl      = model_sim.extra_params.Fx_rl.data;
    Fx_fr      = model_sim.extra_params.Fx_fr.data;
    Fx_fl      = model_sim.extra_params.Fx_fl.data;
    Fy_rr      = model_sim.extra_params.Fy_rr.data;
    Fy_rl      = model_sim.extra_params.Fy_rl.data;
    Fy_fr      = model_sim.extra_params.Fy_fr.data;
    Fy_fl      = model_sim.extra_params.Fy_fl.data;
    Mz_rr      = model_sim.extra_params.Mz_rr.data;
    Mz_rl      = model_sim.extra_params.Mz_rl.data;
    Mz_fr      = model_sim.extra_params.Mz_fr.data;
    Mz_fl      = model_sim.extra_params.Mz_fl.data;
    gamma_rr   = model_sim.extra_params.gamma_rr.data;
    gamma_rl   = model_sim.extra_params.gamma_rl.data;
    gamma_fr   = model_sim.extra_params.gamma_fr.data;
    gamma_fl   = model_sim.extra_params.gamma_fl.data;
    delta_fr   = model_sim.extra_params.delta_fr.data;
    delta_fl   = model_sim.extra_params.delta_fl.data;

    % Chassis side slip angle beta [rad]
    beta = atan(v./u);

    % -----------------
    % Accelerations
    % -----------------
    % Derivatives of u, v [m/s^2]
    dot_u = diff(u)/Ts;
    dot_v = diff(v)/Ts;
    % Total longitudinal and lateral accelerations
    Ax = dot_u(1:end) - Omega(2:end).*v(2:end);
    Ay = dot_v(1:end) + Omega(2:end).*u(2:end);
    % Ax low-pass filtered signal (zero-phase digital low-pass filtering)
    Wn_filter = 0.01;
    [b_butt,a_butt] = butter(4,Wn_filter,'low');
    Ax_filt = filtfilt(b_butt,a_butt,Ax);  
    dot_u_filt = filtfilt(b_butt,a_butt,dot_u);  
    % Steady state lateral acceleration
    % Ay_ss = Omega(round(length(Omega)*0.9)) * u(round(length(u)*0.9));
    % Longitudinal jerk [m/s^3]
    % jerk_x = diff(dot_u)/Ts;

    % -----------------
    % Other parameters
    % -----------------
    % Total CoM speed [m/s]
    vG = sqrt(u.^2 + v.^2);
    % Steady state and transient curvature [m]
    rho_ss   = Omega./vG;
    rho_tran = ((dot_v.*u(1:end-1) - dot_u.*v(1:end-1)) ./ ((vG(1:end-1)).^3)) + rho_ss(1:end-1);
    % Desired sinusoidal steering angle for the equivalent single track front wheel
    desired_steer_atWheel = delta_D/tau_D;
    Rho_ss= (desired_steer_atWheel+ (alpha_rr-alpha_fr))/L;

   % ---------------------------------
    %% PARAMETERS
    % ---------------------------------
    %% Lateral load transfer
    % Wr, Wf, m, h_GS=y_CoM , h_rc_f, h_rc_r 
    Fyr= Fy_rl + Fy_rr;
    Fyf= Fy_fl + Fy_fr;
    
    Delta_Fzr_inst= m.*Ay.*(h_rc_r*Lf/(L*Wr));
    Delta_Fzf_inst= m.*Ay.*(h_rc_f*Lr/(L*Wf));

    Ks_r_0 = Ks_r + Karb_r;
    Ks_f_0 = Ks_f + Karb_f;
    Ks_r_High = Ks_r + Karb_r + 10*Karb_f;
    Ks_f_High = Ks_f + 10*Karb_f;

    eps_phi = Ks_f_0/(Ks_f_0+Ks_r_0);

    Delta_Fzr_trans0=m.*Ay.*(h_s*(1-eps_phi)/Wr);
    Delta_Fzf_trans0=m.*Ay.*(h_s*eps_phi/Wf);
    Delta_Fzr_trans1=m.*Ay.*(hGs*Ks_r_0/(Wr*(Ks_f_High+Ks_r_0)));
    Delta_Fzf_trans1=m.*Ay.*(hGs*Ks_f_High/(Wf*(Ks_f_High+Ks_r_0)));
    Delta_Fzr_trans2=m.*Ay.*(hGs*Ks_r_High/(Wr*(Ks_f_0+Ks_r_High)));
    Delta_Fzf_trans2=m.*Ay.*(hGs*Ks_f_0/(Wf*(Ks_f_0+Ks_r_High)));

    Delta_Fzr_tot0=Delta_Fzr_inst+Delta_Fzr_trans0;
    Delta_Fzf_tot0=Delta_Fzf_inst+Delta_Fzf_trans0;
    Delta_Fzr_tot1=Delta_Fzr_inst+Delta_Fzr_trans1;
    Delta_Fzf_tot1=Delta_Fzf_inst+Delta_Fzf_trans1;
    Delta_Fzr_tot2=Delta_Fzr_inst+Delta_Fzr_trans2;
    Delta_Fzf_tot2=Delta_Fzf_inst+Delta_Fzf_trans2;

    %% Normalized axle characteristic
    % Yr=m*Ay*Lr/L 
    % Fyf=m*Ay*Lf/L
    % Fzr0= m*g*Lf/L
    % Fzf0= m*g*Lr/L
    % Y_r=m.*Ay.*Lf/L;
    % Y_f=m.*Ay.*Lr/L;
    Fzr_ss=m*g*Lf/L;
    Fzf_ss=m*g*Lr/L;

    Fzr0 = Fzr_ss + Delta_Fzr_tot0;
    Fzf0 = Fzf_ss + Delta_Fzf_tot0;
    Fzr1 = Fzr_ss + Delta_Fzr_tot1;
    Fzf1 = Fzf_ss + Delta_Fzf_tot1;
    Fzr2 = Fzr_ss + Delta_Fzr_tot2;
    Fzf2 = Fzf_ss + Delta_Fzf_tot2;

    % mi_r= Y_r/Fzr0;
    % mi_f= Y_f/Fzf0;
    mi_r0=Fyr(1:length(Ay),1)./Fzr0;
    mi_f0=Fyf(1:length(Ay),1)./Fzf0;
    mi_r1=Fyr(1:length(Ay),1)./Fzr1;
    mi_f1=Fyf(1:length(Ay),1)./Fzf1;
    mi_r2=Fyr(1:length(Ay),1)./Fzr2;
    mi_f2=Fyf(1:length(Ay),1)./Fzf2;

    %% Handling diagram
    % oc_rhs_t = desired_steer_atWheel(1:end-1)/L - rho_tran*L;
    oc_rhs_ss = desired_steer_atWheel(1:end-1) - Rho_ss(1:end-1)*L;

    hand_start_X = Ay(1)/g;
    hand_start_Y = oc_rhs_ss(1);
    hand_end_X   = Ay(round(length(Ay)*0.2))/g;
    hand_end_Y   = oc_rhs_ss(round(length(Ay)*0.2));

    coeffs = polyfit([hand_start_X,hand_end_X],[hand_start_Y,hand_end_Y],1);

    hand_lin = polyval(coeffs, Ay./g);

    hand_point_2_X = Ay(round(end*0.25))/g;
    hand_point_2_Y = oc_rhs_ss(round(length(Ay)*0.25));
    hand_point_3_X = Ay(round(end*0.50))/g;
    hand_point_3_Y = oc_rhs_ss(round(length(Ay)*0.50));
    hand_point_4_X = Ay(round(end*0.75))/g;
    hand_point_4_Y = oc_rhs_ss(round(length(Ay)*0.75));
    hand_end_X   = Ay(end)/g;
    hand_end_Y   = oc_rhs_ss(end);

    coeffs = polyfit([hand_start_X,hand_point_2_X,hand_point_3_X,hand_point_4_X,hand_end_X], ...
                     [hand_start_Y,hand_point_2_Y,hand_point_3_Y,hand_point_4_Y,hand_end_Y], 3);

    hand_fit = polyval(coeffs, Ay./g);

    %% K_us
    Kus_theor = (-m/L^2)*(Lf/Ks_r - Lr/Ks_f);
    Kus_pract = deg2rad(coeffs(3));

    fprintf('Theoretical Kus = %f \n' , Kus_theor);
    fprintf('Practical Kus   = %f \n' , Kus_pract);

    % ---------------------------------
    %% PLOTS
    % ---------------------------------

    % ---------------------------------
    %% Plot vehicle inputs
    % ---------------------------------
    figure('Name','Inputs','NumberTitle','off'), clf   
    % --- pedal --- %
    ax(1) = subplot(211);
    hold on
    plot(time_sim,ped_0,'LineWidth',2)
    grid on
    title('pedal $p_0$ [-]')
    xlim([0 time_sim(end)])
    % --- delta_0 --- %
    ax(2) = subplot(212);
    plot(time_sim,delta_D,'LineWidth',2)
    grid on
    title('steering angle $\delta_D$ [deg]')
    xlim([0 time_sim(end)])
    
    if test_type==1
        exportgraphics(gcf, 'Plots\Speed_Ramp_Inputs.png')
    else
        exportgraphics(gcf, 'Plots\Steer_Ramp_Inputs.png')
    end
    
    % ---------------------------------
    %% Plot vehicle motion
    % ---------------------------------
    figure('Name','veh motion','NumberTitle','off'), clf   
    % --- u --- %
    ax(1) = subplot(221);
    plot(time_sim,u*3.6,'LineWidth',2)
    grid on
    title('$u$ [km/h]')
    xlim([0 time_sim(end)])
    % --- v --- %
    ax(2) = subplot(222);
    plot(time_sim,v,'LineWidth',2)
    grid on
    title('$v$ [m/s]')
    xlim([0 time_sim(end)])
    % --- Omega --- %
    ax(3) = subplot(223);
    plot(time_sim,Omega,'LineWidth',2)
    grid on
    title('$\Omega$ [rad/s]')
    xlim([0 time_sim(end)])

    if test_type==1
        exportgraphics(gcf, 'Plots\Vehicle_Motion.png')
    end

    % ---------------------------------
    %% Plot steering angles
    % ---------------------------------
    figure('Name','steer','NumberTitle','off'), clf   
    % --- delta_0 --- %
    ax(1) = subplot(221);
    plot(time_sim,delta_D,'LineWidth',2)
    grid on
    title('$\delta_0$ [deg]')
    xlim([0 time_sim(end)])
    % --- delta_fr --- %
    ax(2) = subplot(222);
    plot(time_sim,delta_fr,'LineWidth',2)
    grid on
    title('$\delta_{fr}$ [deg]')
    xlim([0 time_sim(end)])
    % --- delta_fl --- %
    ax(3) = subplot(223);
    hold on
    plot(time_sim,delta_fl,'LineWidth',2)
    grid on
    title('$\delta_{fl}$ [deg]')
    xlim([0 time_sim(end)])
    % --- comparison --- %
    ax(4) = subplot(224);
    hold on
    plot(time_sim,delta_D/tau_D,'LineWidth',2)
    plot(time_sim,delta_fr,'LineWidth',2)
    plot(time_sim,delta_fl,'LineWidth',2)
    grid on
    legend('$\delta_D/\tau_D$','$\delta_{fr}$','$\delta_{fl}$','location','best')
    xlim([0 time_sim(end)])

    if test_type==1
        exportgraphics(gcf, 'Plots\Steering_Angles.png')
    end

    % -------------------------------
    %% Plot lateral tire slips and lateral forces
    % -------------------------------
    figure('Name','Lateral slips & forces','NumberTitle','off'), clf
    % --- alpha_rr --- %
    ax(1) = subplot(331);
    plot(time_sim,alpha_rr,'LineWidth',2)
    grid on
    title('$\alpha_{rr}$ [deg]')
    xlim([0 time_sim(end)])
    % --- alpha_rl --- %
    ax(2) = subplot(332);
    plot(time_sim,alpha_rl,'LineWidth',2)
    grid on
    title('$\alpha_{rl}$ [deg]')
    xlim([0 time_sim(end)])
    % --- alpha_fr --- %
    ax(3) = subplot(333);
    plot(time_sim,alpha_fr,'LineWidth',2)
    grid on
    title('$\alpha_{fr}$ [deg]')
    xlim([0 time_sim(end)])
    % --- alpha_fl --- %
    ax(4) = subplot(334);
    plot(time_sim,alpha_fl,'LineWidth',2)
    grid on
    title('$\alpha_{fl}$ [deg]')
    xlim([0 time_sim(end)])
    % --- Fy_rr --- %
    ax(5) = subplot(335);
    plot(time_sim,Fy_rr,'LineWidth',2)
    grid on
    title('$Fy_{rr}$ [N]')
    xlim([0 time_sim(end)])
    % --- Fy_rl --- %
    ax(6) = subplot(336);
    plot(time_sim,Fy_rl,'LineWidth',2)
    grid on
    title('$Fy_{rl}$ [Nm]')
    xlim([0 time_sim(end)])
    % --- Fy_fr --- %
    ax(7) = subplot(337);
    plot(time_sim,Fy_fr,'LineWidth',2)
    grid on
    title('$Fy_{fr}$ [N]')
    xlim([0 time_sim(end)])
    % --- Fy_fl --- %
    ax(8) = subplot(338);
    plot(time_sim,Fy_fl,'LineWidth',2)
    grid on
    title('$Fy_{fl}$ [N]')
    xlim([0 time_sim(end)])

    % linkaxes(ax,'x')
    clear ax

    if test_type==1
        exportgraphics(gcf, 'Plots\Lateral_Slips_and_Forces.png')
    end

    % ---------------------------------
    %% Plot longitudinal tire slips and longitudinal forces
    % ---------------------------------
    figure('Name','Long slips & forces','NumberTitle','off'), clf
    % --- kappa_rr --- %
    ax(1) = subplot(331);
    plot(time_sim,kappa_rr,'LineWidth',2)
    grid on
    title('$\kappa_{rr}$ [-]')
    xlim([0 time_sim(end)])
    % --- kappa_rl --- %
    ax(2) = subplot(332);
    plot(time_sim,kappa_rl,'LineWidth',2)
    grid on
    title('$\kappa_{rl}$ [-]')
    xlim([0 time_sim(end)])
    % --- kappa_fr --- %
    ax(3) = subplot(333);
    plot(time_sim,kappa_fr,'LineWidth',2)
    grid on
    title('$\kappa_{fr}$ [-]')
    xlim([0 time_sim(end)])
    % --- kappa_fl --- %
    ax(4) = subplot(334);
    plot(time_sim,kappa_fl,'LineWidth',2)
    grid on
    title('$\kappa_{fl}$ [-]')
    xlim([0 time_sim(end)])
    % --- Fx_rr --- %
    ax(5) = subplot(335);
    plot(time_sim,Fx_rr,'LineWidth',2)
    grid on
    title('$Fx_{rr}$ [N]')
    xlim([0 time_sim(end)])
    % --- Fx_rl --- %
    ax(6) = subplot(336);
    plot(time_sim,Fx_rl,'LineWidth',2)
    grid on
    title('$Fx_{rl}$ [N]')
    xlim([0 time_sim(end)])
    % --- Fx_fr --- %
    ax(7) = subplot(337);
    plot(time_sim,Fx_fr,'LineWidth',2)
    grid on
    title('$Fx_{fr}$ [N]')
    xlim([0 time_sim(end)])
    % --- Fx_fl --- %
    ax(8) = subplot(338);
    plot(time_sim,Fx_fl,'LineWidth',2)
    grid on
    title('$Fx_{fl}$ [N]')
    xlim([0 time_sim(end)])
    
    % linkaxes(ax,'x')
    clear ax

    if test_type==1
        exportgraphics(gcf, 'Plots\Longitudinal_Slips_and_Forces.png')
    end

    % ---------------------------------
    %% Plot wheel torques and wheel rates
    % ---------------------------------
    figure('Name','Wheel rates & torques','NumberTitle','off'), clf
    % --- omega_rr --- %
    ax(1) = subplot(331);
    plot(time_sim,omega_rr,'LineWidth',2)
    grid on
    title('$\omega_{rr}$ [rad/s]')
    xlim([0 time_sim(end)])
    % --- omega_rl --- %
    ax(2) = subplot(332);
    plot(time_sim,omega_rl,'LineWidth',2)
    grid on
    title('$\omega_{rl}$ [rad/s]')
    xlim([0 time_sim(end)])
    % --- omega_fr --- %
    ax(3) = subplot(333);
    plot(time_sim,omega_fr,'LineWidth',2)
    grid on
    title('$\omega_{fr}$ [rad/s]')
    xlim([0 time_sim(end)])
    % --- omega_fl --- %
    ax(4) = subplot(334);
    plot(time_sim,omega_fl,'LineWidth',2)
    grid on
    title('$\omega_{fl}$ [rad/s]')
    xlim([0 time_sim(end)])
    % --- Tw_rr --- %
    ax(5) = subplot(335);
    plot(time_sim,Tw_rr,'LineWidth',2)
    grid on
    title('$Tw_{rr}$ [Nm]')
    xlim([0 time_sim(end)])
    % --- Tw_rl --- %
    ax(6) = subplot(336);
    plot(time_sim,Tw_rl,'LineWidth',2)
    grid on
    title('$Tw_{rl}$ [Nm]')
    xlim([0 time_sim(end)])
    % --- Tw_fr --- %
    ax(7) = subplot(337);
    plot(time_sim,Tw_fr,'LineWidth',2)
    grid on
    title('$Tw_{fr}$ [Nm]')
    xlim([0 time_sim(end)])
    % --- Tw_fl --- %
    ax(8) = subplot(338);
    plot(time_sim,Tw_fl,'LineWidth',2)
    grid on
    title('$Tw_{fl}$ [Nm]')
    xlim([0 time_sim(end)])

    % linkaxes(ax,'x')
    clear ax
    
    if test_type==1
        exportgraphics(gcf, 'Plots\Wheel_Rates_and_Torques.png')
    end

    % ---------------------------------
    %% Plot vertical tire loads and self-aligning torques
    % ---------------------------------
    figure('Name','Vert loads & aligning torques','NumberTitle','off'), clf
    % --- Fz_rr --- %
    ax(1) = subplot(331);
    plot(time_sim,Fz_rr,'LineWidth',2)
    grid on
    title('$Fz_{rr}$ [N]')
    xlim([0 time_sim(end)])
    % --- Fz_rl --- %
    ax(2) = subplot(332);
    plot(time_sim,Fz_rl,'LineWidth',2)
    grid on
    title('$Fz_{rl}$ [N]')
    xlim([0 time_sim(end)])
    % --- Fz_fr --- %
    ax(3) = subplot(333);
    plot(time_sim,Fz_fr,'LineWidth',2)
    grid on
    title('$Fz_{fr}$ [N]')
    xlim([0 time_sim(end)])
    % --- Fz_fl --- %
    ax(4) = subplot(334);
    plot(time_sim,Fz_fl,'LineWidth',2)
    grid on
    title('$Fz_{fl}$ [N]')
    xlim([0 time_sim(end)])
    % --- Mz_rr --- %
    ax(5) = subplot(335);
    plot(time_sim,Mz_rr,'LineWidth',2)
    grid on
    title('$Mz_{rr}$ [Nm]')
    xlim([0 time_sim(end)])
    % --- Mz_rl --- %
    ax(6) = subplot(336);
    plot(time_sim,Mz_rl,'LineWidth',2)
    grid on
    title('$Mz_{rl}$ [Nm]')
    xlim([0 time_sim(end)])
    % --- Mz_fr --- %
    ax(7) = subplot(337);
    plot(time_sim,Mz_fr,'LineWidth',2)
    grid on
    title('$Mz_{fr}$ [Nm]')
    xlim([0 time_sim(end)])
    % --- Mz_fl --- %
    ax(8) = subplot(338);
    plot(time_sim,Mz_fl,'LineWidth',2)
    grid on
    title('$Mz_{fl}$ [Nm]')
    xlim([0 time_sim(end)])

    % linkaxes(ax,'x')
    clear ax

    if test_type==1
        exportgraphics(gcf, 'Plots\Vertical_Loads_and_Self-Aligning_Torques.png')
    end

    % ---------------------------------
    %% Plot wheel camber
    % ---------------------------------
    figure('Name','Wheel camber','NumberTitle','off'), clf
    % --- gamma_rr --- %
    ax(1) = subplot(221);
    plot(time_sim,gamma_rr,'LineWidth',2)
    grid on
    title('$\gamma_{rr}$ [deg]')
    xlim([0 time_sim(end)])
    % --- gamma_rl --- %
    ax(2) = subplot(222);
    plot(time_sim,gamma_rl,'LineWidth',2)
    grid on
    title('$\gamma_{rl}$ [deg]')
    xlim([0 time_sim(end)])
    % --- gamma_fr --- %
    ax(3) = subplot(223);
    plot(time_sim,gamma_fr,'LineWidth',2)
    grid on
    title('$\gamma_{fr}$ [deg]')
    xlim([0 time_sim(end)])
    % --- gamma_fl --- %
    ax(4) = subplot(224);
    plot(time_sim,gamma_fl,'LineWidth',2)
    grid on
    title('$\gamma_{fl}$ [deg]')
    xlim([0 time_sim(end)])

    % linkaxes(ax,'x')
    clear ax

    if test_type==1
        exportgraphics(gcf, 'Plots\Wheel_Camber.png')
    end

    % ---------------------------------
    %% Plot accelerations, chassis side slip angle and curvature
    % ---------------------------------
    figure('Name','Pars extra','NumberTitle','off'), clf
    % --- ax --- %
    ax(1) = subplot(221);
    plot(time_sim(2:end),dot_u - Omega(2:end).*v(2:end),'LineWidth',2)
    hold on
    plot(time_sim(2:end),diff(u)/Ts,'--g','LineWidth',2)
    plot(time_sim(2:end),Ax_filt,'-.b','LineWidth',1)
    plot(time_sim(2:end),dot_u_filt,'-.r','LineWidth',1)
    grid on
    title('$a_{x}$ $[m/s^2]$')
    legend('$\dot{u}-\Omega v$','$\dot{u}$','filt $\dot{u}-\Omega v$','filt $\dot{u}$','Location','northeast')
    xlim([0 time_sim(end)])
    % --- ay --- %
    ax(2) = subplot(222);
    plot(time_sim(2:end),dot_v + Omega(2:end).*u(2:end),'LineWidth',2)
    hold on
    plot(time_sim(2:end),Omega(2:end).*u(2:end),'--g','LineWidth',1)
    grid on
    title('$a_{y}$ $[m/s^2]$')
    legend('$\dot{v}+\Omega u$','$\Omega u$','Location','best')
    xlim([0 time_sim(end)])
    % --- beta --- %
    ax(3) = subplot(223);
    plot(time_sim,rad2deg(beta),'LineWidth',2)
    grid on
    title('$\beta$ [deg]')
    xlim([0 time_sim(end)])
    % --- rho --- %
    ax(4) = subplot(224);
    plot(time_sim,rho_ss,'LineWidth',2)
    hold on
    plot(time_sim(1:end-1),rho_tran,'--g','LineWidth',1)
    grid on
    title('$\rho$ [$m^{-1}$]')
    legend('$\rho_{ss}$','$\rho_{transient}$','Location','best')
    xlim([0 time_sim(end)])

    % linkaxes(ax,'x')
    clear ax

    if test_type==1
        exportgraphics(gcf, 'Plots\Accelerations,_Side_Slip_Angle_and_Curvature.png')
    end

    % ---------------------------------
    %% Plot vehicle pose x,y,psi
    % ---------------------------------
    figure('Name','Pose','NumberTitle','off'), clf 
    % --- x --- %
    ax(1) = subplot(221);
    plot(time_sim,x_CoM,'LineWidth',2)
    grid on
    title('$x$ [m]')
    xlim([0 time_sim(end)])
    % --- y --- %
    ax(2) = subplot(222);
    plot(time_sim,y_CoM,'LineWidth',2)
    grid on
    title('$y$ [m]')
    xlim([0 time_sim(end)])
    % --- psi --- %
    ax(3) = subplot(223);
    plot(time_sim,rad2deg(psi),'LineWidth',2)
    grid on
    title('$\psi$ [deg]')
    xlim([0 time_sim(end)])

    % linkaxes(ax,'x')
    clear ax

    if test_type==1
        exportgraphics(gcf, 'Plots\Vehicle_Pose.png')
    end

    % -------------------------------
    %% Plot G-G diagram from simulation data
    % -------------------------------
    figure('Name','G-G plot','NumberTitle','off'), clf
    axis equal
    hold on
    plot3(Ay,Ax_filt,u(1:end-1),'Color',color('purple'),'LineWidth',3)
    xlabel('$a_y$ [m/s$^2$]')
    ylabel('$a_x$ [m/s$^2$]')
    zlabel('$u$ [m/s]')
    title('G-G diagram from simulation data','FontSize',18)
    grid on

    if test_type==1
        exportgraphics(gcf, 'Plots\G-G_Diagram.png')
    end

    % -------------------------------
    %% Plot vehicle path
    % -------------------------------
    N = length(time_sim);
    figure('Name','Real Vehicle Path','NumberTitle','off'), clf
    set(gca,'fontsize',16)
    hold on
    axis equal
    xlabel('x-coord [m]')
    ylabel('y-coord [m]')
    title('Real Vehicle Path','FontSize',18)
    plot(x_CoM,y_CoM,'Color',color('gold'),'LineWidth',2)
    for i = 1:floor(N/20):N
        rot_mat = [cos(psi(i)) -sin(psi(i)) ; sin(psi(i)) cos(psi(i))];
        pos_rr = rot_mat*[-Lr -Wr/2]';
        pos_rl = rot_mat*[-Lr +Wr/2]';
        pos_fr = rot_mat*[+Lf -Wf/2]';
        pos_fl = rot_mat*[+Lf +Wf/2]';
        pos = [pos_rr pos_rl pos_fl pos_fr];
        p = patch(x_CoM(i) + pos(1,:),y_CoM(i) + pos(2,:),'blue');
        quiver(x_CoM(i), y_CoM(i), u(i)*cos(psi(i)), u(i)*sin(psi(i)), 'color', [1,0,0]);
        quiver(x_CoM(i), y_CoM(i), -v(i)*sin(psi(i)), v(i)*cos(psi(i)), 'color', [0.23,0.37,0.17]);
    end
    grid on
    hold off
    
    if test_type==1
        exportgraphics(gcf, 'Plots\Speed_Test_Vehicle_Path.png')
    else
        exportgraphics(gcf, 'Plots\Steer_Test_Vehicle_Path.png')
    end

    %% Plot Load Transfer 
    t_s=time_sim(1:length(Delta_Fzr_tot0),1);
    figure ('Name', 'Lat. load transfer','NumberTitle','off')
    hold on
    plot(t_s,Delta_Fzr_tot0, 'DisplayName', 'Rear', 'Color','r');
    plot(t_s,Delta_Fzf_tot0, 'DisplayName', 'Front','Color','b');
    xlim([0 time_sim(end)])
    legend('Location','eastoutside')
    xlabel("Time [s]")
    ylabel("Lateral load transfer [N]")
    title("Lateral load transfer")
    grid on

    if test_type==1
        exportgraphics(gcf, 'Plots\Lateral_Load_Transfer_Speed.png')
    end
    if test_type==2
        exportgraphics(gcf, 'Plots\Lateral_Load_Transfer_Steer.png')
    end

    %% Plot normalized axle characteristic
    if test_type==2
        figure ('Name', 'Normalized axle characteristic','NumberTitle','off')
        hold on
        plot(mi_r0(1:length(Ay),1),Ay./g, 'DisplayName', 'Rear Norm. axle char', 'Color','r');
        plot(mi_f0(1:length(Ay),1),Ay./g, 'DisplayName', 'Front Norm. axle char','Color','b');
        plot(mi_r1(1:length(Ay),1),Ay./g, 'DisplayName', 'Rear Norm. axle char High Front', 'Color','r', 'LineStyle', '--');
        plot(mi_f1(1:length(Ay),1),Ay./g, 'DisplayName', 'Front Norm. axle char High Front','Color','b', 'LineStyle', '--');
        plot(mi_r2(1:length(Ay),1),Ay./g, 'DisplayName', 'Rear Norm. axle char High Rear', 'Color','r', 'LineStyle', ':');
        plot(mi_f2(1:length(Ay),1),Ay./g, 'DisplayName', 'Front Norm. axle char High Rear','Color','b', 'LineStyle', ':');
        legend('Location','best')
        ylabel("Fyr,Fyf")
        xlabel("$\alpha_r,\alpha_f$")
        title("Normalized axle characteristic")
        grid on

        exportgraphics(gcf, 'Plots\Normalized_Axle_Characteristics.png')
    end

    %% Rho_tot_=Rho_tot(1:length(Ay),1);
    if test_type==1
        figure ('Name', 'Curvature vs lateral acc','NumberTitle','off')
        hold on
        plot(Ay(dataset_start:length(Ay),1)./g,rho_ss(dataset_start:length(Ay),1), ...
            'DisplayName', 'Steering angle vs lateral acc', 'Color','b');
        ylabel("$\rho$")
        xlabel("$A_y/g$")
        title("Curvature vs lateral acc")
        grid on

        exportgraphics(gcf, 'Plots\Curvature_against_Lateral_Acceleration.png')
    
        figure ('Name', 'Curvature vs time','NumberTitle','off')
        hold on
        plot(t_s(dataset_start:length(t_s),1),rho_ss(dataset_start:length(t_s),1), ...
            'DisplayName', 'Curvature vs time', 'Color','b');
        ylabel("$\rho$")
        xlabel("$t$")
        title("Curvature")
        grid on

        exportgraphics(gcf, 'Plots\Curvature_against_Time.png')
    end


    %% Handling diagram

    figure ('Name', 'Handling curve','NumberTitle','off')
    hold on
    plot(Ay(1:length(Ay),1)./g, oc_rhs_ss(1:length(Ay),1), ...
        'DisplayName', 'Handling Curve ($\delta _H / \tau - \rho_{ss} * L$)', ...
        'LineWidth', 2);
    plot(Ay(1:length(Ay),1)./g, hand_lin, ...
        'DisplayName', 'Linearized Handling Curve');
    plot(Ay(1:length(Ay),1)./g, hand_fit, ...
        'DisplayName', 'Fitted Handling Curve', ...
        'Color', 'k', 'LineStyle', '--', 'LineWidth', 2);
    legend('Location','best')
    ylabel("$\delta _H / \tau - \rho L$")
    xlabel("$A_y/g$")
    title("Handling diagram")
    grid on

    if test_type==1
        exportgraphics(gcf, 'Plots\Handling_Diagram.png')
    end

    %% Yaw rate
    if test_type==1
        figure ('Name', 'Yaw rate','NumberTitle','off')
        hold on
        plot(u*3.6, Omega./delta_D, 'DisplayName', 'Yaw rate', 'Color','b');
        xline(u(end)*3.6+2, 'DisplayName', 'Critical Speed $u_{cr}$', 'LineStyle', '--', 'Color', 'k')
        legend('Location','best')
        ylabel("$\Omega / \delta$")
        xlabel("$u$")
        title("Yaw rate")
        grid on

        exportgraphics(gcf, 'Plots\Yaw_Rate.png')
    end

    %% Body slip gain
    if test_type==2
        figure ('Name', 'Body slip gain','NumberTitle','off')
        hold on
        plot(u(10000:length(Ay),1)*3.6, ...
            rad2deg(beta(10000:length(Ay),1))./desired_steer_atWheel(10000:length(Ay),1), ...
            'DisplayName', 'Body slip gain', 'Color','b');
        xline(u(end)*3.6+2, 'DisplayName', 'Critical Speed $u_{cr}$', 'LineStyle', '--', 'Color', 'k')
        legend('Location','best')
        ylabel("$\beta / \delta$")
        xlabel("$u$")
        title("Body slip gain")
        grid on

        exportgraphics(gcf, 'Plots\Body_Slip_Gain.png')
    
        figure ('Name', 'Body slip gain vs \alpha _y','NumberTitle','off')
        hold on
        plot(Ay(1:end,1),rad2deg(beta(1:length(Ay),1)), 'DisplayName', 'Body slip gain', 'Color','b');
        legend('Location','best')
        ylabel("$\beta$")
        xlabel("$Ay$")
        title("Body slip gain vs $\alpha_y$ ")
        grid on

        exportgraphics(gcf, 'Plots\Body_Slip_Gain_vs_Lateral_Acceleration.png')
    end
end
    
