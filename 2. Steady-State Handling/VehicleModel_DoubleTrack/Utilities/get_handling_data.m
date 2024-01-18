function [handling_data] = get_handling_data(model_sim, vehicle_data)

    % Lf = vehicle_data.vehicle.Lf;
    % Lr = vehicle_data.vehicle.Lr;
    % L  = vehicle_data.vehicle.L;
    % m  = vehicle_data.vehicle.m;
    g  = vehicle_data.vehicle.g;

    % time_sim = model_sim.states.u.time;
    u        = model_sim.states.u.data;
    Omega    = model_sim.states.Omega.data;
    Fz_rr    = model_sim.states.Fz_rr.data;
    Fz_rl    = model_sim.states.Fz_rl.data;
    Fz_fr    = model_sim.states.Fz_fr.data;
    Fz_fl    = model_sim.states.Fz_fl.data;
    alpha_rr = model_sim.states.alpha_rr.data;
    alpha_rl = model_sim.states.alpha_rl.data;
    alpha_fr = model_sim.states.alpha_fr.data;
    alpha_fl = model_sim.states.alpha_fl.data;
    % Fy_rr    = model_sim.states.Fy_rr.data;
    % Fy_rl    = model_sim.states.Fy_rl.data;
    % Fy_fr    = model_sim.states.Fy_fr.data;
    % Fy_fl    = model_sim.states.Fy_fl.data;

    Ay_ss   = Omega.*u;
    dFz_f   = 0.5.*(Fz_fr-Fz_fl);
    dFz_r   = 0.5.*(Fz_rr-Fz_rl);
    alpha_f = 0.5.*deg2rad(alpha_fr+alpha_fl);
    alpha_r = 0.5.*deg2rad(alpha_rr+alpha_rl);
    Dalpha  = alpha_r - alpha_f;

    handling_data.Ay     = Ay_ss./g;
    handling_data.Dalpha = Dalpha;
    handling_data.alpha_r = alpha_r;
    handling_data.alpha_f = alpha_f;
    handling_data.dFz_f  = dFz_f;
    handling_data.dFz_r  = dFz_r;

end