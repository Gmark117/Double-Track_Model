function camber_effect(model_sim_camber, vehicle_data, Ts, camber)

    % ---------------------------------
    % Load vehicle data
    % ---------------------------------
    % Lf = vehicle_data.vehicle.Lf;  % [m] Distance between vehicle CoG and front wheels axle
    % Lr = vehicle_data.vehicle.Lr;  % [m] Distance between vehicle CoG and front wheels axle
    % L  = vehicle_data.vehicle.L;   % [m] Vehicle length
    % Wf = vehicle_data.vehicle.Wf;  % [m] Width of front wheels axle 
    % Wr = vehicle_data.vehicle.Wr;  % [m] Width of rear wheels axle                   
    % m  = vehicle_data.vehicle.m;   % [kg] Vehicle Mass
    g  = vehicle_data.vehicle.g;   % [m/s^2] Gravitational acceleration
    % tau_D = vehicle_data.steering_system.tau_D;  % [-] steering system ratio (pinion-rack)
    % h_rc_f = vehicle_data.front_suspension.h_rc_f; % roll center height front
    % h_rc_r = vehicle_data.rear_suspension.h_rc_r; % roll center height rear
    % hGs = vehicle_data.vehicle.hGs;
    % h_r = h_rc_r + (h_rc_f - h_rc_r)*Lr/L;
    % h_s = hGs - h_r;
    % Ks_r=vehicle_data.rear_suspension.Ks_r;
    % Karb_r = vehicle_data.rear_suspension.Karb_r;
    % Ks_f=vehicle_data.front_suspension.Ks_f;
    % Karb_f = vehicle_data.front_suspension.Karb_f;

    % ---------------------------------
    %% Extract data from simulink model
    % ---------------------------------
    % for i=camber
    %     if i<0
    %         fieldname = strcat('camber_min', num2str(abs(i)));
    %     else
    %         fieldname = strcat('camber_', num2str(i));
    %     end
    % 
    %     time_sim = model_sim_camber.(fieldname).states.u.time;
    %     dt = time_sim(2)-time_sim(1);
    % end
    
    % -----------------
    %% Inputs
    % -----------------
    % for i=camber
    %     if i<0
    %         fieldname = strcat('camber_min', num2str(abs(i)));
    %     else
    %         fieldname = strcat('camber_', num2str(i));
    %     end
    % 
    %     ped_0      = model_sim_camber.(fieldname).inputs.ped_0.data;
    %     delta_D    = model_sim_camber.(fieldname).inputs.delta_D.data;
    % end
    % 

    % -----------------
    %% States
    % -----------------
    for i=camber
        if i<0
            fieldname = strcat('camber_min', num2str(abs(i)));
        else
            fieldname = strcat('camber_', num2str(i));
        end

        camber_states.(fieldname) = model_sim_camber.(fieldname).states;
    end
    
    % -----------------
    %% Extra Params
    % -----------------
    % for i=camber
    %     if i<0
    %         fieldname = strcat('camber_min', num2str(abs(i)));
    %     else
    %         fieldname = strcat('camber_', num2str(i));
    %     end
    % 
    %     camber_extraParams.(fieldname) = model_sim_camber.(fieldname).extra_params;
    % end

    % -----------------
    %% Post-Processing
    % -----------------
    % desired_steer_atWheel = delta_D/tau_D;

    for i=camber
        if i<0
            fieldname = strcat('camber_min', num2str(abs(i)));
        else
            fieldname = strcat('camber_', num2str(i));
        end

        dot_v.(fieldname) = diff(camber_states.(fieldname).v.data)./Ts;
        % dot_u.(fieldname) = diff(camber_states.(fieldname).u)/Ts;

        Ay_ss.(fieldname) = camber_states.(fieldname).Omega.data(2:end) .* camber_states.(fieldname).u.data(2:end);
        Ay.(fieldname) = dot_v.(fieldname) + Ay_ss.(fieldname);
        
        alpha_r = 0.5*(camber_states.(fieldname).alpha_rr + camber_states.(fieldname).alpha_rl);
        alpha_f = 0.5*(camber_states.(fieldname).alpha_fr + camber_states.(fieldname).alpha_fl);
        Dalpha.(fieldname) = alpha_r - alpha_f;
    end

    % ---------------------------------
    %% Plots
    % ---------------------------------
    figure('Name', 'Camber Effect', 'NumberTitle', 'off'), clf
    for i=camber
        hold on
        grid on

        if i<0
            fieldname = strcat('camber_min', num2str(abs(i)));
        else
            fieldname = strcat('camber_', num2str(i));
        end
        
        
        plot(Ay.(fieldname)/g,  - Dalpha.(fieldname).data(1:length(Ay.(fieldname)),1) ...
                                + Dalpha.(fieldname).data(30000) ...
                                - Dalpha.camber_0.data(30000), ...
            'DisplayName', ['$\gamma = $', num2str(i)])
    end
    
    xlabel('$\frac{a_y}{g}$')
    ylabel('$-\Delta\alpha$')
    legend()
    title('Handling Diagram varying $\gamma$')

    exportgraphics(gcf, 'Plots\Camber_Angle_Effect.png')

end