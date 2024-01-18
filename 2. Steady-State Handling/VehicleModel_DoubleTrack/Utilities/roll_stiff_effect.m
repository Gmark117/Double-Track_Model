function roll_stiff_effect(model_sim_roll_stiff, vehicle_data, Ts, roll_stiff)

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
    % for i=roll_stiff
    %     if i<0
    %         fieldname = strcat('roll_stiff_min', num2str(abs(i)));
    %     else
    %         fieldname = strcat('roll_stiff_', num2str(i));
    %     end
    % 
    %     time_sim = model_sim_roll_stiff.(fieldname).states.u.time;
    %     dt = time_sim(2)-time_sim(1);
    % end
    
    % -----------------
    %% Inputs
    % -----------------
    % for i=roll_stiff
    %     if i<0
    %         fieldname = strcat('roll_stiff_min', num2str(abs(i)));
    %     else
    %         fieldname = strcat('roll_stiff_', num2str(i));
    %     end
    % 
    %     ped_0      = model_sim_roll_stiff.(fieldname).inputs.ped_0.data;
    %     delta_D    = model_sim_roll_stiff.(fieldname).inputs.delta_D.data;
    % end
    % 

    % -----------------
    %% States
    % -----------------
    for i=roll_stiff
        fieldname = strcat('roll_stiff_', num2str(i*10));

        roll_stiff_states.(fieldname) = model_sim_roll_stiff.(fieldname).states;
    end
    
    % -----------------
    %% Extra Params
    % -----------------
    % for i=roll_stiff
    %     if i<0
    %         fieldname = strcat('roll_stiff_min', num2str(abs(i)));
    %     else
    %         fieldname = strcat('roll_stiff_', num2str(i));
    %     end
    % 
    %     roll_stiff_extraParams.(fieldname) = model_sim_roll_stiff.(fieldname).extra_params;
    % end

    % -----------------
    %% Post-Processing
    % -----------------
    % desired_steer_atWheel = delta_D/tau_D;

    for i=roll_stiff
        fieldname = strcat('roll_stiff_', num2str(i*10));

        dot_v.(fieldname) = diff(roll_stiff_states.(fieldname).v.data)./Ts;
        % dot_u.(fieldname) = diff(roll_stiff_states.(fieldname).u)/Ts;

        Ay.(fieldname) = dot_v.(fieldname) + roll_stiff_states.(fieldname).Omega.data(2:end).*roll_stiff_states.(fieldname).u.data(2:end);

        Dalpha.(fieldname) = roll_stiff_states.(fieldname).alpha_rr - roll_stiff_states.(fieldname).alpha_fr;
    end

    % ---------------------------------
    %% Plots
    % ---------------------------------
    figure('Name', 'Roll Stiffness Effect', 'NumberTitle', 'off'), clf
    for i=roll_stiff
        hold on
        grid on

        fieldname = strcat('roll_stiff_', num2str(i*10));
        
        plot(Ay.(fieldname)/g, -Dalpha.(fieldname).data(1:length(Ay.(fieldname)),1), 'DisplayName', ['$\epsilon_\phi = $', num2str(i)])
    end
    
    xlabel('$\frac{a_y}{g}$')
    ylabel('$-\Delta\alpha$')
    legend()
    title('Handling Diagram varying $\epsilon_\phi$')

    exportgraphics(gcf, 'Plots\Roll_Stiffness_Effect.png')

end