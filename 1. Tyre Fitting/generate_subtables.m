function [tables] = generate_subtables(tyre_data, vec_samples, name)
    
    to_rad = pi/180;
    to_deg = 180/pi;
    
    % Constant Inclination Angle GAMMA
    GAMMA_tol = 0.05*to_rad;
    idx.GAMMA_0 = 0.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 0.0*to_rad+GAMMA_tol;
    idx.GAMMA_1 = 1.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 1.0*to_rad+GAMMA_tol;
    idx.GAMMA_2 = 2.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 2.0*to_rad+GAMMA_tol;
    idx.GAMMA_3 = 3.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 3.0*to_rad+GAMMA_tol;
    idx.GAMMA_4 = 4.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 4.0*to_rad+GAMMA_tol;
    idx.GAMMA_5 = 5.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 5.0*to_rad+GAMMA_tol;
    GAMMA_0  = tyre_data( idx.GAMMA_0, : );
    GAMMA_1  = tyre_data( idx.GAMMA_1, : );
    GAMMA_2  = tyre_data( idx.GAMMA_2, : );
    GAMMA_3  = tyre_data( idx.GAMMA_3, : );
    GAMMA_4  = tyre_data( idx.GAMMA_4, : );
    GAMMA_5  = tyre_data( idx.GAMMA_5, : );

    constant_GAMMA = struct('GAMMA_0', GAMMA_0, ...
                            'GAMMA_1', GAMMA_1, ...
                            'GAMMA_2', GAMMA_2, ...
                            'GAMMA_3', GAMMA_3, ...
                            'GAMMA_4', GAMMA_4, ...
                            'GAMMA_5', GAMMA_5);
    
    % Constant Vertical Load FZ
    FZ_tol = 100;
    idx.FZ_220  = 220-FZ_tol < tyre_data.FZ & tyre_data.FZ < 220+FZ_tol;
    idx.FZ_440  = 440-FZ_tol < tyre_data.FZ & tyre_data.FZ < 440+FZ_tol;
    idx.FZ_700  = 700-FZ_tol < tyre_data.FZ & tyre_data.FZ < 700+FZ_tol;
    idx.FZ_900  = 900-FZ_tol < tyre_data.FZ & tyre_data.FZ < 900+FZ_tol;
    idx.FZ_1120 = 1120-FZ_tol < tyre_data.FZ & tyre_data.FZ < 1120+FZ_tol;
    FZ_220  = tyre_data( idx.FZ_220, : );
    FZ_440  = tyre_data( idx.FZ_440, : );
    FZ_700  = tyre_data( idx.FZ_700, : );
    FZ_900  = tyre_data( idx.FZ_900, : );
    FZ_1120 = tyre_data( idx.FZ_1120, : );

    constant_FZ = struct('FZ_220', FZ_220, ...
                         'FZ_440', FZ_440, ...
                         'FZ_700', FZ_700, ...
                         'FZ_900', FZ_900, ...
                         'FZ_1120', FZ_1120);
    
    % Constant Side Slip Angle ALPHA
    SA_tol = 0.5*to_rad;
    idx.SA_0    =  0-SA_tol          < tyre_data.SA & tyre_data.SA < 0+SA_tol;
    idx.SA_3neg = -(3*to_rad+SA_tol) < tyre_data.SA & tyre_data.SA < -3*to_rad+SA_tol;
    idx.SA_6neg = -(6*to_rad+SA_tol) < tyre_data.SA & tyre_data.SA < -6*to_rad+SA_tol;
    SA_0     = tyre_data( idx.SA_0, : );
    SA_3neg  = tyre_data( idx.SA_3neg, : );
    SA_6neg  = tyre_data( idx.SA_6neg, : );

    constant_ALPHA = struct('SA_0', SA_0, ...
                           'SA_3neg', SA_3neg, ...
                           'SA_6neg', SA_6neg);
    
    % Constant Slip Angle KAPPA
    SL_tol = 0.05;
    idx.SL_0      =    0-SL_tol  < tyre_data.SL & tyre_data.SL <     0+SL_tol;
    idx.SL_02pos  =  0.2-SL_tol  < tyre_data.SL & tyre_data.SL <   0.2+SL_tol;
    idx.SL_02neg  = -0.2-SL_tol  < tyre_data.SL & tyre_data.SL <  -0.2+SL_tol;
    idx.SL_01pos  =  0.1-SL_tol  < tyre_data.SL & tyre_data.SL <   0.1+SL_tol;
    idx.SL_01neg  = -0.1-SL_tol  < tyre_data.SL & tyre_data.SL <  -0.1+SL_tol;
    SL_0      = tyre_data( idx.SL_0, : );
    SL_02pos  = tyre_data( idx.SL_02pos, : );
    SL_02neg  = tyre_data( idx.SL_02neg, : );
    SL_01pos  = tyre_data( idx.SL_01pos, : );
    SL_01neg  = tyre_data( idx.SL_01neg, : );

    constant_KAPPA = struct('SL_0',SL_0, ...
                            'SL_02pos', SL_02pos, ...
                            'SL_02neg', SL_02neg, ...
                            'SL_01pos', SL_01pos, ...
                            'SL_01neg', SL_01neg);

    tables = struct('GAMMA', constant_GAMMA, ...
                    'FZ', constant_FZ, ...
                    'ALPHA', constant_ALPHA, ...
                    'KAPPA', constant_KAPPA);

    %% Plots

    figure('Name', name)
    tiledlayout(4,1)
    
    ax_list(1) = nexttile;
    plot(tyre_data.IA*to_deg)
    hold on
    plot(vec_samples(idx.GAMMA_0),GAMMA_0.IA*to_deg,'.');
    plot(vec_samples(idx.GAMMA_1),GAMMA_1.IA*to_deg,'.');
    plot(vec_samples(idx.GAMMA_2),GAMMA_2.IA*to_deg,'.');
    plot(vec_samples(idx.GAMMA_3),GAMMA_3.IA*to_deg,'.');
    plot(vec_samples(idx.GAMMA_4),GAMMA_4.IA*to_deg,'.');
    plot(vec_samples(idx.GAMMA_5),GAMMA_5.IA*to_deg,'.');
    title('Camber angle')
    ylabel('[deg]')
    
    ax_list(2) = nexttile;
    plot(tyre_data.FZ)
    hold on
    plot(vec_samples(idx.FZ_220),FZ_220.FZ,'.');
    plot(vec_samples(idx.FZ_440),FZ_440.FZ,'.');
    plot(vec_samples(idx.FZ_700),FZ_700.FZ,'.');
    plot(vec_samples(idx.FZ_900),FZ_900.FZ,'.');
    plot(vec_samples(idx.FZ_1120),FZ_1120.FZ,'.');
    title('Vertical force')
    ylabel('[N]')
    
    ax_list(3) = nexttile;
    plot(tyre_data.SA*to_deg)
    hold on
    plot(vec_samples(idx.SA_0),   SA_0.SA*to_deg,'.');
    plot(vec_samples(idx.SA_3neg),SA_3neg.SA*to_deg,'.');
    plot(vec_samples(idx.SA_6neg),SA_6neg.SA*to_deg,'.');
    title('Side slip')
    ylabel('[rad]')

    ax_list(4) = nexttile;
    plot(tyre_data.SL)
    hold on
    plot(vec_samples(idx.SL_0),SL_0.SL,'.');
    plot(vec_samples(idx.SL_02pos),SL_02pos.SL,'.');
    plot(vec_samples(idx.SL_02neg),SL_02neg.SL,'.');
    plot(vec_samples(idx.SL_01pos),SL_01pos.SL,'.');
    plot(vec_samples(idx.SL_01neg),SL_01neg.SL,'.');
    title('Longitudinal Slip')
    xlabel('Samples [-]')
    ylabel('[rad]')

    sgtitle(name)

end