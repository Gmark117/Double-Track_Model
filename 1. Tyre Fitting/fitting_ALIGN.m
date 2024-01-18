%% Fitting with Fz = Fz$_700$ and $\gamma$ = 0, $\alpha$ = 0

fitting_params = {'qBz1';'qBz9';'qBz10';'qCz1';'qDz1'; 'qDz2' ;'qDz3'; 'qDz4';'qDz6'; 'qEz1';'qEz4';'qHz1'}; 

% Tables Intersection
[Tdata0, ~] = intersect_table_data(Tsubs_LAT.KAPPA.SL_0, ...
                                   Tsubs_LAT.GAMMA.GAMMA_0, ...
                                   Tsubs_LAT.FZ.FZ_700);

plot_selected_data(Tdata0, 'Aligning Moment Selected Data')

% Initialization
KAPPA           =  Tdata0.SL;
ALPHA           =  Tdata0.SA;
GAMMA           =  Tdata0.IA;
% FZ              =  Tdata0.FZ;
MZ              = -Tdata0.MZ;
% FY              =  Tdata0.FY;
FZ0             = mean(Tsubs_LAT.FZ.FZ_700.FZ);
tyre_coeffs.FZ0 = 700;
SA_vec          = -0.3:0.001:0.3;
zeros_ALPHA     = zeros(size(ALPHA));
ones_ALPHA      = ones(size(ALPHA));
% zeros_SA        = zeros(size(SA_vec));
ones_SA         = ones(size(SA_vec));

% Initial Guess
%    [qBz1  qBz9  qBz10   qCz1    qDz1   qDz2     qDz3  qDz4  qDz6  qEz1   qEz4     qHz1]
% P0 = [    5,    10,     0,     1,  0.01,   -0.001,   0.1,      0,     0,     -2,   0.1,    0 ];
% P0 = [   9,    13,     1,   1.2, -0.13,      0.4, -0.22,      0,     0,     -1,  -0.4,-0.003 ]; %R=0.931
P0 = [  10,   13,  1.28,  1.28,  -0.16,  -0.2,  -0.032,    0,    0,    0,  0.28,  -0.005];
% P0 = [  30,    13,   1.4,     1, -0.13,     0.18, -0.06,      0,  0.11,    1.6,    -1,-0.005 ];
% lb = [  10,     5,    -1,     0,    -5,       -1,    -2,     -1,    -1,     -7,    -2,    -1 ];
% ub = [  40,    20,     3,     2,     0,        2,     2,      1,     1,      5,     0,     1 ];
lb = [   7,   10,   0.5,   0.5,     -1,    -1,      -1,   -1,   -1,   -2,   -10,      -1];
ub = [  15,   16,     2,     2,      0,     1,       1,    1,    1,    0,     3,       0];
% lb=[];
% ub=[];
% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised

% First Guess
for i = 1:length(fitting_params)
    tyre_coeffs.(fitting_params{i}) = P0(i);
end 
MZ0_guess = MF96_MZ0_vec(zeros_ALPHA, ...
                         ALPHA, ...
                         zeros_ALPHA, ...
                         tyre_coeffs.FZ0*ones_ALPHA, ...
                         tyre_coeffs);

% Optimization
[P_fz_nom, ~] = fmincon(@(P)resid_pure_Mz(fitting_params, ...
                                          P, ...
                                          MZ, ...
                                          KAPPA, ...
                                          ALPHA, ...
                                          zeros_ALPHA, ...
                                          FZ0*ones_ALPHA, ...
                                          tyre_coeffs), ...
                        P0, [], [], [], [], lb, ub);

% Optimal Parameters Assignment
for i = 1:length(fitting_params)
%     tyre_coeffs.(fitting_params{i}) = P_fz_nom(i);
    tyre_coeffs.(fitting_params{i}) = P_fz_nom(i);
end

% Force and Residuals Computation with Optimal Parameters
MZ0_fz_nom = MF96_MZ0_vec(zeros_ALPHA, ...
                          SA_vec, ...
                          zeros_ALPHA, ...
                          FZ0.*ones_SA, ...
                          tyre_coeffs);

% Plotting
figure('Name', 'First Guess Aligning Moment Nominal FZ')
plot(ALPHA,MZ,'.')
hold on
plot(ALPHA,MZ0_guess,'-', 'LineWidth',2)
title('First Guess Aligning Moment Nominal FZ')

% -------------------------------------------------------------------------

figure('Name', 'Fitted Aligning Moment Nominal FZ')
plot(ALPHA,MZ,'o')
hold on
plot(SA_vec,MZ0_fz_nom,'-','LineWidth',2)
xlabel('$\alpha$ [-]')
ylabel('$M_{z0}$ [Nm]')
title('Fitted Aligning Moment Nominal FZ')

% Conclusions
res_fz_nom = resid_pure_Mz(fitting_params, ...
                           P_fz_nom, ...
                           MZ, ...
                           KAPPA, ...
                           ALPHA, ...
                           GAMMA, ...
                           FZ0*ones_ALPHA, ...
                           tyre_coeffs);
r_squared = 1-res_fz_nom;
disp(['Lower boundary:     ', num2str(lb      , "%.3f     ")])
disp(['Upper boundary:     ', num2str(ub      , "%.3f     ")])
disp(['Optimal solution:   ', num2str(P_fz_nom, "%.3f      ")])
fprintf('Self-Aligning Moment R-squared FZ_nom = %6.3f\n', r_squared);

%% Fitting with Variable Fz and $\gamma$ = 0, $\alpha$ = 0
fitting_params = {'qBz2';'qBz3';'qDz7';'qEz2';'qEz3'};

% Tables Intersection
[TdataDFz, ~] = intersect_table_data(Tsubs_LAT.KAPPA.SL_0, ...
                                     Tsubs_LAT.GAMMA.GAMMA_0);

% Initialization
KAPPA           = TdataDFz.SL;
ALPHA           = TdataDFz.SA;
GAMMA           = TdataDFz.IA;
FZ              = TdataDFz.FZ;
% FY              = TdataDFz.FY;
MZ              = -TdataDFz.MZ;
FZ0             = mean(Tsubs_LAT.FZ.FZ_220.FZ);
FZ1             = mean(Tsubs_LAT.FZ.FZ_440.FZ);
FZ2             = mean(Tsubs_LAT.FZ.FZ_700.FZ);
FZ3             = mean(Tsubs_LAT.FZ.FZ_900.FZ);
FZ4             = mean(Tsubs_LAT.FZ.FZ_1120.FZ);
tyre_coeffs.FZ0 = 700;
SA_vec          = -0.3:0.001:0.3;
zeros_ALPHA     = zeros(size(ALPHA));
ones_ALPHA      = ones(size(ALPHA));
zeros_SA        = zeros(size(SA_vec));
ones_SA         = ones(size(SA_vec));

% Initial Guess
%    [ qBz2  qBz3  qDz7  qEz2    qEz3] 
% P0 = [ -1.8, 5,   0, 5.5,  -7];
P0 = [  1.9,  0.2,    0,    0,  -0.01];
lb = [   -5,   -5,   -1,   -1,     -1];
ub = [    5,    5,    1,    1,      1];
% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised

% First Guess
for i = 1:length(fitting_params)
    tyre_coeffs.(fitting_params{i}) = P0(i);
end 
MZ0_guess = MF96_MZ0_vec(zeros_ALPHA, ...
                         ALPHA, ...
                         zeros_ALPHA, ...
                         tyre_coeffs.FZ0*ones_ALPHA, ...
                         tyre_coeffs);

% Optimization
[P_DFz, ~] = fmincon(@(P)resid_pure_Mz(fitting_params, ...
                                       P, ...
                                       MZ, ...
                                       KAPPA, ...
                                       ALPHA, ...
                                       GAMMA, ...
                                       FZ, ...
                                       tyre_coeffs), ...
                     P0, [], [], [], [], lb, ub);

% Optimal Parameters Assignment
for i = 1:length(fitting_params)
    tyre_coeffs.(fitting_params{i}) = P_DFz(i);
end

% Force and Residuals Computation with Optimal Parameters
MZ0_DFz0 = MF96_MZ0_vec(zeros_SA, ...
                        SA_vec, ...
                        zeros_SA, ...
                        FZ0.*ones_SA, ...
                        tyre_coeffs);
MZ0_DFz1 = MF96_MZ0_vec(zeros_SA, ...
                        SA_vec, ...
                        zeros_SA, ...
                        FZ1.*ones_SA, ...
                        tyre_coeffs);
MZ0_DFz2 = MF96_MZ0_vec(zeros_SA, ...
                        SA_vec, ...
                        zeros_SA, ...
                        FZ2.*ones_SA, ...
                        tyre_coeffs);
MZ0_DFz3 = MF96_MZ0_vec(zeros_SA, ...
                        SA_vec, ...
                        zeros_SA, ...
                        FZ3.*ones_SA, ...
                        tyre_coeffs);
MZ0_DFz4 = MF96_MZ0_vec(zeros_SA, ...
                        SA_vec, ...
                        zeros_SA, ...
                        FZ4.*ones_SA, ...
                        tyre_coeffs);

% Plotting
figure('Name', 'First Guess Aligning Moment Variable FZ')
plot(ALPHA,MZ,'.')
hold on
plot(ALPHA,MZ0_guess,'-', 'LineWidth',2)
title('First Guess Aligning Moment Variable FZ')

% -------------------------------------------------------------------------

figure('Name', 'Fitted Aligning Moment Variable FZ')
plot(ALPHA,MZ,'o')
hold on
plot(SA_vec,MZ0_DFz0, 'DisplayName', '220','LineWidth',2)
plot(SA_vec,MZ0_DFz1, 'DisplayName', '440','LineWidth',2)
plot(SA_vec,MZ0_DFz2, 'DisplayName', '700','LineWidth',2)
plot(SA_vec,MZ0_DFz3, 'DisplayName', '900','LineWidth',2)
plot(SA_vec,MZ0_DFz4, 'DisplayName', '1120','LineWidth',2)
legend()
xlabel('$\alpha$ [-]')
ylabel('$M_{Z0}$ [Nm]')
title('Fitted Aligning Moment Variable FZ')

% Conclusions
res_DFz = resid_pure_Mz(fitting_params, ...
                        P_DFz, ...
                        MZ, ...
                        KAPPA, ...
                        ALPHA, ...
                        GAMMA, ...
                        FZ, ...
                        tyre_coeffs);
r_squared = 1-res_DFz;
disp(['Lower boundary:     ', num2str(lb   , "%.3f     ")])
disp(['Upper boundary:     ', num2str(ub   , "%.3f      ")])
disp(['Optimal solution:   ', num2str(P_DFz, "%.3f      ")])
fprintf('Self-Aligning Moment R-squared FZ_var = %6.3f\n', r_squared);

%% Fitting with Variable $\gamma$
fitting_params = {'qBz4','qBz5','qDz3','qDz4','qDz8','qDz9','qEz5','qHz3','qHz4'};

% Tables Intersection
[TdataGamma, ~] = intersect_table_data(Tsubs_LAT.KAPPA.SL_0, ...
                                       Tsubs_LAT.FZ.FZ_700);

% Initialization
KAPPA           = TdataGamma.SL;
ALPHA           = TdataGamma.SA;
GAMMA           = TdataGamma.IA;
FZ              = TdataGamma.FZ;
FY              = TdataGamma.FY;
MZ              = -TdataGamma.MZ;
FZ0             = mean(Tsubs_LAT.FZ.FZ_700.FZ);
tyre_coeffs.FZ0 = 700;
SA_vec          = -0.3:0.001:0.3;
% zeros_ALPHA     = zeros(size(ALPHA));
ones_ALPHA      = ones(size(KAPPA));
zeros_SA        = zeros(size(SA_vec));
ones_SA         = ones(size(SA_vec));

% Initial Guess
%    [qBz4  qBz5  qDz3  qDz4  qDz8  qDz9  qEz5  qHz3   qHz4] 
% P0 = [-0.1, 0.01,   0.05,-0.005, 0.1,  -0.1,   0.3]; % R=0,926
% P0 = [-1.4,    1,   -1,   -1,  1.5,   -1,    1,  0.5,  0.07];
P0 = [-1.4,    0.9,   -2,   -1.8,  1.5,   -1,    1,  0.65,  0.07];
lb = [  -3,    0,   -2.5,   -2.5,    0,  -5,   -1,   -1,    -1];
ub = [   3,    2,      2,      2,    2,   5,    2,    1,     1];
% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised

% First Guess
for i = 1:length(fitting_params)
    tyre_coeffs.(fitting_params{i}) = P0(i);
end 
MZ0_guess = MF96_MZ0_vec(KAPPA, ...
                         ALPHA, ...
                         GAMMA, ...
                         tyre_coeffs.FZ0*ones_ALPHA, ...
                         tyre_coeffs);

% Optimization
[P_varGamma, ~] = fmincon(@(P)resid_pure_Mz(fitting_params, ...
                                            P, ...
                                            MZ, ...
                                            KAPPA, ...
                                            ALPHA, ...
                                            GAMMA, ...
                                            FZ0*ones_ALPHA, ...
                                            tyre_coeffs), ...
                          P0, [], [], [], [], lb, ub);

% Optimal Parameters Assignment
for i = 1:length(fitting_params)
    tyre_coeffs.(fitting_params{i}) = P_varGamma(i);
end

% Force and Residuals Computation with Optimal Parameters
MZ0_varGamma = MF96_MZ0_vec(zeros_SA, ...
                            SA_vec, ...
                            zeros_SA, ...
                            FZ0.*ones_SA, ...
                            tyre_coeffs);

% Plotting
figure('Name', 'First Guess Aligning Moment Variable Gamma')
plot(ALPHA,MZ,'-')
hold on
plot(ALPHA,MZ0_guess,'-', 'LineWidth',2)
title('First Guess Aligning Moment Variable Gamma')

% -------------------------------------------------------------------------

figure('Name', 'Fitted Aligning Moment Variable Gamma')
tiledlayout(3,2)
hold on

xlabel('$\alpha$ [-]')
ylabel('$M_{z0}$ [N]')

GAMMA_0 = Tsubs_LAT.GAMMA.GAMMA_0;
GAMMA_1 = Tsubs_LAT.GAMMA.GAMMA_1;
GAMMA_2 = Tsubs_LAT.GAMMA.GAMMA_2;
GAMMA_3 = Tsubs_LAT.GAMMA.GAMMA_3;
GAMMA_4 = Tsubs_LAT.GAMMA.GAMMA_4;

camber_TXT    =  [       0,       1,       2,       3,       4];
camber_ANGLES =  { GAMMA_0, GAMMA_1, GAMMA_2, GAMMA_3, GAMMA_4};

for i = 1:length(camber_ANGLES)
    [TDataSub, ~] = intersect_table_data(TdataGamma, ...
                                         camber_ANGLES{i});
    KAPPA = TDataSub.SL;
    ALPHA = TDataSub.SA;
    GAMMA = TDataSub.IA;
    FZ = TDataSub.FZ;
    FY = TDataSub.FY;   
    
    [fy0_fit] = MF96_FY0_vec(KAPPA, ...
                             ALPHA, ...
                             GAMMA, ...
                             FZ, ...
                             tyre_coeffs);

    nexttile
    hold on
    grid on

    plot(ALPHA, FY, 'o',...
        'DisplayName', strcat('Raw $\gamma$ = ',num2str(camber_TXT(i)), ' deg'));
    plot(ALPHA, fy0_fit, '-', ...
        'DisplayName', strcat('Fitted $\gamma$ = ',num2str(camber_TXT(i)), ' deg'), ...
        'LineWidth',2);
    
    legend('location','southeast')
end

sgtitle('Fitted Aligning Moment Variable Gamma')

% Conclusions
res_varGamma = resid_pure_Mz(fitting_params, ...
                             P_varGamma, ...
                             MZ, ...
                             KAPPA, ...
                             ALPHA, ...
                             GAMMA, ...
                             FZ, ...
                             tyre_coeffs);
r_squared = 1-res_varGamma;
disp(['Lower boundary:     ', num2str(lb        , "%.3f     ")])
disp(['Upper boundary:     ', num2str(ub        , "%.3f      ")])
disp(['Optimal solution:   ', num2str(P_varGamma, "%.3f      ")])
fprintf('Self-Aligning Moment R-squared Gamma_var = %6.3f\n', r_squared);

%% Save tyre data structure to mat file

save(['tyre_' dataset_LAT,'.mat'],'tyre_coeffs');