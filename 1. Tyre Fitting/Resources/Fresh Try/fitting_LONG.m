%% Fitting with Fz = Fz$_220$ and $\gamma$ = 0, $\alpha$ = 0
fitting_params = {'pCx1';'pDx1';'pEx1';'pEx4';'pHx1';'pKx1';'pVx1'};

% Tables Intersection
[Tdata0, ~] = intersect_table_data(Tsubs_LONG.ALPHA.SA_0, ...
                                   Tsubs_LONG.GAMMA.GAMMA_0, ...
                                   Tsubs_LONG.FZ.FZ_220);

plot_selected_data(Tdata0, 'Longitudinal Selected Data')

% Initialization
KAPPA           = Tdata0.SL;
ALPHA           = Tdata0.SA;
GAMMA           = Tdata0.IA;
FZ              = Tdata0.FZ;
FX              = Tdata0.FX;
FZ0             = mean(FZ);
tyre_coeffs.FZ0 = 220;
SL_vec          = -0.3:0.001:0.3;
zeros_KAPPA     = zeros(size(KAPPA));
ones_KAPPA      = ones(size(KAPPA));
zeros_SL        = zeros(size(SL_vec));
ones_SL         = ones(size(SL_vec));

% Initial Guess
%    [pCx1 pDx1 pEx1 pEx4 pHx1 pKx1 pVx1] 
P0 = [   1,   2,   1,   0,   0,  75,   0];
lb = [   1, 0.1,   0,   0, -10,   0, -10];
ub = [   2,   4,   1,   1,  10, 100,  10];
% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised

% Fitting Settings



% First Guess
for i = 1:length(fitting_params)
    tyre_coeffs.(fitting_params{i}) = P0(i);
end 
FX0_guess = MF96_FX0_vec(KAPPA, ...
                         zeros_KAPPA, ...
                         zeros_KAPPA, ...
                         tyre_coeffs.FZ0*ones_KAPPA, ...
                         tyre_coeffs);

% Optimization
[P_fz_nom, ~] = fmincon(@(P)resid_pure_Fx(fitting_params, ...
                                          P, ...
                                          FX, ...
                                          KAPPA, ...
                                          ALPHA, ...
                                          GAMMA, ...
                                          FZ0*ones_KAPPA, ...
                                          tyre_coeffs), ...
                        P0, [], [], [], [], lb, ub);

% Optimal Parameters Assignment
for i = 1:length(fitting_params)
    tyre_coeffs.(fitting_params{i}) = P_fz_nom(i);
end

% Force and Residuals Computation with Optimal Parameters
FX0_fz_nom = MF96_FX0_vec(SL_vec, ...
                          zeros_SL, ...
                          zeros_SL, ...
                          FZ0.*ones_SL, ...
                          tyre_coeffs);

% Plotting
figure('Name', 'First Guess Longitudinal')
plot(KAPPA,FX,'.')
hold on
plot(KAPPA,FX0_guess,'-')

% -------------------------------------------------------------------------

figure('Name','Fx0(Fz0)')
plot(KAPPA,FX,'o')
hold on
plot(SL_vec,FX0_fz_nom,'-','LineWidth',2)
xlabel('$\kappa$ [-]')
ylabel('$F_{x0}$ [N]')

% Conclusions
res_fz_nom = resid_pure_Fx(fitting_params, ...
                           P_fz_nom, ...
                           FX, ...
                           KAPPA, ...
                           ALPHA, ...
                           GAMMA, ...
                           FZ0*ones_KAPPA, ...
                           tyre_coeffs);
r_squared = 1-res_fz_nom;
disp(['Lower boundary:     ', num2str(lb      , "%.3f     ")])
disp(['Upper boundary:     ', num2str(ub      , "%.3f     ")])
disp(['Optimal solution:   ', num2str(P_fz_nom, "%.3f      ")])
fprintf('R-squared = %6.3f\n', r_squared);

%% Fitting with Variable Fz and $\gamma$ = 0, $\alpha$ = 0
fitting_params = {'pDx2';'pEx2';'pEx3';'pHx2';'pKx2';'pKx3';'pVx2'};

% Tables Intersection
[TdataDFz, ~] = intersect_table_data(Tsubs_LONG.ALPHA.SA_0, ...
                                   Tsubs_LONG.GAMMA.GAMMA_0);

% Initialization
KAPPA           = TdataDFz.SL;
ALPHA           = TdataDFz.SA;
GAMMA           = TdataDFz.IA;
FZ              = TdataDFz.FZ;
FX              = TdataDFz.FX;
FZ0             = mean(Tsubs_LONG.FZ.FZ_220.FZ);
FZ1             = mean(Tsubs_LONG.FZ.FZ_700.FZ);
FZ2             = mean(Tsubs_LONG.FZ.FZ_900.FZ);
FZ3             = mean(Tsubs_LONG.FZ.FZ_1120.FZ);
tyre_coeffs.FZ0 = 220;
SL_vec          = -0.3:0.001:0.3;
zeros_KAPPA     = zeros(size(KAPPA));
ones_KAPPA      = ones(size(KAPPA));
zeros_SL        = zeros(size(SL_vec));
ones_SL         = ones(size(SL_vec));

% Initial Guess
%    [pDx1 pEx2 pEx3 pHx2  pKx2 pKx3 pVx2] 
P0 = [-0.1,   0,   0,   0,    0,   0,   0];
lb = [  -1,  -1,  -1,  -1, -100, -10,  -1];
ub = [   0,   1,   1,   1,  100,  10,   1];
% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised

% Fitting Settings



% First Guess
for i = 1:length(fitting_params)
    tyre_coeffs.(fitting_params{i}) = P0(i);
end 
FX0_guess = MF96_FX0_vec(KAPPA, ...
                         zeros_KAPPA, ...
                         zeros_KAPPA, ...
                         tyre_coeffs.FZ0*ones_KAPPA, ...
                         tyre_coeffs);

% Optimization
[P_DFz, ~] = fmincon(@(P)resid_pure_Fx(fitting_params, ...
                                       P, ...
                                       FX, ...
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
FX0_DFz0 = MF96_FX0_vec(SL_vec, ...
                        zeros_SL, ...
                        zeros_SL, ...
                        FZ0.*ones_SL, ...
                        tyre_coeffs);
FX0_DFz1 = MF96_FX0_vec(SL_vec, ...
                        zeros_SL, ...
                        zeros_SL, ...
                        FZ1.*ones_SL, ...
                        tyre_coeffs);
FX0_DFz2 = MF96_FX0_vec(SL_vec, ...
                        zeros_SL, ...
                        zeros_SL, ...
                        FZ2.*ones_SL, ...
                        tyre_coeffs);
FX0_DFz3 = MF96_FX0_vec(SL_vec, ...
                        zeros_SL, ...
                        zeros_SL, ...
                        FZ3.*ones_SL, ...
                        tyre_coeffs);

% Cornering Stiffness
[kappa__x, Bx, Cx, Dx, Ex, SVx] =MF96_FX0_coeffs(0, 0, 0, FZ0, tyre_coeffs);
Calfa_vec0_0 = magic_formula_stiffness(kappa__x, Bx, Cx, Dx, Ex, SVx);
[kappa__x, Bx, Cx, Dx, Ex, SVx] =MF96_FX0_coeffs(0, 0, 0, FZ1, tyre_coeffs);
Calfa_vec1_0 = magic_formula_stiffness(kappa__x, Bx, Cx, Dx, Ex, SVx);
[kappa__x, Bx, Cx, Dx, Ex, SVx] =MF96_FX0_coeffs(0, 0, 0, FZ2, tyre_coeffs);
Calfa_vec2_0 = magic_formula_stiffness(kappa__x, Bx, Cx, Dx, Ex, SVx);
[kappa__x, Bx, Cx, Dx, Ex, SVx] =MF96_FX0_coeffs(0, 0, 0, FZ3, tyre_coeffs);
Calfa_vec3_0 = magic_formula_stiffness(kappa__x, Bx, Cx, Dx, Ex, SVx);

Calfa_vec0 = MF96_CorneringStiffness_LONG(SL_vec, zeros_SL ,zeros_SL, FZ0*ones_SL,tyre_coeffs);
Calfa_vec1 = MF96_CorneringStiffness_LONG(SL_vec, zeros_SL ,zeros_SL, FZ1*ones_SL,tyre_coeffs);
Calfa_vec2 = MF96_CorneringStiffness_LONG(SL_vec, zeros_SL ,zeros_SL, FZ2*ones_SL,tyre_coeffs);
Calfa_vec3 = MF96_CorneringStiffness_LONG(SL_vec, zeros_SL ,zeros_SL, FZ3*ones_SL,tyre_coeffs);


% Plotting
figure('Name', 'Variable FZ Raws')
plot(KAPPA,FX,'.')
hold on
plot(KAPPA,FX0_guess,'-')

% -------------------------------------------------------------------------

figure('Name','Fx0(Fz0)')
plot(KAPPA,FX,'o')
hold on
plot(SL_vec,FX0_DFz0,'-','LineWidth',2)
plot(SL_vec,FX0_DFz1,'-','LineWidth',2)
plot(SL_vec,FX0_DFz2,'-','LineWidth',2)
plot(SL_vec,FX0_DFz3,'-','LineWidth',2)

xlabel('$\kappa$ [-]')
ylabel('$F_{x0}$ [N]')

% -------------------------------------------------------------------------

figure('Name','C_alpha')
subplot(2,1,1)
hold on
plot(FZ0,Calfa_vec0_0,'+','LineWidth',2)
plot(FZ1,Calfa_vec1_0,'+','LineWidth',2)
plot(FZ2,Calfa_vec2_0,'+','LineWidth',2)
plot(FZ3,Calfa_vec3_0,'+','LineWidth',2)
legend({'$Fz_{220}$','$Fz_{700}$','$Fz_{900}$','$Fz_{1120}$'})

subplot(2,1,2)
hold on
plot(SL_vec,Calfa_vec0,'-','LineWidth',2)
plot(SL_vec,Calfa_vec1,'-','LineWidth',2)
plot(SL_vec,Calfa_vec2,'-','LineWidth',2)
plot(SL_vec,Calfa_vec3,'-','LineWidth',2)
legend({'$Fz_{220}$','$Fz_{700}$','$Fz_{900}$','$Fz_{1120}$'})

% Conclusions
res_DFz = resid_pure_Fx(fitting_params, ...
                        P_DFz, ...
                        FX, ...
                        KAPPA, ...
                        ALPHA, ...
                        GAMMA, ...
                        FZ, ...
                        tyre_coeffs);
r_squared = 1-res_DFz;
disp(['Lower boundary:     ', num2str(lb   , "%.3f     ")])
disp(['Upper boundary:     ', num2str(ub   , "%.3f      ")])
disp(['Optimal solution:   ', num2str(P_DFz, "%.3f      ")])
fprintf('R-squared = %6.3f\n', r_squared);

%% Fitting with Variable $\gamma$
fitting_params = {'pDx3'};

% Tables Intersection
[TdataGamma, ~] = intersect_table_data(Tsubs_LONG.ALPHA.SA_0, ...
                                       Tsubs_LONG.FZ.FZ_220);

% Initialization
KAPPA           = TdataGamma.SL;
ALPHA           = TdataGamma.SA;
GAMMA           = TdataGamma.IA;
FZ              = TdataGamma.FZ;
FX              = TdataGamma.FX;
FZ0             = mean(Tsubs_LONG.FZ.FZ_220.FZ);
tyre_coeffs.FZ0 = 220;
SL_vec          = -0.3:0.001:0.3;
% zeros_KAPPA     = zeros(size(KAPPA));
ones_KAPPA      = ones(size(KAPPA));
zeros_SL        = zeros(size(SL_vec));
ones_SL         = ones(size(SL_vec));

% Initial Guess
%  [pDx3] 
P0 = 2.5;
lb =  -5;
ub =   5;
% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised

% Fitting Settings



% First Guess
for i = 1:length(fitting_params)
    tyre_coeffs.(fitting_params{i}) = P0(i);
end 
FX0_guess = MF96_FX0_vec(KAPPA, ...
                         ALPHA, ...
                         GAMMA, ...
                         tyre_coeffs.FZ0*ones_KAPPA, ...
                         tyre_coeffs);

% Optimization
[P_varGamma, ~] = fmincon(@(P)resid_pure_Fx(fitting_params, ...
                                            P, ...
                                            FX, ...
                                            KAPPA, ...
                                            ALPHA, ...
                                            GAMMA, ...
                                            FZ0*ones_KAPPA, ...
                                            tyre_coeffs), ...
                          P0, [], [], [], [], lb, ub);

% Optimal Parameters Assignment
for i = 1:length(fitting_params)
    tyre_coeffs.(fitting_params{i}) = P_varGamma(i);
end

% Force and Residuals Computation with Optimal Parameters
FX0_varGamma = MF96_FX0_vec(SL_vec, ...
                        zeros_SL, ...
                        zeros_SL, ...
                        FZ0.*ones_SL, ...
                        tyre_coeffs);

% Plotting
figure('Name', 'Variable Gamma Raws')
plot(KAPPA,FX,'-')
hold on
plot(KAPPA,FX0_guess,'-')

% -------------------------------------------------------------------------

figure('Name', 'Fx0 vs Gamma')
plot(KAPPA,FX,'-')
hold on
plot(SL_vec,FX0_varGamma,'-','LineWidth',2)

xlabel('$\kappa$ [-]')
ylabel('$F_{x0}$ [N]')

% Conclusions
res_varGamma = resid_pure_Fx(fitting_params, ...
                             P_varGamma, ...
                             FX, ...
                             KAPPA, ...
                             ALPHA, ...
                             GAMMA, ...
                             FZ, ...
                             tyre_coeffs);
r_squared = 1-res_varGamma;
disp(['Lower boundary:     ', num2str(lb        , "%.3f     ")])
disp(['Upper boundary:     ', num2str(ub        , "%.3f      ")])
disp(['Optimal solution:   ', num2str(P_varGamma, "%.3f      ")])
fprintf('R-squared = %6.3f\n', r_squared);

%% Fitting Combined Case
fitting_params = {'rBx1'; 'rBx2'; 'rCx1'; 'rHx1'};

% Tables Intersection
TdataFX = intersect_table_data(Tsubs_LONG.GAMMA.GAMMA_0, Tsubs_LONG.FZ.FZ_220);
% TdataFX = tyre_data_LONG;

% Initialization
KAPPA           = TdataFX.SL;
ALPHA           = TdataFX.SA;
GAMMA           = TdataFX.IA;
FZ              = TdataFX.FZ;
FX              = TdataFX.FX;
FZ0             = mean(Tsubs_LONG.FZ.FZ_220.FZ);
tyre_coeffs.FZ0 = 220;
SL_vec          = -0.3:0.001:0.3;
zeros_KAPPA     = zeros(size(KAPPA));
ones_KAPPA      = ones(size(KAPPA));
zeros_SL        = zeros(size(SL_vec));
ones_SL         = ones(size(SL_vec));

% Initial Guess
%    [rBx1  rBx2  rCx1  rHx1]
P0 = [   3,    1,    1,  -10];
% lb = [ -20,  -20,  -20,  -20];
% ub = [  20,   20,   20,   20];
lb = [   0,    0,    0,  -15];
ub = [   5,    2,    2,   -5];
% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised

% Fitting Settings



% First Guess
for i = 1:length(fitting_params)
    tyre_coeffs.(fitting_params{i}) = P0(i);
end 
FX_guess = MF96_FX_vec(KAPPA, ...
                       ALPHA, ...
                       GAMMA, ...
                       FZ0*ones_KAPPA, ...
                       tyre_coeffs);

% Optimization
[P_comb, ~] = fmincon(@(P)resid_Fx_Comb(fitting_params, ...
                                        P, ...
                                        FX, ...
                                        KAPPA, ...
                                        ALPHA, ...
                                        GAMMA, ...
                                        FZ, ...
                                        tyre_coeffs), ...
                      P0, [], [], [], [], lb, ub);

% Optimal Parameters Assignment
for i = 1:length(fitting_params)
    tyre_coeffs.(fitting_params{i}) = P_comb(i);
end

% Force and Residuals Computation with Optimal Parameters
FX_comb = MF96_FX_vec(KAPPA, ...
                      ALPHA, ...
                      GAMMA, ...
                      FZ, ...
                      tyre_coeffs);

% Plotting
figure('Name', 'Combined Raws')
plot(KAPPA,FX,'-')
hold on
plot(KAPPA,FX_guess,'-', LineWidth=5)
plot(KAPPA,FX_comb, '-')

% -------------------------------------------------------------------------

figure('Name', 'Combined longitudinal slip at variable slip angles')
hold on

ALPHA_LEVELS  = {'SA_0', 'SA_3neg', 'SA_6neg'};
ALPHA_LVL_VAL = [    0 ,       -3 ,       -6 ];
for i = 1:length(ALPHA_LEVELS)
    [TdataComb, ~] = intersect_table_data(TdataFX, ...
                                          Tsubs_LONG.ALPHA.(ALPHA_LEVELS{i}));
    KAPPA  = TdataComb.SL;
    ALPHA  = TdataComb.SA;
    GAMMA  = TdataComb.IA;
    FZ     = TdataComb.FZ;
    FX     = TdataComb.FX;
    SL_vec = -0.3:(1/length(KAPPA)):0.3;

    [fx_comb] = MF96_FX_vec(SL_vec, ...
                            ALPHA_LVL_VAL(i)*ones_SL*to_rad, ...
                            zeros_SL, ...
                            tyre_coeffs.FZ0*ones_SL, ...
                            tyre_coeffs);

    plot(KAPPA, FX, '+', ...
         'DisplayName', strcat('Raw $\alpha$ = ', num2str(ALPHA_LVL_VAL(i)), ...
         ' deg'));
    plot(SL_vec, fx_comb, ...
         'DisplayName', strcat('Fitted $\alpha$ = ', num2str(ALPHA_LVL_VAL(i)), ...
         ' deg'));
end

title(strcat({'Combined longitudinal slip at variable slip angles', ...
              '$F_{z} = F_{z0}, \gamma = 0$'}));
xlabel('longitudinal slip $\kappa$ [-]')
ylabel('$F_x$ [N]')
legend('location', 'southeast')
grid on

% -------------------------------------------------------------------------

kappa_vals = [0, 0.05, 0.1, 0.15, 0.2];
alpha_vals = linspace(min(tyre_data_LONG.SA), max(tyre_data_LONG.SA), 100);

figure("Name",'Ckappa')
hold on

for i = 1:length(kappa_vals)
    gxa_fit = [];
    for k = 1:length(alpha_vals)
        [Gxa, Gyk, SVyk] = MF96_FXcomb_coeffs(kappa_vals(i), ...
                                              alpha_vals(k), ...
                                              0, ...
                                              FZ0, ...
                                              tyre_coeffs);
        gxa_fit = [gxa_fit, Gxa];
    end
    plot(alpha_vals*to_deg, gxa_fit, ...
         'DisplayName', strcat('$\kappa = $ ', num2str(kappa_vals(i))));
end

title(strcat({'Weighting function G$_{xa}$ as a function of $\alpha$', ...
                     ' $F_z = F_{z0} , \gamma = 0$ '}));
xlabel('side slip angle $\alpha$ (deg)')
ylabel('G$_{xa}$')
legend('location', 'southeast')
grid on

% -------------------------------------------------------------------------

alpha_vals = [-6, -4.5, -3, -1.5, 0];
kappa_vals = linspace(min(tyre_data_LONG.SL), -min(tyre_data_LONG.SL), 100);

figure("Name",'Calpha')
hold on

for i = 1:length(alpha_vals)
    gxa_fit = [];
    for k = 1:length(kappa_vals)
        [Gxa, Gyk, SVyk] = MF96_FXcomb_coeffs(kappa_vals(k), ...
                                              alpha_vals(i), ...
                                              0, ...
                                              tyre_coeffs.FZ0, ...
                                              tyre_coeffs);
        gxa_fit = [gxa_fit, Gxa];
    end
    plot(kappa_vals*to_deg, gxa_fit, ...
         'DisplayName', strcat('$\alpha = $ ', num2str(alpha_vals(i))));
end

title(strcat({'Weighting function G$_{xa}$ as a function of $\kappa$', ...
                     ' $F_z = F_{z0} , \gamma = 0$ '}));
xlabel('slip angle $\kappa$ (deg)')
ylabel('G$_{xa}$')
legend('location', 'southeast')
grid on

% Conclusions
res_FX = resid_Fx_Comb(fitting_params, ...
                       P_comb, ...
                       FX, ...
                       KAPPA, ...
                       ALPHA, ...
                       GAMMA, ...
                       FZ, ...
                       tyre_coeffs);
r_squared = 1-res_FX;
disp(['Lower boundary:     ', num2str(lb    , "%.3f     ")])
disp(['Upper boundary:     ', num2str(ub    , "%.3f      ")])
disp(['Optimal solution:   ', num2str(P_comb, "%.3f      ")])
fprintf('R-squared = %6.3f\n', r_squared);

%% Save tyre data structure to mat file
%
% save(['tyre_' data_set,'.mat'],'tyre_coeffs');