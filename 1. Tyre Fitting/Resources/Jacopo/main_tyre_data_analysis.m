%
% /////////////////////////////////////////////////////////////////////// %
% ////////////////// DYNAMICS OF VEHICLES - TYRE FITTING //////////////// %
% /////////////////////////////////////////////////////////////////////// %



%% 1 LONGITUDINAL CASE
% /////////////////////////////////////////////////////////////////////// %


%% 1.1 Initialization
clc
clearvars 
close all   

% Set LaTeX as default interpreter for axis labels, ticks and legends
set(0,'defaulttextinterpreter','latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

set(0,'DefaultFigureWindowStyle','docked');
set(0,'defaultAxesFontSize',  16)
set(0,'DefaultLegendFontSize',16)

addpath('tyre_lib/')

to_rad = pi/180;
to_deg = 180/pi;


%% 1.2 Select longitudinal tyre dataset
% Dataset path
data_set_path = 'dataset/';

% Dataset selection and loading
data_set = 'B1464run30';  % braking/traction (pure log. force) + combined

% Tyre geometric data:
% Hoosier	18.0x6.0-10
% 18 diameter in inches
% 6.0 section width in inches
% tread width in inches
diameter = 18*2.56; % Diameter of the tyre, [mm] (x2.56 is the converting
                    % factor from inches to millimeters)
Fz0 = 220;   % Nominal load is given, [N]
R0  = diameter/2/100; % Radius of the tyre, [m] (/100 is to get the value
                        % expressed in meters)

fprintf('Loading dataset ...')
switch data_set
  case 'B1464run23'
  load ([data_set_path, 'Hoosier_B1464run23.mat']); % pure lateral
  cut_start = 27760;
  cut_end   = 54500;
  case 'B1464run30'
  load ([data_set_path, 'Hoosier_B1464run30.mat']); % pure longitudinal
  cut_start = 19028;
  cut_end   = 37643;
  case 'B1464run58'
  load ([data_set_path, 'Goodyear_B1464run58.mat']); % pure longitudinal
  cut_start = 19028;
  cut_end   = 37643;
  case 'B1464run13'
  load ([data_set_path, 'Goodyear_B1464run13.mat']); % pure lateral
  cut_start = 27760;
  cut_end   = 54500;
  case 'B1464run51'
  load ([data_set_path, 'Continental_B1464run51.mat']); % pure longitudinal
  cut_start = 19028;
  cut_end   = 54500;
  case 'B1464run8'
  load ([data_set_path, 'Continental_B1464run30.mat']); % pure lateral
  cut_start = 27760;
  cut_end   = 37643;
  otherwise 
  error('Not found dataset: `%s`\n', data_set) ;
end

% select dataset portion
smpl_range = cut_start:cut_end;

fprintf('completed!\n')


%% 1.3 Plot raw data
figure
tiledlayout(6,1)

ax_list(1) = nexttile; y_range = [min(min(-FZ),0) round(max(-FZ)*1.1)];
plot(-FZ)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Vertical force')
xlabel('Samples [-]')
ylabel('[N]')

ax_list(2) = nexttile; y_range = [min(min(IA),0) round(max(IA)*1.1)];
plot(IA)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Camber angle')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list(3) = nexttile; y_range = [min(min(SA),0) round(max(SA)*1.1)];
plot(SA)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Side slip')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list(4) = nexttile; y_range = [min(min(SL),0) round(max(SL)*1.1)];
plot(SL)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Longitudinal slip')
xlabel('Samples [-]')
ylabel('[-]')

ax_list(5) = nexttile; y_range = [min(min(P),0) round(max(P)*1.1)];
plot(P)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Tyre pressure')
xlabel('Samples [-]')
ylabel('[psi]')

ax_list(6) = nexttile;  y_range = [min(min(TSTC),0) round(max(TSTC)*1.1)];
plot(TSTC,'DisplayName','Center')
hold on
plot(TSTI,'DisplayName','Internal')
plot(TSTO,'DisplayName','Outboard')
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Tyre temperatures')
xlabel('Samples [-]')
ylabel('[degC]')

linkaxes(ax_list,'x')


%% 1.4 Select some specific data
% Cut crappy data and select only 12 psi data (so just one set of pressure)

vec_samples = 1:1:length(smpl_range);
tyre_data = table(); % create empty table

% store raw data in table
tyre_data.SL =  SL(smpl_range);
tyre_data.SA =  SA(smpl_range)*to_rad;
tyre_data.FZ =  -FZ(smpl_range);    % "-" to get the correct standard
tyre_data.FX =  FX(smpl_range);
tyre_data.FY =  -FY(smpl_range);    % "-" to get the correct standard        
tyre_data.MZ =  -MZ(smpl_range);    % "-" to get the correct standard
tyre_data.IA =  IA(smpl_range)*to_rad;

% Extract points at constant inclination angle
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

% Extract points at constant vertical load
% Test data done at: 
%  - 50lbf  ( 50*0.453592*9.81 =  223N )
%  - 150lbf (150*0.453592*9.81 =  667N )
%  - 200lbf (200*0.453592*9.81 =  890N )
%  - 250lbf (250*0.453592*9.81 = 1120N )

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

% The slip angle is varied continuously between -4 and +12° and then
% between -12° and +4° for the pure slip case

% The slip angle is varied step wise for longitudinal slip tests
% 0° , - 3° , -6 °
SA_tol = 0.5*to_rad;
idx.SA_0    =  0-SA_tol          < tyre_data.SA & tyre_data.SA < 0+SA_tol;
idx.SA_3neg = -(3*to_rad+SA_tol) < tyre_data.SA & tyre_data.SA < -3*to_rad+SA_tol;
idx.SA_6neg = -(6*to_rad+SA_tol) < tyre_data.SA & tyre_data.SA < -6*to_rad+SA_tol;
SA_0     = tyre_data( idx.SA_0, : );
SA_3neg  = tyre_data( idx.SA_3neg, : );
SA_6neg  = tyre_data( idx.SA_6neg, : );

figure()
tiledlayout(3,1)

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
xlabel('Samples [-]')
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
xlabel('Samples [-]')
ylabel('[N]')

ax_list(3) = nexttile;
plot(tyre_data.SA*to_deg)
hold on
plot(vec_samples(idx.SA_0),   SA_0.SA*to_deg,'.');
plot(vec_samples(idx.SA_3neg),SA_3neg.SA*to_deg,'.');
plot(vec_samples(idx.SA_6neg),SA_6neg.SA*to_deg,'.');
title('Slide slip')
xlabel('Samples [-]')
ylabel('[rad]')


%% 1.5 Intersecting tables to obtain specific sub-datasets
[TData0, ~] = intersect_table_data(SA_0, GAMMA_0, FZ_220); 


%% 1.6 Plotting selected data for longitudinal dataset
figure('Name','Selected-data')
plot_selected_data(TData0);


%% 1.7 Fitting 
% initialise tyre data
tyre_coeffs = initialise_tyre_data(R0, Fz0);


%% 1.8 Longitudinal slip
% Fitting with Fz = Fz_nom = 220N and camber = 0, alpha = 0, VX = 10

FZ0 = mean(TData0.FZ);
zeros_vec = zeros(size(TData0.SL));
ones_vec  = ones(size(TData0.SL));
FX0_guess = MF96_FX0_vec(TData0.SL,zeros_vec,zeros_vec,tyre_coeffs.FZ0*ones_vec,tyre_coeffs);

% Plotting check guess 
figure()
plot(TData0.SL,TData0.FX,'.')
hold on
grid on
plot(TData0.SL,FX0_guess,'-')
title('Plot raw data and initial guess')
xlabel('$\kappa$ [-]')
ylabel('$F_{x}$ [N]')
legend({'Raw data','Fit guess'})

% Guess values for parameters to be optimised
% [pCx1 pDx1 pEx1 pEx4  pHx1  pKx1  pVx1 
P0 = [ 1, 2, 1, 0, 0, 1, 0]; 

% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised
% 1< pCx1 < 2 
% 0< pEx1 < 1 
% Lower bounds
lb = [1, 0.1, 0, 0, -10, 0, -10];
% Upper bounds
ub = [2, 4, 1, 1, 10, 100, 10];

KAPPA_vec = TData0.SL;
FX_vec    = TData0.FX;
SL_vec = -0.3:0.001:0.3;

% LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 
[P_fz_nom,fval,exitflag] = fmincon(@(P)resid_pure_Fx(P,FX_vec,KAPPA_vec,0,FZ0,tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

% Update tyre data with new optimal values                             
tyre_coeffs.pCx1 = P_fz_nom(1);
tyre_coeffs.pDx1 = P_fz_nom(2);  
tyre_coeffs.pEx1 = P_fz_nom(3);
tyre_coeffs.pEx4 = P_fz_nom(4);
tyre_coeffs.pHx1 = P_fz_nom(5); 
tyre_coeffs.pKx1 = P_fz_nom(6);
tyre_coeffs.pVx1 = P_fz_nom(7);

% Residual from optimal parameters
res_FX0_nom_vec = resid_pure_Fx(P_fz_nom,FX_vec,KAPPA_vec,0,FZ0,tyre_coeffs);

% Vector for plot
FX0_fz_nom_vec = MF96_FX0_vec(SL_vec,zeros(size(SL_vec)),zeros(size(SL_vec)), ...
                              FZ0.*ones(size(SL_vec)),tyre_coeffs);

% Plotting longitudinal slip fit for fixed Fz, no camber, no side slip
figure('Name','Fx0(Fz0)')
plot(TData0.SL,TData0.FX,'o')
hold on
grid on
plot(SL_vec,FX0_fz_nom_vec,'-','LineWidth',2)
xlabel('$\kappa$ [-]')
ylabel('$F_{x0}$ [N]')
title('Pure Long. Slip fit at vertical load $F_z = 220$ N, zero camber, zero side slip, Press. = 12 psi.')

% R-squared is 
% 1-SSE/SST
% SSE/SST = res_Fx0_nom_vec
% Printing the R-squared
fprintf('Printing R-squared for Pure Long. slip, no camber, no side slip :\n');
fprintf('R-squared = %6.3f\n',1-res_FX0_nom_vec);


%% 1.9 Longitudinal slip with variable vertical load
% Fit coefficient with variable load
% extract data with variable load
[TDataDFz, ~] = intersect_table_data(SA_0, GAMMA_0);

zeros_vec = zeros(size(TDataDFz.SL));
ones_vec  = ones(size(TDataDFz.SL));

% Guess fit
FX0_guess = MF96_FX0_vec(TDataDFz.SL,zeros_vec,zeros_vec,tyre_coeffs.FZ0*ones_vec,tyre_coeffs);

% Check guess 
figure()
plot(TDataDFz.SL,TDataDFz.FX,'.')
hold on
grid on
plot(TDataDFz.SL,FX0_guess,'-')
title('Plot raw data and initial guess')
xlabel('$\kappa$ [-]')
ylabel('$F_{x}$ [N]')
legend({'Raw Data','Initial guess fit'})

% Guess values for parameters to be optimised
% Parameters to optimise: [pDx2 pEx2 pEx3 pHx2  pKx2  pKx3  pVx2] 
P0 = [ 0, 0, 0, 0, 0, 0, 0]; 

% Limits for parameters to be optimised
% Lower bounds
lb = [-0.3,-0.37,-0.1,-0.1,-0.1,-0.1,-0.1];
% Upper bounds
ub = [0.1,0.1,0.11,0.1,0.1,0.16,0.1];

KAPPA_vec = TDataDFz.SL;
FX_vec    = TDataDFz.FX;
FZ_vec    = TDataDFz.FZ;
SL_vec = -0.3:0.001:0.3;

% Optimization of the parameters through least square approach 
[P_dfz,~,~] = fmincon(@(P)resid_pure_Fx_varFz(P,FX_vec,KAPPA_vec,0,FZ_vec,tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

% Change tyre data with new optimal values                             
tyre_coeffs.pDx2 = P_dfz(1); 
tyre_coeffs.pEx2 = P_dfz(2);  
tyre_coeffs.pEx3 = P_dfz(3);
tyre_coeffs.pHx2 = P_dfz(4);
tyre_coeffs.pKx2 = P_dfz(5); 
tyre_coeffs.pKx3 = P_dfz(6);
tyre_coeffs.pVx2 = P_dfz(7);

% Residuals from optimized paramters
res_FX0_dfz_vec = resid_pure_Fx_varFz(P_dfz,FX_vec,KAPPA_vec,0,FZ_vec,tyre_coeffs);

% Vectors for plot
tmp_zeros = zeros(size(SL_vec));
tmp_ones = ones(size(SL_vec));

FX0_fz_var_vec1 = MF96_FX0_vec(SL_vec,tmp_zeros,tmp_zeros,mean(FZ_220.FZ)*tmp_ones,tyre_coeffs);
FX0_fz_var_vec2 = MF96_FX0_vec(SL_vec,tmp_zeros,tmp_zeros,mean(FZ_700.FZ)*tmp_ones,tyre_coeffs);
FX0_fz_var_vec3 = MF96_FX0_vec(SL_vec,tmp_zeros,tmp_zeros,mean(FZ_900.FZ)*tmp_ones,tyre_coeffs);
FX0_fz_var_vec4 = MF96_FX0_vec(SL_vec,tmp_zeros,tmp_zeros,mean(FZ_1120.FZ)*tmp_ones,tyre_coeffs);

% Plotting longitudinal force fit for different vertical loads
figure('Name','Fx0(Fz0)')
plot(TDataDFz.SL,TDataDFz.FX,'o')
hold on
grid on
plot(SL_vec,FX0_fz_var_vec1,'-','LineWidth',2)
plot(SL_vec,FX0_fz_var_vec2,'-','LineWidth',2)
plot(SL_vec,FX0_fz_var_vec3,'-','LineWidth',2)
plot(SL_vec,FX0_fz_var_vec4,'-','LineWidth',2)
title('Pure Long. Slip fit at different vertical load $F_z$.')
xlabel('$\kappa$ [-]')
ylabel('$F_{x0}$ [N]')
legend({'Data','$Fz_{220}$','$Fz_{700}$','$Fz_{900}$','$Fz_{1120}$'})

% Generating arameters and vectors needed for the plot
[kappa__x, Bx, Cx, Dx, Ex, SVx] = MF96_FX0_coeffs(0, 0, 0, mean(FZ_220.FZ), tyre_coeffs);
Calfa_vec1_0 = magic_formula_stiffness(kappa__x, Bx, Cx, Dx, Ex, SVx);
[kappa__x, Bx, Cx, Dx, Ex, SVx] = MF96_FX0_coeffs(0, 0, 0, mean(FZ_700.FZ), tyre_coeffs);
Calfa_vec2_0 = magic_formula_stiffness(kappa__x, Bx, Cx, Dx, Ex, SVx);
[kappa__x, Bx, Cx, Dx, Ex, SVx] = MF96_FX0_coeffs(0, 0, 0, mean(FZ_900.FZ), tyre_coeffs);
Calfa_vec3_0 = magic_formula_stiffness(kappa__x, Bx, Cx, Dx, Ex, SVx);
[kappa__x, Bx, Cx, Dx, Ex, SVx] = MF96_FX0_coeffs(0, 0, 0, mean(FZ_1120.FZ), tyre_coeffs);
Calfa_vec4_0 = magic_formula_stiffness(kappa__x, Bx, Cx, Dx, Ex, SVx);

Calfa_vec1 = MF96_CorneringStiffness(SL_vec,tmp_zeros, tmp_zeros, mean(FZ_220.FZ)*tmp_ones, tyre_coeffs);
Calfa_vec2 = MF96_CorneringStiffness(SL_vec,tmp_zeros, tmp_zeros, mean(FZ_700.FZ)*tmp_ones, tyre_coeffs);
Calfa_vec3 = MF96_CorneringStiffness(SL_vec,tmp_zeros, tmp_zeros, mean(FZ_900.FZ)*tmp_ones, tyre_coeffs);
Calfa_vec4 = MF96_CorneringStiffness(SL_vec,tmp_zeros, tmp_zeros, mean(FZ_1120.FZ)*tmp_ones, tyre_coeffs);

% Plotting tyre stiffness
figure('Name','C_\alpha')
subplot(2,1,1)
hold on
grid on
plot(mean(FZ_220.FZ),Calfa_vec1_0,'+','LineWidth',2)
plot(mean(FZ_700.FZ),Calfa_vec2_0,'+','LineWidth',2)
plot(mean(FZ_900.FZ),Calfa_vec3_0,'+','LineWidth',2)
plot(mean(FZ_1120.FZ),Calfa_vec4_0,'+','LineWidth',2)
legend({'$Fz_{220}$','$Fz_{700}$','$Fz_{900}$','$Fz_{1120}$'})
xlabel('$F_Z$ [N]')
ylabel('$C_{\alpha}$')

subplot(2,1,2)
hold on
grid on
plot(SL_vec,Calfa_vec1,'-','LineWidth',2)
plot(SL_vec,Calfa_vec2,'-','LineWidth',2)
plot(SL_vec,Calfa_vec3,'-','LineWidth',2)
plot(SL_vec,Calfa_vec4,'-','LineWidth',2)
legend({'$Fz_{220}$','$Fz_{700}$','$Fz_{900}$','$Fz_{1120}$'})
xlabel('$\kappa$ [-]')
ylabel('$C_{\alpha}$')

% Printing the R-squared
fprintf('Printing R-squared for Pure Long. slip, variation of Load:\n');
fprintf('R-squared = %6.3f\n',1-res_FX0_dfz_vec);


%% 1.10 Longitudinal slip with variable camber angle
% Fit coefficient with variable camber

% Extract data with variable camber
[TDataGamma, ~] = intersect_table_data(SA_0, FZ_220);

% Guess values for parameters to be optimised
% Parameter: [pDx3]
P0 = [0]; 

% Limits for parameters to be optimised
% Lower bounds
lb = [0];
% Upper bounds
ub = [19];

zeros_vec = zeros(size(TDataGamma.SL));
ones_vec  = ones(size(TDataGamma.SL));
KAPPA_vec = TDataGamma.SL;
GAMMA_vec = TDataGamma.IA; 
FX_vec    = TDataGamma.FX;

figure()
hold on
grid on
plot(KAPPA_vec,FX_vec);
xlabel('$\kappa$ [-]')
ylabel('$F_{x}$ [N]')
title('Raw Data')

% Optimization of the parameters through least square approach 
[P_varGamma,~,~] = fmincon(@(P)resid_pure_Fx_varGamma(P,FX_vec, KAPPA_vec,GAMMA_vec,tyre_coeffs.FZ0, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

% Change tyre data with new optimal values                             
tyre_coeffs.pDx3 = P_varGamma(1);

% Vector for plot
FX0_varGamma_vec = MF96_FX0_vec(KAPPA_vec,zeros_vec , GAMMA_vec, tyre_coeffs.FZ0*ones_vec,tyre_coeffs);

% Plotting the pure longitudinal slip for variable camber
figure('Name','Fx0 vs Gamma')
plot(KAPPA_vec,TDataGamma.FX,'o')
hold on
grid on
plot(KAPPA_vec,FX0_varGamma_vec,'-')
xlabel('$\kappa$ [-]')
ylabel('$F_{x0}$ [N]')
legend({'Raw data','$F_{x0}$'})
title('Pure longitudinal force fit for variable $\gamma$')

% Calculate the residuals with the optimal solution found above
res_Fx0_varGamma  = resid_pure_Fx_varGamma(P_varGamma,FX_vec, KAPPA_vec,GAMMA_vec,tyre_coeffs.FZ0, tyre_coeffs);

% Printing the R-squared
fprintf('Printing R-squared for Pure Long. slip, variation of camber:\n');
fprintf('R-squared = %6.3f\n',1-res_Fx0_varGamma);

% Generating paramters of interest
[kappa__x, Bx, Cx, Dx, Ex, SVx] = MF96_FX0_coeffs(0, 0, GAMMA_vec(3), tyre_coeffs.FZ0, tyre_coeffs);

% Printing parameters of interest 
fprintf('Bx      = %6.3f\n',Bx);
fprintf('Cx      = %6.3f\n',Cx);
fprintf('mux      = %6.3f\n',Dx/tyre_coeffs.FZ0);
fprintf('Ex      = %6.3f\n',Ex);
fprintf('SVx     = %6.3f\n',SVx);
fprintf('kappa_x = %6.3f\n',kappa__x);
fprintf('Kx      = %6.3f\n',Bx*Cx*Dx/tyre_coeffs.FZ0);


%% 1.12 Combined longitudinal slip
% Combined Slip, longitudinal force with variable side slip, zero camber

% ADDING MISSING PARAMETER FROM DATASET, INITIALISE AT 0
tyre_coeffs.rBx1 = 0;

% extract data with variable side slip
[TDataCombX, ~] = intersect_table_data(GAMMA_0, FZ_220);

% Plotting raw data 
figure()
plot(TDataCombX.SL,TDataCombX.FX,'.')
hold on
grid on
xlabel('$\kappa$ [-]')
ylabel('$F_{x}$ [N]')
title('Raw combined longitudinal force ')

% Guess values for parameters to be optimised
% Parameters list: [rBx1, rBx2, rCx1, rHx1]
P0 = [ 10, 9, 0.01, 0.01] ;

% Limits for parameters to be optimised
% Lower bounds 
lb = [ -20, -20, -20, -20];
% Upper bounds
ub = [ 20, 20, 20, 20];

% Useful vectors
KAPPA_vec = TDataCombX.SL;
ALPHA_vec = TDataCombX.SA;
FX_vec    = TDataCombX.FX;
SL_vec = -0.3:0.001:0.3;
FZ0 = mean(TDataCombX.FZ);

% Optimization of the parameters through least square approach
[P_comb_x,~,~] = fmincon(@(P)resid_comb_Fx(P,FX_vec, KAPPA_vec, ALPHA_vec,0,FZ0, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

% Change tyre data with new optimal values                             
tyre_coeffs.rBx1 = P_comb_x(1);
tyre_coeffs.rBx2 = P_comb_x(2);  
tyre_coeffs.rCx1 = P_comb_x(3);
tyre_coeffs.rHx1 = P_comb_x(4);

% Residuals from optimal parameters obtained
res_FX_comb_vec = resid_comb_Fx(P_comb_x,FX_vec,KAPPA_vec,ALPHA_vec,0 , FZ0,tyre_coeffs);

% Generating vectors for the plot
tmp_zeros = zeros(size(SL_vec));
tmp_ones = ones(size(SL_vec));
SA_vec_3 = -(3*to_rad)*tmp_ones;
SA_vec_6 = -(6*to_rad)*tmp_ones;

FX_comb_vec1 = MF96_COMB_FX_vec(SL_vec,tmp_zeros ,tmp_zeros, tyre_coeffs.FZ0.*tmp_ones,tyre_coeffs);
FX_comb_vec2 = MF96_COMB_FX_vec(SL_vec,SA_vec_3 ,tmp_zeros, tyre_coeffs.FZ0.*tmp_ones,tyre_coeffs);
FX_comb_vec3 = MF96_COMB_FX_vec(SL_vec,SA_vec_6 ,tmp_zeros, tyre_coeffs.FZ0.*tmp_ones,tyre_coeffs);

% Plotting Combined longitudinal slip fit
figure()
plot(TDataCombX.SL,TDataCombX.FX,'o')
hold on
grid on
plot(SL_vec,FX_comb_vec1,'-','LineWidth',2)
plot(SL_vec,FX_comb_vec2,'-','LineWidth',2)
plot(SL_vec,FX_comb_vec3,'-','LineWidth',2)
title('Combined Long. Slip fit for $F_z = 220$N, zero camber, different side slip.')
xlabel('$\kappa$ [-]')
ylabel('$F_{x}$ [N]')
legend({'Raw Data','$\alpha_{0}$','$\alpha_{-3}$','$\alpha_{-6}$'})

% Printing the R-squared
fprintf('Printing R-squared for Combined Long. slip:\n');
fprintf('R-squared = %6.3f\n',1-res_FX_comb_vec);


%% 1.13 Weighting Function as function of kappa
% Pre allocating memory
Gxa_vec_K_1 = zeros(size(SL_vec));
Gxa_vec_K_2 = zeros(size(SL_vec));
Gxa_vec_K_3 = zeros(size(SL_vec));

% Generating vectors for plot
for i = 1:length(SL_vec)
[Gxa_vec_K_1(i), ~, ~] = MF96_FXFYCOMB_coeffs(SL_vec(i),tmp_zeros(i),0,FZ0,tyre_coeffs);
end
for i = 1:length(SL_vec)
[Gxa_vec_K_2(i), ~, ~] = MF96_FXFYCOMB_coeffs(SL_vec(i),SA_vec_3(i),0,FZ0,tyre_coeffs);
end
for i = 1:length(SL_vec)
[Gxa_vec_K_3(i), ~, ~] = MF96_FXFYCOMB_coeffs(SL_vec(i),SA_vec_6(i),0,FZ0,tyre_coeffs);
end

% Plotting longitudinal tyre stiffness as function of Kappa 
figure('Name','Gxa(K)')
hold on
grid on
plot(SL_vec,Gxa_vec_K_1,'-','LineWidth',2)
plot(SL_vec,Gxa_vec_K_2,'-','LineWidth',2)
plot(SL_vec,Gxa_vec_K_3,'-','LineWidth',2)
title('Weigthting function Gxa as a function of the Long. slip')
xlabel('$\kappa$ [-]')
ylabel('$Gxa$ [-]')
legend({'$\alpha_0$','$\alpha_{-3}$','$\alpha_{-6}$'})

%% 1.14 Weighting Function as function of alpha
ALPHA_vec = (-30:0.1:30);
ALPHA_vec = ALPHA_vec*to_rad;

% Pre allocating memory
Gxa_vec_A_1 = zeros(size(ALPHA_vec));
Gxa_vec_A_2 = zeros(size(ALPHA_vec));
Gxa_vec_A_3 = zeros(size(ALPHA_vec));
Gxa_vec_A_4 = zeros(size(ALPHA_vec)); 
Gxa_vec_A_5 = zeros(size(ALPHA_vec));
Gxa_vec_A_6 = zeros(size(ALPHA_vec));

% Generating vectors for the stiffness plot
for i = 1:length(ALPHA_vec)
[Gxa_vec_A_1(i), ~, ~] = MF96_FXFYCOMB_coeffs(0,ALPHA_vec(i),0,FZ0,tyre_coeffs);
end
for i = 1:length(ALPHA_vec)
[Gxa_vec_A_2(i), ~, ~] = MF96_FXFYCOMB_coeffs(0.1,ALPHA_vec(i),0,FZ0,tyre_coeffs);
end
for i = 1:length(ALPHA_vec)
[Gxa_vec_A_3(i), ~, ~] = MF96_FXFYCOMB_coeffs(0.2,ALPHA_vec(i),0,FZ0,tyre_coeffs);
end
for i = 1:length(ALPHA_vec)
[Gxa_vec_A_4(i), ~, ~] = MF96_FXFYCOMB_coeffs(0.5,ALPHA_vec(i),0,FZ0,tyre_coeffs);
end
for i = 1:length(ALPHA_vec)
[Gxa_vec_A_5(i), ~, ~] = MF96_FXFYCOMB_coeffs(0.8,ALPHA_vec(i),0,FZ0,tyre_coeffs);
end
for i = 1:length(ALPHA_vec)
[Gxa_vec_A_6(i), ~, ~] = MF96_FXFYCOMB_coeffs(1,ALPHA_vec(i),0,FZ0,tyre_coeffs);
end

% Plotting Longitudinal tyre stiffness as a function of side slip
figure('Name','Gxa(\alpha)')
hold on
grid on
plot(ALPHA_vec,Gxa_vec_A_1,'-','LineWidth',2)
plot(ALPHA_vec,Gxa_vec_A_2,'-','LineWidth',2)
plot(ALPHA_vec,Gxa_vec_A_3,'-','LineWidth',2)
plot(ALPHA_vec,Gxa_vec_A_4,'-','LineWidth',2)
plot(ALPHA_vec,Gxa_vec_A_5,'-','LineWidth',2)
plot(ALPHA_vec,Gxa_vec_A_6,'-','LineWidth',2)
title('weigthting function Gxa as a function of the side slip')
xlabel('$\alpha$ [rad]')
ylabel('$Gxa$ [-]')
legend({'$\kappa = 0 $','$\kappa = 0.1 $','$\kappa = 0.2 $','$\kappa = 0.5 $','$\kappa = 0.8 $','$\kappa = 1 $'})


%% 1.15 Save tyre data structure to mat file

fprintf('Saving dataset ...\n');
save(['tyre_' data_set,'.mat'],'tyre_coeffs');
fprintf('Dataset saved!\n');



%% CLEAN UP 
% /////////////////////////////////////////////////////////////////////// %
% //////////////// CLEANING VARIABLES FOR LATERAL DATASET /////////////// %
% /////////////////////////////////////////////////////////////////////// %

% Cleaning variables to avoid mistakes during code writing
% clearvars;
% clc;
% close all;

% Set LaTeX as default interpreter for axis labels, ticks and legends
% set(0,'defaulttextinterpreter','latex')
% set(groot, 'defaultAxesTickLabelInterpreter','latex');
% set(groot, 'defaultLegendInterpreter','latex');
% 
% set(0,'DefaultFigureWindowStyle','docked');
% set(0,'defaultAxesFontSize',  16)
% set(0,'DefaultLegendFontSize',16)
% 
% addpath('tyre_lib/')

% Setting up useful parameters
% to_rad = pi/180;
% to_deg = 180/pi;
% data_set_path = 'dataset/';
% diameter = 18*2.56; %
% Fz0 = 220;   % [N] nominal load is given
% R0  = diameter/2/100; % [m] get from nominal load R0 (m) *** TO BE CHANGED ***



%% 2 LATERAL CASE
% /////////////////////////////////////////////////////////////////////// %


%% 2.1 Loading pure lateral dataset
fprintf('Loading dataset ...')
data_set = 'B1464run23'; % pure lateral forces
switch data_set
  case 'B1464run23'
  load ([data_set_path, 'Hoosier_B1464run23.mat']); % pure lateral
  cut_start = 31331;    % cutting useless initial data
  cut_end   = 54500; 
  case 'B1464run30'
  load ([data_set_path, 'Hoosier_B1464run30.mat']); % pure longitudinal
  cut_start = 19028;
  cut_end   = 37643;
  case 'B1464run58'
  load ([data_set_path, 'Goodyear_B1464run58.mat']); % pure longitudinal
  cut_start = 19028;
  cut_end   = 37643;
  case 'B1464run13'
  load ([data_set_path, 'Goodyear_B1464run13.mat']); % pure lateral
  cut_start = 27760;
  cut_end   = 54500;
  case 'B1464run51'
  load ([data_set_path, 'Continental_B1464run51.mat']); % pure longitudinal
  cut_start = 19028;
  cut_end   = 54500;
  case 'B1464run8'
  load ([data_set_path, 'Continental_B1464run30.mat']); % pure lateral
  cut_start = 27760;
  cut_end   = 37643;
  otherwise 
  error('Not found dataset: `%s`\n', data_set) ;
  
end

% select sample range
smpl_range = cut_start:cut_end;

fprintf('completed!\n')


%% 2.2 Plot raw data of pure lateral dataset
figure
tiledlayout(6,1)

ax_list(1) = nexttile; y_range = [min(min(-FZ),0) round(max(-FZ)*1.1)];
plot(-FZ)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Vertical force')
xlabel('Samples [-]')
ylabel('[N]')

ax_list(2) = nexttile; y_range = [min(min(IA),0) round(max(IA)*1.1)];
plot(IA)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Camber angle')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list(3) = nexttile; y_range = [min(min(SA),0) round(max(SA)*1.1)];
plot(SA)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Side slip')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list(4) = nexttile; y_range = [min(min(SL),0) round(max(SL)*1.1)];
plot(SL)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Longitudinal slip')
xlabel('Samples [-]')
ylabel('[-]')

ax_list(5) = nexttile; y_range = [min(min(P),0) round(max(P)*1.1)];
plot(P)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Tyre pressure')
xlabel('Samples [-]')
ylabel('[psi]')

ax_list(6) = nexttile;  y_range = [min(min(TSTC),0) round(max(TSTC)*1.1)];
plot(TSTC,'DisplayName','Center')
hold on
plot(TSTI,'DisplayName','Internal')
plot(TSTO,'DisplayName','Outboard')
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Tyre temperatures')
xlabel('Samples [-]')
ylabel('[degC]')

linkaxes(ax_list,'x')


%% 2.3 Select specific data of lateral dataset
% Cut undesirable data and select only 12 psi data
vec_samples = 1:1:length(smpl_range);
tyre_data_lat = table(); % create empty table

% store raw data in table
tyre_data_lat.SL =  SL(smpl_range);
tyre_data_lat.SA =  SA(smpl_range)*to_rad;
tyre_data_lat.FZ = -FZ(smpl_range); 
tyre_data_lat.FX =  FX(smpl_range);
tyre_data_lat.FY =  -FY(smpl_range);
tyre_data_lat.MZ =  -MZ(smpl_range);
tyre_data_lat.IA =  IA(smpl_range)*to_rad;

% Extract points at constant inclination angle
GAMMA_tol = 0.05*to_rad;
idy.GAMMA_0_LAT = 0.0*to_rad-GAMMA_tol < tyre_data_lat.IA & tyre_data_lat.IA < 0.0*to_rad+GAMMA_tol;
idy.GAMMA_1_LAT = 1.0*to_rad-GAMMA_tol < tyre_data_lat.IA & tyre_data_lat.IA < 1.0*to_rad+GAMMA_tol;
idy.GAMMA_2_LAT = 2.0*to_rad-GAMMA_tol < tyre_data_lat.IA & tyre_data_lat.IA < 2.0*to_rad+GAMMA_tol;
idy.GAMMA_3_LAT = 3.0*to_rad-GAMMA_tol < tyre_data_lat.IA & tyre_data_lat.IA < 3.0*to_rad+GAMMA_tol;
idy.GAMMA_4_LAT = 4.0*to_rad-GAMMA_tol < tyre_data_lat.IA & tyre_data_lat.IA < 4.0*to_rad+GAMMA_tol;
idy.GAMMA_5_LAT = 5.0*to_rad-GAMMA_tol < tyre_data_lat.IA & tyre_data_lat.IA < 5.0*to_rad+GAMMA_tol;
GAMMA_0_LAT  = tyre_data_lat( idy.GAMMA_0_LAT, : );
GAMMA_1_LAT  = tyre_data_lat( idy.GAMMA_1_LAT, : );
GAMMA_2_LAT  = tyre_data_lat( idy.GAMMA_2_LAT, : );
GAMMA_3_LAT  = tyre_data_lat( idy.GAMMA_3_LAT, : );
GAMMA_4_LAT  = tyre_data_lat( idy.GAMMA_4_LAT, : );
GAMMA_5_LAT  = tyre_data_lat( idy.GAMMA_5_LAT, : );

% Extract points at constant vertical load
% Test data done at: 
%  - 50lbf  ( 50*0.453592*9.81 =  223N )
%  - 150lbf (150*0.453592*9.81 =  667N )
%  - 200lbf (200*0.453592*9.81 =  890N )
%  - 250lbf (250*0.453592*9.81 = 1120N )

FZ_tol = 100;
idy.FZ_220_LAT  = 220-FZ_tol < tyre_data_lat.FZ & tyre_data_lat.FZ < 220+FZ_tol;
idy.FZ_440_LAT  = 440-FZ_tol < tyre_data_lat.FZ & tyre_data_lat.FZ < 440+FZ_tol;
idy.FZ_700_LAT  = 700-FZ_tol < tyre_data_lat.FZ & tyre_data_lat.FZ < 700+FZ_tol;
idy.FZ_900_LAT  = 900-FZ_tol < tyre_data_lat.FZ & tyre_data_lat.FZ < 900+FZ_tol;
idy.FZ_1120_LAT = 1120-FZ_tol < tyre_data_lat.FZ & tyre_data_lat.FZ < 1120+FZ_tol;
FZ_220_LAT  = tyre_data_lat( idy.FZ_220_LAT, : );
FZ_440_LAT  = tyre_data_lat( idy.FZ_440_LAT, : );
FZ_700_LAT  = tyre_data_lat( idy.FZ_700_LAT, : );
FZ_900_LAT  = tyre_data_lat( idy.FZ_900_LAT, : );
FZ_1120_LAT = tyre_data_lat( idy.FZ_1120_LAT, : );

% The slip angle is varied continuously between -4 and +12° and then
% between -12° and +4° for the pure slip case

% The slip angle is varied step wise for lateral slip tests
% 0° , - 3° , -6 °
SL_tol = 0.5*to_rad;
idy.SL_0    =  0-SL_tol          < tyre_data_lat.SL & tyre_data_lat.SL < 0+SL_tol;
idy.SL_3neg = -(3*to_rad+SL_tol) < tyre_data_lat.SL & tyre_data_lat.SL < -3*to_rad+SL_tol;
idy.SL_6neg = -(6*to_rad+SL_tol) < tyre_data_lat.SL & tyre_data_lat.SL < -6*to_rad+SL_tol;
SL_0     = tyre_data_lat( idy.SL_0, : );
SL_3neg  = tyre_data_lat( idy.SL_3neg, : );
SL_6neg  = tyre_data_lat( idy.SL_6neg, : );

figure()
tiledlayout(3,1)

ax_list(1) = nexttile;
plot(tyre_data_lat.IA*to_deg)
hold on
plot(vec_samples(idy.GAMMA_0_LAT),GAMMA_0_LAT.IA*to_deg,'.');
plot(vec_samples(idy.GAMMA_1_LAT),GAMMA_1_LAT.IA*to_deg,'.');
plot(vec_samples(idy.GAMMA_2_LAT),GAMMA_2_LAT.IA*to_deg,'.');
plot(vec_samples(idy.GAMMA_3_LAT),GAMMA_3_LAT.IA*to_deg,'.');
plot(vec_samples(idy.GAMMA_4_LAT),GAMMA_4_LAT.IA*to_deg,'.');
plot(vec_samples(idy.GAMMA_5_LAT),GAMMA_5_LAT.IA*to_deg,'.');
title('Camber angle')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list(2) = nexttile;
plot(tyre_data_lat.FZ)
hold on
plot(vec_samples(idy.FZ_220_LAT),FZ_220_LAT.FZ,'.');
plot(vec_samples(idy.FZ_440_LAT),FZ_440_LAT.FZ,'.');
plot(vec_samples(idy.FZ_700_LAT),FZ_700_LAT.FZ,'.');
plot(vec_samples(idy.FZ_900_LAT),FZ_900_LAT.FZ,'.');
plot(vec_samples(idy.FZ_1120_LAT),FZ_1120_LAT.FZ,'.');
title('Vertical force')
xlabel('Samples [-]')
ylabel('[N]')

ax_list(3) = nexttile;
plot(tyre_data_lat.SL*to_deg)
hold on
plot(vec_samples(idy.SL_0),   SL_0.SL*to_deg,'.');
plot(vec_samples(idy.SL_3neg),SL_3neg.SL*to_deg,'.');
plot(vec_samples(idy.SL_6neg),SL_6neg.SL*to_deg,'.');
title('Long. slip')
xlabel('Samples [-]')
ylabel('[rad]')


%% 2.4 Intersecting tables to obtain specific sub-datasets
[TData0Y, ~] = intersect_table_data(SL_0, GAMMA_0_LAT, FZ_1120_LAT);
% After many attempts, we've decided to change the initial dataset,
% swipping from FZ_220 to FZ_1120 and we get important improvements.


%% 2.5 Plotting selected data for lateral dataset
figure('Name','Selected-data')
plot_selected_data(TData0Y);


%% 2.6 Fitting
% Initializing tyre data
tyre_coeffs_lat = initialise_tyre_data(R0, 1120);


%% 2.7 Lateral slip
% Fitting with Fz = 1120N and camber = 0, kappa = 0, VX = 10

% Plotting raw data 
figure()
plot(TData0Y.SA,TData0Y.FY,'.')
hold on
grid on
xlabel('$\alpha$ [rad]')
ylabel('$F_y [N]$')
title('Raw data')

% Guess values for parameters to be optimised
% Parameters list: [pCy1 pDy1 pEy1 pHy1  pKy1  pKy2  pVy1] 
P0 = [ 1, 2, 1, 0, 0, 1, 0];        

% Parameters limits
% Lower bound
lb = [0,0,0,-1,0,0,0];
% Upper bounds
ub = [1,15,2,0,25,2,1];

% Useful vectors
ALPHA_vec = TData0Y.SA;
FY_vec    = TData0Y.FY;
SA_vec = -0.3:0.001:0.3;
FZ0 = mean(TData0Y.FZ);

% Optimization of the parameters through least square approach
[P_fz_nom_Y,~,~] = fmincon(@(P)resid_pure_Fy(P,FY_vec, ALPHA_vec,0,FZ0, tyre_coeffs_lat),...
                               P0,[],[],[],[],lb,ub);

% Updating tyre data with new optimal values   
tyre_coeffs_lat.pCy1 = P_fz_nom_Y(1);
tyre_coeffs_lat.pDy1 = P_fz_nom_Y(2);  
tyre_coeffs_lat.pEy1 = P_fz_nom_Y(3);
tyre_coeffs_lat.pHy1 = P_fz_nom_Y(4);
tyre_coeffs_lat.pKy1 = P_fz_nom_Y(5); 
tyre_coeffs_lat.pKy2 = P_fz_nom_Y(6);
tyre_coeffs_lat.pVy1 = P_fz_nom_Y(7);

% Residuals from optimal parameters
res_FY0_nom_vec = resid_pure_Fy(P_fz_nom_Y,FY_vec, ALPHA_vec,0,FZ0, tyre_coeffs_lat);

% Generating vector for plot
FY0_fz_nom_vec = MF96_FY0_vec(zeros(size(SA_vec)), SA_vec, zeros(size(SA_vec)), ...
                              FZ0.*ones(size(SA_vec)),tyre_coeffs_lat);

% Plotting pure lateral slip for fixed Fz
figure('Name','Fy0(Fz0)')
plot(TData0Y.SA,TData0Y.FY,'o')
hold on
grid on
plot(SA_vec,FY0_fz_nom_vec,'-','LineWidth',2)
xlabel('$\alpha$ [rad]')
ylabel('$F_{y0}$ [N]')
title('Pure Lat. Slip at vertical load $F_z = 1120$ N, no camber, zero $\kappa$, Press. $12$ psi')

% Printing the R-squared
fprintf('Printing R-squared for Pure Lat. slip, no camber, no Long. slip :\n');
fprintf('R-squared = %6.3f\n',1-res_FY0_nom_vec);


%% 2.8 Lateral slip with variable vertical load
% Extracting anf fitting coefficients with variable load
[TDataDFzY, ~] = intersect_table_data(SL_0, GAMMA_0_LAT);

% Plotting raw data 
figure()
plot(TDataDFzY.SA,TDataDFzY.FY,'.')
hold on
grid on
xlabel('$\alpha$ [rad]')
ylabel('$F_y [N]$')
title('Raw Data')

% Initial values of parameters to optimize
% Parameters to optimize: [pDy2  pEy2 pHy2 pVy2]  
P0 = [2,2,1,1]; %[1.5,2,1,1]; 

% Limits for parameters to be optimised
% Lower bound
lb = [-4,-2,-5,-2]; 
% Upper bound
ub = [4,2,5,2]; 

% Generating and extracting usefull vectors from table
ALPHA_vec = TDataDFzY.SA;
FY_vec    = TDataDFzY.FY;
FZ_vec    = TDataDFzY.FZ;
SA_vec = -0.3:0.001:0.3;

% Optimization of the parameters through least square approach
[P_dfz_Y,~,~] = fmincon(@(P)resid_pure_Fy_varFz(P,FY_vec, ALPHA_vec,0,FZ_vec, tyre_coeffs_lat),...
                               P0,[],[],[],[],lb,ub);

% Change tyre data with new optimal values    
tyre_coeffs_lat.pDy2 = P_dfz_Y(1);
tyre_coeffs_lat.pEy2 = P_dfz_Y(2);  
tyre_coeffs_lat.pHy2 = P_dfz_Y(3);
tyre_coeffs_lat.pVy2 = P_dfz_Y(4);

% Residuals from optimized parameters
res_FY0_dfz_vec = resid_pure_Fy_varFz(P_dfz_Y,FY_vec,SA_vec,0 , FZ_vec,tyre_coeffs_lat);

% Generating vectors for plot
tmp_zeros = zeros(size(SA_vec));
tmp_ones = ones(size(SA_vec));

FY0_fz_var_vec1 = MF96_FY0_vec(tmp_zeros,SA_vec ,tmp_zeros, mean(FZ_220_LAT.FZ)*tmp_ones,tyre_coeffs_lat);
FY0_fz_var_vec2 = MF96_FY0_vec(tmp_zeros,SA_vec ,tmp_zeros, mean(FZ_440_LAT.FZ)*tmp_ones,tyre_coeffs_lat);
FY0_fz_var_vec3 = MF96_FY0_vec(tmp_zeros,SA_vec ,tmp_zeros, mean(FZ_700_LAT.FZ)*tmp_ones,tyre_coeffs_lat);
FY0_fz_var_vec4 = MF96_FY0_vec(tmp_zeros,SA_vec ,tmp_zeros, mean(FZ_900_LAT.FZ)*tmp_ones,tyre_coeffs_lat);
FY0_fz_var_vec5 = MF96_FY0_vec(tmp_zeros,SA_vec ,tmp_zeros, mean(FZ_1120_LAT.FZ)*tmp_ones,tyre_coeffs_lat);

% Plotting pure lateral slip for different Fz
figure('Name','Fy0(Fz0)')
plot(TDataDFzY.SA,TDataDFzY.FY,'o')
hold on
grid on
plot(SA_vec,FY0_fz_var_vec1,'-','LineWidth',2)
plot(SA_vec,FY0_fz_var_vec2,'-','LineWidth',2)
plot(SA_vec,FY0_fz_var_vec3,'-','LineWidth',2)
plot(SA_vec,FY0_fz_var_vec4,'-','LineWidth',2)
plot(SA_vec,FY0_fz_var_vec5,'-','LineWidth',2)
legend({'Data','$Fz_{220}$','$Fz_{440}$','$Fz_{700}$','$Fz_{900}$','$Fz_{1120}$'})
xlabel('$\alpha$ [deg]')
ylabel('$F_{y0}$ [N]')
title('Pure Lat. Slip for different vertical load $F_z$')

% Generating vectors for tyre stiffness plot
[alpha__y, By, Cy, Dy, Ey, SVy, ~, ~, ~] =MF96_FY0_coeffs(0, 0, 0, mean(FZ_220_LAT.FZ), tyre_coeffs_lat);
Calfa_vec1_0 = magic_formula_stiffness(alpha__y, By, Cy, Dy, Ey, SVy);
[alpha__y, By, Cy, Dy, Ey, SVy, ~, ~, ~] =MF96_FY0_coeffs(0, 0, 0, mean(FZ_440_LAT.FZ), tyre_coeffs_lat);
Calfa_vec2_0 = magic_formula_stiffness(alpha__y, By, Cy, Dy, Ey, SVy);
[alpha__y, By, Cy, Dy, Ey, SVy, ~, ~, ~] =MF96_FY0_coeffs(0, 0, 0, mean(FZ_700_LAT.FZ), tyre_coeffs_lat);
Calfa_vec3_0 = magic_formula_stiffness(alpha__y, By, Cy, Dy, Ey, SVy);
[alpha__y, By, Cy, Dy, Ey, SVy, ~, ~, ~] =MF96_FY0_coeffs(0, 0, 0, mean(FZ_900_LAT.FZ), tyre_coeffs_lat);
Calfa_vec4_0 = magic_formula_stiffness(alpha__y, By, Cy, Dy, Ey, SVy);
[alpha__y, By, Cy, Dy, Ey, SVy, ~, ~, ~] =MF96_FY0_coeffs(0, 0, 0, mean(FZ_1120_LAT.FZ), tyre_coeffs_lat);
Calfa_vec5_0 = magic_formula_stiffness(alpha__y, By, Cy, Dy, Ey, SVy);

Calfa_vec1 = MF96_CorneringStiffness_lat(tmp_zeros,SA_vec ,tmp_zeros, mean(FZ_220_LAT.FZ)*tmp_ones,tyre_coeffs_lat);
Calfa_vec2 = MF96_CorneringStiffness_lat(tmp_zeros,SA_vec ,tmp_zeros, mean(FZ_440_LAT.FZ)*tmp_ones,tyre_coeffs_lat);
Calfa_vec3 = MF96_CorneringStiffness_lat(tmp_zeros,SA_vec ,tmp_zeros, mean(FZ_700_LAT.FZ)*tmp_ones,tyre_coeffs_lat);
Calfa_vec4 = MF96_CorneringStiffness_lat(tmp_zeros,SA_vec ,tmp_zeros, mean(FZ_900_LAT.FZ)*tmp_ones,tyre_coeffs_lat);
Calfa_vec5 = MF96_CorneringStiffness_lat(tmp_zeros,SA_vec ,tmp_zeros, mean(FZ_1120_LAT.FZ)*tmp_ones,tyre_coeffs_lat);

% Plotting cornering tyre stiffness
figure('Name','C_{F_\alpha}')
subplot(2,1,1)
hold on
grid on
plot(mean(FZ_220_LAT.FZ),Calfa_vec1_0,'+','LineWidth',2)
plot(mean(FZ_440_LAT.FZ),Calfa_vec2_0,'+','LineWidth',2)
plot(mean(FZ_700_LAT.FZ),Calfa_vec3_0,'+','LineWidth',2)
plot(mean(FZ_900_LAT.FZ),Calfa_vec4_0,'+','LineWidth',2)
plot(mean(FZ_1120_LAT.FZ),Calfa_vec5_0,'+','LineWidth',2)
legend({'$Fz_{220}$','$Fz_{440}$','$Fz_{700}$','$Fz_{900}$','$Fz_{1120}$'})

subplot(2,1,2)
hold on
grid on
plot(SA_vec,Calfa_vec1,'-','LineWidth',2)
plot(SA_vec,Calfa_vec2,'-','LineWidth',2)
plot(SA_vec,Calfa_vec3,'-','LineWidth',2)
plot(SA_vec,Calfa_vec4,'-','LineWidth',2)
plot(SA_vec,Calfa_vec5,'-','LineWidth',2)
legend({'$Fz_{220}$','$Fz_{440}$','$Fz_{700}$','$Fz_{900}$','$Fz_{1120}$'})

% Printing the R-squared
fprintf('Printing R-squared for Pure Lat. slip, variable loads :\n');
fprintf('R-squared = %6.3f\n',1-res_FY0_dfz_vec);


%% 2.9 Lateral slip with variable camber angle 
% Fit coefficient with variable camber for lateral dataset
% Extract data with variable camber
[TDataGammaY, ~] = intersect_table_data( SL_0, FZ_1120_LAT );

% Guess values for parameters to be optimised
% Parameters to optimize: [pDy3  pEy3  pEy4  pHy4  pKy3  pVy3   pVy4]
P0 = [1, 1, 1, 1, 0, 0, 0]; 

% Limits for parameters to be optimised
% Lower bound
lb = [-3,0,-2,0,0,-4,-1];
% Upper bound
ub = [3,2,2,2,2,1,1];

% Generating and extracting vectors needed
ALPHA_vec = TDataGammaY.SA;
GAMMA_vec = TDataGammaY.IA;
FY_vec = TDataGammaY.FY;

% Plotting raw data
figure()
plot(ALPHA_vec,FY_vec);
hold on
xlabel('$\alpha$ [rad]')
ylabel('$F_y$ [N]')
grid on
title('Raw Data')

% Optimization of the parameters through least square approach
[P_varGamma_Y,~,~] = fmincon(@(P)resid_pure_Fy_varGamma(P,FY_vec, ALPHA_vec,GAMMA_vec,tyre_coeffs_lat.FZ0, tyre_coeffs_lat),...
                               P0,[],[],[],[],lb,ub);

% Change tyre data with new optimal values                             
tyre_coeffs_lat.pDy3 = P_varGamma_Y(1);
tyre_coeffs_lat.pEy3 = P_varGamma_Y(2);
tyre_coeffs_lat.pEy4 = P_varGamma_Y(3);
tyre_coeffs_lat.pHy4 = P_varGamma_Y(4);
tyre_coeffs_lat.pKy3 = P_varGamma_Y(5);
tyre_coeffs_lat.pVy3 = P_varGamma_Y(6);
tyre_coeffs_lat.pVy4 = P_varGamma_Y(7);

% Extract data with variable camber for plot
[TDataGammaY0, ~] = intersect_table_data(GAMMA_0_LAT, SL_0, FZ_220_LAT );
[TDataGammaY1, ~] = intersect_table_data(GAMMA_1_LAT, SL_0, FZ_220_LAT );
[TDataGammaY2, ~] = intersect_table_data(GAMMA_2_LAT, SL_0, FZ_220_LAT );
[TDataGammaY3, ~] = intersect_table_data(GAMMA_3_LAT, SL_0, FZ_220_LAT );
[TDataGammaY4, ~] = intersect_table_data(GAMMA_4_LAT, SL_0, FZ_220_LAT );

% Generating vectors for plot
FY0_varGamma_vec_0 = MF96_FY0_vec(zeros(size(TDataGammaY0.IA)), TDataGammaY0.SA ,TDataGammaY0.IA ,tyre_coeffs_lat.FZ0*ones(size(TDataGammaY0.IA)) ,tyre_coeffs_lat );
FY0_varGamma_vec_1 = MF96_FY0_vec(zeros(size(TDataGammaY1.IA)), TDataGammaY1.SA ,TDataGammaY1.IA ,tyre_coeffs_lat.FZ0*ones(size(TDataGammaY1.IA)) ,tyre_coeffs_lat );
FY0_varGamma_vec_2 = MF96_FY0_vec(zeros(size(TDataGammaY2.IA)), TDataGammaY2.SA ,TDataGammaY2.IA ,tyre_coeffs_lat.FZ0*ones(size(TDataGammaY2.IA)) ,tyre_coeffs_lat );
FY0_varGamma_vec_3 = MF96_FY0_vec(zeros(size(TDataGammaY3.IA)), TDataGammaY3.SA ,TDataGammaY3.IA ,tyre_coeffs_lat.FZ0*ones(size(TDataGammaY3.IA)) ,tyre_coeffs_lat );
FY0_varGamma_vec_4 = MF96_FY0_vec(zeros(size(TDataGammaY4.IA)), TDataGammaY4.SA ,TDataGammaY4.IA ,tyre_coeffs_lat.FZ0*ones(size(TDataGammaY4.IA)) ,tyre_coeffs_lat );

% Plotting the pure lateral slip fit for different gamma
figure('Name','Fy0 vs Gamma')
plot(ALPHA_vec,FY_vec,'o')
hold on
grid on
plot(TDataGammaY0.SA,FY0_varGamma_vec_0,'-','LineWidth',1)
plot(TDataGammaY1.SA,FY0_varGamma_vec_1,'-','LineWidth',1)
plot(TDataGammaY2.SA,FY0_varGamma_vec_2,'-','LineWidth',1)
plot(TDataGammaY3.SA,FY0_varGamma_vec_3,'-','LineWidth',1)
plot(TDataGammaY4.SA,FY0_varGamma_vec_4,'-','LineWidth',1)
xlabel('$\alpha$ [-]')
ylabel('$F_{y0}$ [N]')
title('Pure Lat. Slip for different camber angle, vertical load $F_z = 1120$ N')
legend({'Raw data','$\gamma = 0$','$\gamma = 1$','$\gamma = 2$','$\gamma = 3$','$\gamma = 4$'})

% Calculate the residuals with the optimal solution found above
res_Fy0_varGamma  = resid_pure_Fy_varGamma(P_varGamma_Y,FY_vec, ALPHA_vec,GAMMA_vec,tyre_coeffs_lat.FZ0, tyre_coeffs_lat);

% Printing the R-squared
fprintf('Printing R-squared for pure Lat. slip, variable camber angle:\n');
fprintf('R-squared = %6.3f\n',1-res_Fy0_varGamma);

% Obtaining parameters to print
[alpha__y, By, Cy, Dy, Ey, SVy, ~, ~, ~] = MF96_FY0_coeffs(0, 0, GAMMA_vec(3), tyre_coeffs_lat.FZ0, tyre_coeffs_lat);

% Printing parameters of interest
fprintf('By      = %6.3f\n',By);
fprintf('Cy      = %6.3f\n',Cy);
fprintf('muy      = %6.3f\n',Dy/tyre_coeffs_lat.FZ0);
fprintf('Ey      = %6.3f\n',Ey);
fprintf('SVy     = %6.3f\n',SVy);
fprintf('alpha_y = %6.3f\n',alpha__y);
fprintf('Ky      = %6.3f\n',By*Cy*Dy/tyre_coeffs_lat.FZ0);



%% 3 SELF ALIGNING MOMENT
% /////////////////////////////////////////////////////////////////////// %


%% 3.1 Aligning moment with fixed parameters
% Aligning moment Mz with Fz = 1120 N, zero camber and no longitudinal slip
% Adding missing scaling factor to the tyre data, initializing to 1
tyre_coeffs_lat.LMR = 1;
tyre_coeffs_lat.LT = 1;

% Plotting raw data
figure()
plot(TData0Y.SA,TData0Y.MZ,'.')
hold on
grid on
title('Raw Data Plot')
xlabel('$\alpha$ [rad]')
ylabel('$M_z$ [Nm]' )

% Guess values for parameters to be optimised
% Parameters to optimize : [qBz1 qBz10 qBz9 qCz1 qDz1 qDz2 qDz3 qDz4 qDz6 qEz1 qEz4 qHz1]
P0 = [5, 0, 10, 1, 0.01, -0.001, 0.1, 0, 0, -2, 0.1, 0]; %5, 0, 10, 1, 0.01, -0.001, 0.1, 0, 0, -2, 0.1, 0

% Limits for parameters to be optimised
% Lower bound
lb = []; 
% Upper bound
ub = [];

% Extraction and generation of usefull vectors from the data table
ALPHA_vec = TData0Y.SA;
MZ_vec    = TData0Y.MZ; 
SA_vec = -0.3:0.001:0.3;
FZ0 = mean(TData0Y.FZ);

% Optimization of the parameters through least square approach
[P_Mz_fz_nom,~,~] = fmincon(@(P)resid_pure_Mz(P,MZ_vec, ALPHA_vec,0,FZ0, tyre_coeffs_lat),...
                               P0,[],[],[],[],lb,ub);

% Update tyre data with new optimal values   
tyre_coeffs_lat.qBz1 = P_Mz_fz_nom(1);
tyre_coeffs_lat.qBz10 = P_Mz_fz_nom(2);  
tyre_coeffs_lat.qBz9 = P_Mz_fz_nom(3);
tyre_coeffs_lat.qCz1 = P_Mz_fz_nom(4);
tyre_coeffs_lat.qDz1 = P_Mz_fz_nom(5); 
tyre_coeffs_lat.qDz2 = P_Mz_fz_nom(6);
tyre_coeffs_lat.qDz3 = P_Mz_fz_nom(7);
tyre_coeffs_lat.qDz4 = P_Mz_fz_nom(8);
tyre_coeffs_lat.qDz6 = P_Mz_fz_nom(9); 
tyre_coeffs_lat.qEz1 = P_Mz_fz_nom(10);
tyre_coeffs_lat.qEz4 = P_Mz_fz_nom(11);
tyre_coeffs_lat.qHz1 = P_Mz_fz_nom(12);

% Residuals from optimized tyre parameters
res_MZ0_nom_vec = resid_pure_Mz(P_Mz_fz_nom,MZ_vec, ALPHA_vec,0,FZ0, tyre_coeffs_lat);

% Generating vector for plot
MZ0_fz_nom_vec = MF96_MZ0_vec(zeros(size(SA_vec)), SA_vec, zeros(size(SA_vec)), ...
                              FZ0.*ones(size(SA_vec)),tyre_coeffs_lat);

% Plotting the aligning moment fit
figure('Name','Mz0(Fz0)')
plot(TData0Y.SA,TData0Y.MZ,'o')
hold on
grid on
plot(SA_vec,MZ0_fz_nom_vec,'-','LineWidth',2)
xlabel('$\alpha$ [rad]')
ylabel('$M_{z0}$ [N]')
title('Aligning moment $M_z$ at load $F_z = 1120$ N, no camber, zero $\kappa$, Press. $12$ psi')

% Printing the R-squared
fprintf('Printing R-squared for Aligning moment, no camber, no Long. slip :\n');
fprintf('R-squared = %6.3f\n',1-res_MZ0_nom_vec);


%% 3.2 Aligning moment Mz with variable vertical load
% Plotting raw data
figure()
plot(TDataDFzY.SA,TDataDFzY.MZ,'.')
hold on
grid on
title('Raw Data Plot')
xlabel('$\alpha$ [rad]')
ylabel('$M_z$ [Nm]' )

% Initial values for parameters to be optimised
% Parameters to optimize : [qBz2 qBz3 qDz7 qEz2 qEz3]
P0 = [ -1.33, 0.6, -0.002, 0.9, 1];

% Limits for parameters to be optimised
% Lower bound
lb = [ -2, -2, -2, -2, -2];
% Upper bound
ub = [ 2, 2, 2, 2, 2];

% Extraction and generation of usefull vectors from the data table
ALPHA_vec = TDataDFzY.SA;
MZ_vec    = TDataDFzY.MZ;
FZ_vec    = TDataDFzY.FZ;
SA_vec = -0.3:0.001:0.3;

% Optimization of the parameters through least square approach
[P_Mz_dfz,~,~] = fmincon(@(P)resid_pure_Mz_varFz(P,MZ_vec, ALPHA_vec,0,FZ_vec, tyre_coeffs_lat),...
                               P0,[],[],[],[],lb,ub);

% Change tyre data with new optimal values    
tyre_coeffs_lat.qBz2 = P_Mz_dfz(1) ;
tyre_coeffs_lat.qBz3 = P_Mz_dfz(2) ;  
tyre_coeffs_lat.qDz7 = P_Mz_dfz(3) ;
tyre_coeffs_lat.qEz2 = P_Mz_dfz(4) ;
tyre_coeffs_lat.qEz3 = P_Mz_dfz(5) ;

% Residuals from optimized tyre parameters
res_MZ0_dfz_vec = resid_pure_Mz_varFz(P_Mz_dfz,MZ_vec,SA_vec,0 , FZ_vec,tyre_coeffs_lat);

% Generating vectors for plot
tmp_zeros = zeros(size(SA_vec));
tmp_ones = ones(size(SA_vec));

MZ0_fz_var_vec1 = MF96_MZ0_vec(tmp_zeros,SA_vec ,tmp_zeros, FZ_220_LAT.FZ,tyre_coeffs_lat);
MZ0_fz_var_vec2 = MF96_MZ0_vec(tmp_zeros,SA_vec ,tmp_zeros, FZ_440_LAT.FZ,tyre_coeffs_lat);
MZ0_fz_var_vec3 = MF96_MZ0_vec(tmp_zeros,SA_vec ,tmp_zeros, FZ_700_LAT.FZ,tyre_coeffs_lat);
MZ0_fz_var_vec4 = MF96_MZ0_vec(tmp_zeros,SA_vec ,tmp_zeros, FZ_900_LAT.FZ,tyre_coeffs_lat);
MZ0_fz_var_vec5 = MF96_MZ0_vec(tmp_zeros,SA_vec ,tmp_zeros, FZ_1120_LAT.FZ,tyre_coeffs_lat);

% Plotting the aligning moment fit for different loads
figure('Name','Mz0(Fz0)')
plot(TDataDFzY.SA,TDataDFzY.MZ,'o')
hold on
grid on
plot(SA_vec,MZ0_fz_var_vec1,'-','LineWidth',2)
plot(SA_vec,MZ0_fz_var_vec2,'-','LineWidth',2)
plot(SA_vec,MZ0_fz_var_vec3,'-','LineWidth',2)
plot(SA_vec,MZ0_fz_var_vec4,'-','LineWidth',2)
plot(SA_vec,MZ0_fz_var_vec5,'-','LineWidth',2)
legend({'Data','$Fz_{220}$','$Fz_{440}$','$Fz_{700}$','$Fz_{900}$','$Fz_{1120}$'})
xlabel('$\alpha$ [rad]')
ylabel('$M_{z0}$ [Nm]')
title('Aligning moment for different vertical load $F_z$')

% Printing the R-squared
fprintf('Printing R-squared for alinging moment, variable loads :\n');
fprintf('R-squared = %6.3f\n',1-res_MZ0_dfz_vec);


%% 3.3 Aligning moment Mz with variable camber angle
% Extraction and generation of usefull vectors from the data table
ALPHA_vec = TDataGammaY.SA;
GAMMA_vec = TDataGammaY.IA; 
MZ_vec    = TDataGammaY.MZ;

zeros_vec = zeros(size(TDataGammaY.SA));
ones_vec  = ones(size(TDataGammaY.SA));

% Plotting raw data
figure()
plot(ALPHA_vec,MZ_vec);
hold on
xlabel('$\alpha$ [rad]')
ylabel('$Mz_y$ [Nm]')
grid on
title('Raw Data')

% Initial values for parameters to be optimised
% Parameters to optimize : [qBz4  qBz5  qDz8  qDz9  qEz5  qHz3   qHz4]
P0 = [-0.1, 0.01, 0.05, -0.005, 0.1, -0.1, 0.3]; 

% Limits for parameters to be optimised
% [qBz4  qBz5  qDz8  qDz9  qEz5  qHz3   qHz4]
% Lower bound
lb = [ -1, -1, -1, -1, -1, -1, -1];
% Upper bound
ub = [ 1, 1, 1, 1, 1, 1, 1];

% Optimization of the parameters through least square approach
[P_Mz_varGamma,~,~] = fmincon(@(P)resid_pure_Mz_varGamma(P,MZ_vec, ALPHA_vec,GAMMA_vec,tyre_coeffs_lat.FZ0, tyre_coeffs_lat),...
                               P0,[],[],[],[],lb,ub);

% Change tyre data with new optimal values                             
tyre_coeffs_lat.qBz4 = P_Mz_varGamma(1);
tyre_coeffs_lat.qBz5 = P_Mz_varGamma(2);
tyre_coeffs_lat.qDz8 = P_Mz_varGamma(3);
tyre_coeffs_lat.qDz9 = P_Mz_varGamma(4);
tyre_coeffs_lat.qEz5 = P_Mz_varGamma(5);
tyre_coeffs_lat.qHz3 = P_Mz_varGamma(6);
tyre_coeffs_lat.qHz4 = P_Mz_varGamma(7);

% Generating vectors for plot
MZ0_varGamma_vec = MF96_MZ0_vec(zeros_vec, ALPHA_vec ,GAMMA_vec ,tyre_coeffs_lat.FZ0*ones_vec ,tyre_coeffs_lat );

% Plotting the aligning moment fit for different camber angles
figure('Name','Mz0 vs Gamma')
plot(ALPHA_vec,TDataGammaY.MZ,'o')
hold on
grid on
plot(ALPHA_vec,MZ0_varGamma_vec,'-')
xlabel('$\alpha$ [-]')
ylabel('$M_{z0}$ [Nm]')
title('Aligning moment for different camber angle, vertical load $F_z = 1120$ N')

% Residuals from optimized tyre parameters
res_Mz0_varGamma  = resid_pure_Mz_varGamma(P_Mz_varGamma,MZ_vec, ALPHA_vec,GAMMA_vec,tyre_coeffs_lat.FZ0, tyre_coeffs_lat);

% Printing the R-squared
fprintf('Printing R-squared for aligning moment, variable camber angle:\n');
fprintf('R-squared = %6.3f\n',1-res_Mz0_varGamma);

% Printing paramters of Mz
[Br, Bt, Ct, Dr, Dt, Et, ~, ~, ~, ~, alpha__r, alpha__t, ~] = MF96_MZ0_coeffs(0, 0, GAMMA_vec(3), tyre_coeffs_lat.FZ0, tyre_coeffs_lat);
fprintf('Br      = %6.3f\n',Br);
fprintf('Bt      = %6.3f\n',Bt);
fprintf('Ct      = %6.3f\n',Ct);
fprintf('Dr      = %6.3f\n',Dr);
fprintf('Dt      = %6.3f\n',Dt);
fprintf('Et      = %6.3f\n',Et);
fprintf('alpha_r     = %6.3f\n',alpha__r);
fprintf('alpha_t = %6.3f\n',alpha__t);



%% 4 COMBINED LATERAL SLIP
% /////////////////////////////////////////////////////////////////////// %


%% 4.1 Combined slip
% ADDING MISSING PARAMETER FROM DATASET, INITIALIZE AT 0
tyre_coeffs_lat.rBx1 = 1;

% Intersecting needed data into a table
[TDataCombY, ~] = intersect_table_data(GAMMA_0, FZ_700);

% -----------------------------------------  COMBINED LATERAL SLIP -----------------------------------------

% Plotting the raw data 
figure()
plot(TDataCombY.SL,TDataCombY.FY,'.')
hold on
grid on
xlabel('$\kappa$ [-]')
ylabel('$F_y$ [N]')
title('Raw combined lateral force')

% Initial values for parameters to be optimised
% Parameters to optimize : [rBy1, rBy2, rBy3, rCy1, rHy1, rVy1, rVy4, rVy5, rVy6]
P0 = [14,   13,   -0.49,    0.9,    0.03,   -0.27,    3.76,   -0.09,   28.38]; 
%-14, -13, +0.49, -0.9, -0.03, +0.27, -3.76, +0.09, -28.38
% 14,   13,   -0.49,    0.9,    0.03,   -0.27,    3.76,   -0.09,   28.38

% Limits for parameters to be optimised
% Upper bound
ub = [20,20,20,20,20,20,20,20,20];
% Lower bound
lb = [-5,-5,-5,-5,-5,-5,-5,-5,-5];

% Extraction and generation of usefull vectors from the data table
FZ0 = mean(TDataCombY.FZ);
KAPPA_vec = TDataCombY.SL;          
ALPHA_vec = TDataCombY.SA;
FY_vec    = TDataCombY.FY;

% Optimization of the parameters through least square approach
[P_comb_y,~,~] = fmincon(@(P)resid_comb_Fy(P,FY_vec, KAPPA_vec, ALPHA_vec,0,FZ0, tyre_coeffs_lat),...
                               P0,[],[],[],[],lb,ub);

% Changing tyre data with new optimal values                             
tyre_coeffs_lat.rBy1 = P_comb_y(1) ; 
tyre_coeffs_lat.rBy2 = P_comb_y(2) ;  
tyre_coeffs_lat.rBy3 = P_comb_y(3) ;
tyre_coeffs_lat.rCy1 = P_comb_y(4) ;
tyre_coeffs_lat.rHy1 = P_comb_y(5) ; 
tyre_coeffs_lat.rVy1 = P_comb_y(6) ;  
tyre_coeffs_lat.rVy4 = P_comb_y(7) ;
tyre_coeffs_lat.rVy5 = P_comb_y(8) ;
tyre_coeffs_lat.rVy6 = P_comb_y(9) ; 

% Residuals from optimized tyre parameters
res_FY_comb_vec = resid_comb_Fy(P_comb_y,FY_vec,KAPPA_vec,ALPHA_vec,0 , FZ0,tyre_coeffs_lat);

% Generating vectors for plot
[Comb_y_0, ~] = intersect_table_data(SA_0, GAMMA_0, FZ_700);
[Comb_y_1, ~] = intersect_table_data(SA_3neg, GAMMA_0, FZ_700);
[Comb_y_2, ~] = intersect_table_data(SA_6neg, GAMMA_0, FZ_700);

FY_comb_vec0 = MF96_COMB_FY_vec(Comb_y_0.SL,Comb_y_0.SA,zeros(size(Comb_y_0.SA)), FZ0*ones(size(Comb_y_0.SA)), tyre_coeffs_lat);
FY_comb_vec3 = MF96_COMB_FY_vec(Comb_y_1.SL,Comb_y_1.SA,zeros(size(Comb_y_1.SA)), FZ0*ones(size(Comb_y_1.SA)), tyre_coeffs_lat);
FY_comb_vec6 = MF96_COMB_FY_vec(Comb_y_2.SL,Comb_y_2.SA,zeros(size(Comb_y_2.SA)), FZ0*ones(size(Comb_y_2.SA)), tyre_coeffs_lat);

% % Plotting the combined lateral slip fit
figure()
hold on
grid on
plot(Comb_y_0.SL,Comb_y_0.FY,'or')
plot(Comb_y_1.SL,Comb_y_1.FY,'og')
plot(Comb_y_2.SL,Comb_y_2.FY,'ob')
plot(Comb_y_0.SL, FY_comb_vec0,'-r','LineWidth',2)
plot(Comb_y_1.SL, FY_comb_vec3,'-g','LineWidth',2)
plot(Comb_y_2.SL, FY_comb_vec6,'-b','LineWidth',2)
title('Combined Lateral Slip fit for $F_z = 700$N, zero camber, different side slip.')
xlabel('$\kappa$ [-]')
ylabel('$F_{y}$ [N]')
legend({'Raw data at $0^{\circ}$','Raw data at $3^{\circ}$','Raw data at $6^{\circ}$','$\alpha = 0^{\circ}$','$\alpha = 3^{\circ}$','$\alpha = 6^{\circ}$'})

% Printing the R-squared
fprintf('Printing R-squared for Combined Lat. slip:\n');
fprintf('R-squared = %6.3f\n',1-res_FY_comb_vec);


%% 4.2 Weighting Function as function of kappa 
% Initializing vectors
SL_vec = -0.3:0.001:0.3;

Gyk_vec_K_3 = zeros(size(SL_vec));
Gyk_vec_K_6 = zeros(size(SL_vec));
Gyk_vec_K_9 = zeros(size(SL_vec));
Gyk_vec_K_12 = zeros(size(SL_vec));
Gyk_vec_K_15 = zeros(size(SL_vec));
Gyk_vec_K_18 = zeros(size(SL_vec));

SA_3_vec = -3*to_rad*ones(size(SL_vec));
SA_6_vec = -6*to_rad*ones(size(SL_vec));
SA_9_vec = -9*to_rad*ones(size(SL_vec));
SA_12_vec = -12*to_rad*ones(size(SL_vec));
SA_15_vec = -15*to_rad*ones(size(SL_vec));
SA_18_vec = -18*to_rad*ones(size(SL_vec));

% Generating weighting function vectors
for i = 1:length(SA_3_vec)
[~, Gyk_vec_K_3(i), ~] = MF96_FXFYCOMB_coeffs(SL_vec(i),SA_3_vec(i),0,FZ0,tyre_coeffs_lat);
end
for i = 1:length(SA_6_vec)
[~, Gyk_vec_K_6(i), ~] = MF96_FXFYCOMB_coeffs(SL_vec(i),SA_6_vec(i),0,FZ0,tyre_coeffs_lat);
end
for i = 1:length(SA_9_vec)
[~, Gyk_vec_K_9(i), ~] = MF96_FXFYCOMB_coeffs(SL_vec(i),SA_9_vec(i),0,FZ0,tyre_coeffs_lat);
end
for i = 1:length(SA_12_vec)
[~, Gyk_vec_K_12(i), ~] = MF96_FXFYCOMB_coeffs(SL_vec(i),SA_12_vec(i),0,FZ0,tyre_coeffs_lat);
end
for i = 1:length(SA_15_vec)
[~, Gyk_vec_K_15(i), ~] = MF96_FXFYCOMB_coeffs(SL_vec(i),SA_15_vec(i),0,FZ0,tyre_coeffs_lat);
end
for i = 1:length(SA_18_vec)
[~, Gyk_vec_K_18(i), ~] = MF96_FXFYCOMB_coeffs(SL_vec(i),SA_18_vec(i),0,FZ0,tyre_coeffs_lat);
end

% Plotting the wighting function as a function of longitudinal slip
figure('Name','Gyk($\kappa$)')
hold on
grid on
plot(SL_vec,Gyk_vec_K_3,'-','LineWidth',2)
plot(SL_vec,Gyk_vec_K_6,'-','LineWidth',2)
plot(SL_vec,Gyk_vec_K_9,'-','LineWidth',2)
plot(SL_vec,Gyk_vec_K_12,'-','LineWidth',2)
plot(SL_vec,Gyk_vec_K_15,'-','LineWidth',2)
plot(SL_vec,Gyk_vec_K_18,'-','LineWidth',2)
title('Weighting function Gyk as a function of the Long. slip for Fz = 700 N')
xlabel('$\kappa$ [-]')
ylabel('$Gyk$ [-]')
legend({'$\alpha = -3^{\circ}$', '$\alpha = -6^{\circ}$', '$\alpha = -9^{\circ}$','$\alpha = -12^{\circ}$', '$\alpha = -15^{\circ}$' , '$\alpha = -18^{\circ}$',})


%% 4.3 Save tyre data structure to mat file
fprintf('Saving dataset ...\n');
save(['tyre_' data_set,'.mat'],'tyre_coeffs');
fprintf('Dataset saved!\n');
