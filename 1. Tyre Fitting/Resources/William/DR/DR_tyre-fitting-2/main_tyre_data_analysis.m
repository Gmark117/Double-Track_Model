%% Initialisation
clc
clearvars 
close all   

% Set LaTeX as default interpreter for axis labels, ticks and legends
set(0,'defaulttextinterpreter','latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

set(0,'DefaultFigureWindowStyle','docked');
set(0,'defaultAxesFontSize',20)
set(0,'DefaultLegendFontSize',20)

%% load data

addpath('tyre_lib/')

load 'dataset/front_longitudinal_dataset.mat'; 
table = front_longitudinal_dataset;

table = table(table.FZ ~= 0,:);

%% Longitudinal test

kappa_vec = linspace( -1.0, +1.0, 200 ).^3;                  % [-]
alpha_vec = [ 0.0, 2.5, 5.0, 7.5, 10.0 ]*pi/180;             % [rad]
fz_vec    = [ 100, 1000, 2500, 4500, 6500, 8500, 10000 ];    % [N]
gamma_vec = [ 0.0, 2.5, 5.0 ]*pi/180;                        % [rad]
Vx_vec    = [ 5.0, 20.0, 45.0, 90.0 ];                       % [m/s]

%% Plot raw data

plot_raw_data

%% Select some specific data

% Extract points at constant side slip

ALPHA_tol = 0.005;
ALPHA_00  = table( alpha_vec(1)-ALPHA_tol < table.ALPHA & table.ALPHA < alpha_vec(1)+ALPHA_tol, : );
ALPHA_25  = table( alpha_vec(2)-ALPHA_tol < table.ALPHA & table.ALPHA < alpha_vec(2)+ALPHA_tol, : );
ALPHA_50  = table( alpha_vec(3)-ALPHA_tol < table.ALPHA & table.ALPHA < alpha_vec(3)+ALPHA_tol, : );
ALPHA_75  = table( alpha_vec(4)-ALPHA_tol < table.ALPHA & table.ALPHA < alpha_vec(4)+ALPHA_tol, : );
ALPHA_100 = table( alpha_vec(5)-ALPHA_tol < table.ALPHA & table.ALPHA < alpha_vec(5)+ALPHA_tol, : );

% Extract points at constant vertical load
FZ_tol   = 250.0;
FZ_1  = table(  fz_vec(1)-FZ_tol < table.FZ & table.FZ < fz_vec(1)+FZ_tol, : );
FZ_2  = table(  fz_vec(2)-FZ_tol < table.FZ & table.FZ < fz_vec(2)+FZ_tol, : );
FZ_3  = table(  fz_vec(3)-FZ_tol < table.FZ & table.FZ < fz_vec(3)+FZ_tol, : );
FZ_4  = table(  fz_vec(4)-FZ_tol < table.FZ & table.FZ < fz_vec(4)+FZ_tol, : );
FZ_5  = table(  fz_vec(5)-FZ_tol < table.FZ & table.FZ < fz_vec(5)+FZ_tol, : );
FZ_6  = table(  fz_vec(6)-FZ_tol < table.FZ & table.FZ < fz_vec(6)+FZ_tol, : );
FZ_7  = table(  fz_vec(7)-FZ_tol < table.FZ & table.FZ < fz_vec(7)+FZ_tol, : );

% Extract points at constant inclination angle
GAMMA_tol = 0.1*pi/180;
GAMMA_00  = table( gamma_vec(1)*pi/180-GAMMA_tol < table.GAMMA & table.GAMMA < gamma_vec(1)*pi/180+GAMMA_tol, : );
GAMMA_25  = table( gamma_vec(2)*pi/180-GAMMA_tol < table.GAMMA & table.GAMMA < gamma_vec(2)*pi/180+GAMMA_tol, : );
GAMMA_50  = table( gamma_vec(3)*pi/180-GAMMA_tol < table.GAMMA & table.GAMMA < gamma_vec(3)*pi/180+GAMMA_tol, : );

% Extract points at constant road speed
VX_tol = 0.5;
VX_05  = table( Vx_vec(1)-VX_tol < table.VX & table.VX < Vx_vec(1)+VX_tol, : );
VX_20  = table( Vx_vec(2)-VX_tol < table.VX & table.VX < Vx_vec(2)+VX_tol, : );
VX_45  = table( Vx_vec(3)-VX_tol < table.VX & table.VX < Vx_vec(3)+VX_tol, : );
VX_90  = table( Vx_vec(4)-VX_tol < table.VX & table.VX < Vx_vec(4)+VX_tol, : );

%% Intersect tables to obtain specific sub-datasets

%[tableSelectedData, ~] = intersectTableData( KAPPA_00, GAMMA_00, VX_10 );
[TDataSub, ~] = intersectTableData( ALPHA_00, GAMMA_00, VX_05, FZ_5 );

%[tableSelectedData, ~] = intersectTableData( KAPPA_00, GAMMA_00, VX_10, FZ_6500 );

%%
%plot_selected_data
plot_selected_data(TDataSub)

%% FITTING 

% Tyre structure data initialization
tyre_data.Fz0             = fz_vec(5); % Fz0            
tyre_data.pCx1            = 1; % pCx1           
tyre_data.pDx1            = 1; % pDx1           
tyre_data.pDx2            = 0; % pDx2           
tyre_data.pDx3            = 0; % pDx3           
tyre_data.pEx1            = 0; % pEx1           
tyre_data.pEx2            = 0; % pEx2           
tyre_data.pEx3            = 0; % pEx3           
tyre_data.pEx4            = 0; % pEx4           
tyre_data.pHx1            = 0; % pHx1           
tyre_data.pHx2            = 0; % pHx2           
tyre_data.pKx1            = 1; % pKx1           
tyre_data.pKx2            = 0; % pKx2           
tyre_data.pKx3            = 0; % pKx3           
tyre_data.pVx1            = 0; % pVx1           
tyre_data.pVx2            = 0; % pVx2           
tyre_data.pCy1            = 1; % pCy1           
tyre_data.pDy1            = 1; % pDy1           
tyre_data.pDy2            = 0; % pDy2           
tyre_data.pDy3            = 0; % pDy3           
tyre_data.pEy1            = 1; % pEy1           
tyre_data.pEy2            = 0; % pEy2           
tyre_data.pEy3            = 0; % pEy3           
tyre_data.pEy4            = 0; % pEy4           
tyre_data.pHy1            = 1; % pHy1           
tyre_data.pHy2            = 0; % pHy2           
tyre_data.pHy3            = 0; % pHy3           
tyre_data.pKy1            = 1; % pKy1           
tyre_data.pKy2            = 1; % pKy2           
tyre_data.pKy3            = 0; % pKy3           
tyre_data.pVy1            = 0; % pVy1           
tyre_data.pVy2            = 0; % pVy2           
tyre_data.pVy3            = 0; % pVy3           
tyre_data.pVy4            = 0; % pVy4           
tyre_data.qBz10           = 0; % qBz10          
tyre_data.qBz2            = 0; % qBz2           
tyre_data.qBz3            = 0; % qBz3           
tyre_data.qBz4            = 0; % qBz4           
tyre_data.qBz5            = 0; % qBz5           
tyre_data.qBz9            = 0; % qBz9           
tyre_data.qCz1            = 0; % qCz1           
tyre_data.qDz1            = 0; % qDz1           
tyre_data.qDz2            = 0; % qDz2           
tyre_data.qDz3            = 0; % qDz3           
tyre_data.qDz4            = 0; % qDz4           
tyre_data.qDz6            = 0; % qDz6           
tyre_data.qDz7            = 0; % qDz7           
tyre_data.qDz8            = 0; % qDz8           
tyre_data.qDz9            = 0; % qDz9           
tyre_data.qEz1            = 0; % qEz1           
tyre_data.qEz2            = 0; % qEz2           
tyre_data.qEz3            = 0; % qEz3           
tyre_data.qEz4            = 0; % qEz4           
tyre_data.qEz5            = 0; % qEz5           
tyre_data.qHz1            = 0; % qHz1           
tyre_data.qHz2            = 0; % qHz2           
tyre_data.qHz3            = 0; % qHz3           
tyre_data.qHz4            = 0; % qHz4           
tyre_data.rBx2            = 0; % rBx2           
tyre_data.rBy1            = 0; % rBy1           
tyre_data.rBy2            = 0; % rBy2           
tyre_data.rBy3            = 0; % rBy3           
tyre_data.rCx1            = 0; % rCx1           
tyre_data.rCy1            = 0; % rCy1           
tyre_data.rHx1            = 0; % rHx1           
tyre_data.rHy1            = 0; % rHy1           
tyre_data.rVy1            = 0; % rVy1           
tyre_data.rVy2            = 0; % rVy2           
tyre_data.rVy3            = 0; % rVy3           
tyre_data.rVy4            = 0; % rVy4           
tyre_data.rVy5            = 0; % rVy5           
tyre_data.rVy6            = 0; % rVy6           
tyre_data.lambda__Cx      = 1; % lambda__Cx     
tyre_data.lambda__Cy      = 1; % lambda__Cy     
tyre_data.lambda__Ex      = 1; % lambda__Ex     
tyre_data.lambda__Ey      = 1; % lambda__Ey     
tyre_data.lambda__Fz0     = 1; % lambda__Fz0    
tyre_data.lambda__Hx      = 1; % lambda__Hx     
tyre_data.lambda__Hy      = 1; % lambda__Hy     
tyre_data.lambda__Kxk     = 1; % lambda__Kxk    
tyre_data.lambda__Kya     = 1; % lambda__Kya    
tyre_data.lambda__Vx      = 1; % lambda__Vx     
tyre_data.lambda__Vy      = 1; % lambda__Vy     
tyre_data.lambda__Vyk     = 1; % lambda__Vyk    
tyre_data.lambda__gamma__y = 1; % lambda__gamma__y
tyre_data.lambda__mu__x   = 1; % lambda__mu__x  
tyre_data.lambda__mu__y   = 1; % lambda__mu__y  
tyre_data.lambda__xa      = 1; % lambda__xa     
tyre_data.lambda__yk      = 1; % lambda__yk     

%% Fitting with load = load_vec(5) = 6500 N and gamma = 0, alpha = 0, Vx = 5
% ------------------
% long slip

% Fit the coeffs {pCx1, pDx1, pEx1, pEx4, pHx1, pKx1, pVx1}

zeros_vec = zeros(size(TDataSub.KAPPA));
ones_vec  = ones(size(TDataSub.KAPPA));
FX0_guess = MF96_FX0_vec(TDataSub.KAPPA, zeros_vec, zeros_vec, tyre_data.Fz0.*ones_vec, tyre_data);

figure()
plot(TDataSub.KAPPA,TDataSub.FX,'o')
hold on
plot(TDataSub.KAPPA,FX0_guess,'x')

% Guess values for parameters to be optimised
% P0 = [1,1,0,0,10,0,0];
% P0 = [1.83872818168040,2.06497792957922,0.833608891995960,0.00216175725649305,1.76119619738288e-06,67.5054374726643,0.000160367103930572];
P0 = [2,2,1,0.1,0.1,70,0.1];

% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised
% 1< pCx1 < 2 
% 0< pEx1 < 1 
lb = [1, 0,   0, 0,  0,  0,  0];
ub = [2, 1e6, 1, 1,1e1,1e2,1e2];

KAPPA_vec = TDataSub.KAPPA;
FX_vec    = TDataSub.FX;

% LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem

% [P_fz_nom,fval,exitflag]
[P_fz_nom,~,~] = fmincon(@(P)resid_pure_Fx(P, FX_vec, KAPPA_vec, 0, fz_vec(5), tyre_data),...
                         P0,[],[],[],[],lb,ub);
                     
% P_fz_nom =
% [1.83872818168040,2.06497792957922,0.833608891995960,0.00216175725649305,1.76119619738288e-06,67.5054374726643,0.000160367103930572];

% Change tyre data with new optimal values                             
tyre_data.pCx1 = P_fz_nom(1) ; % 1
tyre_data.pDx1 = P_fz_nom(2) ;  
tyre_data.pEx1 = P_fz_nom(3) ;
tyre_data.pEx4 = P_fz_nom(4) ;
tyre_data.pHx1 = P_fz_nom(5) ; 
tyre_data.pKx1 = P_fz_nom(6) ;
tyre_data.pVx1 = P_fz_nom(7) ;


FX0_fz_nom_vec = MF96_FX0_vec(TDataSub.KAPPA, zeros_vec, zeros_vec, tyre_data.Fz0.*ones_vec, tyre_data);

figure()
plot(TDataSub.KAPPA,TDataSub.FX,'o')
hold on
plot(TDataSub.KAPPA,FX0_fz_nom_vec,'-')

% Calculate the residuals with the optimal solution found above
res_Fx0_nom  = resid_pure_Fx(P_fz_nom, FX_vec, KAPPA_vec, 0, fz_vec(5), tyre_data);

% R-squared is 
% 1-SSE/SST
% SSE/SST = res_Fx0_nom

% SSE is the sum of squared error,  SST is the sum of squared total
fprintf('R-squared = %6.3f\n',1-res_Fx0_nom);

[kappa__x, Bx, Cx, Dx, Ex, SVx] = MF96_FX0_coeffs(0, 0, 0, tyre_data.Fz0, tyre_data);

fprintf('Bx      = %6.3f\n',Bx);
fprintf('Cx      = %6.3f\n',Cx);
fprintf('Ex      = %6.3f\n',Ex);
fprintf('SVx     = %6.3f\n',SVx);
fprintf('kappa_x = %6.3f\n',kappa__x);
fprintf('Kx      = %6.3f\n',Bx*Cx*Dx/tyre_data.Fz0);

%% Fitting with Fz = var and gamma = 0, alpha = 0, VX = 5
% ------------------
% extract data with variable load
[TDataSubFzAll, ~] = intersectTableData( ALPHA_00, GAMMA_00, VX_05 );

% Fit the coeffs { pDx2, pEx2, pEx3, pKx2,pKx3, pHx2, pVx2}

% Guess values for parameters to be optimised
P0 = [0,0,0,0,0,0,0]; 
lb = [];
ub = [];

KAPPA_vec = TDataSubFzAll.KAPPA;
FX_vec    = TDataSubFzAll.FX;
FZ_vec    = TDataSubFzAll.FZ;
zeros_vec = zeros(size(FX_vec));
ones_vec  = ones(size(FX_vec));

resid_pure_Fx_varFz(P0, FX_vec, KAPPA_vec, 0, FZ_vec, tyre_data)

% LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 
[P_fz_nom,~,~] = fmincon(@(P)resid_pure_Fx_varFz(P, FX_vec, KAPPA_vec, 0, FZ_vec, tyre_data),...
                               P0,[],[],[],[],lb,ub);

% Change tyre data with new optimal values                             
tyre_data.pDx2 = P_fz_nom(1) ; % 1
tyre_data.pEx2 = P_fz_nom(2) ;
tyre_data.pEx3 = P_fz_nom(3) ;
tyre_data.pHx2 = P_fz_nom(4) ; 
tyre_data.pKx2 = P_fz_nom(5) ;
tyre_data.pKx3 = P_fz_nom(6) ;
tyre_data.pVx2 = P_fz_nom(7) ;

FX0_fz_var_vec = MF96_FX0_vec(TDataSubFzAll.KAPPA, zeros_vec, zeros_vec, FZ_vec, tyre_data);

figure()
plot(TDataSubFzAll.KAPPA,TDataSubFzAll.FX,'o')
hold on
plot(TDataSubFzAll.KAPPA,FX0_fz_var_vec,'-')

% Calculate the residuals with the optimal solution found above
res_Fx0_fz_var  = resid_pure_Fx_varFz(P_fz_nom, FX_vec, KAPPA_vec, 0, FZ_vec, tyre_data);

% R-squared is 
% 1-SSE/SST
% SSE/SST = res_Fx0_nom

% SSE is the sum of squared error,  SST is the sum of squared total
fprintf('R-squared = %6.3f\n',1-res_Fx0_fz_var);

[kappa__x, Bx, Cx, Dx, Ex, SVx] = MF96_FX0_coeffs(0, 0, 0, fz_vec(5).*ones(size(tyre_data.Fz0)), tyre_data);
% 
fprintf('Bx      = %6.3f\n',Bx);
fprintf('Cx      = %6.3f\n',Cx);
fprintf('mux      = %6.3f\n',Dx/fz_vec(3));
fprintf('Ex      = %6.3f\n',Ex);
fprintf('SVx     = %6.3f\n',SVx);
fprintf('kappa_x = %6.3f\n',kappa__x);
fprintf('Kx      = %6.3f\n',Bx*Cx*Dx/tyre_data.Fz0);

%% Fitting with gamma = var and alpha = 0, VX = 5, load = 6500 N
% ------------------
% extract data with variable load
[TDataSubFzAll, ~] = intersectTableData( ALPHA_00, FZ_5, VX_05 );

% Fit the coeff { pDx3 }

% Guess values for parameters to be optimised
P0 = [0]; 
lb = [];
ub = [];

KAPPA_vec = TDataSubFzAll.KAPPA;
FX_vec    = TDataSubFzAll.FX;
GAMMA_vec    = TDataSubFzAll.GAMMA;
zeros_vec = zeros(size(FX_vec));
ones_vec  = ones(size(FX_vec));

% LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 
[P_fz_nom,fval,exitflag] = fmincon(@(P)resid_pure_Fx_varGamma(P, FX_vec, KAPPA_vec, GAMMA_vec, fz_vec(5), tyre_data),...
                               P0,[],[],[],[],lb,ub);

% Change tyre data with new optimal values                             
tyre_data.pDx3 = P_fz_nom(1) ;

FX0_gamma_var_vec = MF96_FX0_vec(TDataSubFzAll.KAPPA, zeros_vec, GAMMA_vec, ones_vec*fz_vec(5), tyre_data);

figure()
plot(TDataSubFzAll.KAPPA,TDataSubFzAll.FX,'o')
hold on
plot(TDataSubFzAll.KAPPA,FX0_gamma_var_vec,'-')

% Calculate the residuals with the optimal solution found above
res_Fx0_gamma_var = resid_pure_Fx_varGamma(P_fz_nom, FX_vec, KAPPA_vec, GAMMA_vec, fz_vec(5), tyre_data);

% R-squared is 
% 1-SSE/SST
% SSE/SST = res_Fx0_nom

% SSE is the sum of squared error,  SST is the sum of squared total
fprintf('R-squared = %6.3f\n',1-res_Fx0_gamma_var);

[kappa__x, Bx, Cx, Dx, Ex, SVx] = MF96_FX0_coeffs(0, 0, 0, fz_vec(5).*ones(size(tyre_data.Fz0)), tyre_data);
% 
fprintf('Bx      = %6.3f\n',Bx);
fprintf('Cx      = %6.3f\n',Cx);
fprintf('mux      = %6.3f\n',Dx/fz_vec(3));
fprintf('Ex      = %6.3f\n',Ex);
fprintf('SVx     = %6.3f\n',SVx);
fprintf('kappa_x = %6.3f\n',kappa__x);
fprintf('Kx      = %6.3f\n',Bx*Cx*Dx/tyre_data.Fz0);

%% Lateral test

load 'dataset/front_lateral_dataset.mat'; 
table = front_lateral_dataset;

table = table(table.FZ ~= 0,:);

clear alpha_vec kappa_vec ...
      ALPHA_00 ALPHA_25 ALPHA_50 ALPHA_75 ALPHA_100 ...
      
alpha_vec = linspace( -15.0, +15.0, 200 )*pi/180;
kappa_vec = [ 0.00, 0.05, 0.10, 0.15, 0.20 ];

%% Selecting specific data

% Extract points at constant slip ratio

KAPPA_tol = 0.005;
KAPPA_00  = table( kappa_vec(1)-KAPPA_tol < table.KAPPA & table.KAPPA < kappa_vec(1)+KAPPA_tol, : );
KAPPA_05  = table( kappa_vec(2)-KAPPA_tol < table.KAPPA & table.KAPPA < kappa_vec(2)+KAPPA_tol, : );
KAPPA_10  = table( kappa_vec(3)-KAPPA_tol < table.KAPPA & table.KAPPA < kappa_vec(3)+KAPPA_tol, : );
KAPPA_15  = table( kappa_vec(4)-KAPPA_tol < table.KAPPA & table.KAPPA < kappa_vec(4)+KAPPA_tol, : );
KAPPA_20 =  table( kappa_vec(5)-KAPPA_tol < table.KAPPA & table.KAPPA < kappa_vec(5)+KAPPA_tol, : );

% Extract points at constant vertical load
FZ_tol   = 25.0;
FZ_1  = table(  fz_vec(1)-FZ_tol < table.FZ & table.FZ < fz_vec(1)+FZ_tol, : );
FZ_2  = table(  fz_vec(2)-FZ_tol < table.FZ & table.FZ < fz_vec(2)+FZ_tol, : );
FZ_3  = table(  fz_vec(3)-FZ_tol < table.FZ & table.FZ < fz_vec(3)+FZ_tol, : );
FZ_4  = table(  fz_vec(4)-FZ_tol < table.FZ & table.FZ < fz_vec(4)+FZ_tol, : );
FZ_5  = table(  fz_vec(5)-FZ_tol < table.FZ & table.FZ < fz_vec(5)+FZ_tol, : );
FZ_6  = table(  fz_vec(6)-FZ_tol < table.FZ & table.FZ < fz_vec(6)+FZ_tol, : );
FZ_7  = table(  fz_vec(7)-FZ_tol < table.FZ & table.FZ < fz_vec(7)+FZ_tol, : );

% Extract points at constant inclination angle
GAMMA_tol = 0.1*pi/180;
GAMMA_00  = table( gamma_vec(1)*pi/180-GAMMA_tol < table.GAMMA & table.GAMMA < gamma_vec(1)*pi/180+GAMMA_tol, : );
GAMMA_25  = table( gamma_vec(2)*pi/180-GAMMA_tol < table.GAMMA & table.GAMMA < gamma_vec(2)*pi/180+GAMMA_tol, : );
GAMMA_50  = table( gamma_vec(3)*pi/180-GAMMA_tol < table.GAMMA & table.GAMMA < gamma_vec(3)*pi/180+GAMMA_tol, : );

% Extract points at constant road speed
VX_tol = 0.5;
VX_05  = table( Vx_vec(1)-VX_tol < table.VX & table.VX < Vx_vec(1)+VX_tol, : );
VX_20  = table( Vx_vec(2)-VX_tol < table.VX & table.VX < Vx_vec(2)+VX_tol, : );
VX_45  = table( Vx_vec(3)-VX_tol < table.VX & table.VX < Vx_vec(3)+VX_tol, : );
VX_90  = table( Vx_vec(4)-VX_tol < table.VX & table.VX < Vx_vec(4)+VX_tol, : );

%% Fitting with load = load_vec(5) = 6500 N and gamma = 0, kappa = 0, Vx = 5
% ------------------
% side slip

[TDataSub, ~] = intersectTableData( KAPPA_00, GAMMA_00, VX_05, FZ_5 );
plot_selected_data(TDataSub)

% Fit the coeffs {pCx1, pDx1, pEx1, pEx4, pHx1, pKx1, pVx1}
% Fit the coeffs {pCy1, pDy1, pEy1, pHy1, pKy1, pKy2, pVy1}

zeros_vec = zeros(size(TDataSub.ALPHA));
ones_vec  = ones(size(TDataSub.ALPHA));
FY0_guess = MF96_FY0_vec(zeros_vec, TDataSub.ALPHA, zeros_vec, tyre_data.Fz0.*ones_vec, tyre_data);

figure()
plot(TDataSub.ALPHA,(-1)*TDataSub.FY,'o')
hold on
plot(TDataSub.ALPHA,FY0_guess,'x')

P0 = [2,2,1,0.1,70,0.1,0.05];

lb = [1, 0.1, 0, 0.001, 0.1,  0,  0.001];
ub = [2, 1e6, 1, 1, 1e2, 1e2, 0.01];

ALPHA_vec = TDataSub.ALPHA;
FY_vec    = (-1)*TDataSub.FY;

[P_fz_nom,~,~] = fmincon(@(P)resid_pure_Fy(P, FY_vec, ALPHA_vec, 0, fz_vec(5), tyre_data),...
                         P0,[],[],[],[],lb,ub);
                     
% Change tyre data with new optimal values                             
tyre_data.pCy1 = P_fz_nom(1) ;
tyre_data.pDy1 = P_fz_nom(2) ;  
tyre_data.pEy1 = P_fz_nom(3) ;
tyre_data.pHy1 = P_fz_nom(4) ;
tyre_data.pKy1 = P_fz_nom(5) ; 
tyre_data.pKy2 = P_fz_nom(6) ;
tyre_data.pVy1 = P_fz_nom(7) ;

FY0_fz_nom_vec = MF96_FY0_vec(zeros_vec, TDataSub.ALPHA, zeros_vec, tyre_data.Fz0.*ones_vec, tyre_data);

figure()
plot(TDataSub.ALPHA,(-1)*TDataSub.FY,'o')
hold on
plot(TDataSub.ALPHA,FY0_fz_nom_vec,'-')

% Calculate the residuals with the optimal solution found above
res_Fy0_nom  = resid_pure_Fy(P_fz_nom, FY_vec, ALPHA_vec, 0, fz_vec(5), tyre_data);

fprintf('R-squared = %6.3f\n',1-res_Fy0_nom);

[alpha__y, By, Cy, Dy, Ey, SVy] = MF96_FY0_coeffs(0, 0, 0, tyre_data.Fz0, tyre_data);

fprintf('By      = %6.3f\n',By);
fprintf('Cy      = %6.3f\n',Cy);
fprintf('Ey      = %6.3f\n',Ey);
fprintf('SVy     = %6.3f\n',SVy);
fprintf('alpha_y = %6.3f\n',alpha__y);
fprintf('Ky      = %6.3f\n',By*Cy*Dy/tyre_data.Fz0);

%% Fitting with Fz = var and gamma = 0, kappa = 0, VX = 5
% ------------------
% extract data with variable load
[TDataSubFzAll, ~] = intersectTableData( KAPPA_00, GAMMA_00, VX_05 );

% Fit the coeffs { pDy2, pEy2, pHy2, pVy2 }

% Guess values for parameters to be optimised
P0 = [0,0,0,0];
lb = [];
ub = [];

ALPHA_vec = TDataSubFzAll.ALPHA;
FY_vec    = (-1)*TDataSubFzAll.FY;
FZ_vec    = TDataSubFzAll.FZ;
zeros_vec = zeros(size(FY_vec));
ones_vec  = ones(size(FY_vec));

[P_fz_nom,~,~] = fmincon(@(P)resid_pure_Fy_varFz(P, FY_vec, ALPHA_vec, 0, FZ_vec, tyre_data),...
                               P0,[],[],[],[],lb,ub);

% Change tyre data with new optimal values                             
tyre_data.pDy2 = P_fz_nom(1) ;
tyre_data.pEy2 = P_fz_nom(2) ;
tyre_data.pHy2 = P_fz_nom(3) ; 
tyre_data.pVy2 = P_fz_nom(4) ;

FY0_fz_var_vec = MF96_FY0_vec(zeros_vec, TDataSubFzAll.ALPHA, zeros_vec, FZ_vec, tyre_data);

figure()
plot(TDataSubFzAll.ALPHA,(-1)*TDataSubFzAll.FY,'o')
hold on
plot(TDataSubFzAll.ALPHA,FY0_fz_var_vec,'-')

% Calculate the residuals with the optimal solution found above
res_Fy0_fz_var  = resid_pure_Fy_varFz(P_fz_nom, FY_vec, ALPHA_vec, 0, FZ_vec, tyre_data);
fprintf('R-squared = %6.3f\n',1-res_Fy0_fz_var);

[alpha__y, By, Cy, Dy, Ey, SVy] = MF96_FY0_coeffs(0, 0, 0, fz_vec(5).*ones(size(tyre_data.Fz0)), tyre_data);
% 
fprintf('By      = %6.3f\n',By);
fprintf('Cy      = %6.3f\n',Cy);
fprintf('muy     = %6.3f\n',Dy/fz_vec(5));
fprintf('Ey      = %6.3f\n',Ey);
fprintf('SVy     = %6.3f\n',SVy);
fprintf('alpha_y = %6.3f\n',alpha__y);
fprintf('Ky      = %6.3f\n',By*Cy*Dy/tyre_data.Fz0);

%% Fitting with gamma = var and kappa = 0, VX = 5, load = 6500 N
% ------------------
% extract data with variable load
[TDataSubFzAll, ~] = intersectTableData( KAPPA_00, FZ_5, VX_05 );

% Fit the coeff { pDy3, pEy3, pEy4, pHy3, pKy3, pVy3 }

% Guess values for parameters to be optimised
P0 = [0,0,0,0,0,0]; 
lb = [];
ub = [];

ALPHA_vec = TDataSubFzAll.ALPHA;
FY_vec    = (-1)*TDataSubFzAll.FY;
GAMMA_vec = TDataSubFzAll.GAMMA;
zeros_vec = zeros(size(FY_vec));
ones_vec  = ones(size(FY_vec));

[P_fz_nom,fval,exitflag] = fmincon(@(P)resid_pure_Fy_varGamma(P, FY_vec, ALPHA_vec, GAMMA_vec, fz_vec(5), tyre_data),...
                               P0,[],[],[],[],lb,ub);

% Change tyre data with new optimal values
tyre_data.pDy3 = P_fz_nom(1);
tyre_data.pEy3 = P_fz_nom(2);
tyre_data.pEy4 = P_fz_nom(3);
tyre_data.pHy3 = P_fz_nom(4);
tyre_data.pKy3 = P_fz_nom(5);
tyre_data.pVy3 = P_fz_nom(6);

FY0_gamma_var_vec = MF96_FY0_vec(zeros_vec, TDataSubFzAll.ALPHA, GAMMA_vec, ones_vec*fz_vec(5), tyre_data);

figure()
plot(TDataSubFzAll.ALPHA,(-1)*TDataSubFzAll.FY,'o')
hold on
plot(TDataSubFzAll.ALPHA,FY0_gamma_var_vec,'-')

% Calculate the residuals with the optimal solution found above
res_Fy0_gamma_var = resid_pure_Fy_varGamma(P_fz_nom, FY_vec, ALPHA_vec, GAMMA_vec, fz_vec(5), tyre_data);

% R-squared is 
% 1-SSE/SST
% SSE/SST = res_Fx0_nom

% SSE is the sum of squared error,  SST is the sum of squared total
fprintf('R-squared = %6.3f\n',1-res_Fy0_gamma_var);

[alpha__y, By, Cy, Dy, Ey, SVy] = MF96_FY0_coeffs(0, 0, 0, fz_vec(5).*ones(size(tyre_data.Fz0)), tyre_data);
% 
fprintf('By      = %6.3f\n',By);
fprintf('Cy      = %6.3f\n',Cy);
fprintf('muy     = %6.3f\n',Dy/fz_vec(5));
fprintf('Ey      = %6.3f\n',Ey);
fprintf('SVy     = %6.3f\n',SVy);
fprintf('alpha_y = %6.3f\n',alpha__y);
fprintf('Ky      = %6.3f\n',By*Cy*Dy/tyre_data.Fz0);

%% Fitting with load = var and kappa = 0, VX = 5, gamma = 2.5
% ------------------
% extract data with variable load
[TDataSubFzAll, ~] = intersectTableData( KAPPA_00, GAMMA_25, VX_05 );

% Fit the coeff { pVy4 }

% Guess values for parameters to be optimised
P0 = [0]; 
lb = [];
ub = [];

ALPHA_vec = TDataSubFzAll.ALPHA;
FY_vec    = (-1)*TDataSubFzAll.FY;
FZ_vec    = TDataSubFzAll.FZ;
zeros_vec = zeros(size(FY_vec));
ones_vec  = ones(size(FY_vec));

[P_fz_nom,fval,exitflag] = fmincon(@(P)resid_pure_Fy_varFz_gamma(P, FY_vec, ALPHA_vec, gamma_vec(2), FZ_vec, tyre_data),...
                               P0,[],[],[],[],lb,ub);

% Change tyre data with new optimal values                             
tyre_data.pVy4 = P_fz_nom(1) ;

FY0_fz_var_vec_gamma = MF96_FY0_vec(zeros_vec, TDataSubFzAll.ALPHA, ones_vec*gamma_vec(2), FZ_vec, tyre_data);

figure()
plot(TDataSubFzAll.ALPHA,(-1)*TDataSubFzAll.FY,'o')
hold on
plot(TDataSubFzAll.ALPHA,FY0_fz_var_vec_gamma,'-')

% Calculate the residuals with the optimal solution found above
res_Fy0_fz_var_gamma = resid_pure_Fy_varFz_gamma(P_fz_nom, FY_vec, ALPHA_vec, gamma_vec(2), FZ_vec, tyre_data);

fprintf('R-squared = %6.3f\n',1-res_Fy0_fz_var_gamma);

[alpha__y, By, Cy, Dy, Ey, SVy] = MF96_FY0_coeffs(0, 0, gamma_vec(2), fz_vec(5).*ones(size(tyre_data.Fz0)), tyre_data);
% 
fprintf('By      = %6.3f\n',By);
fprintf('Cy      = %6.3f\n',Cy);
fprintf('muy     = %6.3f\n',Dy/fz_vec(5));
fprintf('Ey      = %6.3f\n',Ey);
fprintf('SVy     = %6.3f\n',SVy);
fprintf('alpha_y = %6.3f\n',alpha__y);
fprintf('Ky      = %6.3f\n',By*Cy*Dy/tyre_data.Fz0);