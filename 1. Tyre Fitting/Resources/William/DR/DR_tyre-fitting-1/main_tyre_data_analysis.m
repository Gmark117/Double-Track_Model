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

kappa_vec = linspace( -1.0, +1.0, 200 ).^3;                  % [-]
alpha_vec = [ 0.0, 2.5, 5.0, 7.5, 10.0 ]*pi/180;             % [rad]
load_vec  = [ 100, 1000, 2500, 4500, 6500, 8500, 10000 ];    % [N]
gamma_vec = [ 0.0, 2.5, 5.0 ]*pi/180;                        % [rad]
Vx_vec    = [ 5.0, 20.0, 45.0, 90.0 ];                       % [m/s]

%% Plot raw data

plot_raw_data

%% Select some specific data

% Extract points at constant slip ratio
if length(kappa_vec) == 5
    KAPPA_tol = 0.005;
    KAPPA_00  = table( kappa_vec(1)-KAPPA_tol < table.KAPPA & table.KAPPA < kappa_vec(1)+KAPPA_tol, : );
    KAPPA_05  = table( kappa_vec(2)-KAPPA_tol < table.KAPPA & table.KAPPA < kappa_vec(2)+KAPPA_tol, : );
    KAPPA_10  = table( kappa_vec(3)-KAPPA_tol < table.KAPPA & table.KAPPA < kappa_vec(3)+KAPPA_tol, : );
    KAPPA_15  = table( kappa_vec(4)-KAPPA_tol < table.KAPPA & table.KAPPA < kappa_vec(4)+KAPPA_tol, : );
    KAPPA_20 = table( kappa_vec(5)-KAPPA_tol < table.KAPPA & table.KAPPA < kappa_vec(5)+KAPPA_tol, : );
end

% Extract points at constant side slip
if length(alpha_vec) == 5
    ALPHA_tol = 0.005;
    ALPHA_00  = table( alpha_vec(1)-ALPHA_tol < table.ALPHA & table.ALPHA < alpha_vec(1)+ALPHA_tol, : );
    ALPHA_25  = table( alpha_vec(2)-ALPHA_tol < table.ALPHA & table.ALPHA < alpha_vec(2)+ALPHA_tol, : );
    ALPHA_50  = table( alpha_vec(3)-ALPHA_tol < table.ALPHA & table.ALPHA < alpha_vec(3)+ALPHA_tol, : );
    ALPHA_75  = table( alpha_vec(4)-ALPHA_tol < table.ALPHA & table.ALPHA < alpha_vec(4)+ALPHA_tol, : );
    ALPHA_100 = table( alpha_vec(5)-ALPHA_tol < table.ALPHA & table.ALPHA < alpha_vec(5)+ALPHA_tol, : );
end

% Extract points at constant vertical load
LOAD_tol   = 25.0;
LOAD_1  = table(  load_vec(1)-LOAD_tol < table.FZ & table.FZ < load_vec(1)+LOAD_tol, : );
LOAD_2  = table(  load_vec(2)-LOAD_tol < table.FZ & table.FZ < load_vec(2)+LOAD_tol, : );
LOAD_3  = table(  load_vec(3)-LOAD_tol < table.FZ & table.FZ < load_vec(3)+LOAD_tol, : );
LOAD_4  = table(  load_vec(4)-LOAD_tol < table.FZ & table.FZ < load_vec(4)+LOAD_tol, : );
LOAD_5  = table(  load_vec(5)-LOAD_tol < table.FZ & table.FZ < load_vec(5)+LOAD_tol, : );
LOAD_6  = table(  load_vec(6)-LOAD_tol < table.FZ & table.FZ < load_vec(6)+LOAD_tol, : );
LOAD_7  = table(  load_vec(7)-LOAD_tol < table.FZ & table.FZ < load_vec(7)+LOAD_tol, : );

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
[TDataSub, ~] = intersectTableData( ALPHA_00, GAMMA_00, VX_05, LOAD_5 );

%[tableSelectedData, ~] = intersectTableData( KAPPA_00, GAMMA_00, VX_10, FZ_6500 );



%%
%plot_selected_data
plot_selected_data(TDataSub)

%% FITTING 

% Tyre structure data initialization
tyre_data.Fz0             = load_vec(5); % Fz0            
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
tyre_data.pCy1            = 0; % pCy1           
tyre_data.pDy1            = 0; % pDy1           
tyre_data.pDy2            = 0; % pDy2           
tyre_data.pDy3            = 0; % pDy3           
tyre_data.pEy1            = 0; % pEy1           
tyre_data.pEy2            = 0; % pEy2           
tyre_data.pEy3            = 0; % pEy3           
tyre_data.pEy4            = 0; % pEy4           
tyre_data.pHy1            = 0; % pHy1           
tyre_data.pHy2            = 0; % pHy2           
tyre_data.pHy3            = 0; % pHy3           
tyre_data.pKy1            = 0; % pKy1           
tyre_data.pKy2            = 0; % pKy2           
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

% Fit the coeffs {pCx1, pDx1, pEx1, pEx4, pKx1, pHx1, pVx1}

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
[P_fz_nom,~,~] = fmincon(@(P)resid_pure_Fx(P, FX_vec, KAPPA_vec, 0, load_vec(5), tyre_data),...
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
res_Fx0_nom  = resid_pure_Fx(P_fz_nom, FX_vec, KAPPA_vec, 0, load_vec(5), tyre_data);

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

%% Fitting with load = var and gamma = 0, alpha = 0, VX = 5
% ------------------
% extract data with variable load
[TDataSubFzAll, ~] = intersectTableData( ALPHA_00, GAMMA_00, VX_05 );

% Fit the coeffs { pDx2, pDx3, pEx2, pEx3, pKx2,pKx3, pHx2, pVx2}

% zeros_vec = zeros(size(TDataSub.KAPPA));
% ones_vec  = ones(size(TDataSub.KAPPA));
% FX0_guess = MF96_FX0_vec(TDataSub.KAPPA,zeros_vec , zeros_vec, tyre_data.Fz0.*ones_vec, tyre_data);
% 
% figure()
% plot(TDataSub.KAPPA,TDataSub.FX,'o')
% hold on
% plot(TDataSub.KAPPA,FX0_guess,'x')

%
% Guess
% X0 = [1.5,2.22,0,0,0.34,0,0,-0.05543,-0.00087,0,46.06,0,0,0.00295,0];

% Guess values for parameters to be optimised
P0 = [0,0,0,0,0,0,0,0];

% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised
% 1< pCx1 < 2 
% 0< pEx1 < 1 
%lb = [0, 0,  0, 0,  0,  0,  0];
%ub = [2, 1e6,1, 1,1e1,1e2,1e2];
lb = [];
ub = [];

KAPPA_vec = TDataSubFzAll.KAPPA;
FX_vec    = TDataSubFzAll.FX;
FZ_vec    = TDataSubFzAll.FZ;
zeros_vec = zeros(size(FX_vec));
ones_vec  = ones(size(FX_vec));

% LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 
[P_fz_nom,fval,exitflag] = fmincon(@(P)resid_pure_Fx_varFz(P, FX_vec, KAPPA_vec, 0, FZ_vec, tyre_data),...
                               P0,[],[],[],[],lb,ub);

% Change tyre data with new optimal values                             
tyre_data.pDx2 = P_fz_nom(1) ; % 1
tyre_data.pDx3 = P_fz_nom(2) ;  
tyre_data.pEx2 = P_fz_nom(3) ;
tyre_data.pEx3 = P_fz_nom(4) ;
tyre_data.pHx2 = P_fz_nom(5) ; 
tyre_data.pKx2 = P_fz_nom(6) ;
tyre_data.pKx3 = P_fz_nom(7) ;
tyre_data.pVx2 = P_fz_nom(8) ;

FX0_fz_nom_vec = MF96_FX0_vec(TDataSubFzAll.KAPPA, zeros_vec, zeros_vec, FZ_vec, tyre_data);

figure()
plot(TDataSubFzAll.KAPPA,TDataSubFzAll.FX,'o')
hold on
plot(TDataSubFzAll.KAPPA,FX0_fz_nom_vec,'-')

% Calculate the residuals with the optimal solution found above
res_Fx0_nom  = resid_pure_Fx_varFz(P_fz_nom, FX_vec, KAPPA_vec, 0, FZ_vec, tyre_data);

% R-squared is 
% 1-SSE/SST
% SSE/SST = res_Fx0_nom

% SSE is the sum of squared error,  SST is the sum of squared total
fprintf('R-squared = %6.3f\n',1-res_Fx0_nom);

[kappa__x, Bx, Cx, Dx, Ex, SVx] = MF96_FX0_coeffs(0, 0, 0, load_vec(3).*ones(size(tyre_data.Fz0)), tyre_data);
% 
fprintf('Bx      = %6.3f\n',Bx);
fprintf('Cx      = %6.3f\n',Cx);
fprintf('mux      = %6.3f\n',Dx/load_vec(3));
fprintf('Ex      = %6.3f\n',Ex);
fprintf('SVx     = %6.3f\n',SVx);
fprintf('kappa_x = %6.3f\n',kappa__x);
fprintf('Kx      = %6.3f\n',Bx*Cx*Dx/tyre_data.Fz0);

%% Fitting for change of Fz



