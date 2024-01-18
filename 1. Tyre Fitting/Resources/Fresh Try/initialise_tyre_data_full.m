function ty_data = initialise_tyre_data_full(R0, Fz0)


% Tyre structure data initialization

ty_data.FZ0             = Fz0; % Fz0  % Normal load
ty_data.R0              = R0; % R0  % nominal radius
ty_data.pCx1            = 1.44094240439010;
ty_data.pDx1            = 2.35256495098336;
ty_data.pDx2            = -0.810191533417750;
ty_data.pDx3            = 10.8761875333065;
ty_data.pEx1            = -0.00184214062456847;
ty_data.pEx2            = 0.00296190301240869;
ty_data.pEx3            = 0.00361681692585595;
ty_data.pEx4            = 73.9065661693525;
ty_data.pHx1            = 0.000851251478635586;
ty_data.pHx2            = 0.000906494654404505;
ty_data.pKx1            = 54.0397564627930;
ty_data.pKx2            = -0.0101430342064303;
ty_data.pKx3            = 0.548143673871934;
ty_data.pVx1            = -0.0868210803989277;
ty_data.pVx2            = -0.0261373983233012;

ty_data.Fz01            = 0; % Fz01
ty_data.pCy1            = 1.54027309068454;
ty_data.pDy1            = 2.43404156744141;
ty_data.pDy2            = -0.274390567365254;
ty_data.pDy3            = 5.16349361588200;
ty_data.pEy1            = 0.346228681636638;
ty_data.pEy2            = -0.278086233494965;
ty_data.pEy3            = 0.883205306970767;
ty_data.pEy4            = -4.20062483943437;
ty_data.pHy1            = -0.00489082598162627;
ty_data.pHy2            = -0.00117998967899828;
ty_data.pHy3            = -0.0313517996563904;
ty_data.pKy1            = -29.6801903969105;
ty_data.pKy2            = -1.01712562049631;
ty_data.pKy3            = 1.34159743861844;
ty_data.pVy1            = 0.0476749940408608;
ty_data.pVy2            = -0.0401680086940341;
ty_data.pVy3            = -2.92370451981899;
ty_data.pVy4            = -2.81666715709586;

ty_data.qBz1            = 7.35347522395674;
ty_data.qBz2            = -4.71595003290210;
ty_data.qBz3            = -6.17019496981862;
ty_data.qBz4            = -0.281057344706274;
ty_data.qBz5            = -0.281558360232466;
ty_data.qBz9            = 0;
ty_data.qBz10           = 0;
ty_data.qCz1            = 1.67097884314219;
ty_data.qDz1            = 0.197223329858040;
ty_data.qDz2            = -0.0549111816154103;
ty_data.qDz3            = -3.31843325657878;
ty_data.qDz4            = 42.5363302221465;
ty_data.qDz6            = 0.00307532725022005;
ty_data.qDz7            = 0.0171847774709159;
ty_data.qDz8            = -1.82587609907207;
ty_data.qDz9            = -0.796929452914326;
ty_data.qEz1            = 0.515207229031359;
ty_data.qEz2            = -1.69617162475273;
ty_data.qEz3            = -3.39063846665826;
ty_data.qEz4            = 0.164755141111871;
ty_data.qEz5            = -9.68793770573582;
ty_data.qHz1            = -0.0106709285011952;
ty_data.qHz2            = -0.00767395160043653;
ty_data.qHz3            = 0.720767969228803;
ty_data.qHz4            = 0.0934620620316319;

ty_data.rBx1            = 23.3086620367818;
ty_data.rBx2            = 19.4145696771092; % rBx2
ty_data.rCx1            = 0.926013692178950;
ty_data.rEx1            = 0;
ty_data.rEx2            = 0;
ty_data.rHx1            = -0.00125774506013463;

ty_data.rBy1            = 14.1655485377604;
ty_data.rBy2            = 13.2919779719267;
ty_data.rBy3            = -0.499338655070305;
ty_data.rBy4            = 0;
ty_data.rCy1            = 0.979136167579569; % rCy1
ty_data.rEy1            = 0;
ty_data.rEy2            = 0;
ty_data.rHy1            = 0.0298165069376798; % rHy1
ty_data.rVy1            = -0.237740939096051;
ty_data.rVy2            = 0;
ty_data.rVy3            = 0;
ty_data.rVy4            = 3.76811371158387;
ty_data.rVy5            = -0.0906119916069004;
ty_data.rVy6            = 28.3775771964863;




% scaling factor
ty_data.LCX             = 1; % LCX
ty_data.LCY             = 1; % LCY
ty_data.LEX             = 1; % LEX
ty_data.LEY             = 1; % LEY
ty_data.LFZ0            = 1; % LFZ0
ty_data.LGAMMAY         = 1; % LGAMMAY
ty_data.LHX             = 1; % LHX
ty_data.LHY             = 1; % LHY
ty_data.LKA             = 1; % LKA
ty_data.LKXK            = 1; % LKXK
ty_data.LMUX            = 1; % LMUX
ty_data.LMUY            = 1; % LMUY
ty_data.LVX             = 1; % LVX
ty_data.LVY             = 1; % LVY
ty_data.LVYK            = 1; % LVYK
ty_data.LXA             = 1; % LXA
ty_data.LKY             = 1; % LYK
ty_data.LT              = 1; % LT
ty_data.LMR             = 1; % LMR

end