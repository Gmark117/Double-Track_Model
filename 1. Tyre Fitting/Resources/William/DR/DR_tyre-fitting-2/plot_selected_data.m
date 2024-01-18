function plot_selected_data(tab)

  color     = '#0072BD';
  linewidth = 1.1;
  linestyle = 'none';
  marker    = '*';

  % Plot selected data
  figure;
  tiledlayout(3,2);

  % FX ---------
  nexttile
  hold on;
  grid on;
  grid minor;
  xlabel('$\kappa$ (-)');
  ylabel('$F_x$ (N)');
  plot(tab.KAPPA, tab.FX, 'Color', color, 'LineStyle', linestyle, ...
       'LineWidth', linewidth, 'Marker', marker);
  hold on;

  nexttile
  hold on;
  grid on;
  grid minor;
  xlabel('$\alpha$ (deg)');
  ylabel('$F_x$ (N)');
  plot(tab.ALPHA, tab.FX, 'Color', color, 'LineStyle', linestyle, ...
       'LineWidth', linewidth, 'Marker', marker);
  hold on;

  % FY ---------
  nexttile
  hold on;
  grid on;
  grid minor;
  xlabel('$\kappa$ (-)');
  ylabel('$F_y$ (N)');
  plot(tab.KAPPA, tab.FY, 'Color', color, 'LineStyle', linestyle, ...
       'LineWidth', linewidth, 'Marker', marker);
  hold on;

  nexttile
  hold on;
  grid on;
  grid minor;
  xlabel('$\alpha$ (deg)');
  ylabel('$F_y$ (N)');
  plot(tab.ALPHA, tab.FY, 'Color', color, 'LineStyle', linestyle, ...
       'LineWidth', linewidth, 'Marker', marker);
  hold on;

  % MZ ---------
  nexttile
  hold on;
  grid on;
  grid minor;
  xlabel('$\kappa$ (-)');
  ylabel('$M_z$ (Nm)');
  plot(tab.KAPPA, tab.MZ, 'Color', color, 'LineStyle', linestyle, ...
       'LineWidth', linewidth, 'Marker', marker);
  hold on;

  nexttile
  hold on;
  grid on;
  grid minor;
  xlabel('$\alpha$ (deg)');
  ylabel('$M_z$ (Nm)');
  plot(tab.ALPHA, tab.MZ, 'Color', color, 'LineStyle', linestyle, ...
       'LineWidth', linewidth, 'Marker', marker);
  hold on;

  hold off;
end