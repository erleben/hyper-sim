% Copyright 2011, Kenny Erleben
clear all;
close all;
clc;

meshes = {
  'bar_T1098_V300.mat',...
  'bar_T1197_V325.mat',...
  'bar_T2293_V576.mat',...
  'bar_T7737_V1701.mat',...
  'bar_T8150_V1782.mat',...
  'bar_T11741_V2500.mat',...
  'bar_T20423_V4176.mat',...
  'bar_T27977_V5577.mat'...
  };

for m=1:length(meshes)
  
  load( meshes{m} );
  [ Qrr, Qrl, Qtheta, Qvl ] = compute_quality_measures(T, X, Y, Z);
  fh = figure(1);
  clf;
  set(gca,'FontSize',18);
  hist(Qrr, 40);
  ylabel('#Elements','FontSize',18);
  xlabel('$Q_{rr} = 3 \frac{R_{in}}{R_{out}}$','FontSize',18,'interpreter','latex');
  grid on
  filename = strcat( 'Qrr_m', num2str(m));
  print(fh, '-depsc2', filename);
  
  fh = figure(2);
  clf;
  set(gca,'FontSize',18);
  hist(Qrl, 40);
  ylabel('#Elements','FontSize',18);
  xlabel('$Q_{rl} = 2  \sqrt{6}  \frac{R_{in}}{L_{max}}$','FontSize',18,'interpreter','latex');
  grid on
  filename = strcat( 'Qrl_m', num2str(m));
  print(fh, '-depsc2', filename);
  
  fh = figure(3);
  clf;
  set(gca,'FontSize',18);
  hist(Qtheta, 40);
  ylabel('#Elements','FontSize',18);
  xlabel('$Q_{\theta} = \frac{9 \sqrt{2}}{8}  V  S_{min}$','FontSize',18,'interpreter','latex');
  grid on
  filename = strcat( 'Qtheta_m', num2str(m));
  print(fh, '-depsc2', filename);
  
  fh = figure(4);
  clf;
  set(gca,'FontSize',18);
  hist(Qvl, 40);
  ylabel('#Elements','FontSize',18);
  xlabel('$Q_{vl} = 12  \frac{{\left(3V\right)}^{\frac{2}{3}}}{L_2}$','FontSize',18,'interpreter','latex');
  grid on
  filename = strcat( 'Qvl_m', num2str(m));
  print(fh, '-depsc2', filename);
  
end
