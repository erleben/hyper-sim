k_1e5 = load('./results/exp_pull_1_k=100000_pdata.mat', 'pdata');
k_5e5 = load('./results/exp_pull_1_k=500000_pdata.mat', 'pdata');
k_1e6 = load('./results/exp_pull_1_k=1000000_pdata.mat', 'pdata');
k_2_5e6 = load('./results/exp_pull_1_k=2500000_pdata.mat', 'pdata');
k_5e6 = load('./results/exp_pull_1_k=5000000_pdata.mat', 'pdata');
k_1e7 = load('./results/exp_pull_1_k=10000000_pdata.mat', 'pdata');

pdatas = [k_1e5, k_5e5, k_1e6, k_2_5e6, k_5e6, k_1e7];

number_cable_points = 140;

% Min z positon of mesh
time_min_z = {length(pdatas(:, 1))};
min_z = {length(pdatas(:, 1))};
% Min z positon of cable
time_min_c = {length(pdatas(:, 1))};
min_c = {length(pdatas(:, 1))};
%CP(number_cable_points, 3)
% Cable lengths
cable_lengths = {length(pdatas(:, 1))};
for n = 1:length(pdatas)
   %pdata_name = pdata_names(n, :);
   p = pdatas(n);
   %--- Profiling data members --------------------------------------------
   MIN_Z = p.pdata.MIN_Z;
   CP = p.pdata.CABLE_POINTS;
   CP(end, :)
   T = floor(length(CP(:, 1)) / number_cable_points);
   time = [];
   c_pos  = [];
   t_real = 0;
   C_lengths = [];
   for t = 1:T
       % Computing min cable position
       c_pos = [c_pos; CP(t*number_cable_points, 3)];
       % Computing cable length 
       cl = CP((1 + (t-1)*number_cable_points):t*number_cable_points, :);
       lbar = cl(2:end, :) - cl(1:end-1, :);
       l = sum(vecnorm(lbar, 2, 2));
       C_lengths = [C_lengths; l];
       % Computing time
       time = [time; t_real];
       t_real = t_real + (2./T);
   end
   % Add cable lengths and mid point positions to global arrays
   min_z{n} = MIN_Z;
   min_c{n} = c_pos;
   cable_lengths{n} = C_lengths;
   time_min_c{n} = time;
end
time = [];
for t=1:floor(1.5 * 30)
    time = [time; t*1/30;];
end
for n = 1:length(pdatas)
    time_min_z{n} = time;
end

% Plots
% -------- Plot: Min z-position of mesh nodes --------------
fh = figure('Visible','off');
clf;
% Plot settings
x0=10;y0=10;width=850;height=600;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',18);
title('z-position of last cable point over time [s]','FontSize',18);
ylabel('Position [z]','FontSize',18);
xlabel('Time [s]','FontSize',18);
hold on;
% Plotting for all k's
styles = ['-', '--', ':', '--', ':'];
for s=1:length(pdatas)
    style = '-';%styles(s);
    c_pos = min_c{s};
    time = time_min_c{s};
    [~] = plot(time, c_pos(:, 1), style,'LineWidth',2);
    hold on;
end
% Add legend
legend('k=1e5', 'k=5e5', 'k=1e6', 'k=2.5e6','k=5e6', 'k=1e7', ...
       'Interpreter','latex',  'Location' , 'northwest');
% Save plot
hold off;
axis tight;
filename = './plots/exp_pull_1_min_cable_z_position.eps';
print(fh, '-depsc2', filename);
close(fh);

% -------- Plot: cable lengths --
fh = figure('Visible','off');
clf;
% Plot settings
x0=10;y0=10;width=850;height=600;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',18);
title('Cable length [mm] as a function of time [s]','FontSize',18);
ylabel('Cable length [mm]','FontSize',18);
xlabel('Time [s]','FontSize',18);
hold on;
% Plotting for all k's
styles = ['-', '--', ':', '--', ':'];
for s=1:length(pdatas)
    style = '-';%styles(s);
    cl = cable_lengths{s};
    time = time_min_c{s};
    [~] = plot(time, cl, style,'LineWidth',2);
    hold on;
end
% Add legend
legend('k=1e5', 'k=5e5', 'k=1e6', 'k=2.5e6','k=5e6', 'k=1e7', ...
    'Interpreter','latex',  'Location' , 'northwest');
% Save plot
hold off;
axis tight;
filename = './plots/exp_pull_1_cable_lengths.eps';
print(fh, '-depsc2', filename);
close(fh);

% -------- Plot: z-position of last cable point --
fh = figure('Visible','off');
clf;
% Plot settings
x0=10;y0=10;width=850;height=600;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',18);
title('Min z-position of mesh nodes over time [s]','FontSize',18);
ylabel('Position [z]','FontSize',18);
xlabel('Time [s]','FontSize',18);
hold on;
% Plotting for all k's
styles = ['-', '--', ':', '--', ':'];
for s=1:length(pdatas)
    style = '-';%styles(s);
    z_pos = min_z{s};
    time = time_min_z{s};
    [~] = plot(time, z_pos(:, 1), style,'LineWidth',2);
    hold on;
end
% Add legend
legend('k=1e5', 'k=5e5', 'k=1e6', 'k=2.5e6','k=5e6', 'k=1e7', ...
    'Interpreter','latex',  'Location' , 'northwest');
% Save plot
hold off;
axis tight;
filename = './plots/exp_pull_1_min_mesh_z_position.eps';
print(fh, '-depsc2', filename);
close(fh);

