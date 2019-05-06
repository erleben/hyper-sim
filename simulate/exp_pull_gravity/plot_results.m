pdata_names = {...
    './results/exp_pull_gravity_k=0_pdata.mat'; ...
    './results/exp_pull_gravity_k=100000_pdata.mat'; ...
    './results/exp_pull_gravity_k=10000000_pdata.mat'; ...
};


number_cable_points = 131;

stiffnesses = [0, 1e5, 1e7];

% Cable lengths
time_CLs = {length(pdata_names(:, 1))};
CLs = {length(pdata_names(:, 1))};
% Mid Point positions
time_MidPoints = {length(pdata_names(:, 1))};
MidPoints = {length(pdata_names(:, 1))};
min_zs = {length(pdata_names(:, 1))};
% Computing cable length and position of middle cable point over time
for n = 1:length(pdata_names)
   pdata_name = pdata_names(n, :);
   pdata = load(pdata_names{1});
   pdata = pdata.pdata;
   %--- Profiling data members --------------------------------------------
   CP = pdata.CABLE_POINTS;
   T = floor(length(CP(:, 1)) / number_cable_points);
   time = [];
   cable_lengths = [];
   midPoints  = [];
   t_real = 0;
   for t = 1:T
       % Computing cable length
       cl = CP((1 + (t-1)*number_cable_points):t*number_cable_points, :);
       lbar = cl(2:end, :) - cl(1:end-1, :);
       l = sum(vecnorm(lbar, 2, 2));
       cable_lengths = [cable_lengths; l];
       % Computing midpoint
       midPoints = [midPoints ; CP((1 + (t-1)*number_cable_points)+floor(number_cable_points/2), :)];
       % Computing time
       time = [time; t_real];
       t_real = t_real + (4./T);
   end
   % Add cable lengths and mid point positions to global arrays
   CLs{n} = cable_lengths;
   MidPoints{n} = midPoints;
   min_zs{n} = pdata.MIN_Z;
   time_CLs{n} = time;
   time_MidPoints{n} = time;
end
time_min_zs = {length(pdata_names(:, 1))};
time = [];
for t=1:floor(4 * 30 + 1)
    time = [time; t*1/30;];
end
for n = 1:length(pdata_names(:, 1))
    time_min_zs{n} = time;
end

% Plots
% -------- Plot: Cable Length [mm] as a function of time [s] --------------
fh = figure('Visible','off');
clf;
% Plot settings
x0=10;y0=10;width=850;height=600;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',18);
title('Cable Length [mm] as a function of time [s]','FontSize',18);
ylabel('Cable Length [mm]','FontSize',18);
xlabel('Time [s]','FontSize',18);
hold on;

% Plotting for all k's
styles = {'-', '--', ':'};
for s=1:length(pdata_names)
    style = styles{s};
    cable_lengths = CLs{s};
    time = time_CLs{s};
    [~] = plot(time, cable_lengths(:, 1), style,'LineWidth',2);
    hold on;
end
% Add legend
legend('k=0', 'k=1e5', 'k=1e7', ...
       'Interpreter','latex',  'Location' , 'northeast');
% Save plot
hold off;
axis tight;
filename = './plots/exp_gravity_cable_lengths.eps';
print(fh, '-depsc2', filename);
close(fh);


% -------- Plot: Middle Cable Position [z] as a function of time [s] --
fh = figure('Visible','off');
clf;
% Plot settings
x0=10;y0=10;width=850;height=600;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',18);
title('Middle Cable Position [z] as a function of time [s]','FontSize',18);
ylabel('Position [z]','FontSize',18);
xlabel('Time [s]','FontSize',18);
hold on;

% Plotting for all k's
styles = {'-'; '--'; ':'};
for s=1:length(pdata_names)
    style = styles{s};
    midPoints = MidPoints{s};
    time = time_MidPoints{s};
    [~] = plot(time, midPoints(:, 3), style, 'LineWidth',2);
    hold on
    grid on
end
% Add legend
legend('k=0', 'k=1e5', 'k=1e7', ...
       'Interpreter','latex',  'Location' , 'southeast');
% Save plot
hold off;
axis tight;
filename = './plots/exp_gravity_mid_z_position.eps';
print(fh, '-depsc2', filename);
close(fh);

% -------- Plot: Minimum Position [z] of the mesh as a function of time [s] --
fh = figure('Visible','off');
clf;
% Plot settings
x0=10;y0=10;width=850;height=600;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',18);
title('Minimum Position [z] of the mesh as a function of time [s] ','FontSize',18);
ylabel('Position [z]','FontSize',18);
xlabel('Time [s]','FontSize',18);
hold on;

% Plotting for all k's
styles = {'-'; '--'; ':'};
for s=1:length(pdata_names)
    style = styles{s};
    midPoints = min_zs{s};
    time = time_min_zs{s};
    size(midPoints)
    size(time)
    [~] = plot(time, midPoints(:, 1), style, 'LineWidth',2);
    %[~] = plot(midPoints(:, 1), style, 'LineWidth',2);
    hold on
    grid on
end
% Add legend
legend('k=0', 'k=1e5', 'k=1e7', ...
       'Interpreter','latex',  'Location' , 'southeast');
% Save plot
hold off;
axis tight;
filename = './plots/exp_gravity_min_point_position.eps';
print(fh, '-depsc2', filename);
close(fh);