beta_0_5 = load('./results/exp_pressure_uniform_beta=-500_pdata.mat', 'pdata');
beta_1_0 = load('./results/exp_pressure_uniform_beta=1000_pdata.mat', 'pdata');


pdatas = [beta_0_5, beta_1_0];
betas = [-500, 1000];
time_simulated = 1;
% COM of meshes
time_coms = {length(pdatas(:, 1))};
coms = {length(pdatas(:, 1))};
maxxs = {length(pdatas(:, 1))};
maxys = {length(pdatas(:, 1))};
maxzs = {length(pdatas(:, 1))};
minxs = {length(pdatas(:, 1))};
minys = {length(pdatas(:, 1))};
minzs = {length(pdatas(:, 1))};
for n = 1:length(pdatas)
   %pdata_name = pdata_names(n, :);
   p = pdatas(n);
   
   %--- Profiling data members --------------------------------------------
   COM = p.pdata.COM;
   MAX_X = p.pdata.MAX_X;
   MAX_Y = p.pdata.MAX_Y;
   MAX_Z = p.pdata.MAX_Z;
   MIN_X = p.pdata.MIN_X;
   MIN_Y = p.pdata.MIN_Y;
   MIN_Z = p.pdata.MIN_Z;

   T = length(COM(:, 1));
   size(COM)
   time = [];
   t_real = 0;
   for t = 1:T
       % Computing time
       time = [time; t_real];
       t_real = t_real + (time_simulated./T);
   end
   coms{n} = COM;
   maxxs{n} = MAX_X;
   maxys{n} = MAX_Y;
   maxzs{n} = MAX_Z;
   minxs{n} = MIN_X;
   minys{n} = MIN_Y;
   minzs{n} = MIN_Z;
   time_coms{n} = time;
end

% Plots
% -------- Plot: COM of mesh nodes over time ------------------------------
for s = 1:length(pdatas)
    fh = figure('Visible','off');
    clf;
    % Plot settings
    x0=10;y0=10;width=850;height=600;
    set(gcf,'position',[x0,y0,width,height])
    set(gca,'FontSize',18);
    title('COM as a function of time [s]','FontSize',18);
    ylabel('Position','FontSize',18);
    xlabel('Time [s]','FontSize',18);
    hold on;
    % Plotting for all k's
    styles = ['-', '--', ':', '--', ':'];
    pos = coms{s};
    time = time_coms{s};
    xs = pos(:, 1);
    ys = pos(:, 2);
    zs = pos(:, 3);
    [~] = plot(time, xs, '-','LineWidth',2);
    hold on;
    [~] = plot(time, ys, '--','LineWidth',2);
    hold on;
    [~] = plot(time, zs, ':','LineWidth',2);
    hold on;
    % Add legend
    legend('Center of mass [x]', 'Center of mass [y]', 'Center of mass [z]', ...
           'Interpreter','latex',  'Location' , 'northwest');
    % Save plot
    hold off;
    axis tight;
    b = betas(s);
    filename = strcat('./plots/com_beta_', string(b), '.eps');
    print(fh, '-depsc2', filename);
    close(fh);
end

for s = 1:length(pdatas)

    % Plots
    % -------- Plot: Max position of mesh nodes over time ------------------------------
    fh = figure('Visible','off');
    clf;
    % Plot settings
    x0=10;y0=10;width=850;height=600;
    set(gcf,'position',[x0,y0,width,height])
    set(gca,'FontSize',18);
    title('Max-x, Max-y, Max-z positions over time [s]','FontSize',18);
    ylabel('Position','FontSize',18);
    xlabel('Time [s]','FontSize',18);
    hold on;
    % Plotting for all k's
    styles = ['-', '--', ':', '--', ':'];
    pos = coms{s};
    time = time_coms{s};
    xs = maxxs{s};
    ys = maxys{s};
    zs = maxzs{s};
    [~] = plot(time, xs, '-','LineWidth',2);
    hold on;
    [~] = plot(time, ys, '--','LineWidth',2);
    hold on;
    [~] = plot(time, zs, ':','LineWidth',2);
    hold on;
    % Add legend
    legend('Max positions [x]', 'Max positions  [y]', 'Max positions  [z]', ...
           'Interpreter','latex',  'Location' , 'northwest');
    % Save plot
    hold off;
    axis tight;
    b = betas(s);
    filename = strcat('./plots/max_beta_', string(b), '.eps');
    print(fh, '-depsc2', filename);
    close(fh);
end

% Plots
% -------- Plot: Min position of mesh nodes over time ------------------------------
for s = 1:length(pdatas)
    fh = figure('Visible','off');
    clf;
    % Plot settings
    x0=10;y0=10;width=850;height=600;
    set(gcf,'position',[x0,y0,width,height])
    set(gca,'FontSize',18);
    title('Min-x, Min-y, Min-z positions over time [s]','FontSize',18);
    ylabel('Position','FontSize',18);
    xlabel('Time [s]','FontSize',18);
    hold on;
    % Plotting for all k's

    pos = coms{s};
    time = time_coms{s};
    xs = minxs{s};
    ys = minys{s};
    zs = minzs{s};
    [~] = plot(time, xs, '-','LineWidth',2);
    hold on;
    [~] = plot(time, ys, '--','LineWidth',2);
    hold on;
    [~] = plot(time, zs, ':','LineWidth',2);
    hold on;
    % Add legend
    legend('Min positions [x]', 'Min positions  [y]', 'Min positions  [z]', ...
           'Interpreter','latex',  'Location' , 'northwest');
    % Save plot
    hold off;
    axis tight;
    b = betas(s);
    filename = strcat('./plots/min_beta_', string(b), '.eps');
    print(fh, '-depsc2', filename);
    close(fh);

end
