function postprocess( profile_info, profile_data)

%--- Profiling data members -----------------------------------------------
H  = profile_data.H;   % Time step sizes
W  = profile_data.W;   % Wall clock times of time steps
V  = profile_data.V;   % Volume monitoring
KE = profile_data.KE;  % Kinetic energy monitoring
P  = profile_data.P;   % Test point monitoring
C  = profile_data.C;   % Convergence rate monitoring

cumH = zeros(size(H));
cumH(1) = H(1);
for h=2:length(H)
  cumH(h) = cumH(h-1) + H(h);
end

cumW = zeros(size(W));
cumW(1) = W(1);
for h=2:length(W)
  cumW(h) = cumW(h-1) + W(h);
end

fh = figure('Visible','off');
clf;
set(gca,'FontSize',18);
h1 = plot(cumW./cumH, '-','LineWidth',2,'Color',[0.7 0.1, 0.1]);
grid on
hold on
title('Real Time Ratio','FontSize',18);
xlabel('Iteration','FontSize',18);
ylabel('Wall clock / Simulated  [%]','FontSize',18);
hold off;
axis tight;
filename = strcat(  profile_info.output_path,  profile_info.filename_prefix, 'ratio'  );
print(fh, '-depsc2', filename);

fh = figure('Visible','off');
clf;
set(gca,'FontSize',18);
h1 = plot(W, '-','LineWidth',2,'Color',[0.7 0.1, 0.1]);
grid on
hold on
h2 = plot(sort(W), '-','LineWidth',1.75,'Color',[0.1 0.1, 0.7]);
title('Wall clock time of time step','FontSize',18);
xlabel('Iteration','FontSize',18);
ylabel('Value [s]','FontSize',18);
legend([h1, h2],'Unsorted','Sorted');
hold off;
axis tight;
filename = strcat(  profile_info.output_path,  profile_info.filename_prefix, 'wall_clock'  );
print(fh, '-depsc2', filename);

Gmin = V(:,1);
Gmax = V(:,2);
Gtot = V(:,3);

fh = figure('Visible','off');
clf;
set(gca,'FontSize',18);
h1 = plot(Gtot, '-','LineWidth',2,'Color',[0.7 0.1, 0.1]);
grid on
hold on
h2 = plot(Gmin, '-','LineWidth',2,'Color',[0.1 0.7, 0.1]);
h3 = plot(Gmax, '-','LineWidth',2,'Color',[0.1 0.1, 0.7]);
title('Volume Gain','FontSize',18);
xlabel('Iteration','FontSize',18);
ylabel('Procentage [%]','FontSize',18);
legend([h1, h2, h3],'Total','Min Element','Max Element','Location','SouthWest');
hold off;
axis tight;
filename = strcat(  profile_info.output_path,  profile_info.filename_prefix, 'volume'  );
print(fh, '-depsc2', filename);

fh = figure('Visible','off');
clf;
set(gca,'FontSize',18);
plot(KE, '-','LineWidth',2,'Color',[0.7 0.1, 0.1]);
grid on
hold on
title('Kinetric Energy','FontSize',18);
xlabel('Iteration','FontSize',18);
ylabel('Value [J]','FontSize',18);
hold off;
axis tight;
filename = strcat(  profile_info.output_path,  profile_info.filename_prefix, 'energy'  );
print(fh, '-depsc2', filename);

fh = figure('Visible','off');
clf;
set(gca,'FontSize',18);
h1 = plot(P(:,1), '-','LineWidth',2,'Color',[0.7 0.1, 0.1]);
grid on
hold on
h2 = plot(P(:,2), '.','LineWidth',2,'Color',[0.1 0.7, 0.1]);
h3 = plot(P(:,3), '--','LineWidth',2,'Color',[0.1 0.1, 0.7]);
title('Test Point','FontSize',18);
xlabel('Iteration','FontSize',18);
ylabel('Value [m]','FontSize',18);
legend([h1, h2, h3],'X','Y','Z');
hold off;
axis tight;
filename = strcat(  profile_info.output_path,  profile_info.filename_prefix, 'coord'  );
print(fh, '-depsc2', filename);


stats = zeros(3,4);

stats(1,1) = mean(H);
stats(1,2) = min(H);
stats(1,3) = max(H);
stats(1,4) = std(H);

stats(2,1) = mean(W);
stats(2,2) = min(W);
stats(2,3) = max(W);
stats(2,4) = std(W);

stats(3,1) = mean(V(:,3));
stats(3,2) = min(V(:,3));
stats(3,3) = max(V(:,3));
stats(3,4) = std(V(:,3));


filename = strcat(  profile_info.output_path,  profile_info.filename_prefix, 'stats.tex'  );
fid = fopen( filename, 'w');

fprintf(fid, '\\hline \n');
fprintf(fid, '{} & Mean & Min & Max & Std \\\\ \n');
fprintf(fid, '\\hline \n');
fprintf(fid, 'Time step size & %6.4f & %6.4f & %6.4f & %6.4f \\\\ \n', stats(1,1), stats(1,2), stats(1,3), stats(1,4) );
fprintf(fid, 'Wall clock     & %6.4f & %6.4f & %6.4f & %6.4f \\\\ \n', stats(2,1), stats(2,2), stats(2,3), stats(2,4) );
fprintf(fid, 'Volume         & %6.4f & %6.4f & %6.4f & %6.4f \\\\ \n', stats(3,1), stats(3,2), stats(3,3), stats(3,4) );
fprintf(fid, '\\hline \n');
fclose(fid);


if ~isempty(C)
  fh = figure('Visible','off');
  clf;
  set(gca,'FontSize',18);
  c = C{1};
  color = rand(1,3)*0.6 + 0.1;
  semilogy(c, '-','LineWidth',2,'Color', color);
  grid on
  hold on
  for i=2:length(C)
    c = C{i};
    color = rand(1,3)*0.6 + 0.1;
    semilogy(c, '-','LineWidth',2,'Color', color);
  end
  title('Convergence Rate','FontSize',18);
  xlabel('Iteration','FontSize',18);
  ylabel('Error','FontSize',18);
  hold off;
  axis tight;
  filename = strcat(  profile_info.output_path,  profile_info.filename_prefix, 'convergence'  );
  print(fh, '-depsc2', filename);
end

end
