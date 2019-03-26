% Copyright 2011, Kenny Erleben
clear all;
close all;
clc;
m=10;
load('model1.10.mat');
X = n(:, 1);
Y = n(:, 2);
Z = n(:, 3);
T = e;
[ Qrr, Qrl, Qtheta, Qvl ] = compute_quality_measures(double(T)+1.0, X, Y, Z);
fh = figure(1);
clf;
set(gca,'FontSize',18);
hist(Qrr, 40);
ylabel('value','FontSize',18);
xlabel('Q_{rr}','FontSize',18);
filename = strcat( 'Qrr_m', num2str(m));
print(fh, '-depsc2', filename);

fh = figure(2);
clf;
set(gca,'FontSize',18);
hist(Qrl, 40);
ylabel('value','FontSize',18);
xlabel('Q_{rl}','FontSize',18);
filename = strcat( 'Qrl_m', num2str(m));
print(fh, '-depsc2', filename);

fh = figure(3);
clf;
set(gca,'FontSize',18);
hist(Qtheta, 40);
ylabel('value','FontSize',18);
xlabel('Q_{\theta}','FontSize',18);
filename = strcat( 'Qtheta_m', num2str(m));
print(fh, '-depsc2', filename);

fh = figure(4);
clf;
set(gca,'FontSize',18);
hist(Qvl, 40);
ylabel('value','FontSize',18);
xlabel('Q_{vl}','FontSize',18);
filename = strcat( 'Qvl_m', num2str(m));
print(fh, '-depsc2', filename);
