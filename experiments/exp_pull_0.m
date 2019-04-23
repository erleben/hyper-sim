%-- Experiment 0: Pulling forces on cantilever beam without gravity
%-- The experiment will be performed with a static value alpha = 0.9, which
%-- means that the beam should contract to 90% of its length and then stay
%-- in approximately that state.
%-- The tetrahedron vertices which contains one of the end points will be
%-- affected by boundary conditions which means that only the right side
%-- (positive) of the beam should move.
%
%-- We create cable points at -59 and 59, and therefore expect the maximum
%-- x to be at 48+/-1 and the minimum x to be at -60 (as before simulation
%-- begins). We also need to show that it doesn't get (much) lower than 49.

clear all;
close all;
clc;
% Add libraries
addpath('../fem');
addpath('../meshing');
addpath('../simulate');
addpath('../simulate/pull');
% Use a FEM method for simulating cable forces
method = fem_method();
% The cable that will be embedded in the beam mesh
cable = [-59, 0.0, 0.0; 0.0, 0.0, 0.0; 59, 0.0, 0.0];
% Create the scene for the simulation
scene  = pull_create_scene(1, '../', cable, false);
% Create the parameters for 'ecoflex-00-50' silicone, and use an adaptive
% step in the simulation
params = create_params('ecoflex-00-50','adaptive', 0.9, 0.001, false, true, 50, 1e3, false); 
params.E = 5000;
params.k = 1e5;
% Create profile info
pinfo  = create_profile_info();
pinfo.output_path = '../output/exp_pull_0/';
% Data generated by the simulation
pdata = simulate( params, method, scene, pinfo); 
save('exp_pull_0_pdata.mat', 'pdata');
save('exp_pull_0_pinfo.mat', 'pinfo');
% Process the data
postprocess(pinfo,pdata)