% Copyright 2011, Kenny Erleben
clear all;
close all;
clc;

addpath('../fem');
addpath('../cor');
addpath('../fvm');
addpath('./twist');
addpath('./twist2');
addpath('./twist3');
addpath('./bend');
addpath('./bend2');
addpath('./bend3');
addpath('./squeeze');
addpath('./stretch');

% method = fem_method();
% method = cor_method();
method = fvm_method();

% scene  = twist_create_scene();   % Constant traction
% scene  = twist2_create_scene();  % Increasing traction
% scene  = twist3_create_scene();  % No traction
% scene  = bend_create_scene();    % Constant traction
 scene  = bend2_create_scene();   % Increasing traction
% scene  = bend3_create_scene();   % No traction
% scene  = squeeze_create_scene();
% scene  = stretch_create_scene();

params = create_params();
%params = create_params([],'adaptive');
%--- COR method is based on Hookean material but FVM and FEM uses 
%--- a St. Venant-Kirchoff material. This gives COR a more stiff behavior/
%--- appearance. Thus, for animation one would most like use unrealistic
%--- small values of Yound modulos for the COR method.
%
% params.E = 5000;
%
params.T  = 2.0;

pinfo  = create_profile_info();
% pinfo.debug_level = 1;


pdata = simulate( params, method, scene, pinfo ); 

postprocess(pinfo,pdata);