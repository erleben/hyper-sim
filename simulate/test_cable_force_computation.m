addpath('./exp_pull_alpha_0.95');
%scene = pull_create_scene();
%--- Load mesh data -------------------------------------------------------
load('../meshing/model1.12_upside_down.mat');

%T = double(T); %-- This is needed due to T being exported from python
%-- Embed cable C inside
%if isfield(scene, 'cable')
%  C = scene.cable; 
%end
%--- Setup computational mesh and scene -----------------------------------
%method = fem_method();
%mesh   = method.create_mesh(T,X,Y,Z);
