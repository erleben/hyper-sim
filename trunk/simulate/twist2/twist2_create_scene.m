function [ scene ] = twist2_create_scene(mesh_no, absolute_path)
% Copyright 2011, Kenny Erleben

if nargin<1 || isempty(mesh_no)
  mesh_no = 1;
end
mesh_no = min( max( 1, mesh_no ), 8);

if nargin<2 || isempty(absolute_path)
  absolute_path = '../';
end

switch(mesh_no)
 case 1
  meshfile = '/meshing/bar_T1098_V300.mat';
 case 2
  meshfile = '/meshing/bar_T1197_V325';
 case 3
  meshfile = '/meshing/bar_T2293_V576';
 case 4
  meshfile = '/meshing/bar_T7737_V1701';
 case 5
  meshfile = '/meshing/bar_T8150_V1782';
  case 6
   meshfile = '/meshing/bar_T11741_V2500';
 case 7
  meshfile = '/meshing/bar_T20423_V4176';
 case 8
  meshfile = '/meshing/bar_T27977_V5577';
end

meshfile = strcat(absolute_path, meshfile);

cls = @twist2_create_surface_traction_info;
cbc = @twist2_create_boundary_conditions;

scene = struct(...
  'meshfile', meshfile,...
  'create_surface_traction_info', cls,...
  'create_boundary_conditions', cbc...
  );

end