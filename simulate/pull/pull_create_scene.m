function [ scene ] = pull_create_scene(mesh_no, absolute_path, cable)

if nargin<1 || isempty(mesh_no)
  mesh_no = 1;
end
mesh_no = min( max( 1, mesh_no ), 8);

if nargin<2 || isempty(absolute_path)
  absolute_path = '../';
end

meshfile = '/meshing/model1.10.mat';
meshfile = strcat(absolute_path, meshfile);

cls = @pull_create_surface_traction_info;
cbc = @pull_create_boundary_conditions;
cbf = @pull_create_cable_forces;

scene = struct(...
  'meshfile', meshfile,...
  'create_surface_traction_info', cls,...
  'create_boundary_conditions', cbc,...
  'create_cable_forces', cbf, ...
  'cable', cable ...
  );

end