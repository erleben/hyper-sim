function [ scene ] = pull_gravity_create_scene(mesh_no, absolute_path, cable, robot_name)

if nargin<1 || isempty(mesh_no)
  mesh_no = 1;
end
mesh_no = min( max( 1, mesh_no ), 8);

if nargin<2 || isempty(absolute_path)
  absolute_path = '../';
end


if nargin<4
    meshfile = '/meshing/robot.mat';
else
    meshfile = strcat('/meshing/', robot_name);
end
meshfile = strcat(absolute_path, meshfile);


cbc = @pull_gravity_create_boundary_conditions;
cls = @pull_gravity_create_surface_traction_info;
cbf = @pull_gravity_create_cable_forces;

scene = struct(...
  'meshfile', meshfile,...
  'create_surface_traction_info', cls,...
  'create_boundary_conditions', cbc,...
  'create_cable_forces', cbf, ...
  'cable', cable ...
  );

end