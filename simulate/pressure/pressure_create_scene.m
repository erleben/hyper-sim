function [ scene ] = pressure_create_scene(mesh_no, absolute_path, robot_name, initial_pressure, beta, inflation)

if nargin<1 || isempty(mesh_no)
  mesh_no = 1;
end
mesh_no = min( max( 1, mesh_no ), 8);

if nargin<2 || isempty(absolute_path)
  absolute_path = '../';
end


if nargin<5
    meshfile = '/meshing/robot.mat';
else
    meshfile = strcat('/meshing/', robot_name);
end
meshfile = strcat(absolute_path, meshfile);


cbc = @pressure_create_boundary_conditions;
cls = @pressure_create_surface_traction_info;
cpf = @create_pressure_forces;

pressure = struct(...
    'pressure', initial_pressure, ...
    'beta', beta, ...
    'inflation', inflation ...
    );


scene = struct(...
  'meshfile', meshfile,...
  'create_surface_traction_info', cls,...
  'create_boundary_conditions', cbc,...
  'create_pressure_forces', cpf, ...
  'pressure', pressure ...
  );

end