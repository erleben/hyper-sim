function [ state ] = cor_create_state( mesh, params )
% Copyright 2011, Kenny Erleben

cntV = length(mesh.x0);        % Number of vertices

K   = zeros(cntV*3,cntV*3);    % Allocate global stiffness matrix
fu  = zeros(cntV*3,1);         % Allocate global force offset vector ( - K p0 term)

vx = zeros(size(mesh.x0));     % Velocity field
vy = zeros(size(mesh.y0));
vz = zeros(size(mesh.z0));

fx = zeros(size(mesh.x0));     % External force field
fy = zeros(size(mesh.y0));
fz = zeros(size(mesh.z0));

x  = mesh.x0;                  % Spatial coordintes
y  = mesh.y0;
z  = mesh.z0;

%--- Do material mesh pre-computations ------------------------------------
Me    = cor_compute_mass_elements(mesh, params);
Ke    = cor_compute_stiffness_elements(mesh, params);
%Ce    = cor_compute_rayleigh_damping_elements(mesh, params);
Ce    = cor_compute_damping_elements(mesh, params);
fue   = cor_compute_offset_force_elements(Ke, mesh);
M     = cor_assemble_global_matrix( mesh, Me );
C     = cor_assemble_global_matrix( mesh, Ce );

state = struct('x',  x, 'y', y, 'z', z,...
  'vx', vx, 'vy', vy, 'vz', vz,...
  'fx', fx, 'fy', fy, 'fz', fz,...
  'K', K, 'M', M, 'C', C,...
  'fue', fue,...
  'Ke', Ke,...
  'fu', fu...
  );

end