function [ state ] = fvm_create_state( mesh, params )
% Copyright 2011, Kenny Erleben

E = length(mesh.T(:,1));   % Number of tetrahedrons
V = length(mesh.x0);       % Number of vertices

M   = zeros(V,1);          % Allocate global mass matrix

vx = zeros(V,1);           % Velocity field
vy = zeros(V,1);
vz = zeros(V,1);

fx = zeros(V,1);           % External force field
fy = zeros(V,1);
fz = zeros(V,1);

fex = zeros(V,1);           % Elastic forces
fey = zeros(V,1);
fez = zeros(V,1);

x  = mesh.x0;              % Spatial coordintes (initialized to be equal to material coordiantes)
y  = mesh.y0;
z  = mesh.z0;

for e=1:E
  
  % Get tetrahedron indices
  i = mesh.T(e,1);
  j = mesh.T(e,2);
  k = mesh.T(e,3);
  m = mesh.T(e,4);
  
  % Get vertex coordinates
  Pi = [ x(i),  y(i),  z(i) ]';
  Pj = [ x(j),  y(j),  z(j) ]';
  Pk = [ x(k),  y(k),  z(k) ]';
  Pm = [ x(m),  y(m),  z(m) ]';
      
  % Compute the tetrahedron element volumes.
  %
  % Using triple scalar product the tetrahedron volumes can be computed as
  %
  %   V   =  \frac{1}{6}  (P_m - P_i)  \cdot   ( P_j - P_i)  \times  ( P_k - P_i)
  %
  V = dot( (Pm - Pi) , cross( (Pj - Pi), (Pk - Pi) ) ) ./ 6.0;
  
  M(i) = M(i) + params.rho*V/4;
  M(j) = M(j) + params.rho*V/4;
  M(k) = M(k) + params.rho*V/4;
  M(m) = M(m) + params.rho*V/4;
  
end

C = (params.c/params.rho)  .* M;  % Global viscous lumped damping matrix

state = struct(   'x',  x, 'y', y, 'z', z,...
  'vx', vx, 'vy', vy, 'vz', vz,...
  'fx', fx, 'fy', fy, 'fz', fz,...  
  'M', M,...
  'C', C,...
  'fex', fex,...
  'fey', fey,...
  'fez', fez...  
  );

end