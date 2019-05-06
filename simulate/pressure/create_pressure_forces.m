function [ pressure_info ] = create_pressure_forces( time, state, mesh )
%PULL_CREATE_CABLE_FORCES -- Computes cable forces and nodes to apply them to.
%
% INPUT:
%
%   time         - The current simulated time.
%   state        - The current state.
%   mesh         - The initial  mesh.
%
% OUTPUT:
%
% pressure_info  - Upon return this struct will hold the cable force to apply
%               at every node, as well as the weighting of the vertices
%               used to interpolate the cable node positions.
%
% Author: Max Kragballe, 2019

%--- Get triangle and pressure info ---------------------------------------
F    = state.F;
if isempty(F)
  return
end

x    = state.x;
y    = state.y;
z    = state.z;

%--- Get current pressure -------------------------------------------------
pressure_struct = state.pressure;
time_to_infalte = pressure_struct.inflation;
beta = max(pressure_struct.beta, time_to_inflate - time*0.1); %-- Adjusting beta always takes time_to_inflate seconds
pressure = beta * pressure_struct.pressure;

%--- Get triangle indices -------------------------------------------------
i = F(:,1);
j = F(:,2);
k = F(:,3);

%--- Get vertex coordinates -----------------------------------------------
Pi = [ x(i), y(i), z(i) ]';
Pj = [ x(j), y(j), z(j) ]';
Pk = [ x(k), y(k), z(k) ]';

%--- Compute face areas ---------------------------------------------------
Avec =  cross( (Pj - Pi), (Pk - Pi) )  ./ 2.0 ;
A    =  sum( Avec.*Avec, 1).^(0.5);

% Compute pressure forces for each face 
Fp = zeros(length(F(:, 1)), 3);
for f=1:length(F(:,1))
  % Normal of face f
  i = F(f, 1);
  j = F(f, 2);
  k = F(f, 3);
  
  Pi = [x(i), y(i), z(i)];
  Pj = [x(j), y(j), z(j)];
  Pk = [x(k), y(k), z(k)];
  
  u = Pj - Pi;
  v = Pk - Pi;
  n = [u(2) * v(3) - u(3) * v(2); ...
       u(3) * v(1) - u(1) * v(3); ...
       u(1) * v(2) - u(2) * v(1)  ...
      ];
  % F = nP / A
  Fp(f,:) = n .* pressure /  A(f);
end

%--- Bundle all info into one structure ----------------------------------
pressure_info = struct(...
  'Fp', Fp, ...
  'F', F ...
  );

end