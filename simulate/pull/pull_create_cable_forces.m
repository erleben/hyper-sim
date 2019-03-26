function [ cable_info ] = pull_create_cable_forces( time, state, mesh )
%PULL_CREATE_CABLE_FORCES -- Computes cable forces and nodes to apply them to.
%
% INPUT:
%
%   time      - The current simulated time.
%   state     - The current state.
%   mesh      - The initial  mesh.
%
% OUTPUT:
%
% cable_info  - Upon return this struct will hold the cable force to apply
%               at every node, as well as the weighting of the vertices
%               used to interpolate the cable node positions.
%
% Author: Max Kragballe, 2019

%-- Get the current cable at time (t)
cable = state.cable;
T = mesh.T;  %-- Mesh tetrahedralization
X = mesh.x0; %-- Mesh x coordinates
Y = mesh.y0; %-- Mesh y coordinates
Z = mesh.z0; %-- Mesh z coordinates

%--- Barycentric weights where
%       Wij := The weight of node j in terms of the position of via point i
W = cable.W;
%--- Cable forces at each point on the cable (in terms of (x,y,z) forces)
Fp = zeros(length(cable.W(:, 1)), 3);
%-- Compute the cable forces given by
%       -k max(l - a*l0, 0)
V = [X(:), Y(:), Z(:)];                                %-- Positions of all vertices
P = W * V;                                             %-- Re-interpolated points of the cable
ls = vecnorm(P(2:end, :) - P(1:end-1, :), 2, 2);       %-- Length of the individual cables
l  = sum(ls);                                          %-- Length of the entire cable
l0 = sum(cable.L0s);                                   %-- Length of cable at time t=0
alpha = 0.5;                                           %-- Control parameter
Fc = -cable.k * max(l - alpha * l0, 0);                %-- The cable forces (identical across the cable)

%--- Compute directions to apply forces
ds = zeros(length(P(:, 1))-1, 3);
for d = 1:length(ds(:, 1))
   ds(d, :) = (P(d+1, :) - P(d)) / ls(d);
end

%--- Add forces based on directions to Fp
for p = 1:length(Fp(:, 1))
    if p == 1
        Fp(p, :) = Fc .* ds(p, :);                  %-- Start attachment point
    else
    if p == length(Fp(:, 1))
        Fp(p, :) = -Fc .* ds(p-1, :);               %-- End attachment point
    else
        Fp(p, :) = Fc .* (ds(p, :) - ds(p-1, :));   %-- Regular via point
    end
    end
end

%--- Bundle all info into one structure ----------------------------------
cable_info = struct(...
  'Fp', Fp, ...
  'W', W ...
  );

end