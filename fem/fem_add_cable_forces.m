function [ state ] = fem_add_cable_forces( state, cable_info )
% Copyright 2011, Kenny Erleben
% Modified by Max Kragballe, 2019

%--- Get cable info -------------------------------------------------------
Fp = cable_info.Fp;     %-- Cable forces at the via points
W = cable_info.W;       %-- Barycentric weighting of via point positions

%--- Compute the nodal forces by solving
%       Fp = W Fv =>
%       W' Fp = W'W Fv =>
%       Fv = (W'W)^-1 W' Fp
% Fv = inv(W' * W) * W' * Fp; %-- Not possible as (W' * W) is not invertible...
%--- Compute the nodal forces by solving
%       W Fv = Fp =>
%       Fv = W \ Fp
Fv = W \ Fp;
% Add the forces to the correct nodes
state.fx(:) = state.fx(:) + Fv(:, 1);
state.fy(:) = state.fy(:) + Fv(:, 2);
state.fz(:) = state.fz(:) + Fv(:, 3);
end