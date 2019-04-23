function [ state ] = fem_add_cable_forces( state, cable_info )
% Copyright 2011, Kenny Erleben
% Modified by Max Kragballe, 2019

%--- Get cable info -------------------------------------------------------
Fp = cable_info.Fp;     %-- Cable forces at the via points
W = cable_info.W;       %-- Barycentric weighting of via point positions
indices = cable_info.indices;
%--- Compute the nodal forces by solving
%       Fp = W Fv =>
%       W' Fp = W'W Fv =>
%       Fv = (W'W)^-1 W' Fp
Fv = pinv(W' * W) * W' * Fp;
fprintf("Cable forces")
max(Fv)
% Add the forces to the correct nodes
state.fx(indices,:) = state.fx(indices,:) + Fv(:, 1);
state.fy(indices,:) = state.fy(indices,:) + Fv(:, 2);
state.fz(indices,:) = state.fz(indices,:) + Fv(:, 3);
end