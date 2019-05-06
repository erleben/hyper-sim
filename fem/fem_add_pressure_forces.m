function [ state ] = fem_add_pressure_forces( state, pressure_info )
% Copyright 2011, Kenny Erleben
% Modified by Max Kragballe, 2019

%--- Get pressure info ----------------------------------------------------
Fp = pressure_info.Fp;     %-- Pressure forces at each triangle of the surface
F = pressure_info.F;
%--- Compute the nodal forces by distributing them uniformly across
%--- vertices
if isempty(F)
    state.fx(:) = state.fx(:) + 0;
    state.fy(:) = state.fy(:) + 0;
    state.fz(:) = state.fz(:) + 0;
    return
end
for f = 1:length(F(:, 1))
  force = Fp(f, :) ./ length(F(f, :));
  for v = 1:length(F(f, :))
     % Add the forces to the correct nodes
     state.fx(v) = state.fx(v) + force(1);
     state.fy(v) = state.fy(v) + force(2);
     state.fz(v) = state.fz(v) + force(3);
  end
end
end