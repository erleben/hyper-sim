function [ BC ] = pressure_create_boundary_conditions( time, state, mesh )
%TWIST_CREATE_BOUNDARY_CONDITIONS -- Creates a fixed position boundary condition.
%
% INPUT:
%
%   time   - The current simulated time.
%   state  - The current state.
%   mesh   - The initial  mesh.
%
% OUTPUT:
%
% BC    - Upon return this struct will hold indices and values for the
%         wanted boundary conditions.
%
% Copyright 2011, Kenny Erleben
%
% Modified by Max Kragballe, 2019

% No boundary conditions
idx = [];
values = [];

BC = struct( 'idx', idx,...
  'values', values...
  );

end