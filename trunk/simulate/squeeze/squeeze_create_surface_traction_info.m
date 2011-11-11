function [ traction_info ] = squeeze_create_surface_traction_info( time, state, mesh )
% Copyright 2011, Kenny Erleben

F = [];
tx = [];
ty = [];
tz = [];

%--- Bundle all info into one structure ----------------------------------
traction_info = struct(...
  'F', F,...
  'tx', tx,...
  'ty', ty,...
  'tz', tz...
  );

end