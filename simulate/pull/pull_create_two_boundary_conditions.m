function [ BC ] = pull_create_boundary_conditions( time, state, mesh )
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

x = mesh.x0;
y = mesh.y0;
z = mesh.z0;

% Get mesh dimensions
%max_x            = max( x );
%min_x            = min( x );
%max_y            = max( y );
%min_y            = min( y );
%max_z            = max( z );
%min_z            = min( z );
delta            = 0.01;
max_x            = max(x);
max_x = 65;
min_x            = min(x);
min_x = -65;
max_y            = max( y );
min_y            = min( y );
max_z            = max( z );
min_z            = min( z );



% Creates boundary conditions at one end of the finger-like robot
box = [...
  min_x-5.1, min_y-delta,...
  min_z-delta, min_x+5.1,...
  max_y+delta, max_z+delta...
  ];

x1_min = box(1) .* ones( size( x ) );
y_min = box(2) .* ones( size( y ) );
z_min = box(3) .* ones( size( z ) );

x1_max = box(4) .* ones( size( x ) );
y_max = box(5) .* ones( size( y ) );
z_max = box(6) .* ones( size( z ) );

x2_min = (max_x-5.1) .* ones( size( x ) );
x2_max = (max_x+5.1) .* ones( size( x ) );

I = find( (((x1_min<x)&(x<x1_max))|((x2_min<x) & x<x2_max))&(y_min<y)&(y<y_max)&(z_min<z)&(z<z_max) );
xoffset = 0;
yoffset = length(x);
zoffset = length(y)+length(x);

idx    = [I+xoffset; I+yoffset; I+zoffset];
values = [x(I); y(I); z(I) ];

BC = struct( 'idx', idx,...
  'values', values...
  );

end