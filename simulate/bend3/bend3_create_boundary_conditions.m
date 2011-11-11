function [ BC ] = bend3_create_boundary_conditions( time, state, mesh )
%BEND3_CREATE_BOUNDARY_CONDITIONS -- Creates a fixed position boundary condition.
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

x = mesh.x0;
y = mesh.y0;
z = mesh.z0;

max_x            = max( x );
min_x            = min( x );
max_y            = max( y );
min_y            = min( y );
max_z            = max( z );
min_z            = min( z );
delta            = 0.01;

box = [...
  min_x-delta, min_y-delta,...
  min_z-delta, min_x+delta,...
  max_y+delta, max_z+delta...
  ];
x_min = box(1) .* ones( size( x ) );
y_min = box(2) .* ones( size( y ) );
z_min = box(3) .* ones( size( z ) );
x_max = box(4) .* ones( size( x ) );
y_max = box(5) .* ones( size( y ) );
z_max = box(6) .* ones( size( z ) );

Il = find( (x_min<x)&(x<x_max)&(y_min<y)&(y<y_max)&(z_min<z)&(z<z_max) );

box = [...
  max_x-delta, min_y-delta, min_z-delta,...
  max_x+delta, max_y+delta, max_z+delta...
  ];
x_min = box(1) .* ones( size( x ) );
y_min = box(2) .* ones( size( y ) );
z_min = box(3) .* ones( size( z ) );
x_max = box(4) .* ones( size( x ) );
y_max = box(5) .* ones( size( y ) );
z_max = box(6) .* ones( size( z ) );

Ir = find( (x_min<x)&(x<x_max)&(y_min<y)&(y<y_max)&(z_min<z)&(z<z_max) );

xoffset = 0;
yoffset = length(x);
zoffset = 2*length(x);

angle = (pi/2)*min(time,1);
disp  = max(time,1)*0;

%xr = x(Ir);
%yr =  cos(angle).*y(Ir) + sin(angle).*z(Ir);
%zr = -sin(angle).*y(Ir) + cos(angle).*z(Ir);

xr =  cos(angle).*x(Ir) + sin(angle).*z(Ir);
yr = y(Ir) ;
zr = -sin(angle).*x(Ir) + cos(angle).*z(Ir) + disp;


xl = x(Il);
yl = y(Il);
zl = z(Il);



idx    = [...
          Il+xoffset; Ir+xoffset;...
          Il+yoffset; Ir+yoffset;...
          Il+zoffset; Ir+zoffset...
          ];
values = [ xl; xr;...
           yl; yr;...
           zl; zr...
           ];

BC = struct( 'idx', idx,...
  'values', values...
  );

end