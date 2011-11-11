function meshplot( mesh, state, bcon, debug_level)
% Copyright 2011, Kenny Erleben

%--- Visualize spatial mesh -----------------------------------------------
if debug_level<=1
  color = [0.5 0.5 0.75];
  T     = mesh.T;  
  TR    = TriRep(T, state.x, state.y, state.z);
  faces = freeBoundary(TR);
  h     = trimesh(faces, state.x, state.y, state.z);
%trisurf(tri, XF(:,1),XF(:,2),XF(:,3), 'FaceColor','cyan', 'FaceAlpha',
%0.8);
  set(h,'facecolor',color,'edgecolor','k');
end

if debug_level>0
  %--- Visualize boundary conditions --------------------------------------
  color = [0.75 0.5 0.5];
  x = state.x;
  y = state.y;
  z = state.z;
  s = 0.2;
  
  [mx, my, mz] = sphere;
  dx = mx*s;
  dy = my*s;
  dz = mz*s;
  
  idx = bcon.idx( bcon.idx <= length(state.x) );
  for i=1:length(idx)
    cx = x(idx(i));
    cy = y(idx(i));
    cz = z(idx(i));
    h = surf( dx + cx, dy + cy, dz + cz ); % Cool 3D drawing
    set(h,'facecolor',color,'edgecolor','k');
  end
  %--- Visualize Velocities -----------------------------------------------
  vx = state.vx;
  vy = state.vy;
  vz = state.vz;
  X = [x'; (vx+x)' ];
  Y = [y'; (vy+y)' ];
  Z = [z'; (vz+z)' ];
  plot3( X,Y,Z,'-r');
  %--- Visualize External Forces ------------------------------------------
  fx = state.fx;
  fy = state.fy;
  fz = state.fz;
  X = [x'; (fx+x)' ];
  Y = [y'; (fy+y)' ];
  Z = [z'; (fz+z)' ];
  plot3( X,Y,Z,'-g');
end
end
