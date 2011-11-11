function [ traction_info ] = twist2_create_surface_traction_info( time, state, mesh )
% Copyright 2011, Kenny Erleben

%--- Find the triangle surfaces where we want to apply the surface --------
%--- traction   -----------------------------------------------------------
trep  = TriRep(mesh.T, mesh.x0, mesh.y0, mesh.z0);
ff    = freeBoundary(trep);
X     = mesh.x0;

x_limit = max( mesh.x0 ) - 0.01;

F = [];

for f=1:length(ff(:,1))
  
  i = ff(f,1);
  j = ff(f,2);
  k = ff(f,3);
  
  tst = (X(i) > x_limit) && (X(j) > x_limit) && (X(k) > x_limit);
  if tst  
    F = [F;  [i j k] ];
  end  
end

%--- We apply a constant traction onto all nodes on the surface -----------
magnitude =  4*min(time,1)*5000;

tx = zeros(  size(mesh.x0) );
ty = zeros(  size(mesh.y0) );
tz = zeros(  size(mesh.z0) );

mask = zeros(  size(mesh.z0) );
mask( F(:,1) ) = 1;
mask( F(:,2) ) = 1;
mask( F(:,3) ) = 1;

idx = find(  mask == 1 );
ty(idx) =   state.z(idx)*magnitude; 
tz(idx) =  -state.y(idx)*magnitude; 

%--- Bundle all info into one structure ----------------------------------
traction_info = struct(...
  'F', F,...
  'tx', tx,...
  'ty', ty,...
  'tz', tz...
  );

end