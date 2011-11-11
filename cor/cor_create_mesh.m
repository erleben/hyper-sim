function [ mesh ] = cor_create_mesh(T, X, Y, Z)
% Copyright 2011, Kenny Erleben

cntT  = length(T(:,1));            % Number of tetrahedrons

invE0 = repmat(zeros(3,3),1,cntT); % Inverse material edge vector matrices

x0 = X;                            % Material coordinates
y0 = Y;
z0 = Z;

% Get tetrahedron indices
i = T(:,1);
j = T(:,2);
k = T(:,3);
m = T(:,4);

% Get vertex coordinates 
Pi = [ x0(i),  y0(i),  z0(i) ]';
Pj = [ x0(j),  y0(j),  z0(j) ]';
Pk = [ x0(k),  y0(k),  z0(k) ]';
Pm = [ x0(m),  y0(m),  z0(m) ]';

% Compute the tetrahedron element volumes.
%
% Using triple scalar product the tetrahedron volumes can be computed as
%
%   V   =  \frac{1}{6}  (P_m - P_i)  \cdot   ( P_j - P_i)  \times  ( P_k - P_i)
%
V = dot( (Pm - Pi) , cross( (Pj - Pi), (Pk - Pi) ) ) ./ 6.0;

for e=1:cntT
  
  % Get tetrahedron indices
  i = T(e,1);
  j = T(e,2);
  k = T(e,3);
  m = T(e,4);
  
  % Get material vertex coordinates
  P0i = [ x0(i);  y0(i);  z0(i) ];
  P0j = [ x0(j);  y0(j);  z0(j) ];
  P0k = [ x0(k);  y0(k);  z0(k) ];
  P0m = [ x0(m);  y0(m);  z0(m) ];
  
  % Define edge matrices
  E0 = [P0j-P0i,  P0m-P0i,  P0k-P0i ];
  
  invE0( 1:3, (3*(e-1)+1) :(3*e)) = inv(E0);
  
end


mesh = struct( 'x0', x0, 'y0', y0, 'z0', z0,...
  'T', T,...
  'V', V,...
  'invE0', invE0...
  );

end