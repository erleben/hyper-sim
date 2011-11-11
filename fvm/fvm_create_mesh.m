function [ mesh ] = fvm_create_mesh(T, X, Y, Z)
% Copyright 2011, Kenny Erleben

E  = length(T(:,1));  % Number of tetrahedrons

x0 = X;               % Material coordinates
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

invD0 = repmat(zeros(3,3),1,E);    % Inverse material edge matrix
G     = repmat(zeros(3,4),1,E);    % Material unit normals multiplied by face area; Gi = Ni*Ai

for e=1:E
  
  % Get tetrahedron indices
  i = T(e,1);
  j = T(e,2);
  k = T(e,3);
  m = T(e,4);
  
  % Get vertex coordinates
  Pi = [ x0(i),  y0(i),  z0(i) ]';
  Pj = [ x0(j),  y0(j),  z0(j) ]';
  Pk = [ x0(k),  y0(k),  z0(k) ]';
  Pm = [ x0(m),  y0(m),  z0(m) ]';
    
  % Define edge matrices
  D0 = [Pj-Pi,  Pk-Pi,  Pm-Pi ];
  invD0( 1:3, (3*(e-1)+1) :(3*e)) = inv(D0);
    
  % Compute outward unit normals multiplied by face area
  NAi = cross( (Pk - Pj), (Pm - Pj) )/2;
  NAj = cross( (Pi - Pk), (Pm - Pk) )/2;
  NAk = cross( (Pj - Pi), (Pm - Pi) )/2;
  NAm = cross( (Pk - Pi), (Pj - Pi) )/2;
  
  % Pre-compute material 'force' directions
  Gi = -(NAj + NAk + NAm)/3;
  Gj = -(NAi + NAk + NAm)/3;
  Gk = -(NAi + NAj + NAm)/3;
  Gm = -(NAi + NAj + NAk)/3;
  
  G(1:3, (4*(e-1)+1) :(4*e)) = [Gi Gj Gk Gm];
  
end

mesh = struct( 'x0', x0, 'y0', y0, 'z0', z0,...
  'V', V,...
  'T', T,...
  'G', G,...
  'invD0', invD0...
  );

end