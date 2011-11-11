function [ Ke ] = cor_compute_stiffness_elements(mesh, params)
% Copyright 2011, Kenny Erleben

cntT = length(mesh.T(:,1));
Ke   = zeros(12,cntT*12);       % Array of element stiffness matrices
E    = params.E;                % Young modulus
nu   = params.nu;               % Poisson ratio

% We use a linear material -- Hookean material law
%
% That means
%
%   \vec \sigma = \mat D \vec \epsilon
%
% where sigma is the Cauchy stress tensor and epsilon the small displacment
% strain tensor (Cauchy strain tensor).
%
% These are symmetric second order tensors which we have re-written the elements as
% vectors
%
% \mat D is the elasticity matrix and depends on the Lamé coefficients
%
D  =  E/((1+nu)*(1-2*nu)) * [ 
  (1-nu) nu     nu     0          0          0;
  nu     (1-nu) nu     0          0          0;
  nu     nu     (1-nu) 0          0          0;
  0      0      0      (1-2*nu)/2 0          0;
  0      0      0      0          (1-2*nu)/2 0 ;
  0      0      0      0          0          (1-2*nu)/2;
  ];


% We have the Cauchy equation
%
%   \vec r = \vec b + \nabla \cdot \sigma
%
% Here we only care about the last term, multiplying by a virtual
% displacment \delta u and taking the volume integral yields
%
%     int_v (\nabla \cdot \sigma) \cdot \delta u dV
%
% Knowing the dispacment field u one may find the strain as
%
%     \vec \epsilon = \mat S \vec u
%
% Where S is an appropriate differrential operator-


% Now compute element stiffness matrices
for e=1:cntT
  
  % Get tetrahedron indices
  i = mesh.T(e,1);
  j = mesh.T(e,2);
  k = mesh.T(e,3);
  m = mesh.T(e,4);
  
  % Get vertex coordinates
  Pi = [ mesh.x0(i);  mesh.y0(i);  mesh.z0(i) ];
  Pj = [ mesh.x0(j);  mesh.y0(j);  mesh.z0(j) ];
  Pk = [ mesh.x0(k);  mesh.y0(k);  mesh.z0(k) ];
  Pm = [ mesh.x0(m);  mesh.y0(m);  mesh.z0(m) ];
  
  % Compute spatial gradients of the barycentric coordinates
  %
  %   Using triple scalar product the tetrahedron volumes can be computed as
  %
  %   V   =  \frac{1}{6}  (P_m - P_i)  \cdot   ( P_j - P_i)  \times  ( P_k - P_i)
  %
  %   V_i =  \frac{1}{6}  (P   - P_j)  \cdot   ( P_m - P_j)  \times  ( P_k - P_j)
  %   V_j =  \frac{1}{6}  (P   - P_i)  \cdot   ( P_k - P_i)  \times  ( P_m - P_i)
  %   V_k =  \frac{1}{6}  (P   - P_i)  \cdot   ( P_m - P_i)  \times  ( P_j - P_i)
  %   V_m =  \frac{1}{6}  (P   - P_i)  \cdot   ( P_j - P_i)  \times  ( P_k - P_i)
  %
  %   And the bary centric coordinates are then defined as the weighted volumes
  %
  %   w_i = \frac{V_i}{V}
  %   w_j = \frac{V_j}{V}
  %   w_k = \frac{V_k}{V}
  %   w_m = \frac{V_m}{V}
  %
  
  dwi =  cross( Pm - Pj, Pk - Pj ) / mesh.V(e);
  dwj =  cross( Pk - Pi, Pm - Pi ) / mesh.V(e);
  dwk =  cross( Pm - Pi, Pj - Pi ) / mesh.V(e);
  dwm =  cross( Pj - Pi, Pk - Pi ) / mesh.V(e);
  
  % Compute the B matrix
  bi =  [
    dwi(1)   0         0;
    0        dwi(2)    0;
    0        0         dwi(3);
    dwi(2)   dwi(1)    0;
    dwi(3)   0         dwi(1);
    0        dwi(3)    dwi(2);
    ];
  
  bj =  [
    dwj(1)   0         0;
    0        dwj(2)    0;
    0        0         dwj(3);
    dwj(2)   dwj(1)    0;
    dwj(3)   0         dwj(1);
    0        dwj(3)    dwj(2);
    ];
  
  bk =  [
    dwk(1)   0         0;
    0        dwk(2)    0;
    0        0         dwk(3);
    dwk(2)   dwk(1)    0;
    dwk(3)   0         dwk(1);
    0        dwk(3)    dwk(2);
    ];
  
  bm =  [
    dwm(1)   0         0;
    0        dwm(2)    0;
    0        0         dwm(3);
    dwm(2)   dwm(1)    0;
    dwm(3)   0         dwm(1);
    0        dwm(3)    dwm(2);
    ];
  
  B  =  [ bi bj bk bm];
  
  % Compute element stiffness matrix and store it in Ke array
  
  tmp = (B' * D * B) * mesh.V(e);
        
  Ke( 1:12, (12*(e-1)+1) :(12*e)) =  tmp;
  
end

end