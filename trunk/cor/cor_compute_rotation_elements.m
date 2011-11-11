function [ Re ] = cor_compute_rotation_elements( mesh, state, params )
% Copyright 2011, Kenny Erleben

cntT  = length(mesh.T(:,1));
Re    = repmat(eye(3,3),1,cntT); % Allocate array of element rotations (stiffness warping)

% From continuum mecahnics theory we know that a unique polar
% decomposition of the deformation gradient exist
%
%  F = R U
%
% The right Cauchy strain tensor then becomes by definition
%
%  C = F'*F = U*U
%
% From this it is clear that C is a symmetric matrix and that it
% is PD is F is non-singular.
%
% This means that an eigenvalue decomposition of C exist.. Now let D be
% the eigenvalues of C and let N be the matrix of eigenvectors then
%
%  C = N * D * N'
%
% and from this we have U
%
%  U = N * sqrt{D} * N'
%
% Knowing U we can now compute R as
%
%  R  = F inv(U)
%
% This of cause only works if det(F)>0
%
% Notes:
%
% Muller and Gross 2004 computes the deformation gradient and then
% applies polar decomposition method directly on the deformation
% gradient.
%
% Irving et. al. 2004 SCA paper has an SVD approach, and adds
% inversion handling.
%
% Doing a SVD of F results in
%
%  F = W S V^T =  W (V^T V) S V^T ) = (w V^T) (V S V^T)
%
% Here (W V^) = R is the orthogonal term we are looking for.
%
% Georgi and Westermann 2008 vriphys paper has a quaternion based
% energy minimization technique that could be used as an alternative
%
%
% Schmedding and Teschner 2008 CGI paper improves inversion handling but
% uses the eigenvalue decomposition strategy outlined above.
%


if params.warp
  
  for e=1:cntT
    
    % Get tetrahedron indices
    i = mesh.T(e,1);
    j = mesh.T(e,2);
    k = mesh.T(e,3);
    m = mesh.T(e,4);
    
    % Get spatial vertex coordinates
    Pi = [ state.x(i);  state.y(i);  state.z(i) ];
    Pj = [ state.x(j);  state.y(j);  state.z(j) ];
    Pk = [ state.x(k);  state.y(k);  state.z(k) ];
    Pm = [ state.x(m);  state.y(m);  state.z(m) ];
    
    % Define edge matrix
    E  = [  Pj-Pi,    Pm-Pi,    Pk-Pi ];
    
    % Compute deformation gradient
    %
    %    E = F * E0
    %
    % note inv(E0) is precomputed
    invE0 = mesh.invE0(:,(3*(e-1)+1) :(3*e));
    
    Fe = E*invE0;
    
    Ce = Fe' *Fe;
    [N, D] = eig(Ce);
    U = N*sqrt(D)*N';
    A = Fe * inv(U);  % Note this can be done faster as U^{-1} = N D^{-1/2} N^T
    
    % Orthonormalization, computational inexpensive and fast, but not very
    % ``correct''.
    %A = orth(Fe);
    
    Re( 1:3, (3*(e-1)+1) :(3*e)) = A;
    
  end
  
end