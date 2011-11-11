function [ state ] = fvm_compute_elastic_forces( mesh, state, params )
% Copyright 2011, Kenny Erleben

cntT = length(mesh.T(:,1));

% Clear elastic force values from any previous time-steps
state.fex = zeros(size(state.fex));
state.fey = zeros(size(state.fey));
state.fez = zeros(size(state.fez));

% Convert Young Modulus and Poisson Ratio into Lame Coefficients
lambda =  (params.E*params.nu) / ((1+params.nu)*(1-2*params.nu));
mu     =  params.E/(2*(1+params.nu));

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
  
  % Create spatial edge matrix
  D  = [  Pj-Pi,    Pk-Pi,    Pm-Pi ];
  
  % Compute deformation gradient
  %
  %    D = F * D0
  %
  % note inv(D0) is precomputed
  invD0 = mesh.invD0(:,(3*(e-1)+1) :(3*e));
  Fe = D*invD0;
  
  % Compute Green strain
  Ee = (Fe'*Fe - eye(3,3))/2;
  
  % Use a Saint Venant Kirchhoff constitutive law to compute second
  % Piola--Kirchhoff stress tensor
  Se = lambda*trace(Ee)*eye(3,3) + 2*mu*Ee;
  
  % Compute first piola kirchhoff stress tensor
  Pe = Fe*Se;
    
  % Compute elastic forces from e'th element
  Ge  = mesh.G(:,(4*(e-1)+1) :(4*e));
  fei = Pe*( Ge(:,1) );
  fej = Pe*( Ge(:,2) );
  fek = Pe*( Ge(:,3) );
  fem = Pe*( Ge(:,4) );
  
  
  % Accumulate elastic force contributions to nodes
  state.fex(i) = state.fex(i) + fei(1);
  state.fey(i) = state.fey(i) + fei(2);
  state.fez(i) = state.fez(i) + fei(3);
  
  state.fex(j) = state.fex(j) + fej(1);
  state.fey(j) = state.fey(j) + fej(2);
  state.fez(j) = state.fez(j) + fej(3);
  
  state.fex(k) = state.fex(k) + fek(1);
  state.fey(k) = state.fey(k) + fek(2);
  state.fez(k) = state.fez(k) + fek(3);
  
  state.fex(m) = state.fex(m) + fem(1);
  state.fey(m) = state.fey(m) + fem(2);
  state.fez(m) = state.fez(m) + fem(3);
  
end

end