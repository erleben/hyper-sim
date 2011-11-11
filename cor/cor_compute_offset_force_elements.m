function [ fue ] = cor_compute_offset_force_elements(Ke, mesh)
% Copyright 2011, Kenny Erleben

E = length(mesh.T(:,1));

fue   = zeros(12,E);          % Allocate array of element force offset vectors

for e=1:E
  
  % Get tetrahedron indices
  i = mesh.T(e,1);
  j = mesh.T(e,2);
  k = mesh.T(e,3);
  m = mesh.T(e,4);
  
  % Create element material coordinate vector
  p0 = [
    mesh.x0(i), mesh.y0(i), mesh.z0(i), ...
    mesh.x0(j), mesh.y0(j), mesh.z0(j), ...
    mesh.x0(k), mesh.y0(k), mesh.z0(k), ...
    mesh.x0(m), mesh.y0(m), mesh.z0(m), ...
    ]';
  
  % Store value
  fue(:,e) = Ke(:,(12*(e-1)+1):(12*e)) * p0;
  
end

end