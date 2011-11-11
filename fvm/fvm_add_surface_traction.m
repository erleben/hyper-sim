function [ state ] = fvm_add_surface_traction( state, traction_info )
% Copyright 2011, Kenny Erleben

%--- Get triangle and traction info ---------------------------------------
F    = traction_info.F;
if isempty(F)
  return
end

tx   = traction_info.tx;
ty   = traction_info.ty;
tz   = traction_info.tz;

x    = state.x;
y    = state.y;
z    = state.z;

%--- Get triangle indices -------------------------------------------------
i = F(:,1);
j = F(:,2);
k = F(:,3);

%--- Get vertex coordinates -----------------------------------------------
Pi = [ x(i), y(i), z(i) ]';
Pj = [ x(j), y(j), z(j) ]';
Pk = [ x(k), y(k), z(k) ]';

%--- Compute face areas ---------------------------------------------------
Avec =  cross( (Pj - Pi), (Pk - Pi) )  ./ 2.0 ;
A    =  sum( Avec.*Avec, 1).^(0.5);

%--- Solving integral of products of shape functions for isoparametric ----
%--- linear triangle element ----------------------------------------------
%
% The corners of the 'sub' control area wrt. the i'th vertex is 
%
%  p1 = [  1    0    0 ];    % 12/12      0      0
%  p2 = [3/4  1/4    0 ];    %  9/12     3/12    0
%  p3 = [3/4    0  1/4 ];    %  9/12      0     4/12
%  p4 = [1/3  1/3  1/3 ];    %  4/12     4/12   3/12
%  
%  where we used barycentric coordinates.
%  
%  The mid point of the control area is therefore
%  
%  p = (p1 + p2 + p3 + p4)/4;     %   (34 / 48     7/48   7/48) 
%  
%  Using a mid point rule we have
%  
%  lf(i) = 1/3 A (34/48 Ti + 7/48 Tj + 7/48 Tk )  
%    
% We find the load (nodal traction distribution)  matrix
%
%  L = 1 / 144  *  [ 34*I_3x3    7*I_3x3   7*I_3x3;
%                     7*I_3x3  34*I_3x3   7*I_3x3;
%                     7*I_3x3   7*I_3x3  34*I_3x3; ]
%
L = [  34*eye(3,3),  7*eye(3,3),  7*eye(3,3);...
        7*eye(3,3), 34*eye(3,3),  7*eye(3,3);...
        7*eye(3,3),  7*eye(3,3), 34*eye(3,3)...
     ] ./ 144;
   
%--- The spatial load force -----------------------------------------------
%
%     lf(f) = A(f) * L * [Ti; Tj; Tk]
%
%  where T is the nodal surface traction.
for f=1:length(F(:,1))
  
  i = F(f,1); 
  j = F(f,2); 
  k = F(f,3);
  
  Ti = [ tx(i); ty(i); tz(i) ];
  Tj = [ tx(j); ty(j); tz(j) ];
  Tk = [ tx(k); ty(k); tz(k) ];
  
  T = [ Ti; Tj; Tk ];
  
  LF = L * T .* A(f);
  
  state.fx(i) = state.fx(i) + LF(1);
  state.fy(i) = state.fy(i) + LF(2);
  state.fz(i) = state.fz(i) + LF(3);
  
  state.fx(j) = state.fx(j) + LF(4);
  state.fy(j) = state.fy(j) + LF(5);
  state.fz(j) = state.fz(j) + LF(6);
  
  state.fx(k) = state.fx(k) + LF(7);
  state.fy(k) = state.fy(k) + LF(8);
  state.fz(k) = state.fz(k) + LF(9);
  
end

end