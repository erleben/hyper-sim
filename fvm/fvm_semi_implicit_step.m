function [ state, conv ] = fvm_semi_implicit_step(dt, state, bcon, ~)
% Copyright 2011, Kenny Erleben

%--- Assemble all vector and matrices needed ------------------------------
fe   = [ state.fex; state.fey;  state.fez ]; % Elastic forces
f    = [ state.fx;  state.fy;   state.fz  ]; % External forces
p    = [ state.x;   state.y;    state.z   ]; % Current spatial position
v    = [ state.vx;  state.vy;   state.vz  ]; % Current spatial velocity
w    = [1./state.M; 1./state.M; 1./state.M]; % Lumped inverted masses
C    = [state.C;    state.C;    state.C   ]; % Lumped viscous damping
fd   = - C .* v;                             % Damping forces

%--- Get information about boundary conditions ----------------------------
V      = length( state.x );      % Number of vertices in mesh
idx    = bcon.idx;               % Get indices of boundary conditions
values = bcon.values;            % Get values of the boundary conditions
free   = setdiff( 1:3*V, idx );  % Get indices of non-boundary conditions

%--- Apply boundary conditions --------------------------------------------
v(idx) = 0;
p(idx) = values;

%--- Do velocity update ---------------------------------------------------
a       = w(free).*(f(free) + fe(free) + fd(free) );

v(free) = v(free)  +  dt * a;

%--- Do position update ---------------------------------------------------
p(free)    = p(free) + dt*v(free);

%--- Store the updated values in state structure --------------------------
state.vx = v(1:V);
state.vy = v(V+1:2*V);
state.vz = v(2*V+1:end);
state.x  = p(1:V);
state.y  = p(V+1:2*V);
state.z  = p(2*V+1:end);

conv = [];   % There is no convergence rate to profile in the FVM method so we always return an empty array

end