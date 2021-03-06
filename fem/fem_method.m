function [ method ] = fem_method( )
% Copyright 2011, Kenny Erleben

create_mesh            = @fem_create_mesh;
create_state           = @fem_create_state;
implicit_step          = @fem_implicit_step;
compute_elastic_forces = @fem_compute_elastic_forces;
semi_implicit_step     = @fem_semi_implicit_step;
clear_forces           = @fem_clear_forces;
add_surface_traction   = @fem_add_surface_traction;
compute_kinetic_energy = @fem_compute_kinetic_energy;

method = struct(...
  'create_mesh', create_mesh,...
  'create_state', create_state,...
  'implicit_step', implicit_step,...
  'compute_elastic_forces', compute_elastic_forces,...
  'semi_implicit_step', semi_implicit_step,...
  'clear_forces', clear_forces,...
  'add_surface_traction', add_surface_traction,...
  'compute_kinetic_energy', compute_kinetic_energy...
  );

end