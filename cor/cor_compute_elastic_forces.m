function [ state ] = cor_compute_elastic_forces( mesh, state, params )
% Copyright 2011, Kenny Erleben

Re          = cor_compute_rotation_elements(mesh, state, params);
[Ke, fue]   = cor_compute_stiffness_warp( Re, state.Ke, state.fue, params );
state.K     = sparse( cor_assemble_global_matrix( mesh, Ke ) );
state.fu    = cor_assemble_global_vector( mesh, fue );

end