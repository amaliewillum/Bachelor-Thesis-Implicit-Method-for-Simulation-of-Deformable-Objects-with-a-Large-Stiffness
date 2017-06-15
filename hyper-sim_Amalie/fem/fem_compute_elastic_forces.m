function [ state ] = fem_compute_elastic_forces( mesh, state, params, material)
% Copyright 2011, Kenny Erleben

ke        = fem_compute_elastic_force_elements(mesh, state, params, material);
state.k   = fem_assemble_global_vector(mesh, ke);

end