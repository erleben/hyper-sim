function [state] = fem_add_gravity(forces, state)
% Update the z forces by gravitational pull
state.fz(:) = state.fz(:) + forces(:);
end

