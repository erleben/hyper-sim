close all;
% Set-up a very small example (Two tetrahedra and two via points)
T = [1, 2, 3, 4;
     5, 6, 7, 8];
X = [1, 1, -1,-1, ...
     2, 2, -2, -2]';
Y = [1, -1, 1, -1, ...
     2, -2, 2, -2];
Z = [-1, 1, 1, -1, ...
     -2, 2, 2, -2]; 

Ps = [X(:), Y(:), Z(:)];
  
C = [0, 0, 0;
     0., 0., 0.5];

%-- Correct barycentric coordinates
val = [0.25, 0.25, 0.25, 0.25; ...
       0.1875, 0.3125, 0.3125, 0.1875];

%-- Compute barycentric coordinates and indices into T
cable = cable_embedding( T, X, Y, Z, C);
W =  cable.W;
indices = cable.indices;

%-- Print test result
LogicalStr = {'Failed', 'Success'};
sprintf('Mesh Coupling (Small Test) Result: %s', LogicalStr{isequal(cable.W * Ps(cable.indices, :), C)+1})
sprintf('Mesh Coupling (Small Test) Result: %s', LogicalStr{isequal((Ps(cable.indices, :)' * cable.W')', C)+1})

%-- Testing a larger example (Embedding a cable in a robot finger
%-- Embed cable C inside
C = [...
     59, 0.0, 0.0;...
     -59, 0.0, 0.0];
scene = pull_create_scene(1, '../', C);
load(scene.meshfile);
T = double(T); %-- This is needed due to T being exported from python
 % Visualize the via points as spheres on the mesh
figure;
tetramesh(T, [X, Y, Z])
[x,y,z] = sphere;
hold on;
for c = 1:length(C(:, 1))
   surf(x+C(c, 1), y+C(c,2), z+C(c,3));
end

plot3(C(:, 1), C(:,2), C(:, 3), 'LineWidth', 4, 'Color', 'Red')

xlabel('x'); ylabel('y'); zlabel('z');

cable = cable_embedding( T, X, Y, Z, C);
% Test that positions C can be recomputed:
P = zeros(length(X(:, 1)), 3);
for c = 1:length(X(:,1))
    P(c, :) = [X(c), Y(c), Z(c)];
end
C_reinterpolated = cable.W * P(cable.indices, :);
error = norm(C_reinterpolated - C);
W = cable.W;
sprintf('Total error between true positions and re-interpolated positions: %d', error)