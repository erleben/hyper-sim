function [ cable ] = cable_embedding( T, X, Y, Z, C)
%CABLE_EMBEDDING -- Computes which tetrahedrons the points of C are
%                   contained in.
%
% INPUT:
%
%   T      - The tetrahedrons of a mesh
%   X      - The X coordinates.
%   Y      - The Y coordinates.
%   Z      - The Z coordinates.
%   C      - 3-by-N via point coordinates.
%
% OUTPUT:
%
% cable       - Upon return this struct will hold the indices into T
%               corresponding to the triangles holding the cable points,
%               and the barycentric coordinates of each point.
%
% Author: Max Kragballe, 2019

%% Helper function: Checks if p lies on the same side as the plane defined
%  by i,j,k, and m
    function ss = check_same_side(i,j,k,m,p)
        n = cross(j-i, k-i);
        ss = sign(dot(n, m-i)) == sign(dot(n, p-i));
    end
%% Helper function: Checks if p lies within tetrahedron i,j,k,m
    function inside = inside_tetrahedron(i,j,k,m,p)
       inside = check_same_side(i,j,k,m,p) && ...
                check_same_side(j,k,m,i,p) && ...
                check_same_side(k,m,i,j,p) && ...
                check_same_side(m,i,j,k,p);
    end

%% Helper function: Computes barycentric coordinates
function [ B ] = coord2barycentric( Pi, Pj, Pk, Pm, Pp )
    % Copyright 2011, Kenny Erleben
    % Modified 2019, by Max Kragballe
    M = [Pj - Pi, Pk- Pi, Pm - Pi];
    b = Pp - Pi;
    q = M \ b;
    B = [1-q(1)-q(2)-q(3);q(1);q(2);q(3)];
end
%% The embedding of the cables
%--- Barycentric weights, where
%       Wij := The weight of node j in terms of the position of via point i
W = zeros(length(C(:, 1)), length(X(:, 1)));
indices = [];
for c=1:length(C(:,1))
    Pp = C(c, :)';
    for t=1:length(T(:,1))
        i = T(t,1);
        j = T(t,2);
        k = T(t,3);
        m = T(t,4);
            
        Pi = [X(i); Y(i); Z(i)];
        Pj = [X(j); Y(j); Z(j)];
        Pk = [X(k); Y(k); Z(k)];
        Pm = [X(m); Y(m); Z(m)];
            
        if inside_tetrahedron(Pi, Pj, Pk, Pm, Pp)
            Wct = coord2barycentric(Pi, Pj, Pk, Pm, Pp)';
            W(c, i) = Wct(1);
            W(c, j) = Wct(2);
            W(c, k) = Wct(3);
            W(c, m) = Wct(4);
            
            indices = [indices; i; j; k; m];
            break; %-- A point can only lie within one tetrahedron
        end
    end
end
%-- The nodes of the mesh to apply forces to
indices = unique(indices);
%-- The weights in terms of only
W_p = zeros(length(C(:, 1)), length(indices(:, 1)));
for c=1:length(C(:,1))
    idx = find(W(c, :));    %-- Non-zero elements of the weight matrix
    i = idx(1); j = idx(2);  k = idx(3);  m = idx(4); % The four corner nodes
    W_p(c, find(indices==i)) = W(c, i); %-- Set the corresponding weight corresponding to W(c, i) in W_p 
    W_p(c, find(indices==j)) = W(c, j); %-- Set the corresponding weight corresponding to W(c, j) in W_p 
    W_p(c, find(indices==k)) = W(c, k); %-- Set the corresponding weight corresponding to W(c, k) in W_p 
    W_p(c, find(indices==m)) = W(c, m); %-- Set the corresponding weight corresponding to W(c, m) in W_p 
end
%-- Computing the initial lengths of the cables
L0s = flip(vecnorm(C(2:end, :) - C(1:end-1, :), 2, 2));
%-- The stifness of the cable
k = 1e4;

cable = struct(...
    'W', W_p, ...
    'L0s', L0s, ...
    'k', k,...
    'indices', indices ...
);

end

