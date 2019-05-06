function [f_gravity] = fem_compute_gravity(T, X, Y, Z, rho)
g = 9.807;
f_gravity = zeros(length(X(:, 1)), 1);
for t=1:length(T(:, 1))
   i = T(t, 1);
   j = T(t, 2);
   k = T(t, 3);
   m = T(t, 4);
   
   Pi = [X(i); Y(i); Z(i)];
   Pj = [X(j); Y(j); Z(j)];
   Pk = [X(k); Y(k); Z(k)];
   Pm = [X(m); Y(m); Z(m)];
   
   a = Pi - Pj;
   b = Pk - Pj;
   c = Pm - Pj;
   
   V = 1/6. * dot(a, cross(c, b));
   mass = rho * V / 4;
   
   f_gravity(i) = f_gravity(i) - mass * g;
   f_gravity(j) = f_gravity(j) - mass * g;
   f_gravity(k) = f_gravity(k) - mass * g;
   f_gravity(m) = f_gravity(m) - mass * g;
end
max(f_gravity)
end

