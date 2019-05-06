function [profile_data] = simulate( params, method, scene, profile_info )

%--- Profiling data members -----------------------------------------------
H  = [];           % Time step sizes
W  = [];           % Wall clock times of time steps
V  = [];           % Volume monitoring
KE = [];           % Kinetic energy monitoring
P  = [];           % Test point monitoring
C  = {};           % Convergence rate monitoring
MAX_X = [];        % Max x-position of mesh nodes
MIN_Z = [];        % Min z-position of mesh nodes
MEDIAN_XYZ = [];   % [x,y,z]-positions of mesh node closest to middle
CABLE_POINTS = []; % Cable positions [x,y,z]
COM = [];          % Center of mass
MEAN_POS = [];     % Mean position of mesh

profile = profile_info.record_convergence_rate; % Flag passed to time steppers 
conv    = [];                                   % Convergence rate array returned by time steppers (if used by them)

%--- Load mesh data -------------------------------------------------------
load(scene.meshfile);
T = double(T);
%--- Setup computational mesh and scene -----------------------------------
mesh   = method.create_mesh(T,X,Y,Z);
state  = method.create_state(mesh,params);
%--- Embed cables into the mesh, if this is a 'pull' simulation -----------
if isfield(scene, 'cable')
   cable_info = cable_embedding(T, X, Y, Z, scene.cable, params.alpha);
   cable_info.k = params.k;
   state.cable = cable_info;
end

delta = 1.5;
median_index  = find(X>-delta & X<delta & Y < delta & Y > -delta);
median_index = median_index(1);

if profile_info.record_test_point
  Q = profile_info.test_point;
  [SI, B] = coord2barycentric(Q,T,X,Y,Z);
  if SI==0
    display('error no enclsoing simplex found');
    return
  end
end

%-- If this is a pressure force simulation, add pressure parameters to
%-- state
if isfield(scene, 'pressure')
    state.pressure = scene.pressure;
    state.F = scene.F;
end


%-- Gravitational Computations

f_gravity = method.compute_gravity(T, X, Y, Z, params.rho);

clear T X Y Z;

%--- Get simulation parameters --------------------------------------------
T          = params.T;     % Total time to simulate
fps        = params.fps;   % Frames per second
max_frames = T*fps-1;      % Total number of frames
frame      = 1;            % Frame counter
kk         = 0;            % Number of time-steps done with out time-step changes (only used by adaptive time-stepping scheme)

%--- Start simulation loop ------------------------------------------------
while T > 0
  
  dt   = min(T,1./fps);
  t    = dt;
  
  while t > 0
    
    dh       = params.h_max;
    cur_time = params.T - T + (dt - t);
    
    state             = method.clear_forces( state );
    % Adding gravitational forces if part of simulation
    if params.gravity
       state              = method.add_gravity(f_gravity, state); 
    end
    % Adding cable forces if part of simulation
    if isfield(scene, 'cable')
        cable_force_info  = scene.create_cable_forces( cur_time, state, mesh );
        state             = method.add_cable_forces( state, cable_force_info );
        CABLE_POINTS = [CABLE_POINTS; cable_force_info.cable_positions]; % For plotting
    end
    % Adding pressure forces if part of simulation
    if isfield(scene, 'pressure')
       pressure_force_info = scene.create_pressure_forces( cur_time, state, mesh );
       state             = method.add_pressure_forces( state, pressure_force_info);

    end
    
    BC                = scene.create_boundary_conditions( cur_time, state, mesh );
    
    if profile_info.record_wall_clock
      t_start = tic;
    end
    
    switch lower(params.integration)
      %--------------------------------------------------------------------
      case 'implicit'
        
        [state, conv] = method.implicit_step(dh, mesh, state, params, BC, profile);
        
        %--------------------------------------------------------------------
      case 'adaptive'
        % Simple scheme we compare two steps with step-size half against
        % one full step. If difference is too large then we reduce the
        % full step-size and try again
        
        s1 = method.compute_elastic_forces(mesh, state, params);
        s2 = method.semi_implicit_step( dh, s1, BC, false );
        
        s3 = method.semi_implicit_step( dh/2, s1, BC, false );
        
        s3 = method.compute_elastic_forces(mesh, s3, params);
        s4 = method.semi_implicit_step( dh/2, s3, BC, false );
        
        p2 = [s2.x; s2.y; s2.z];
        p4 = [s4.x; s4.y; s4.z];
        
        max_diff  =  max( abs(p4 - p2) );
        if( max_diff > params.tol )
          dh = 0;
          params.h_max = params.h_max/2;
          kk = 0;
        else
          state = s4;
          kk = kk + 1;
          if(kk >  params.k_max)
            params.h_max = params.h_max*2;
            kk = 0;
          end
        end
        clear s1 s2 s3 s4 p2 p4 max_diff;
        %------------------------------------------------------------------
      case 'fixed'
        
        state = method.compute_elastic_forces(mesh, state, params);
        [state, conv ]= method.semi_implicit_step( dh, state, BC, profile );
        
        %------------------------------------------------------------------
    end
    
    if profile_info.record_wall_clock
      t_end = toc(t_start);
      W = [W; t_end];
    end
    if profile_info.record_time_steps
      H = [H; dh];
    end
    if profile_info.record_volume
      [v_min, v_max, v_tot] = profile_volume(mesh, state );
      V = [V; [v_min, v_max, v_tot] ];
    end
    if profile_info.record_volume
      KE = [KE; method.compute_kinetic_energy( state ) ];
    end
    if profile_info.record_test_point
      
      i = mesh.T( SI, 1 );
      j = mesh.T( SI, 2 );
      k = mesh.T( SI, 3 );
      m = mesh.T( SI, 4 );
      
      Pi= [state.x(i), state.y(i), state.z(i)];
      Pj= [state.x(j), state.y(j), state.z(j)];
      Pk= [state.x(k), state.y(k), state.z(k)];
      Pm= [state.x(m), state.y(m), state.z(m)];
      
      Q = B(1)*Pi + B(2)*Pj + B(3)*Pk + B(4)*Pm;
      
      P = [P; Q];
    end
    if profile_info.record_convergence_rate    
      if ~isempty(conv)
        C{ length(C) + 1 } = conv;
      end
    end
        
    t = t - dh;
    
  end
  
  T = T - dt;
  
  % Save maximum X point
  MAX_X = [MAX_X; max(state.x(:))];
  MIN_Z = [MIN_Z; min(state.z(:))];
  % Save median coordinates
  MEDIAN_XYZ = [ MEDIAN_XYZ; state.x(median_index), state.y(median_index), state.z(median_index)];
  % Compute and save COM and Mean Position of mesh
  total_volume = 0;
  com = [0,0,0];
  mean_pos = [0,0,0];
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
   
   vol = 1/6. * dot(a, cross(c, b));
   total_volume = total_volume + vol;
   mean_pos = mean_pos + mean([Pi'; Pj'; Pk'; Pm']);
   com  = com + mean([Pi'; Pj'; Pk'; Pm']) .* vol; 
  end
  COM = [COM; com / total_volume];
  MEAN_POS = [MEAN_POS; mean_pos / length(T(:, 1))];
  if profile_info.draw_images
    
    %--- Visualize the computed frame ---------------------------------------
    
    %%%% View 1 %%%%
    fh = figure('Visible','off');  % Should be used if run in commandline mode
    %fh = figure(100);
    clf;
    hold on
    meshplot(mesh, state, BC, profile_info.debug_level);
    axis([-100 100 -100 100 -100 100])
    grid on
    %view(0, 0)
    view(3)
    title( ['T = '  num2str( params.T - T ) ' [s]'] );
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('z [m]');
    hold off;
    if profile_info.save_images
      if(frame<=max_frames)
        filename = strcat(  profile_info.output_path,  profile_info.filename_prefix, num2str( frame, '%.5u'), 'view1');
        for f=1:length(profile_info.image_formats)
          format = profile_info.image_formats{f};
          print(fh, format, filename);
        end        
        %frame = frame + 1;
      end
    end
    %%%% View 2 %%%%
    fh = figure('Visible','off');  % Should be used if run in commandline mode
    %fh = figure(100);
    clf;
    hold on
    meshplot(mesh, state, BC, profile_info.debug_level);
    axis([-100 100 -100 100 -100 100])
    max_y = max(state.y);
    barycentricCoordinates = state.cable.W;
    indices = state.cable.indices;
    vertexPositions = [state.x(indices, :), state.y(indices, :), state.z(indices, :)];     %-- Positions of all vertices
    cablePositions = barycentricCoordinates * vertexPositions;                             %-- Re-interpolated points of the cable
    [xs,ys,zs] = sphere;
    % Place spheres representing cable nodes
    for c = 1:length(cablePositions(:, 1))
       %surf(xs+cablePositions(c, 1), ys+cablePositions(c,2), zs+cablePositions(c,3));
       surf(xs+cablePositions(c, 1), ys-max_y, zs+cablePositions(c,3));
    end
    % Plot lines between cable points
    plot3(cablePositions(:, 1), -ones(length(cablePositions(:, 1)), 1)*max_y, cablePositions(:, 3), 'LineWidth', 3, 'Color', 'Yellow')
    grid on
    view(0, 0)
    %view(3)
    title( ['T = '  num2str( params.T - T ) ' [s]'] );
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('z [m]');
    hold off;
    if profile_info.save_images
      if(frame<=max_frames)
        filename = strcat(  profile_info.output_path,  profile_info.filename_prefix, num2str( frame, '%.5u'), 'view2'  );
        for f=1:length(profile_info.image_formats)
          format = profile_info.image_formats{f};
          print(fh, format, filename);
        end        
        frame = frame + 1;
      end
    end
    close(fh);
        
  end
  
end

profile_data = struct(...
  'H',H,...
  'W',W,...
  'V',V,...
  'KE',KE,...
  'P',P,...
  'C',{C},...
  'MAX_X', MAX_X,...
  'MIN_Z', MIN_Z, ...
  'MEDIAN_XYZ', MEDIAN_XYZ, ...
  'CABLE_POINTS', CABLE_POINTS, ...
  'COM', COM, ...
  'MEAN_POS', MEAN_POS ...
  );

end

function [ v_min, v_max, v_tot ] = profile_volume( mesh, state )
% Copyright 2011, Kenny Erleben

%--- Get current spatial coordinates --------------------------------------
x = state.x;
y = state.y;
z = state.z;

%--- Get tetrahedron indices ----------------------------------------------
i = mesh.T(:,1);
j = mesh.T(:,2);
k = mesh.T(:,3);
m = mesh.T(:,4);

%--- Make vertex points ---------------------------------------------------
Pi = [ x(i),  y(i),  z(i) ]';
Pj = [ x(j),  y(j),  z(j) ]';
Pk = [ x(k),  y(k),  z(k) ]';
Pm = [ x(m),  y(m),  z(m) ]';

%--- Compute signed element volumes ---------------------------------------
V = dot( (Pm - Pi) , cross( (Pj - Pi), (Pk - Pi) ) ) ./ 6.0;

%--- Test if there exist any tetrahedron that has been flipped ------------
v_min = min(  V./mesh.V );
v_max = max(  V./mesh.V );
v_tot = sum(V(:))  /  sum(mesh.V(:));

end

function [ SI, B ] = coord2barycentric( Q, T, X, Y, Z )
% Copyright 2011, Kenny Erleben
  min_x = min(X);
  max_x = max(X);
  min_y = min(Y);
  max_y = max(Y);
  min_z = min(Z);
  max_z = max(Z);
  
  if Q(1) < min_x || Q(1) > max_x
    display('error bad test point');
    return
  end
  if Q(2) < min_y || Q(2) > max_y
    display('error bad test point');
    return
  end
  if Q(3) < min_z || Q(3) > max_z
    display('error bad test point');
    return
  end
  
  SI = 0;
  for t = 1:length(T(:,1))
    
    i = T(t,1);
    j = T(t,2);
    k = T(t,3);
    m = T(t,4);

    Pi = [ X(i); Y(i); Z(i) ];
    Pj = [ X(j); Y(j); Z(j) ];
    Pk = [ X(k); Y(k); Z(k) ];
    Pm = [ X(m); Y(m); Z(m) ];
    
    M = [  Pj-Pi, Pk-Pi, Pm-Pi ];
    
    b = (Q'-Pi);
    q = M \ b;
    
    B = [1 - q(1) - q(2) - q(3); q(1); q(2); q(3) ];
    
    if ( min(B) >= 0 ) && ( max(B) <= 1 )
      SI = t;
      break;
    end
  end  
end
