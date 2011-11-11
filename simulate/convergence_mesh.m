% Copyright 2011, Kenny Erleben
clear all;
close all;
clc;

addpath('../fem');
addpath('../cor');
addpath('../fvm');
addpath('./twist');
addpath('./bend');

%--- FVM on TWIST ------------------------------------------------------
clear all;
close all;
clc;
method = fvm_method();
params = create_params();
pinfo  = create_profile_info();
for tst=1:8
  scene  = twist_create_scene(tst);
  pinfo.filename_prefix = strcat( 'fvm_twist_m', num2str(tst), '_' );
  pdata = simulate( params, method, scene, pinfo ); 
  postprocess(pinfo,pdata);
end

%--- FEM on TWIST ------------------------------------------------------
clear all;
close all;
clc;
method = fem_method();
params = create_params();
pinfo  = create_profile_info();
for tst=1:8
  scene  = twist_create_scene(tst);
  pinfo.filename_prefix = strcat( 'fem_twist_m', num2str(tst), '_' );
  pdata = simulate( params, method, scene, pinfo ); 
  postprocess(pinfo,pdata);
end

%--- COR on TWIST ------------------------------------------------------
clear all;
close all;
clc;
method = cor_method();
params = create_params();
params.E  = 5000; %--- COR needs low values to create interesting motion!
pinfo  = create_profile_info();
for tst=1:8
  scene  = twist_create_scene(tst);
  pinfo.filename_prefix = strcat( 'cor_twist_m', num2str(tst), '_' );
  pdata = simulate( params, method, scene, pinfo ); 
  postprocess(pinfo,pdata);
end

%--- FVM on BEND ------------------------------------------------------
clear all;
close all;
clc;
method = fvm_method();
params = create_params();
pinfo  = create_profile_info();
for tst=1:8
 scene  = bend_create_scene(tst);
 pinfo.filename_prefix = strcat( 'fvm_bend_m', num2str(tst), '_' );
 pdata = simulate( params, method, scene, pinfo ); 
 postprocess(pinfo,pdata);
end

%--- FEM on BEND ------------------------------------------------------
clear all;
close all;
clc;
method = fem_method();
params = create_params();
pinfo  = create_profile_info();
for tst=1:8
 scene  = bend_create_scene(tst);
 pinfo.filename_prefix = strcat( 'fem_bend_m', num2str(tst), '_' );
 pdata = simulate( params, method, scene, pinfo ); 
 postprocess(pinfo,pdata);
end

%--- COR on BEND ------------------------------------------------------
clear all;
close all;
clc;
method = cor_method();
params = create_params();
params.E  = 5000; %--- COR needs low values to create interesting motion!
pinfo  = create_profile_info();
for tst=1:8
 scene  = bend_create_scene(tst);
 pinfo.filename_prefix = strcat( 'cor_bend_m', num2str(tst), '_' );
 pdata = simulate( params, method, scene, pinfo ); 
 postprocess(pinfo,pdata);
end

%--- Implicit FEM on TWIST ------------------------------------------------------
clear all;
close all;
clc;
method = fem_method();
params = create_params([],'implicit');
params.h_max = 0.01;
pinfo  = create_profile_info();
for tst=1:5
  scene  = twist_create_scene(tst);
  pinfo.filename_prefix = strcat( 'ifem_twist_m', num2str(tst), '_' );
  pdata = simulate( params, method, scene, pinfo ); 
  postprocess(pinfo,pdata);
end

%--- Implicit FEM on BEND ------------------------------------------------------
clear all;
close all;
clc;
method = fem_method();
params = create_params([],'implicit');
params.h_max = 0.01;
pinfo  = create_profile_info();
for tst=1:5
  scene  = bend_create_scene(tst);
  pinfo.filename_prefix = strcat( 'ifem_bend_m', num2str(tst), '_' );
  pdata = simulate( params, method, scene, pinfo ); 
  postprocess(pinfo,pdata);
end

%--- Lumped FEM on TWIST ------------------------------------------------------
clear all;
close all;
clc;
method = fem_method();
params = create_params();
params.use_lumped = true;
pinfo  = create_profile_info();
for tst=1:6
  scene  = twist_create_scene(tst);
  pinfo.filename_prefix = strcat( 'lfem_twist_m', num2str(tst), '_' );
  pdata = simulate( params, method, scene, pinfo ); 
  postprocess(pinfo,pdata);
end

%--- Lumped FEM on BEND ------------------------------------------------------
clear all;
close all;
clc;
method = fem_method();
params = create_params();
params.use_lumped = true;
pinfo  = create_profile_info();
for tst=1:6
  scene  = bend_create_scene(tst);
  pinfo.filename_prefix = strcat( 'lfem_bend_m', num2str(tst), '_' );
  pdata = simulate( params, method, scene, pinfo ); 
  postprocess(pinfo,pdata);
end
