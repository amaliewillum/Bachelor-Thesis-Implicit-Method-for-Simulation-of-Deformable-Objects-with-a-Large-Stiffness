% Copyright 2011, Kenny Erleben
clear all;
close all;
clc;

addpath('../fem');
addpath('./twist');
addpath('./bend');


%%
%--- FEM on TWIST ------------------------------------------------------
clear all;
close all;
clc;
method = fem_method();
params = create_params();
pinfo  = create_profile_info();
% params.T = 1.0;
for tst=1:8
  scene  = twist_create_scene(tst);
  pinfo.filename_prefix = strcat( 'fem_twist_m', num2str(tst), '_' );
  pdata = simulate( params, method, scene, pinfo ); 
  postprocess(pinfo,pdata);
end

%%
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

%%
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

%%
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
