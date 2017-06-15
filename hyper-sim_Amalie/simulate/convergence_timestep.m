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
h      = [0.01 0.001 0.0001 0.00001];
scene  = twist_create_scene();
method = fem_method();
params = create_params();
pinfo  = create_profile_info();
% params.T = 10;
% tic
for tst=1:length(h)
  params.h_max = h(tst);
  pinfo.filename_prefix = strcat( 'fem_twist_h', num2str(tst), '_' );
  pdata = simulate( params, method, scene, pinfo ); 
  postprocess(pinfo,pdata);
end
% toc
%%
%--- Implicit FEM on TWIST ------------------------------------------------------
clear all;
close all;
clc;
h      = [0.01 0.001 0.0001];
scene  = twist_create_scene();
method = fem_method();
params = create_params([],'implicit');
pinfo  = create_profile_info();
params.T = 10;
tic
for tst=1:length(h)
  params.h_max = h(tst);
  pinfo.filename_prefix = strcat( 'ifem_twist_h', num2str(tst), '_' );
  pdata = simulate( params, method, scene, pinfo ); 
  postprocess(pinfo,pdata);
end
toc
%%
%--- Lumped FEM on TWIST ------------------------------------------------------
clear all;
close all;
clc;
h      = [0.01 0.001 0.0001 0.00001];
scene  = twist_create_scene();
method = fem_method();
params = create_params();
params.use_lumped = true;
pinfo  = create_profile_info();
params.T = 10;
tic
for tst=1:length(h)
  params.h_max = h(tst);
  pinfo.filename_prefix = strcat( 'lfem_twist_h', num2str(tst), '_' );
  pdata = simulate( params, method, scene, pinfo ); 
  postprocess(pinfo,pdata);
end
toc

%%
%--- FEM on BEND ------------------------------------------------------
clear all;
close all;
clc;
h      = [0.01 0.001 0.0001 0.00001];
scene  = bend_create_scene();
method = fem_method();
params = create_params();
pinfo  = create_profile_info();
params.T = 10;
tic
for tst=1:length(h)
  params.h_max = h(tst);
  pinfo.filename_prefix = strcat( 'fem_bend_h', num2str(tst), '_' );
  pdata = simulate( params, method, scene, pinfo ); 
  postprocess(pinfo,pdata);
end
toc
%%
%--- Implicit FEM on BEND ------------------------------------------------------
clear all;
close all;
clc;
h      = [0.01 0.001 0.0001];
scene  = bend_create_scene();
method = fem_method();
params = create_params([],'implicit');
pinfo  = create_profile_info();
params.T = 10;
tic
for tst=1:length(h)
  params.h_max = h(tst);
  pinfo.filename_prefix = strcat( 'ifem_bend_h', num2str(tst), '_' );
  pdata = simulate( params, method, scene, pinfo ); 
  postprocess(pinfo,pdata);
end
toc

%--- Lumped FEM on BEND ------------------------------------------------------
clear all;
close all;
clc;
h      = [0.01 0.001 0.0001 0.00001];
scene  = bend_create_scene();
method = fem_method();
params = create_params();
params.use_lumped = true;
pinfo  = create_profile_info();
params.T = 10;
tic
for tst=1:length(h)
  params.h_max = h(tst);
  pinfo.filename_prefix = strcat( 'lfem_bend_h', num2str(tst), '_' );
  pdata = simulate( params, method, scene, pinfo ); 
  postprocess(pinfo,pdata);
end
toc
