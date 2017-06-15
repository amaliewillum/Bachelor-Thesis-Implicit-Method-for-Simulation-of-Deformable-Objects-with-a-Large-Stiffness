% Copyright 2011, Kenny Erleben
%clear all;
clear all;
close all;
% clc;

addpath('../fem');
addpath('./twist');
addpath('./twist2');
addpath('./twist3');
addpath('./bend');
addpath('./bend2');
addpath('./bend3');
addpath('./squeeze');
addpath('./stretch');


method = fem_method();

mesh_no = 3;
% scene  = twist_create_scene(mesh_no);   % Constant traction  
% scene  = twist2_create_scene(mesh_no);  % Increasing traction
% scene  = twist3_create_scene(mesh_no);  % No traction
scene  = bend_create_scene(mesh_no);    % Constant traction
%  scene  = bend2_create_scene(mesh_no);   % Increasing traction
% scene  = bend3_create_scene(mesh_no);   % No traction
% scene  = squeeze_create_scene(mesh_no);
% scene  = stretch_create_scene(mesh_no);

params = create_params([],'implicit_amalie'); 

% params = create_params([],'adaptive');
%--- COR method is based on Hookean material but FVM and FEM uses 
%--- a St. Venant-Kirchoff material. This gives COR a more stiff behavior/
%--- appearance. Thus, for animation one would most like use unrealistic
%--- small values of Yound modulos for the COR method.


params.E = 1e+2;
params.h_max = 1e-4;

params.T  = 2;

pinfo  = create_profile_info();
% pinfo.debug_level = 1;

material = 'lin';
% material = 'svk';

pdata = simulate( params, material, method, scene, pinfo ); 

postprocess(pinfo,pdata);