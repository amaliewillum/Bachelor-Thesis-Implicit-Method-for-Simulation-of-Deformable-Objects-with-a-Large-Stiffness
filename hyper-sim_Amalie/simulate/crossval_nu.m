%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   CROSSVALIDATION OF NU AND H
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% nu x h: CV_nu_lin_real_squeeze_1st_h301
clear all 
close all

addpath('../fem');

load('CV_nu_lin_real_squeeze_1st_h301.mat')

% Get the field names of the structure.
fields = fieldnames(CV_nu_lin_real, '-full');
% Find out how many there are - for our loop.
numberOfFields = length(fields);
for f = 1 : numberOfFields
    thisField = fields{f};
    commandLine = sprintf('%s = CV_nu_lin_real.%s', thisField, thisField);
    eval(commandLine);
end
% Release temporary variables.
clear('f', 'thisField', 'numberOfFields');
clear('fields', 'commandLine');

delta_thing = [];
h_thing = [];

val = 1e-3;
for i = 1:length(nu)-1
   tmp = abs(delta(:,i) - val);
   [c idx] = min(tmp);
   h_thing = [h_thing; h(idx)];
   delta_thing = [delta_thing; delta(idx,i)];
end
%% Plot CV_E_lin_real_squeeze_1st_h301_delta.png
close all

figure('units','normalized','position',[.22,.3,.42,.5])
p1 = surf(nu,h,delta,log(delta));
hold on 
p2 = plot3([nu(1) nu(end)],[h_thing(1) h_thing(end)],[delta_thing(1) delta_thing(end)],'k-','LineWidth',3);
% p2 = plot3(nu,h_thing,delta_thing,'k-','LineWidth',3);
p3 = plot3(nu(1),h_thing(1),delta_thing(1),'b.','MarkerSize',40);
p4 = plot3(nu(end),h_thing(end),delta_thing(end),'r.','MarkerSize',40);

L = legend([p3,p4],['( h , $\nu$ ) = ( ' num2str(h_thing(1),'%.2e') ' , ' num2str(nu(1),'%.1e') ' )'   ],...
    ['( h , $\nu$ ) = ( ' num2str(h_thing(end),'%.2e') ' , ' num2str(nu(end-1),'%.1e') ' )' ],...
    'Location','NorthEast');
set(L,'Interpreter','latex','FontSize',12)

% set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

% title('Absolute Error','interpreter','latex','fontsize',18)
xlabel('Poissons Ratio $\nu$','interpreter','latex','fontsize',15)
ylabel('Step Size $h$','interpreter','latex','fontsize',15)
zlabel('Error $\delta$','interpreter','latex','fontsize',15)

alpha 0.8
axis([0,.5,h(end),h(1),1e-80,1e+20])
box on

shading interp 
lighting gouraud
camlight('headlight','infinite')
material dull

view(-50,30)
% view(-90,0)

col = makeColorMap([0 .2 .2],[.02 .58 .66],[.57 .89 .95],100);

colormap(col)
fname= '/users/amaliewillum/dropbox/2016_WILLUM/rapport/images';
filename = 'CV_nu_lin_real_1st_delta';


print( '-depsc' ,fullfile(fname, filename))

% print('-dpng' ,'CV_nu_lin_real_squeeze_1st_h301_delta')
% print('-dpng' ,'CV_nu90_lin_real_squeeze_1st_h301_delta')
hold off

%% u x h: CV_nu_lin_e_squeeze_1st_h301
clear all 
close all

addpath('../fem');

load('CV_nu_lin_e_squeeze_1st_h301.mat')

% Get the field names of the structure.
fields = fieldnames(CV_nu_lin_real, '-full');
% Find out how many there are - for our loop.
numberOfFields = length(fields);
for f = 1 : numberOfFields
    thisField = fields{f};
    commandLine = sprintf('%s = CV_nu_lin_real.%s', thisField, thisField);
    eval(commandLine);
end
% Release temporary variables.
clear('f', 'thisField', 'numberOfFields');
clear('fields', 'commandLine');

delta_thing = [];
h_thing = [];

val = 1e-3;
for i = 1:length(nu)-1
   tmp = abs(delta(:,i) - val);
   [c idx] = min(tmp);
   h_thing = [h_thing; h(idx)];
   delta_thing = [delta_thing; delta(idx,i)];
end
%% Plot CV_E_lin_e_squeeze_1st_h301_delta.png
close all

figure('units','normalized','position',[.22,.3,.42,.5])
p1 = surf(nu,h,delta,log(delta));
hold on 
p2 = plot3([nu(1) nu(end)],[h_thing(1) h_thing(end)],[delta_thing(1) delta_thing(end)],'k-','LineWidth',3);
% p2 = plot3(nu,h_thing,delta_thing,'k-','LineWidth',3);
p3 = plot3(nu(1),h_thing(1),delta_thing(1),'b.','MarkerSize',40);
p4 = plot3(nu(end),h_thing(end),delta_thing(end),'r.','MarkerSize',40);

L = legend([p3,p4],['( h , $\nu$ ) = ( ' num2str(h_thing(1),'%.2e') ' , ' num2str(nu(1),'%.1e') ' )'   ],...
    ['( h , $\nu$ ) = ( ' num2str(h_thing(end),'%.2e') ' , ' num2str(nu(end-1),'%.1e') ' )' ],...
    'Location','NorthEast');
set(L,'Interpreter','latex','FontSize',12)

% set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

% title('Absolute Error','interpreter','latex','fontsize',18)
xlabel('Poissons Ratio $\nu$','interpreter','latex','fontsize',15)
ylabel('Step Size $h$','interpreter','latex','fontsize',15)
zlabel('Error $\delta$','interpreter','latex','fontsize',15)

axis([0,.5,h(end),h(1),1e-80,1e+20])
alpha 0.8

box on

shading interp 
lighting gouraud
camlight('headlight','infinite')
material dull

view(-50,30)
% view(-90,0)
col = makeColorMap([0 .2 .2],[.02 .58 .66],[.57 .89 .95],100);

colormap(col)
fname= '/users/amaliewillum/dropbox/2016_WILLUM/rapport/images';
filename = 'CV_nu_lin_e_1st_delta';


print( '-dpng' ,fullfile(fname, filename))

% print('-dpng' ,'CV_nu_lin_e_squeeze_1st_h301_delta')
% print('-dpng' ,'CV_nu90_lin_e_squeeze_1st_h301_delta')
hold off

%% u x h: CV_nu_lin_rand_squeeze_1st_h301
clear all 
close all

addpath('../fem');

load('CV_nu_lin_rand_squeeze_1st_h301.mat')

% Get the field names of the structure.
fields = fieldnames(CV_nu_lin_real, '-full');
% Find out how many there are - for our loop.
numberOfFields = length(fields);
for f = 1 : numberOfFields
    thisField = fields{f};
    commandLine = sprintf('%s = CV_nu_lin_real.%s', thisField, thisField);
    eval(commandLine);
end
% Release temporary variables.
clear('f', 'thisField', 'numberOfFields');
clear('fields', 'commandLine');

delta_thing = [];
h_thing = [];

val = 1e-3;
for i = 1:length(nu)-1
   tmp = abs(delta(:,i) - val);
   [c idx] = min(tmp);
   h_thing = [h_thing; h(idx)];
   delta_thing = [delta_thing; delta(idx,i)];
end
%% Plot CV_E_lin_rand_squeeze_1st_h301_delta.png
close all

figure('units','normalized','position',[.22,.3,.42,.5])
p1 = surf(nu,h,delta,log(delta));
hold on 
p2 = plot3([nu(1) nu(end)],[h_thing(1) h_thing(end)],[delta_thing(1) delta_thing(end)],'k-','LineWidth',3);
% p2 = plot3(nu,h_thing,delta_thing,'k-','LineWidth',3);
p3 = plot3(nu(1),h_thing(1),delta_thing(1),'b.','MarkerSize',40);
p4 = plot3(nu(end),h_thing(end),delta_thing(end),'r.','MarkerSize',40);

L = legend([p3,p4],['( h , $\nu$ ) = ( ' num2str(h_thing(1),'%.2e') ' , ' num2str(nu(1),'%.1e') ' )'   ],...
    ['( h , $\nu$ ) = ( ' num2str(h_thing(end),'%.2e') ' , ' num2str(nu(end-1),'%.1e') ' )' ],...
    'Location','NorthEast');
set(L,'Interpreter','latex','FontSize',12)

% set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

% title('Absolute Error','interpreter','latex','fontsize',18)
xlabel('Poissons Ratio $\nu$','interpreter','latex','fontsize',15)
ylabel('Step Size $h$','interpreter','latex','fontsize',15)
zlabel('Error $\delta$','interpreter','latex','fontsize',15)

alpha 0.8

box on
axis([0,.5,h(end),h(1),1e-80,1e+20])
shading interp 
lighting gouraud
camlight('headlight','infinite')
material dull

view(-50,30)
% view(-90,0)

col = makeColorMap([0 .2 .2],[.02 .58 .66],[.57 .89 .95],100);

colormap(col)
fname= '/users/amaliewillum/dropbox/2016_WILLUM/rapport/images';
filename = 'CV_nu_lin_rand_1st_delta';


print( '-dpng' ,fullfile(fname, filename))

% print('-dpng' ,'CV_nu_lin_rand_squeeze_1st_h301_delta')
% print('-dpng' ,'CV_nu90_lin_rand_squeeze_1st_h301_delta')
hold off

%% u x h: CV_nu_lin_real_squeeze_2nd_h301 
close all

addpath('../fem');

load('CV_nu_lin_real_squeeze_2nd_h301.mat')

% Get the field names of the structure.
fields = fieldnames(CV_nu_lin_real, '-full');
% Find out how many there are - for our loop.
numberOfFields = length(fields);
for f = 1 : numberOfFields
    thisField = fields{f};
    commandLine = sprintf('%s = CV_nu_lin_real.%s', thisField, thisField);
    eval(commandLine);
end
% Release temporary variables.
clear('f', 'thisField', 'numberOfFields');
clear('fields', 'commandLine');

delta_thing = [];
h_thing = [];

val = 1e-3;
for i = 1:length(nu)-1
   tmp = abs(delta(:,i) - val);
   [c idx] = min(tmp);
   h_thing = [h_thing; h(idx)];
   delta_thing = [delta_thing; delta(idx,i)];
end
%% Plot CV_E_lin_real_squeeze_2nd_h301_delta.png
close all

figure('units','normalized','position',[.22,.3,.42,.5])
p1 = surf(nu,h,delta,log(delta));
hold on 
p2 = plot3([nu(1) nu(end)],[h_thing(1) h_thing(end)],[delta_thing(1) delta_thing(end)],'k-','LineWidth',3);
% p2 = plot3(nu,h_thing,delta_thing,'k-','LineWidth',3);
p3 = plot3(nu(1),h_thing(1),delta_thing(1),'b.','MarkerSize',40);
p4 = plot3(nu(end),h_thing(end),delta_thing(end),'r.','MarkerSize',40);

L = legend([p3,p4],['( h , $\nu$ ) = ( ' num2str(h_thing(1),'%.2e') ' , ' num2str(nu(1),'%.1e') ' )'   ],...
    ['( h , $\nu$ ) = ( ' num2str(h_thing(end),'%.2e') ' , ' num2str(nu(end-1),'%.1e') ' )' ],...
    'Location','NorthEast');
set(L,'Interpreter','latex','FontSize',12)

% set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

% title('Absolute Error','interpreter','latex','fontsize',18)
xlabel('Poissons Ratio $\nu$','interpreter','latex','fontsize',15)
ylabel('Step Size $h$','interpreter','latex','fontsize',15)
zlabel('Error $\delta$','interpreter','latex','fontsize',15)

alpha 0.8
axis([0,.5,h(end),h(1),1e-80,1e+20])
box on

shading interp 
lighting gouraud
camlight('headlight','infinite')
material dull

view(-50,30)
% view(-90,0)

col = makeColorMap([0 .2 .2],[.02 .58 .66],[.57 .89 .95],100);

colormap(col)
fname= '/users/amaliewillum/dropbox/2016_WILLUM/rapport/images';
filename = 'CV_nu_lin_real_2nd_delta';


print( '-dpng' ,fullfile(fname, filename))

% print('-dpng' ,'CV_nu_lin_real_squeeze_2nd_h301_delta')
% print('-dpng' ,'CV_nu90_lin_real_squeeze_2nd_h301_delta')
hold off

%% nu x h: CV_nu_svk_real_squeeze_1st_h301
close all

addpath('../fem');

load('CV_nu_svk_real_squeeze_1st_h301.mat')

% Get the field names of the structure.
fields = fieldnames(CV_nu_svk_real, '-full');
% Find out how many there are - for our loop.
numberOfFields = length(fields);
for f = 1 : numberOfFields
    thisField = fields{f};
    commandLine = sprintf('%s = CV_nu_svk_real.%s', thisField, thisField);
    eval(commandLine);
end
% Release temporary variables.
clear('f', 'thisField', 'numberOfFields');
clear('fields', 'commandLine');

delta_thing = [];
h_thing = [];

val = 1e-3;
for i = 1:length(nu)-1
   tmp = abs(delta(:,i) - val);
   [c idx] = min(tmp);
   h_thing = [h_thing; h(idx)];
   delta_thing = [delta_thing; delta(idx,i)];
end
%% Plot CV_E_svk_real_squeeze_1st_h301_delta.png
close all

figure('units','normalized','position',[.22,.3,.42,.5])
p1 = surf(nu,h,delta,log(delta));
hold on 
p2 = plot3([nu(1) nu(end)],[h_thing(1) h_thing(end)],[delta_thing(1) delta_thing(end)],'k-','LineWidth',3);
% p2 = plot3(nu,h_thing,delta_thing,'k-','LineWidth',3);
p3 = plot3(nu(1),h_thing(1),delta_thing(1),'b.','MarkerSize',40);
p4 = plot3(nu(end),h_thing(end),delta_thing(end),'r.','MarkerSize',40);

L = legend([p3,p4],['( h , $\nu$ ) = ( ' num2str(h_thing(1),'%.2e') ' , ' num2str(nu(1),'%.1e') ' )'   ],...
    ['( h , $\nu$ ) = ( ' num2str(h_thing(end),'%.2e') ' , ' num2str(nu(end-1),'%.1e') ' )' ],...
    'Location','NorthEast');
set(L,'Interpreter','latex','FontSize',12)

% set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

% title('Absolute Error','interpreter','latex','fontsize',18)
xlabel('Poissons Ratio $\nu$','interpreter','latex','fontsize',15)
ylabel('Step Size $h$','interpreter','latex','fontsize',15)
zlabel('Error $\delta$','interpreter','latex','fontsize',15)

alpha 0.8
axis([0,.5,h(end),h(1),1e-80,1e+20])
box on

shading interp 
lighting gouraud
camlight('headlight','infinite')
material dull

view(-50,30)
% view(-90,0)
col = makeColorMap([0 .2 .2],[.02 .58 .66],[.57 .89 .95],100);

colormap(col)
fname= '/users/amaliewillum/dropbox/2016_WILLUM/rapport/images';
filename = 'CV_nu_svk_real_1st_delta';


print( '-dpng' ,fullfile(fname, filename))

% print('-dpng' ,'CV_nu_svk_real_squeeze_1st_h301_delta')
% print('-dpng' ,'CV_nu90_svk_real_squeeze_1st_h301_delta')
hold off

%% u x h: CV_nu_svk_e_squeeze_1st_h301
close all

addpath('../fem');

load('CV_nu_svk_e_squeeze_1st_h301.mat')

% Get the field names of the structure.
fields = fieldnames(CV_nu_svk_real, '-full');
% Find out how many there are - for our loop.
numberOfFields = length(fields);
for f = 1 : numberOfFields
    thisField = fields{f};
    commandLine = sprintf('%s = CV_nu_svk_real.%s', thisField, thisField);
    eval(commandLine);
end
% Release temporary variables.
clear('f', 'thisField', 'numberOfFields');
clear('fields', 'commandLine');

delta_thing = [];
h_thing = [];

val = 1e-3;
for i = 1:length(nu)-1
   tmp = abs(delta(:,i) - val);
   [c idx] = min(tmp);
   h_thing = [h_thing; h(idx)];
   delta_thing = [delta_thing; delta(idx,i)];
end
%% Plot CV_E_svk_e_squeeze_1st_h301_delta.png
close all

figure('units','normalized','position',[.22,.3,.42,.5])
p1 = surf(nu,h,delta,log(delta));
hold on 
p2 = plot3([nu(1) nu(end)],[h_thing(1) h_thing(end)],[delta_thing(1) delta_thing(end)],'k-','LineWidth',3);
% p2 = plot3(nu,h_thing,delta_thing,'k-','LineWidth',3);
p3 = plot3(nu(1),h_thing(1),delta_thing(1),'b.','MarkerSize',40);
p4 = plot3(nu(end),h_thing(end),delta_thing(end),'r.','MarkerSize',40);

L = legend([p3,p4],['( h , $\nu$ ) = ( ' num2str(h_thing(1),'%.2e') ' , ' num2str(nu(1),'%.1e') ' )'   ],...
    ['( h , $\nu$ ) = ( ' num2str(h_thing(end),'%.2e') ' , ' num2str(nu(end-1),'%.1e') ' )' ],...
    'Location','NorthEast');
set(L,'Interpreter','latex','FontSize',12)

% set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

% title('Absolute Error','interpreter','latex','fontsize',18)
xlabel('Poissons Ratio $\nu$','interpreter','latex','fontsize',15)
ylabel('Step Size $h$','interpreter','latex','fontsize',15)
zlabel('Error $\delta$','interpreter','latex','fontsize',15)

alpha 0.8
axis([0,.5,h(end),h(1),1e-80,1e+20])
box on

shading interp 
lighting gouraud
camlight('headlight','infinite')
material dull

view(-50,30)
% view(-90,0)

col = makeColorMap([0 .2 .2],[.02 .58 .66],[.57 .89 .95],100);

colormap(col)
fname= '/users/amaliewillum/dropbox/2016_WILLUM/rapport/images';
filename = 'CV_nu_svk_e_1st_delta';


print( '-dpng' ,fullfile(fname, filename))

% print('-dpng' ,'CV_nu_svk_e_squeeze_1st_h301_delta')
% print('-dpng' ,'CV_nu90_svk_e_squeeze_1st_h301_delta')
hold off

%% u x h: CV_nu_svk_rand_squeeze_1st_h301
close all

addpath('../fem');

load('CV_nu_svk_rand_squeeze_1st_h301.mat')

% Get the field names of the structure.
fields = fieldnames(CV_nu_svk_real, '-full');
% Find out how many there are - for our loop.
numberOfFields = length(fields);
for f = 1 : numberOfFields
    thisField = fields{f};
    commandLine = sprintf('%s = CV_nu_svk_real.%s', thisField, thisField);
    eval(commandLine);
end
% Release temporary variables.
clear('f', 'thisField', 'numberOfFields');
clear('fields', 'commandLine');

delta_thing = [];
h_thing = [];

val = 1e-3;
for i = 1:length(nu)-1
   tmp = abs(delta(:,i) - val);
   [c idx] = min(tmp);
   h_thing = [h_thing; h(idx)];
   delta_thing = [delta_thing; delta(idx,i)];
end
%% Plot CV_E_svk_rand_squeeze_1st_h301_delta.png
close all

figure('units','normalized','position',[.22,.3,.42,.5])
p1 = surf(nu,h,delta,log(delta));
hold on 
p2 = plot3([nu(1) nu(end)],[h_thing(1) h_thing(end)],[delta_thing(1) delta_thing(end)],'k-','LineWidth',3);
% p2 = plot3(nu,h_thing,delta_thing,'k-','LineWidth',3);
p3 = plot3(nu(1),h_thing(1),delta_thing(1),'b.','MarkerSize',40);
p4 = plot3(nu(end),h_thing(end),delta_thing(end),'r.','MarkerSize',40);

L = legend([p3,p4],['( h , $\nu$ ) = ( ' num2str(h_thing(1),'%.2e') ' , ' num2str(nu(1),'%.1e') ' )'   ],...
    ['( h , $\nu$ ) = ( ' num2str(h_thing(end),'%.2e') ' , ' num2str(nu(end-1),'%.1e') ' )' ],...
    'Location','NorthEast');
set(L,'Interpreter','latex','FontSize',12)

% set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

% title('Absolute Error','interpreter','latex','fontsize',18)
xlabel('Poissons Ratio $\nu$','interpreter','latex','fontsize',15)
ylabel('Step Size $h$','interpreter','latex','fontsize',15)
zlabel('Error $\delta$','interpreter','latex','fontsize',15)

alpha 0.8
axis([0,.5,h(end),h(1),1e-80,1e+20])
box on

shading interp 
lighting gouraud
camlight('headlight','infinite')
material dull

view(-50,30)
% view(-90,0)

col = makeColorMap([0 .2 .2],[.02 .58 .66],[.57 .89 .95],100);

colormap(col)
fname= '/users/amaliewillum/dropbox/2016_WILLUM/rapport/images';
filename = 'CV_nu_svk_rand_1st_delta';


print( '-dpng' ,fullfile(fname, filename))

% print('-dpng' ,'CV_nu_svk_rand_squeeze_1st_h301_delta')
% print('-dpng' ,'CV_nu90_svk_rand_squeeze_1st_h301_delta')
hold off

%% u x h: CV_nu_svk_real_squeeze_2nd_h301 

close all

addpath('../fem');

load('CV_nu_svk_real_squeeze_2nd_h301.mat')

% Get the field names of the structure.
fields = fieldnames(CV_nu_svk_real, '-full');
% Find out how many there are - for our loop.
numberOfFields = length(fields);
for f = 1 : numberOfFields
    thisField = fields{f};
    commandLine = sprintf('%s = CV_nu_svk_real.%s', thisField, thisField);
    eval(commandLine);
end
% Release temporary variables.
clear('f', 'thisField', 'numberOfFields');
clear('fields', 'commandLine');

delta_thing = [];
h_thing = [];

val = 1e-3;
for i = 1:length(nu)-1
   tmp = abs(delta(:,i) - val);
   [c idx] = min(tmp);
   h_thing = [h_thing; h(idx)];
   delta_thing = [delta_thing; delta(idx,i)];
end
%% Plot CV_E_svk_real_squeeze_2nd_h301_delta.png
close all

figure('units','normalized','position',[.22,.3,.42,.5])
p1 = surf(nu,h,delta,log(delta));
hold on 
p2 = plot3([nu(1) nu(end)],[h_thing(1) h_thing(end)],[delta_thing(1) delta_thing(end)],'k-','LineWidth',3);
% p2 = plot3(nu,h_thing,delta_thing,'k-','LineWidth',3);
p3 = plot3(nu(1),h_thing(1),delta_thing(1),'b.','MarkerSize',40);
p4 = plot3(nu(end),h_thing(end),delta_thing(end),'r.','MarkerSize',40);

L = legend([p3,p4],['( h , $\nu$ ) = ( ' num2str(h_thing(1),'%.2e') ' , ' num2str(nu(1),'%.1e') ' )'   ],...
    ['( h , $\nu$ ) = ( ' num2str(h_thing(end),'%.2e') ' , ' num2str(nu(end-1),'%.1e') ' )' ],...
    'Location','NorthEast');
set(L,'Interpreter','latex','FontSize',12)

% set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

% title('Absolute Error','interpreter','latex','fontsize',18)
xlabel('Poissons Ratio $\nu$','interpreter','latex','fontsize',15)
ylabel('Step Size $h$','interpreter','latex','fontsize',15)
zlabel('Error $\delta$','interpreter','latex','fontsize',15)

alpha 0.8

box on

shading interp 
lighting gouraud
camlight('headlight','infinite')
material dull

view(-50,30)
% view(-90,0)
col = makeColorMap([0 .2 .2],[.02 .58 .66],[.57 .89 .95],100);

colormap(col)
fname= '/users/amaliewillum/dropbox/2016_WILLUM/rapport/images';
filename = 'CV_nu_svk_real_2nd_delta';

axis([0,.5,h(end),h(1),1e-80,1e+20])
print( '-dpng' ,fullfile(fname, filename))

% print('-dpng' ,'CV_nu_svk_real_squeeze_2nd_h301_delta')
% print('-dpng' ,'CV_nu90_svk_real_squeeze_2nd_h301_delta')
hold off





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Sammenligning af modellerne etc.
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all

addpath('../fem');

load('nu_delta.mat')



nr = 1;

h_nr = h(nr);
% nu_nr = nu(nr);

delta1 = lin_e_1st(nr,1:end-1); 
delta2 = lin_rand_1st(nr,1:end-1);
delta3 = lin_real_1st(nr,1:end-1);
delta4 = lin_real_2nd(nr,1:end-1);
delta5 = svk_e_1st(nr,1:end-1); 
delta6 = svk_rand_1st(nr,1:end-1);
delta7 = svk_real_1st(nr,1:end-1);
delta8 = svk_real_2nd(nr,1:end-1);

%%
ColorSet = fliplr(hsv(8));
% ColorSet = hsv(8);

figure('units','normalized','position',[.26,.3,.46,.5])
p3 = plot(nu(1:end-1),delta1,'color',ColorSet(8,:),'linewidth',2);
hold on
p2 = plot(nu(1:end-1),delta2,'color',ColorSet(7,:),'linewidth',2);
p1 = plot(nu(1:end-1),delta3,'color',ColorSet(6,:),'linewidth',2);
p4 = plot(nu(1:end-1),delta4,'color',ColorSet(5,:),'linewidth',2);
p5 = plot(nu(1:end-1),delta5,'color',ColorSet(4,:),'linewidth',2);
p6 = plot(nu(1:end-1),delta6,'color',ColorSet(3,:),'linewidth',2);
p7 = plot(nu(1:end-1),delta7,'color',ColorSet(2,:),'linewidth',2);
p8 = plot(nu(1:end-1),delta8,'color',ColorSet(1,:),'linewidth',2);

plot(nu(1:end), ones(100,1)*1e-1,'k-','linewidth',1.2)

l = legend([p3,p2,p1,p4,p5,p6,p7,p8],'lin$_e$ 1st','lin$_r$ 1st','lin$_{irl}$ 1st','lin$_{irl}$ 2nd','svk$_e$ 1st','svk$_r$ 1st','svk$_{irl}$ 1st','svk$_{irl}$ 2nd','location','southeastoutside');
set(l,'fontsize',14,'interpreter','latex');

title('2D view of Crossvalidation for $\nu$ and $h$ for $h \approx 10^{-2}$','interpreter','latex','fontsize',18)

ylabel('Error $\delta$','interpreter','latex','fontsize',20)
xlabel('Poisson''s Ratio $\nu$','interpreter','latex','fontsize',15)

set(gca, 'YScale', 'log')

YTick = [-5:3];
ys = cellstr(num2str(round(YTick(:)), '10^{%d}'));
set(gca,'YTickLabel',ys)


h1 = axes('Position', [.76 .65 .23 .23], 'Layer','top');
hold on
p1 = plot(nu(1:end-1),delta3,'color',ColorSet(6,:),'linewidth',2);
p4 = plot(nu(1:end-1),delta4,'color',ColorSet(5,:),'linewidth',2);
p7 = plot(nu(1:end-1),delta7,'color',ColorSet(2,:),'linewidth',2);
p8 = plot(nu(1:end-1),delta8,'color',ColorSet(1,:),'linewidth',2);
set(gca, 'YScale', 'log')
box on
annotation('textbox',[.38 .15 .081 .081]);
annotation('textbox',[.76 .65 .23 .23]);
annotation('line',[.38 + .081 .76], [.15 + .081 .65]);
axis(h1,'off',[0.22,0.275,10^(-2.5),10^(-1.9)])
title('Zoom','interpreter','latex')

h2 = axes('Position', [.74 .24 .23 .23], 'Layer','top');
hold on
set(gca, 'YScale', 'log')
box on
annotation('line',[.653 .78], [.315 .45]);
axis(h2,'off',[1,5,10^(-4.5),10^(-3)])
tt = title({'Feasible Area for', 'tol = $10^{-3}$ and','$h \approx 10^{-2}$'});
set(tt,'interpreter','latex')

h1 = axes('Position', [.76 .45 .23 .23], 'Layer','top');
hold on
set(gca, 'YScale', 'log')
box on
annotation('textbox',[.78 .45 .16 .13]);
axis(h1,'off',[1,5,10^(-4.5),10^(-3)])

fname= '/users/amaliewillum/dropbox/2016_WILLUM/rapport/images';
filename = 'delta_nu_contrast';


print( '-depsc' ,fullfile(fname, filename))

hold off
