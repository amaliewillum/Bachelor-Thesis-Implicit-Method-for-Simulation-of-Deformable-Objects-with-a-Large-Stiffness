%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   CROSSVALIDATION OF L AND H
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% L x h: CV_L_lin_real_1st
clear all 
close all

addpath('../fem');

load('CV_L_lin_real_1st.mat')

% Get the field names of the structure.
fields = fieldnames(CV_L_lin_real, '-full');
% Find out how many there are - for our loop.
numberOfFields = length(fields);
for f = 1 : numberOfFields
    thisField = fields{f};
    commandLine = sprintf('%s = CV_L_lin_real.%s', thisField, thisField);
    eval(commandLine);
end
% Release temporary variables.
clear('f', 'thisField', 'numberOfFields');
clear('fields', 'commandLine');

delta_thing = [];
h_thing = [];

L = 1:9;

val = 1e-3;
for i = 1:length(L)-1
   tmp = abs(delta(:,i) - val);
   [c idx] = min(tmp);
   h_thing = [h_thing; h(idx)];
   delta_thing = [delta_thing; delta(idx,i)];
end
%% Plot CV_L_lin_real_1st_delta.png
close all

figure('units','normalized','position',[.22,.3,.42,.5])
p1 = surf(L',h,delta,log(delta));
hold on 
% p2 = plot3([L(1) L(end)],[h_thing(1) h_thing(end)],[delta_thing(1) delta_thing(end)],'k-','LineWidth',5);
p2 = plot3(L(1:end-1)',h_thing,delta_thing,'k-','LineWidth',3);
p3 = plot3(L(1),h_thing(1),delta_thing(1),'b.','MarkerSize',40);
p4 = plot3(L(end-1),h_thing(end),delta_thing(end),'r.','MarkerSize',40);

l = legend([p3,p4],['( h , $L$ ) = ( ' num2str(h_thing(1),'%.2e') ' , ' num2str(L(1),'%.1e') ' )'   ],...
    ['( h , $L$ ) = ( ' num2str(h_thing(end),'%.2e') ' , ' num2str(L(end-1),'%.1e') ' )' ],...
    'Location','NorthEast');
set(l,'Interpreter','latex','FontSize',12)

% set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

% title('Absolute Error','interpreter','latex','fontsize',18)
xlabel('Mesh Size $L$','interpreter','latex','fontsize',15)
ylabel('Step Size $h$','interpreter','latex','fontsize',15)
zlabel('Error $\delta$','interpreter','latex','fontsize',15)

col = makeColorMap([ .07840 .6 .5961],[0 1 .5],[0.76 1 0.76],100);

colormap(col)

alpha 0.8

box on
axis([1,8,h(end),h(1),1e-40,1e+10])

shading interp 
lighting gouraud
camlight('headlight','infinite')
material dull

view(-50,30)
% view(-90,0)

fname= '/users/amaliewillum/dropbox/2016_WILLUM/rapport/images';
filename = 'CV_L_lin_real_1st_delta';


print( '-dpng' ,fullfile(fname, filename))

% print('-dpng' ,'CV_L_lin_real_1st_delta')
% print('-dpng' ,'CV_L_lin_real_1st_delta')
hold off

%% L x h: CV_L_lin_real_2nd
clear all 
close all

addpath('../fem');

load('CV_L_lin_real_2nd.mat')

% Get the field names of the structure.
fields = fieldnames(CV_L_lin_real, '-full');
% Find out how many there are - for our loop.
numberOfFields = length(fields);
for f = 1 : numberOfFields
    thisField = fields{f};
    commandLine = sprintf('%s = CV_L_lin_real.%s', thisField, thisField);
    eval(commandLine);
end
% Release temporary variables.
clear('f', 'thisField', 'numberOfFields');
clear('fields', 'commandLine');

delta_thing = [];
h_thing = [];

L = 1:9;

val = 1e-3;
for i = 1:length(L)-1
   tmp = abs(delta(:,i) - val);
   [c idx] = min(tmp);
   h_thing = [h_thing; h(idx)];
   delta_thing = [delta_thing; delta(idx,i)];
end
%% Plot CV_L_lin_real_2nd_delta.png
close all

figure('units','normalized','position',[.22,.3,.42,.5])
p1 = surf(L',h,delta,log(delta));
hold on 
% p2 = plot3([L(1) L(end)],[h_thing(1) h_thing(end)],[delta_thing(1) delta_thing(end)],'k-','LineWidth',5);
p2 = plot3(L(1:end-1)',h_thing,delta_thing,'k-','LineWidth',3);
p3 = plot3(L(1),h_thing(1),delta_thing(1),'b.','MarkerSize',40);
p4 = plot3(L(end-1),h_thing(end),delta_thing(end),'r.','MarkerSize',40);

l = legend([p3,p4],['( h , $L$ ) = ( ' num2str(h_thing(1),'%.2e') ' , ' num2str(L(1),'%.1e') ' )'   ],...
    ['( h , $L$ ) = ( ' num2str(h_thing(end),'%.2e') ' , ' num2str(L(end-1),'%.1e') ' )' ],...
    'Location','NorthEast');
set(l,'Interpreter','latex','FontSize',12)

% set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

% title('Absolute Error','interpreter','latex','fontsize',18)
xlabel('Mesh Size $L$','interpreter','latex','fontsize',15)
ylabel('Step Size $h$','interpreter','latex','fontsize',15)
zlabel('Error $\delta$','interpreter','latex','fontsize',15)

col = makeColorMap([ .07840 .6 .5961],[0 1 .5],[0.76 1 0.76],100);

colormap(col)

alpha 0.8

box on
axis([1,8,h(end),h(1),1e-40,1e+10])

shading interp 
lighting gouraud
camlight('headlight','infinite')
material dull

view(-50,30)
% view(-90,0)

col = makeColorMap([ .07840 .6 .5961],[0 1 .5],[0.76 1 0.76],100);

colormap(col)
fname= '/users/amaliewillum/dropbox/2016_WILLUM/rapport/images';
filename = 'CV_L_lin_real_2nd_delta';


print( '-dpng' ,fullfile(fname, filename))

% print('-dpng' ,'CV_L_lin_real_2nd_delta')
% print('-dpng' ,'CV_L90_lin_real_2nd_delta')
hold off

%% L x h: CV_L_lin_rand_1st
clear all 
close all

addpath('../fem');

load('CV_L_lin_rand_1st.mat')

% Get the field names of the structure.
fields = fieldnames(CV_L_lin_rand, '-full');
% Find out how many there are - for our loop.
numberOfFields = length(fields);
for f = 1 : numberOfFields
    thisField = fields{f};
    commandLine = sprintf('%s = CV_L_lin_rand.%s', thisField, thisField);
    eval(commandLine);
end
% Release temporary variables.
clear('f', 'thisField', 'numberOfFields');
clear('fields', 'commandLine');

delta_thing = [];
h_thing = [];

L = 1:9;

val = 1e-3;
for i = 1:length(L)-1
   tmp = abs(delta(:,i) - val);
   [c idx] = min(tmp);
   h_thing = [h_thing; h(idx)];
   delta_thing = [delta_thing; delta(idx,i)];
end
%% Plot CV_L_lin_rand_1st_delta.png
close all

figure('units','normalized','position',[.22,.3,.42,.5])
p1 = surf(L',h,delta,log(delta));
hold on 
% p2 = plot3([L(1) L(end)],[h_thing(1) h_thing(end)],[delta_thing(1) delta_thing(end)],'k-','LineWidth',5);
p2 = plot3(L(1:end-1)',h_thing,delta_thing,'k-','LineWidth',3);
p3 = plot3(L(1),h_thing(1),delta_thing(1),'b.','MarkerSize',40);
p4 = plot3(L(end-1),h_thing(end),delta_thing(end),'r.','MarkerSize',40);

l = legend([p3,p4],['( h , $L$ ) = ( ' num2str(h_thing(1),'%.2e') ' , ' num2str(L(1),'%.1e') ' )'   ],...
    ['( h , $L$ ) = ( ' num2str(h_thing(end),'%.2e') ' , ' num2str(L(end-1),'%.1e') ' )' ],...
    'Location','NorthEast');
set(l,'Interpreter','latex','FontSize',12)

% set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

% title('Absolute Error','interpreter','latex','fontsize',18)
xlabel('Mesh Size $L$','interpreter','latex','fontsize',15)
ylabel('Step Size $h$','interpreter','latex','fontsize',15)
zlabel('Error $\delta$','interpreter','latex','fontsize',15)

col = makeColorMap([ .07840 .6 .5961],[0 1 .5],[0.76 1 0.76],100);

colormap(col)

alpha 0.8

box on
axis([1,8,h(end),h(1),1e-40,1e+10])
shading interp 
lighting gouraud
camlight('headlight','infinite')
material dull

view(-50,30)
% view(-90,0)

col = makeColorMap([ .07840 .6 .5961],[0 1 .5],[0.76 1 0.76],100);

colormap(col)
fname= '/users/amaliewillum/dropbox/2016_WILLUM/rapport/images';
filename = 'CV_L_lin_rand_1st_delta';


print( '-dpng' ,fullfile(fname, filename))

% print('-dpng' ,'CV_L_lin_rand_1st_delta')
% print('-dpng' ,'CV_L90_lin_rand_1st_delta')
hold off


%% L x h: CV_L_lin_e_1st
clear all 
close all

addpath('../fem');

load('CV_L_lin_e_1st.mat')

% Get the field names of the structure.
fields = fieldnames(CV_L_lin_e, '-full');
% Find out how many there are - for our loop.
numberOfFields = length(fields);
for f = 1 : numberOfFields
    thisField = fields{f};
    commandLine = sprintf('%s = CV_L_lin_e.%s', thisField, thisField);
    eval(commandLine);
end
% Release temporary variables.
clear('f', 'thisField', 'numberOfFields');
clear('fields', 'commandLine');

delta_thing = [];
h_thing = [];

L = 1:9;

val = 1e-3;
for i = 1:length(L)-1
   tmp = abs(delta(:,i) - val);
   [c idx] = min(tmp);
   h_thing = [h_thing; h(idx)];
   delta_thing = [delta_thing; delta(idx,i)];
end
%% Plot CV_L_lin_e_1st_delta.png
close all

figure('units','normalized','position',[.22,.3,.42,.5])
p1 = surf(L',h,delta,log(delta));
hold on 
% p2 = plot3([L(1) L(end)],[h_thing(1) h_thing(end)],[delta_thing(1) delta_thing(end)],'k-','LineWidth',5);
p2 = plot3(L(1:end-1)',h_thing,delta_thing,'k-','LineWidth',3);
p3 = plot3(L(1),h_thing(1),delta_thing(1),'b.','MarkerSize',40);
p4 = plot3(L(end-1),h_thing(end),delta_thing(end),'r.','MarkerSize',40);

l = legend([p3,p4],['( h , $L$ ) = ( ' num2str(h_thing(1),'%.2e') ' , ' num2str(L(1),'%.1e') ' )'   ],...
    ['( h , $L$ ) = ( ' num2str(h_thing(end),'%.2e') ' , ' num2str(L(end-1),'%.1e') ' )' ],...
    'Location','NorthEast');
set(l,'Interpreter','latex','FontSize',12)

% set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

% title('Absolute Error','interpreter','latex','fontsize',18)
xlabel('Mesh Size $L$','interpreter','latex','fontsize',15)
ylabel('Step Size $h$','interpreter','latex','fontsize',15)
zlabel('Error $\delta$','interpreter','latex','fontsize',15)

col = makeColorMap([ .07840 .6 .5961],[0 1 .5],[0.76 1 0.76],100);

colormap(col)

alpha 0.8

box on
axis([1,8,h(end),h(1),1e-40,1e+10])

shading interp 
lighting gouraud
camlight('headlight','infinite')
material dull

view(-50,30)
% view(-90,0)

col = makeColorMap([ .07840 .6 .5961],[0 1 .5],[0.76 1 0.76],100);

colormap(col)
fname= '/users/amaliewillum/dropbox/2016_WILLUM/rapport/images';
filename = 'CV_L_lin_e_1st_delta';


print( '-dpng' ,fullfile(fname, filename))

% print('-dpng' ,'CV_L_lin_e_1st_delta')
% print('-dpng' ,'CV_L90_lin_e_1st_delta')
hold off


%% L x h: CV_L_lin_real_1st
clear all 
close all

addpath('../fem');

load('CV_L_svk_real_1st.mat')

% Get the field names of the structure.
fields = fieldnames(CV_L_svk_real, '-full');
% Find out how many there are - for our loop.
numberOfFields = length(fields);
for f = 1 : numberOfFields
    thisField = fields{f};
    commandLine = sprintf('%s = CV_L_svk_real.%s', thisField, thisField);
    eval(commandLine);
end
% Release temporary variables.
clear('f', 'thisField', 'numberOfFields');
clear('fields', 'commandLine');

delta_thing = [];
h_thing = [];

L = 1:9;

val = 1e-3;
for i = 1:length(L)-1
   tmp = abs(delta(:,i) - val);
   [c idx] = min(tmp);
   h_thing = [h_thing; h(idx)];
   delta_thing = [delta_thing; delta(idx,i)];
end
%% Plot CV_L_svk_real_1st_delta.png
close all

figure('units','normalized','position',[.22,.3,.42,.5])
p1 = surf(L',h,delta,log(delta));
hold on 
% p2 = plot3([L(1) L(end)],[h_thing(1) h_thing(end)],[delta_thing(1) delta_thing(end)],'k-','LineWidth',5);
p2 = plot3(L(1:end-1)',h_thing,delta_thing,'k-','LineWidth',3);
p3 = plot3(L(1),h_thing(1),delta_thing(1),'b.','MarkerSize',40);
p4 = plot3(L(end-1),h_thing(end),delta_thing(end),'r.','MarkerSize',40);

l = legend([p3,p4],['( h , $L$ ) = ( ' num2str(h_thing(1),'%.2e') ' , ' num2str(L(1),'%.1e') ' )'   ],...
    ['( h , $L$ ) = ( ' num2str(h_thing(end),'%.2e') ' , ' num2str(L(end-1),'%.1e') ' )' ],...
    'Location','NorthEast');
set(l,'Interpreter','latex','FontSize',12)

% set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

% title('Absolute Error','interpreter','latex','fontsize',18)
xlabel('Mesh Size $L$','interpreter','latex','fontsize',15)
ylabel('Step Size $h$','interpreter','latex','fontsize',15)
zlabel('Error $\delta$','interpreter','latex','fontsize',15)

col = makeColorMap([ .07840 .6 .5961],[0 1 .5],[0.76 1 0.76],100);

colormap(col)

alpha 0.8

box on
axis([1,8,h(end),h(1),1e-40,1e+10])
shading interp 
lighting gouraud
camlight('headlight','infinite')
material dull

view(-50,30)
% view(-90,0)

col = makeColorMap([ .07840 .6 .5961],[0 1 .5],[0.76 1 0.76],100);

colormap(col)
fname= '/users/amaliewillum/dropbox/2016_WILLUM/rapport/images';
filename = 'CV_L_svk_real_1st_delta';


print( '-dpng' ,fullfile(fname, filename))

% print('-dpng' ,'CV_L_svk_real_1st_delta')
% print('-dpng' ,'CV_L_svk_real_1st_delta')
hold off


%% L x h: CV_L_lin_real_2nd
clear all 
close all

addpath('../fem');

load('CV_L_svk_real_2nd.mat')

% Get the field names of the structure.
fields = fieldnames(CV_L_svk_real, '-full');
% Find out how many there are - for our loop.
numberOfFields = length(fields);
for f = 1 : numberOfFields
    thisField = fields{f};
    commandLine = sprintf('%s = CV_L_svk_real.%s', thisField, thisField);
    eval(commandLine);
end
% Release temporary variables.
clear('f', 'thisField', 'numberOfFields');
clear('fields', 'commandLine');

delta_thing = [];
h_thing = [];

L = 1:9;

val = 1e-3;
for i = 1:length(L)-1
   tmp = abs(delta(:,i) - val);
   [c idx] = min(tmp);
   h_thing = [h_thing; h(idx)];
   delta_thing = [delta_thing; delta(idx,i)];
end
%% Plot CV_L_svk_real_2nd_delta.png
close all

figure('units','normalized','position',[.22,.3,.42,.5])
p1 = surf(L',h,delta,log(delta));
hold on 
% p2 = plot3([L(1) L(end)],[h_thing(1) h_thing(end)],[delta_thing(1) delta_thing(end)],'k-','LineWidth',5);
p2 = plot3(L(1:end-1)',h_thing,delta_thing,'k-','LineWidth',3);
p3 = plot3(L(1),h_thing(1),delta_thing(1),'b.','MarkerSize',40);
p4 = plot3(L(end-1),h_thing(end),delta_thing(end),'r.','MarkerSize',40);

l = legend([p3,p4],['( h , $L$ ) = ( ' num2str(h_thing(1),'%.2e') ' , ' num2str(L(1),'%.1e') ' )'   ],...
    ['( h , $L$ ) = ( ' num2str(h_thing(end),'%.2e') ' , ' num2str(L(end-1),'%.1e') ' )' ],...
    'Location','NorthEast');
set(l,'Interpreter','latex','FontSize',12)

% set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

% title('Absolute Error','interpreter','latex','fontsize',18)
xlabel('Mesh Size $L$','interpreter','latex','fontsize',15)
ylabel('Step Size $h$','interpreter','latex','fontsize',15)
zlabel('Error $\delta$','interpreter','latex','fontsize',15)

col = makeColorMap([ .07840 .6 .5961],[0 1 .5],[0.76 1 0.76],100);

colormap(col)

alpha 0.8

box on
axis([1,8,h(end),h(1),1e-40,1e+10])

shading interp 
lighting gouraud
camlight('headlight','infinite')
material dull

view(-50,30)
% view(-90,0)

col = makeColorMap([ .07840 .6 .5961],[0 1 .5],[0.76 1 0.76],100);

colormap(col)
fname= '/users/amaliewillum/dropbox/2016_WILLUM/rapport/images';
filename = 'CV_L_svk_real_2nd_delta';


print( '-dpng' ,fullfile(fname, filename))

% print('-dpng' ,'CV_L_svk_real_2nd_delta')
% print('-dpng' ,'CV_L90_svk_real_2nd_delta')
hold off



%% L x h: CV_L_lin_rand_1st
clear all 
close all

addpath('../fem');

load('CV_L_svk_rand_1st.mat')

% Get the field names of the structure.
fields = fieldnames(CV_L_svk_rand, '-full');
% Find out how many there are - for our loop.
numberOfFields = length(fields);
for f = 1 : numberOfFields
    thisField = fields{f};
    commandLine = sprintf('%s = CV_L_svk_rand.%s', thisField, thisField);
    eval(commandLine);
end
% Release temporary variables.
clear('f', 'thisField', 'numberOfFields');
clear('fields', 'commandLine');

delta_thing = [];
h_thing = [];

L = 1:9;

val = 1e-3;
for i = 1:length(L)-1
   tmp = abs(delta(:,i) - val);
   [c idx] = min(tmp);
   h_thing = [h_thing; h(idx)];
   delta_thing = [delta_thing; delta(idx,i)];
end
%% Plot CV_L_svk_rand_1st_delta.png
close all

figure('units','normalized','position',[.22,.3,.42,.5])
p1 = surf(L',h,delta,log(delta));
hold on 
% p2 = plot3([L(1) L(end)],[h_thing(1) h_thing(end)],[delta_thing(1) delta_thing(end)],'k-','LineWidth',5);
p2 = plot3(L(1:end-1)',h_thing,delta_thing,'k-','LineWidth',3);
p3 = plot3(L(1),h_thing(1),delta_thing(1),'b.','MarkerSize',40);
p4 = plot3(L(end-1),h_thing(end),delta_thing(end),'r.','MarkerSize',40);

l = legend([p3,p4],['( h , $L$ ) = ( ' num2str(h_thing(1),'%.2e') ' , ' num2str(L(1),'%.1e') ' )'   ],...
    ['( h , $L$ ) = ( ' num2str(h_thing(end),'%.2e') ' , ' num2str(L(end-1),'%.1e') ' )' ],...
    'Location','NorthEast');
set(l,'Interpreter','latex','FontSize',12)

% set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

% title('Absolute Error','interpreter','latex','fontsize',18)
xlabel('Mesh Size $L$','interpreter','latex','fontsize',15)
ylabel('Step Size $h$','interpreter','latex','fontsize',15)
zlabel('Error $\delta$','interpreter','latex','fontsize',15)

col = makeColorMap([ .07840 .6 .5961],[0 1 .5],[0.76 1 0.76],100);

colormap(col)

alpha 0.8

box on
axis([1,8,h(end),h(1),1e-40,1e+10])

shading interp 
lighting gouraud
camlight('headlight','infinite')
material dull

view(-50,30)
% view(-90,0)

col = makeColorMap([ .07840 .6 .5961],[0 1 .5],[0.76 1 0.76],100);

colormap(col)
fname= '/users/amaliewillum/dropbox/2016_WILLUM/rapport/images';
filename = 'CV_L_svk_rand_1st_delta';


print( '-dpng' ,fullfile(fname, filename))

% print('-dpng' ,'CV_L_svk_rand_1st_delta')
% print('-dpng' ,'CV_L90_svk_rand_1st_delta')
hold off



%% L x h: CV_L_svk_e_1st
clear all 
close all

addpath('../fem');

load('CV_L_svk_e_1st.mat')

% Get the field names of the structure.
fields = fieldnames(CV_L_svk_e, '-full');
% Find out how many there are - for our loop.
numberOfFields = length(fields);
for f = 1 : numberOfFields
    thisField = fields{f};
    commandLine = sprintf('%s = CV_L_svk_e.%s', thisField, thisField);
    eval(commandLine);
end
% Release temporary variables.
clear('f', 'thisField', 'numberOfFields');
clear('fields', 'commandLine');

delta_thing = [];
h_thing = [];

L = 1:9;

val = 1e-3;
for i = 1:length(L)-1
   tmp = abs(delta(:,i) - val);
   [c idx] = min(tmp);
   h_thing = [h_thing; h(idx)];
   delta_thing = [delta_thing; delta(idx,i)];
end
%% Plot CV_L_svk_e_1st_delta.png
close all

figure('units','normalized','position',[.22,.3,.42,.5])
p1 = surf(L',h,delta,log(delta));
hold on 
% p2 = plot3([L(1) L(end)],[h_thing(1) h_thing(end)],[delta_thing(1) delta_thing(end)],'k-','LineWidth',5);
p2 = plot3(L(1:end-1)',h_thing,delta_thing,'k-','LineWidth',3);
p3 = plot3(L(1),h_thing(1),delta_thing(1),'b.','MarkerSize',40);
p4 = plot3(L(end-1),h_thing(end),delta_thing(end),'r.','MarkerSize',40);

l = legend([p3,p4],['( h , $L$ ) = ( ' num2str(h_thing(1),'%.2e') ' , ' num2str(L(1),'%.1e') ' )'   ],...
    ['( h , $L$ ) = ( ' num2str(h_thing(end),'%.2e') ' , ' num2str(L(end-1),'%.1e') ' )' ],...
    'Location','NorthEast');
set(l,'Interpreter','latex','FontSize',12)

% set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

% title('Absolute Error','interpreter','latex','fontsize',18)
xlabel('Mesh Size $L$','interpreter','latex','fontsize',15)
ylabel('Step Size $h$','interpreter','latex','fontsize',15)
zlabel('Error $\delta$','interpreter','latex','fontsize',15)

col = makeColorMap([ .07840 .6 .5961],[0 1 .5],[0.76 1 0.76],100);

colormap(col)

alpha 0.8

box on
axis([1,8,h(end),h(1),1e-40,1e+10])

shading interp 
lighting gouraud
camlight('headlight','infinite')
material dull

view(-50,30)
% view(-90,0)

col = makeColorMap([ .07840 .6 .5961],[0 1 .5],[0.76 1 0.76],100);

colormap(col)
fname= '/users/amaliewillum/dropbox/2016_WILLUM/rapport/images';
filename = 'CV_L_svk_e_1st_delta';


print( '-dpng' ,fullfile(fname, filename))

% print('-dpng' ,'CV_L_svk_e_1st_delta')
% print('-dpng' ,'CV_L90_svk_e_1st_delta')
hold off



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Sammenligning af modellerne etc.
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
close all
clear all

addpath('../fem');

load('L_delta.mat')

L = 1:9;

nr = 10;

h_nr = h(nr);
% nu_nr = nu(nr);

delta1 = lin_e; 
delta2 = lin_r;
delta3 = lin_real_1st;
delta4 = lin_real_2nd;
delta5 = svk_e; 
delta6 = svk_r;
delta7 = svk_real_1st;
delta8 = svk_real_2nd;

ColorSet = fliplr(hsv(8));
close all
figure('units','normalized','position',[.26,.3,.46,.5])
p3 = plot(L,delta1(nr,:),'color',ColorSet(8,:),'linewidth',4);
hold on
p2 = plot(L,delta2(nr,:),'color',ColorSet(7,:),'linewidth',2);
p1 = plot(L,delta3(nr,:),'color',ColorSet(6,:),'linewidth',2);
p4 = plot(L,delta4(nr,:),'color',ColorSet(5,:),'linewidth',2);
p5 = plot(L,delta5(nr,:),'color',ColorSet(4,:),'linewidth',2);
p6 = plot(L,delta6(nr,:),'color',ColorSet(3,:),'linewidth',2);
p7 = plot(L,delta7(nr,:),'color',ColorSet(2,:),'linewidth',2);
p8 = plot(L,delta8(nr,:),'color',ColorSet(1,:),'linewidth',2);


l = legend([p3,p2,p1,p4,p5,p6,p7,p8],'lin$_e$ 1st','lin$_r$ 1st','lin$_{irl}$ 1st','lin$_{irl}$ 2nd','svk$_e$ 1st','svk$_r$ 1st','svk$_{irl}$ 1st','svk$_{irl}$ 2nd','location','southeastoutside');
set(l,'fontsize',14,'interpreter','latex');

title('2D view of Crossvalidation for $L$ and $h$ with $h = 10^{-4}$','interpreter','latex','fontsize',18)

ylabel('Error $\delta$','interpreter','latex','fontsize',20)
xlabel('Mesh Size $L$','interpreter','latex','fontsize',15)

set(gca, 'YScale', 'log')

plot(linspace(1,5.3,10),ones(1,10)*1e-3,'k-')
plot(ones(1,10)*5.3,logspace(-4.2,-3,10),'k-')



axis tight

h11 = axes('Position', [.74 .24 .23 .23], 'Layer','top');
hold on
set(gca, 'YScale', 'log')
box on
annotation('line',[.455 .78], [.1 + .076 .45]);
axis(h11,'off',[1,5,10^(-4.5),10^(-3)])
title({'Feasible Area for', 'tol = $10^{-3}$ and','$h = 10^{-4}$'},'interpreter','latex')

h1 = axes('Position', [.76 .45 .23 .23], 'Layer','top');
hold on
set(gca, 'YScale', 'log')
box on
annotation('textbox',[.78 .45 .16 .13]);
axis(h1,'off',[1,5,10^(-4.5),10^(-3)])

fname= '/users/amaliewillum/dropbox/2016_WILLUM/rapport/images';
filename = 'delta_L_contrast';


print( '-depsc' ,fullfile(fname, filename))
hold off


ColorSet = fliplr(hsv(8));

nr = 3;
close all

figure('units','normalized','position',[.26,.3,.46,.5])
p3 = plot(h,delta1(:,nr),'color',ColorSet(8,:),'linewidth',4);
hold on
% p2 = plot(h,delta2(:,nr),'color',ColorSet(7,:),'linewidth',2);
p1 = plot(h,delta3(:,nr),'color',ColorSet(6,:),'linewidth',2);
p4 = plot(h,delta4(:,nr),'color',ColorSet(5,:),'linewidth',2);
p5 = plot(h,delta5(:,nr),'color',ColorSet(4,:),'linewidth',2);
% p6 = plot(h,delta6(:,nr),'color',ColorSet(3,:),'linewidth',2);
p7 = plot(h,delta7(:,nr),'color',ColorSet(2,:),'linewidth',2);
p8 = plot(h,delta8(:,nr),'color',ColorSet(1,:),'linewidth',2);


plot(logspace(-2,-2.9,10),ones(1,10)*1e-3,'k--')
plot(ones(1,10)*(1.79)*1e-3,logspace(-6,-3,10),'k--')
plot(ones(1,10)*(1.22)*1e-3,logspace(-6,-3,10),'k--')

point = plot([1.79*1e-3 1.22*1e-3],[1e-6,1e-6],'k.','markersize',20)

l = legend([point,p3,p1,p4,p5,p7,p8],'$h \approx 10^{-3}$','lin$_e$ 1st','lin$_{irl}$ 1st','lin$_{irl}$ 2nd','svk$_e$ 1st','svk$_{irl}$ 1st','svk$_{irl}$ 2nd','location','southeastoutside');
set(l,'fontsize',14,'interpreter','latex');

title('2D view of feasible area from Crossvalidation for $L_3$ and $h$','interpreter','latex','fontsize',18)

ylabel('$\delta$','interpreter','latex','fontsize',20)
xlabel('$-h$','interpreter','latex','fontsize',15)



axis([1e-5,10^(-2),1e-6,1e-2])
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
set(gca, 'XDir', 'reverse')

fname= '/users/amaliewillum/dropbox/2016_WILLUM/rapport/images';
filename = 'delta_Lh_contrast';


print( '-depsc' ,fullfile(fname, filename))
hold off
