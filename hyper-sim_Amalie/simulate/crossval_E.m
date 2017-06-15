%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   CROSSVALIDATION OF E AND H - with h = logspace(-1,-30,100) and p_norm
%   and delta t
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% E x h: CV_E_lin_real_squeeze_1st_h301
clear all 
close all

addpath('../fem');

load('CV_E_lin_real_squeeze_1st_h301.mat')

% Get the field names of the structure.
fields = fieldnames(CV_E_lin_real, '-full');
% Find out how many there are - for our loop.
numberOfFields = length(fields);
for f = 1 : numberOfFields
    thisField = fields{f};
    commandLine = sprintf('%s = CV_E_lin_real.%s', thisField, thisField);
    eval(commandLine);
end
% Release temporary variables.
clear('f', 'thisField', 'numberOfFields');
clear('fields', 'commandLine');

delta_real_thing = [];
h_real_thing = [];
delta_real = delta;

val = 1e-3;
for i = 1:length(E)
   tmp = abs(delta(:,i) - val);
   [c idx] = min(tmp);
   h_real_thing = [h_real_thing; h(idx)];
   delta_real_thing = [delta_real_thing; delta(idx,i)];
end
%% Plot CV_E_lin_real_squeeze_1st_h301_delta.png
close all

E_tmp = E(1:1:end);
h_tmp = h(1:1:end);
delta_tmp = delta(1:1:end,1:1:end);

figure('units','normalized','position',[.22,.3,.42,.5])
p1 = surf(E_tmp,h_tmp,delta_tmp,log(delta_tmp));
hold on 
p2 = plot3([E(1) E(end)],[h_thing(1) h_thing(end)],[delta_thing(1) delta_thing(end)],'k-','LineWidth',3);
p3 = plot3(E(1),h_thing(1),delta_thing(1),'b.','MarkerSize',40);
p4 = plot3(E(end),h_thing(end),delta_thing(end),'r.','MarkerSize',40);

L = legend([p3,p4],['( h , E ) = ( ' num2str(h_thing(1),'%.2e') ' , ' num2str(E(1),'%.1e') ' )'   ],...
    ['( h , E ) = ( ' num2str(h_thing(end),'%.2e') ' , ' num2str(E(end),'%.1e') ' )' ],...
    'Location','NorthEast');
set(L,'Interpreter','latex','FontSize',12)

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

% title('Absolute Error','interpreter','latex','fontsize',18)
xlabel('Youngs Modulus $E$','interpreter','latex','fontsize',15)
ylabel('Step Size $h$','interpreter','latex','fontsize',15)
zlabel('Error $\delta$','interpreter','latex','fontsize',15)

% col = makeColorMap([0 .2 .4], [0 .5 1], [.6 .8 1],100); % blå

col = makeColorMap([.4 0 .2],[1 .4 1 ],[1 .8 1],100); % lilla

% col = makeColorMap([0 .6 .298],[0 1 .5],[.6 1 .8],100);

% col = makeColorMap([ .07840 .6 .5961],[0 1 .5],[0.76 1 0.76],100);

colormap(col)
alpha .8
box on
shading interp 
lighting gouraud
camlight('headlight','infinite')
material dull

view(-50,30)

fname= '/users/amaliewillum/dropbox/2016_WILLUM/rapport/images';
filename = 'CV_E_lin_real_1st_delta';


print( '-dpng' ,fullfile(fname, filename))

hold off

%% E x h: CV_E_lin_e_squeeze_1st_h301
clear all 
close all

addpath('../fem');

load('CV_E_lin_e_squeeze_1st_h301.mat')

% Get the field names of the structure.
fields = fieldnames(CV_E_lin_real, '-full');
% Find out how many there are - for our loop.
numberOfFields = length(fields);
for f = 1 : numberOfFields
    thisField = fields{f};
    commandLine = sprintf('%s = CV_E_lin_real.%s', thisField, thisField);
    eval(commandLine);
end
% Release temporary variables.
clear('f', 'thisField', 'numberOfFields');
clear('fields', 'commandLine');

delta_thing = [];
h_thing = [];

val = 1e-3;
for i = 1:length(E)
   tmp = abs(delta(:,i) - val);
   [c idx] = min(tmp);
   h_thing = [h_thing; h(idx)];
   delta_thing = [delta_thing; delta(idx,i)];
end
%% Plot CV_E_lin_e_squeeze_1st_h301_delta.png
close all

figure('units','normalized','position',[.22,.3,.42,.5])
p1 = surf(E,h,delta,log(delta));
hold on 
p2 = plot3([E(1) E(end)],[h_thing(1) h_thing(end)],[delta_thing(1) delta_thing(end)],'k-','LineWidth',3);
p3 = plot3(E(1),h_thing(1),delta_thing(1),'b.','MarkerSize',40);
p4 = plot3(E(end),h_thing(end),delta_thing(end),'r.','MarkerSize',40);

L = legend([p3,p4],['( h , E ) = ( ' num2str(h_thing(1),'%.2e') ' , ' num2str(E(1),'%.1e') ' )'   ],...
    ['( h , E ) = ( ' num2str(h_thing(end),'%.2e') ' , ' num2str(E(end),'%.1e') ' )' ],...
    'Location','NorthEast');
set(L,'Interpreter','latex','FontSize',12)

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

% title('Absolute Error','interpreter','latex','fontsize',18)
xlabel('Youngs Modulus $E$','interpreter','latex','fontsize',15)
ylabel('Step Size $h$','interpreter','latex','fontsize',15)
zlabel('Error $\delta$','interpreter','latex','fontsize',15)

col = makeColorMap([.4 0 .2],[1 .4 1 ],[1 .8 1],100);

colormap(col)

alpha .8
box on

shading interp 
lighting gouraud
camlight headlight infinite
material dull

view(-50,30)

fname= '/users/amaliewillum/dropbox/2016_WILLUM/rapport/images';
filename = 'CV_E_lin_e_1st_delta';


print( '-dpng' ,fullfile(fname, filename))

hold off


%% E x h: CV_E_lin_rand_squeeze_1st_h301
% clear all 
close all

addpath('../fem');

load('CV_E_lin_rand_squeeze_1st_h301.mat')

% Get the field names of the structure.
fields = fieldnames(CV_E_lin_real, '-full');
% Find out how many there are - for our loop.
numberOfFields = length(fields);
for f = 1 : numberOfFields
    thisField = fields{f};
    commandLine = sprintf('%s = CV_E_lin_real.%s', thisField, thisField);
    eval(commandLine);
end
% Release temporary variables.
clear('f', 'thisField', 'numberOfFields');
clear('fields', 'commandLine');

delta_rand_thing = [];
h_rand_thing = [];
delta_rand = delta;
val = 1e-3;
for i = 1:length(E)
   tmp = abs(delta(:,i) - val);
   [c idx] = min(tmp);
   h_rand_thing = [h_rand_thing; h(idx)];
   delta_rand_thing = [delta_rand_thing; delta(idx,i)];
end

% clearvars -except delta
% lin_r = delta;
% 
% save('E_delta.mat','-append','lin_r');
%% Plot CV_E_lin_rand_squeeze_1st_h301_delta.png
close all

figure('units','normalized','position',[.22,.3,.42,.5])
p1 = surf(E,h,delta,log(delta));
hold on 
p2 = plot3([E(1) E(end)],[h_thing(1) h_thing(end)],[delta_thing(1) delta_thing(end)],'k-','LineWidth',3);
% p2 = plot3(E,h_thing,delta_thing,'k-','LineWidth',3);
p3 = plot3(E(1),h_thing(1),delta_thing(1),'b.','MarkerSize',40);
p4 = plot3(E(end),h_thing(end),delta_thing(end),'r.','MarkerSize',40);

L = legend([p3,p4],['( h , E ) = ( ' num2str(h_thing(1),'%.2e') ' , ' num2str(E(1),'%.1e') ' )'   ],...
    ['( h , E ) = ( ' num2str(h_thing(end),'%.2e') ' , ' num2str(E(end),'%.1e') ' )' ],...
    'Location','NorthEast');
set(L,'Interpreter','latex','FontSize',12)

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

% title('Absolute Error','interpreter','latex','fontsize',18)
xlabel('Youngs Modulus $E$','interpreter','latex','fontsize',15)
ylabel('Step Size $h$','interpreter','latex','fontsize',15)
zlabel('Error $\delta$','interpreter','latex','fontsize',15)

% axis([E(1),E(end),h(34),h(1),0,10^8])
col = makeColorMap([.4 0 .2],[1 .4 1 ],[1 .8 1],100);
colormap(col)
% colorbar
% colormap spring
% colormap cool 
% colorbar eastoutside
% shading flat
alpha .8

box on

shading interp 
lighting gouraud
camlight headlight infinite
material dull 

view(-50,30)

fname= '/users/amaliewillum/dropbox/2016_WILLUM/rapport/images';
filename = 'CV_E_lin_rand_1st_delta';


print( '-dpng' ,fullfile(fname, filename))
% print( '-dpng' ,'CV_E_lin_rand_squeeze_1st_h301_delta')
hold off

%% E x h: CV_E_lin_real_squeeze_2nd_h301
clear all 
close all

addpath('../fem');

load('CV_E_lin_real_squeeze_2nd_h301.mat')

% Get the field names of the structure.
fields = fieldnames(CV_E_lin_real, '-full');
% Find out how many there are - for our loop.
numberOfFields = length(fields);
for f = 1 : numberOfFields
    thisField = fields{f};
    commandLine = sprintf('%s = CV_E_lin_real.%s', thisField, thisField);
    eval(commandLine);
end
% Release temporary variables.
clear('f', 'thisField', 'numberOfFields');
clear('fields', 'commandLine');

delta_thing = [];
h_thing = [];

val = 1e-3;
for i = 1:length(E)
   tmp = abs(delta(:,i) - val);
   [c idx] = min(tmp);
   h_thing = [h_thing; h(idx)];
   delta_thing = [delta_thing; delta(idx,i)];
end
%% Plot CV_E_lin_real_squeeze_2nd_h301_delta.png
close all

E_tmp = E(1:1:end);
h_tmp = h(1:1:end);
delta_tmp = delta(1:1:end,1:1:end);

figure('units','normalized','position',[.22,.3,.42,.5])
p1 = surf(E_tmp,h_tmp,delta_tmp,log(delta_tmp));
hold on 
p2 = plot3([E(1) E(end)],[h_thing(1) h_thing(end)],[delta_thing(1) delta_thing(end)],'k-','LineWidth',3);
p3 = plot3(E(1),h_thing(1),delta_thing(1),'b.','MarkerSize',40);
p4 = plot3(E(end),h_thing(end),delta_thing(end),'r.','MarkerSize',40);

L = legend([p3,p4],['( h , E ) = ( ' num2str(h_thing(1),'%.2e') ' , ' num2str(E(1),'%.1e') ' )'   ],...
    ['( h , E ) = ( ' num2str(h_thing(end),'%.2e') ' , ' num2str(E(end),'%.1e') ' )' ],...
    'Location','NorthEast');
set(L,'Interpreter','latex','FontSize',12)

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

% title('Absolute Error','interpreter','latex','fontsize',18)
xlabel('Youngs Modulus $E$','interpreter','latex','fontsize',15)
ylabel('Step Size $h$','interpreter','latex','fontsize',15)
zlabel('Error $\delta$','interpreter','latex','fontsize',15)

% axis([E(1),E(end),h(34),h(1),0,10^8])

col = makeColorMap([.4 0 .2],[1 .4 1 ],[1 .8 1],100);
colormap(col)

alpha 0.8

box on

shading interp 
lighting gouraud
camlight('headlight','infinite')
material dull

view(-50,30)

fname= '/users/amaliewillum/dropbox/2016_WILLUM/rapport/images';
filename = 'CV_E_lin_real_2nd_delta';


print( '-dpng' ,fullfile(fname, filename))

hold off

%% E x h: CV_E_svk_real_squeeze_1st_h301
clear all 
close all

addpath('../fem');

load('CV_E_svk_real_squeeze_1st_h301.mat')

% Get the field names of the structure.
fields = fieldnames(CV_E_svk_real, '-full');
% Find out how many there are - for our loop.
numberOfFields = length(fields);
for f = 1 : numberOfFields
    thisField = fields{f};
    commandLine = sprintf('%s = CV_E_svk_real.%s', thisField, thisField);
    eval(commandLine);
end
% Release temporary variables.
clear('f', 'thisField', 'numberOfFields');
clear('fields', 'commandLine');

delta_thing = [];
h_thing = [];

val = 1e-3;
for i = 1:length(E)
   tmp = abs(delta(:,i) - val);
   [c idx] = min(tmp);
   h_thing = [h_thing; h(idx)];
   delta_thing = [delta_thing; delta(idx,i)];
end
%% Plot CV_E_svk_real_squeeze_1st_h301_delta.png
close all

figure('units','normalized','position',[.22,.3,.42,.5])
p1 = surf(E,h,delta,log(delta));
hold on 
p2 = plot3([E(1) E(end)],[h_thing(1) h_thing(end)],[delta_thing(1) delta_thing(end)],'k-','LineWidth',3);
% p2 = plot3(E,h_thing,delta_thing,'k-','LineWidth',3);
p3 = plot3(E(1),h_thing(1),delta_thing(1),'b.','MarkerSize',40);
p4 = plot3(E(end),h_thing(end),delta_thing(end),'r.','MarkerSize',40);

L = legend([p3,p4],['( h , E ) = ( ' num2str(h_thing(1),'%.2e') ' , ' num2str(E(1),'%.1e') ' )'   ],...
    ['( h , E ) = ( ' num2str(h_thing(end),'%.2e') ' , ' num2str(E(end),'%.1e') ' )' ],...
    'Location','NorthEast');
set(L,'Interpreter','latex','FontSize',12)

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

% title('Absolute Error','interpreter','latex','fontsize',18)
xlabel('Youngs Modulus $E$','interpreter','latex','fontsize',15)
ylabel('Step Size $h$','interpreter','latex','fontsize',15)
zlabel('Error $\delta$','interpreter','latex','fontsize',15)


col = makeColorMap([.4 0 .2],[1 .4 1 ],[1 .8 1],100);
colormap(col)

alpha .8

box on

shading interp 
lighting gouraud
camlight infinite
material dull

view(-50,30)

fname= '/users/amaliewillum/dropbox/2016_WILLUM/rapport/images';
filename = 'CV_E_svk_real_1st_delta';


print( '-dpng' ,fullfile(fname, filename))

hold off


%% E x h: CV_E_svk_e_squeeze_1st_h301
clear all 
close all

addpath('../fem');

load('CV_E_svk_e_squeeze_1st_h301.mat')

% Get the field names of the structure.
fields = fieldnames(CV_E_svk_real, '-full');
% Find out how many there are - for our loop.
numberOfFields = length(fields);
for f = 1 : numberOfFields
    thisField = fields{f};
    commandLine = sprintf('%s = CV_E_svk_real.%s', thisField, thisField);
    eval(commandLine);
end
% Release temporary variables.
clear('f', 'thisField', 'numberOfFields');
clear('fields', 'commandLine');

delta_thing = [];
h_thing = [];

val = 1e-3;
for i = 1:length(E)
   tmp = abs(delta(:,i) - val);
   [c idx] = min(tmp);
   h_thing = [h_thing; h(idx)];
   delta_thing = [delta_thing; delta(idx,i)];
end
%% Plot CV_E_svk_e_squeeze_1st_h301_delta.png
close all

figure('units','normalized','position',[.22,.3,.42,.5])
p1 = surf(E,h,delta,log(delta));
hold on 
p2 = plot3([E(1) E(end)],[h_thing(1) h_thing(end)],[delta_thing(1) delta_thing(end)],'k-','LineWidth',3);
p3 = plot3(E(1),h_thing(1),delta_thing(1),'b.','MarkerSize',40);
p4 = plot3(E(end),h_thing(end),delta_thing(end),'r.','MarkerSize',40);

L = legend([p3,p4],['( h , E ) = ( ' num2str(h_thing(1),'%.2e') ' , ' num2str(E(1),'%.1e') ' )'   ],...
    ['( h , E ) = ( ' num2str(h_thing(end),'%.2e') ' , ' num2str(E(end),'%.1e') ' )' ],...
    'Location','NorthEast');
set(L,'Interpreter','latex','FontSize',12)

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

% title('Absolute Error','interpreter','latex','fontsize',18)
xlabel('Youngs Modulus $E$','interpreter','latex','fontsize',15)
ylabel('Step Size $h$','interpreter','latex','fontsize',15)
zlabel('Error $\delta$','interpreter','latex','fontsize',15)


col = makeColorMap([.4 0 .2],[1 .4 1 ],[1 .8 1],100);
colormap(col)

alpha .8

box on

shading interp 
lighting gouraud
camlight headlight infinite
material dull 

view(-50,30)

fname= '/users/amaliewillum/dropbox/2016_WILLUM/rapport/images';
filename = 'CV_E_svk_e_1st_delta';


print( '-dpng' ,fullfile(fname, filename))

hold off

%% E x h: CV_E_svk_rand_squeeze_1st_h301
clear all 
close all

addpath('../fem');

load('CV_E_svk_rand_squeeze_1st_h301.mat')

% Get the field names of the structure.
fields = fieldnames(CV_E_svk_real, '-full');
% Find out how many there are - for our loop.
numberOfFields = length(fields);
for f = 1 : numberOfFields
    thisField = fields{f};
    commandLine = sprintf('%s = CV_E_svk_real.%s', thisField, thisField);
    eval(commandLine);
end
% Release temporary variables.
clear('f', 'thisField', 'numberOfFields');
clear('fields', 'commandLine');

delta_thing = [];
h_thing = [];

val = 1e-3;
for i = 1:length(E)
   tmp = abs(delta(:,i) - val);
   [c idx] = min(tmp);
   h_thing = [h_thing; h(idx)];
   delta_thing = [delta_thing; delta(idx,i)];
end
%% Plot CV_E_svk_rand_squeeze_1st_h301_delta.png
close all

figure('units','normalized','position',[.22,.3,.42,.5])
p1 = surf(E,h,delta,log(delta));
hold on 
p2 = plot3([E(1) E(end)],[h_thing(1) h_thing(end)],[delta_thing(1) delta_thing(end)],'k-','LineWidth',3);
p3 = plot3(E(1),h_thing(1),delta_thing(1),'b.','MarkerSize',40);
p4 = plot3(E(end),h_thing(end),delta_thing(end),'r.','MarkerSize',40);

L = legend([p3,p4],['( h , E ) = ( ' num2str(h_thing(1),'%.2e') ' , ' num2str(E(1),'%.1e') ' )'   ],...
    ['( h , E ) = ( ' num2str(h_thing(end),'%.2e') ' , ' num2str(E(end),'%.1e') ' )' ],...
    'Location','NorthEast');
set(L,'Interpreter','latex','FontSize',12)

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

% title('Absolute Error','interpreter','latex','fontsize',18)
xlabel('Youngs Modulus $E$','interpreter','latex','fontsize',15)
ylabel('Step Size $h$','interpreter','latex','fontsize',15)
zlabel('Error $\delta$','interpreter','latex','fontsize',15)

col = makeColorMap([.4 0 .2],[1 .4 1 ],[1 .8 1],100);
colormap(col)

alpha .8

box on

shading interp 
lighting gouraud
camlight headlight infinite
material dull

view(-50,30)

fname= '/users/amaliewillum/dropbox/2016_WILLUM/rapport/images';
filename = 'CV_E_svk_rand_1st_delta';


print( '-dpng' ,fullfile(fname, filename))

hold off

%% E x h: CV_E_svk_real_squeeze_2nd_h301
clear all 
close all

addpath('../fem');

load('CV_E_svk_real_squeeze_2nd_h301.mat')

% Get the field names of the structure.
fields = fieldnames(CV_E_svk_real, '-full');
% Find out how many there are - for our loop.
numberOfFields = length(fields);
for f = 1 : numberOfFields
    thisField = fields{f};
    commandLine = sprintf('%s = CV_E_svk_real.%s', thisField, thisField);
    eval(commandLine);
end
% Release temporary variables.
clear('f', 'thisField', 'numberOfFields');
clear('fields', 'commandLine');

delta_thing = [];
h_thing = [];

val = 1e-3;
for i = 1:length(E)
   tmp = abs(delta(:,i) - val);
   [c idx] = min(tmp);
   h_thing = [h_thing; h(idx)];
   delta_thing = [delta_thing; delta(idx,i)];
end
%% Plot CV_E_svk_real_squeeze_2nd_h301_delta.png
close all

figure('units','normalized','position',[.22,.3,.42,.5]) % uden color bar [.22,.3,.42,.5] 
p1 = surf(E,h,delta,log(delta));
hold on 
p2 = plot3([E(1) E(end)],[h_thing(1) h_thing(end)],[delta_thing(1) delta_thing(end)],'k-','LineWidth',3);
% p2 = plot3(E,h_thing,delta_thing,'k-','LineWidth',3);
p3 = plot3(E(1),h_thing(1),delta_thing(1),'b.','MarkerSize',40);
p4 = plot3(E(end),h_thing(end),delta_thing(end),'r.','MarkerSize',40);

L = legend([p3,p4],['( h , E ) = ( ' num2str(h_thing(1),'%.2e') ' , ' num2str(E(1),'%.1e') ' )'   ],...
    ['( h , E ) = ( ' num2str(h_thing(end),'%.2e') ' , ' num2str(E(end),'%.1e') ' )' ],...
    'Location','NorthEast');
set(L,'Interpreter','latex','FontSize',12)

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

% title('Absolute Error','interpreter','latex','fontsize',18)
xlabel('Youngs Modulus $E$','interpreter','latex','fontsize',15)
ylabel('Step Size $h$','interpreter','latex','fontsize',15)
zlabel('Error $\delta$','interpreter','latex','fontsize',15)

col = makeColorMap([.4 0 .2],[1 .4 1 ],[1 .8 1],100); % lilla

colormap(col)
% colorbar

alpha .8
box on

shading interp 
lighting gouraud
camlight infinite
material dull

view(-50,30)

fname= '/users/amaliewillum/dropbox/2016_WILLUM/rapport/images';
filename = 'CV_E_svk_real_2nd_delta';


print( '-dpng' ,fullfile(fname, filename))

hold off


%%
close all

figure('units','normalized','position',[.22,.3,.42,.5]) % uden color bar [.22,.3,.42,.5] 

hold on
p1 = surf(E,h,delta_real,log(delta_real));
% colormap winter
hold on 
% set(p1,'colormap','winter','alpha',.5)
dut = surf(E,h,delta_rand,log(delta_rand),'summer');
% col = makeColorMap([.4 0 .2],[1 .4 1 ],[1 .8 1],100); % lilla

% colormap(col)
% p2 = plot3([E(1) E(end)],[h_thing(1) h_thing(end)],[delta_thing(1) delta_thing(end)],'k-','LineWidth',3);
% p2 = plot3(E,h_thing,delta_thing,'k-','LineWidth',3);
% p3 = plot3(E(1),h_thing(1),delta_thing(1),'b.','MarkerSize',40);
% p4 = plot3(E(end),h_thing(end),delta_thing(end),'r.','MarkerSize',40);

% L = legend([p3,p4],['( h , E ) = ( ' num2str(h_thing(1),'%.2e') ' , ' num2str(E(1),'%.1e') ' )'   ],...
%     ['( h , E ) = ( ' num2str(h_thing(end),'%.2e') ' , ' num2str(E(end),'%.1e') ' )' ],...
%     'Location','NorthEast');
% set(L,'Interpreter','latex','FontSize',12)

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

% title('Absolute Error','interpreter','latex','fontsize',18)
xlabel('Youngs Modulus $E$','interpreter','latex','fontsize',15)
ylabel('Step Size $h$','interpreter','latex','fontsize',15)
zlabel('Error $\delta$','interpreter','latex','fontsize',15)

grid minor
% colorbar

alpha .8
box on

shading interp 
lighting gouraud
camlight infinite
material dull

view(-90,0)

fname= '/users/amaliewillum/dropbox/2016_WILLUM/rapport/images';
filename = 'CV_E_svk_real_2nd_delta';


print( '-dpng' ,fullfile(fname, filename))

hold off








%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Sammenligning af modellerne etc.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
close all
clear all

addpath('../fem');

load('E_delta.mat')

h_nr = 10;

delta1 = lin_e(h_nr,:); 
delta2 = lin_r(h_nr,:);
delta3 = lin_real_1st(h_nr,:);
delta4 = lin_real_2nd(h_nr,:);
delta5 = svk_e(h_nr,:); 
delta6 = svk_r(h_nr,:);
delta7 = svk_real_1st(h_nr,:);
delta8 = svk_real_2nd(h_nr,:);

ColorSet = fliplr(hsv(8));

figure('units','normalized','position',[.26,.3,.46,.5])
p3 = plot(E,delta1,'color',ColorSet(8,:),'linewidth',2);
hold on
p2 = plot(E,delta2,'color',ColorSet(7,:),'linewidth',2);
p1 = plot(E,delta3,'color',ColorSet(6,:),'linewidth',2);
p4 = plot(E,delta4,'color',ColorSet(5,:),'linewidth',2);
p5 = plot(E,delta5,'color',ColorSet(4,:),'linewidth',2);
p6 = plot(E,delta6,'color',ColorSet(3,:),'linewidth',2);
p7 = plot(E,delta7,'color',ColorSet(2,:),'linewidth',2);
p8 = plot(E,delta8,'color',ColorSet(1,:),'linewidth',2);

plot(logspace(2,4.7,20),ones(1,20)*1e-3,'k-')
plot(10^(4.7)*ones(1,20), logspace(-6,-3,20),'k-')
% 
l = legend([p3,p2,p1,p4,p5,p6,p7,p8],'lin$_e$ 1st','lin$_r$ 1st','lin$_{irl}$ 1st','lin$_{irl}$ 2nd','svk$_e$ 1st','svk$_r$ 1st','svk$_{irl}$ 1st','svk$_{irl}$ 2nd','location','southeastoutside');
set(l,'fontsize',14,'interpreter','latex');
% 
title('2D view of Crossvalidation for $E$ and $h$ with $h = 10^{-4}$','interpreter','latex','fontsize',18)
% 
ylabel('Error $\delta$','interpreter','latex','fontsize',20)
xlabel('Youngs Modulus $E$','interpreter','latex','fontsize',15)

axis tight
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
% 
h = axes('Position', [.74 .24 .23 .23], 'Layer','top');
hold on
set(gca, 'YScale', 'log')
box on
annotation('line',[.243 .78], [.1 + .076 .45]);
axis(h,'off',[1,5,10^(-4.5),10^(-3)])
title({'Feasible Area for', 'tol = $10^{-3}$ and','$h = 10^{-4}$'},'interpreter','latex')

h1 = axes('Position', [.76 .45 .23 .23], 'Layer','top');
hold on
set(gca, 'YScale', 'log')
box on
annotation('textbox',[.78 .45 .16 .13]);
axis(h1,'off',[1,5,10^(-4.5),10^(-3)])


h1 = axes('Position', [.77 .65 .23 .23], 'Layer','top');
p3 = plot(E,delta1,'color',ColorSet(8,:),'linewidth',2);
hold on
p2 = plot(E,delta2,'color',ColorSet(7,:),'linewidth',2);
p1 = plot(E,delta3,'color',ColorSet(6,:),'linewidth',2);
p4 = plot(E,delta4,'color',ColorSet(5,:),'linewidth',2);
p5 = plot(E,delta5,'color',ColorSet(4,:),'linewidth',2);
p6 = plot(E,delta6,'color',ColorSet(3,:),'linewidth',2);
p7 = plot(E,delta7,'color',ColorSet(2,:),'linewidth',2);
p8 = plot(E,delta8,'color',ColorSet(1,:),'linewidth',2);
axis tight
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
annotation('textbox',[.77 .65 .23 .23]);
annotation('line',[.243 .77], [.1 + .076 .65]);
axis(h1,'off',[100,320,1e-6,1e-3])
title('Zoom','interpreter','latex')

fname= '/users/amaliewillum/dropbox/2016_WILLUM/rapport\images';
filename = 'delta_E_contrast';


print( '-depsc' ,fullfile(fname, filename))
hold off
