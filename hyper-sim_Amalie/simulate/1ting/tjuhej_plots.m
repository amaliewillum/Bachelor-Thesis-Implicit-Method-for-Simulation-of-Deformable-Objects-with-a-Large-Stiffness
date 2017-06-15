clear all
close all
% clc


dut = load('CV_E_lin1_e.mat');

CV = dut.CV_E_lin1_e;
E = CV.E;
E_best = CV.E_best;
h = CV.h;
h_best = CV.h_best;
delta = CV.delta;
delta_best = CV.delta_best;

%%
% close all 

figure(100)

surf(E,h,delta,log(delta));
hold on 
p2 = plot3(E_best,h_best,delta_best,'r.', 'MarkerSize',40);
% text(E_best,h_best,delta_best,['(',num2str(E_best),',',num2str(h_best),',',num2str(delta_best),')'],'Color','red','FontSize',15)

title({'Crossvalidation of parameters  on 1st order FDM: E, h. '; '$|\delta| = K_{approx} - K_{exact}$'},'Interpreter','LaTex','FontSize',15)

% title('Crossvalidation of parameters  on 1st order FDM: E, h. Point (E,h,$|\delta|$)','Interpreter','LaTex','FontSize',15)
legend( [p2],  ['Point (E,h,|\delta|) : ' '(',num2str(round(E_best,3,'significant')),' , ',num2str(round(h_best,3,'significant')),' , ', num2str(round(delta_best,3,'significant')),')'],'Location','NorthEast')

view(-90,0)
xlabel('E Youngs Modulus','Interpreter','LaTex','FontSize',15)
ylabel('h Step size','Interpreter','LaTex','FontSize',15)
zlabel('$|\delta|$ Error','Interpreter','LaTex','FontSize',15)


% view(-48,20)
colormap spring 
shading interp
camlight headlight
% camlight

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')


%%
close all

[row, col] = find(delta < 1e-1);

x = unique(col);
y = unique(row);
z = delta(min(y):max(y),min(x):max(x));

max(h(y))
max(E(x))


figure(56789)

surf(E(x),h(y),z,log(z))
colormap spring 
% shading interp
% camlight left
% view(-90,0)
% view(-48,20)

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')
xlabel('E Youngs Modulus','Interpreter','LaTex','FontSize',15)
ylabel('h Step size','Interpreter','LaTex','FontSize',15)
zlabel('$|\delta|$ Error','Interpreter','LaTex','FontSize',15)
