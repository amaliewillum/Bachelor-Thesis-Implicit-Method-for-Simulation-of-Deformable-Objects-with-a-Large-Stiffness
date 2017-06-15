close all
clear all


dut = load('CV_E_1storder_random');

first_order_CV_E1_random       = dut.CV_E1_rand;
first_order_E               = first_order_CV_E1_random.E;
first_order_E_best          = first_order_CV_E1_random.E_best;
first_order_E_best1         = first_order_CV_E1_random.E_best1;
first_order_h               = first_order_CV_E1_random.h;
first_order_h_best          = first_order_CV_E1_random.h_best;
first_order_h_best1         = first_order_CV_E1_random.h_best1;
first_order_delta           = first_order_CV_E1_random.delta;
first_order_delta_best      = first_order_CV_E1_random.delta_best;
first_order_delta1          = first_order_CV_E1_random.delta1;
first_order_delta_best1     = first_order_CV_E1_random.delta_best1;

%%
close all

figure(100)

surf(first_order_E,first_order_h,first_order_delta,log(first_order_delta));
hold on 
p2 = plot3(first_order_E_best,first_order_h_best,first_order_delta_best,'r.', 'MarkerSize',40);
% text(E_best,h_best,delta_best,['(',num2str(E_best),',',num2str(h_best),',',num2str(delta_best),')'],'Color','red','FontSize',15)

title({'Crossvalidation of parameters  on 1st order FDM: E, h. '; '$|\delta| = K_{approx} - K_{exact}$'},'Interpreter','LaTex','FontSize',15)

% title('Crossvalidation of parameters  on 1st order FDM: E, h. Point (E,h,$|\delta|$)','Interpreter','LaTex','FontSize',15)
legend( [p2],  ['Point (E,h,|\delta|) : ' '(',num2str(round(first_order_E_best,3,'significant')),' , ',num2str(round(first_order_h_best,3,'significant')),' , ', num2str(round(first_order_delta_best,3,'significant')),')'],'Location','NorthEast')

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

print -depsc 1st_random_E
savefig('1st_random_E.fig')



%%
close all 
clear all

dut = load('CV_nu_1storder_random');

first_order_CV_nu1_rand     = dut.CV_nu1_rand;
first_order_nu              = first_order_CV_nu1_rand.nu;
first_order_nu_best         = first_order_CV_nu1_rand.nu_best;
first_order_nu_best1        = first_order_CV_nu1_rand.nu_best1;
first_order_h               = first_order_CV_nu1_rand.h;
first_order_h_best          = first_order_CV_nu1_rand.h_best;
first_order_h_best1         = first_order_CV_nu1_rand.h_best1;
first_order_delta           = first_order_CV_nu1_rand.delta;
first_order_delta_best      = first_order_CV_nu1_rand.delta_best;
first_order_delta1          = first_order_CV_nu1_rand.delta1;
first_order_delta_best1     = first_order_CV_nu1_rand.delta_best1;

%%
close all 


figure(200)
surf(first_order_nu,first_order_h,first_order_delta,log(first_order_delta));
hold on 
p2 = plot3(first_order_nu_best,first_order_h_best,first_order_delta_best,'r.', 'MarkerSize',40);

title({'Crossvalidation of parameters  on 1st order FDM: $\nu$, h. '; '$|\delta| = K_{approx} - K_{exact}$'},'Interpreter','LaTex','FontSize',15)

legend( [p2],  ['Point (\nu,h,|\delta|) : ' '(',num2str(round(first_order_nu_best,3,'significant')),' , ',num2str(round(first_order_h_best,3,'significant')),' , ', num2str(round(first_order_delta_best,3,'significant')),')'],'Location','NorthEast')


xlabel('$\nu$ Poissons Ratio','Interpreter','LaTex','FontSize',15)
ylabel('h Step size','Interpreter','LaTex','FontSize',15)
zlabel('$|\delta|$ Error','Interpreter','LaTex','FontSize',15)

view(-90,0)
% view(-58,30)
shading interp
camlight right

set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

print -depsc 1st_random_nu
savefig('1st_random_nu.fig')

