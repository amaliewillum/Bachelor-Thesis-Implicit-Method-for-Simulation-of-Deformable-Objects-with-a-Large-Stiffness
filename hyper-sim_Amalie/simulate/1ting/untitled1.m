close all

dut = load('CV_nu_1storder');

first_order_CV_nu1       = dut.CV_nu1;
first_order_nu           = first_order_CV_nu1.nu;
first_order_nu_best      = first_order_CV_nu1.nu_best;
first_order_nu_best1     = first_order_CV_nu1.nu_best1;
first_order_h            = first_order_CV_nu1.h;
first_order_h_best       = first_order_CV_nu1.h_best;
first_order_h_best1      = first_order_CV_nu1.h_best1;
first_order_delta        = first_order_CV_nu1.delta;
first_order_delta_best   = first_order_CV_nu1.delta_best;
first_order_delta1       = first_order_CV_nu1.delta1;
first_order_delta_best1  = first_order_CV_nu1.delta_best1;



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


%%

dut = load('CV_nu_2ndorder');

second_order_CV_nu2       = dut.CV_nu2;
second_order_nu           = second_order_CV_nu2.nu;
second_order_nu_best      = second_order_CV_nu2.nu_best;
second_order_nu_best1     = second_order_CV_nu2.nu_best1;
second_order_h            = second_order_CV_nu2.h;
second_order_h_best       = second_order_CV_nu2.h_best;
second_order_h_best1      = second_order_CV_nu2.h_best1;
second_order_delta        = second_order_CV_nu2.delta;
second_order_delta_best   = second_order_CV_nu2.delta_best;
second_order_delta1       = second_order_CV_nu2.delta1;
second_order_delta_best1  = second_order_CV_nu2.delta_best1;

dut = load('CV_nu_3rdorder');

third_order_CV_nu3       = dut.CV_nu3;
third_order_nu           = third_order_CV_nu3.nu;
third_order_nu_best      = third_order_CV_nu3.nu_best;
third_order_nu_best1     = third_order_CV_nu3.nu_best1;
third_order_h            = third_order_CV_nu3.h;
third_order_h_best       = third_order_CV_nu3.h_best;
third_order_h_best1      = third_order_CV_nu3.h_best1;
third_order_delta        = third_order_CV_nu3.delta;
third_order_delta_best   = third_order_CV_nu3.delta_best;
third_order_delta1       = third_order_CV_nu3.delta1;
third_order_delta_best1  = third_order_CV_nu3.delta_best1;

dut = load('CV_nu_4thorder');

fourth_order_CV_nu4       = dut.CV_nu4;
fourth_order_nu           = fourth_order_CV_nu4.nu;
fourth_order_nu_best      = fourth_order_CV_nu4.nu_best;
fourth_order_nu_best1     = fourth_order_CV_nu4.nu_best1;
fourth_order_h           = fourth_order_CV_nu4.h;
fourth_order_h_best      = fourth_order_CV_nu4.h_best;
fourth_order_h_best1     = fourth_order_CV_nu4.h_best1;
fourth_order_delta       = fourth_order_CV_nu4.delta;
fourth_order_delta_best  = fourth_order_CV_nu4.delta_best;
fourth_order_delta1      = fourth_order_CV_nu4.delta1;
fourth_order_delta_best1 = fourth_order_CV_nu4.delta_best1;



%%  HISTOGRAM plots for E

close all
clc

figure(1)
one = semilogxhist(first_order_delta1(:),'b');
hold on 
two = semilogxhist(second_order_delta1(:),'m');
three = semilogxhist(third_order_delta1(:),'g');
four = semilogxhist(fourth_order_delta1(:),'y');
legend([one,two,three,four],'1st order','2nd order','3rd order','4th order','Location','NorthWest')
title({'Error of FDMs of different orders from crossvalidation of $h$ and $\nu$'},'Interpreter','LaTex','FontSize',15)
xlabel('$|\delta| = \frac{K_{approx}}{K_{exact}}$','FontSize',20,'Interpreter','LaTex')
ylabel('Number of occurences','Interpreter','LaTex')

print -depsc hist_approximations
savefig('hist_approximations.fig')
hold off

%%  SURFACE plots for E 

close all
clc

figh = figure(100);
pos = get(figh,'position');
set(figh,'position',[pos(1:2)/4 pos(3:4)*2])

subplot(2,2,1)
surf(first_order_nu,first_order_h,first_order_delta,log(first_order_delta));
hold on 
p2 = plot3(first_order_nu_best,first_order_h_best,first_order_delta_best,'r.', 'MarkerSize',40);

title({'Crossvalidation of parameters  on 1st order FDM: $\nu$, h. '; '$|\delta| = K_{approx} - K_{exact}$'},'Interpreter','LaTex','FontSize',15)

legend( [p2],  ['Point (\nu,h,|\delta|) : ' '(',num2str(round(first_order_nu_best,3,'significant')),' , ',num2str(round(first_order_h_best,3,'significant')),' , ', num2str(round(first_order_delta_best,3,'significant')),')'],'Location','NorthEast')


xlabel('$\nu$ Poissons Ratio','Interpreter','LaTex','FontSize',15)
ylabel('h Step size','Interpreter','LaTex','FontSize',15)
zlabel('$|\delta|$ Error','Interpreter','LaTex','FontSize',15)


view(-58,30)
shading interp
camlight right

set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')



subplot(2,2,2)
surf(second_order_nu,second_order_h,second_order_delta,log(second_order_delta));
hold on 
p2 = plot3(second_order_nu_best,second_order_h_best,second_order_delta_best,'r.', 'MarkerSize',40);
title({'Crossvalidation of parameters  on 2nd order FDM: $\nu$, h. '; '$|\delta| = K_{approx} - K_{exact}$'},'Interpreter','LaTex','FontSize',15)
legend( [p2],  ['Point (\nu,h,|\delta|) : ' '(',num2str(round(second_order_nu_best,3,'significant')),' , ',num2str(round(second_order_h_best,3,'significant')),' , ', num2str(round(second_order_delta_best,3,'significant')),')'],'Location','NorthEast')


xlabel('$\nu$ Poissons Ratio','Interpreter','LaTex','FontSize',15)
ylabel('h Step size','Interpreter','LaTex','FontSize',15)
zlabel('$|\delta|$ Error','Interpreter','LaTex','FontSize',15)


view(-58,30)
shading interp
camlight right

set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')


subplot(2,2,3)
surf(third_order_nu,third_order_h,third_order_delta,log(third_order_delta));
hold on 
p2 = plot3(third_order_nu_best,third_order_h_best,third_order_delta_best,'r.', 'MarkerSize',40);
title({'Crossvalidation of parameters  on 3rd order FDM: $\nu$, h. '; '$|\delta| = K_{approx} - K_{exact}$'},'Interpreter','LaTex','FontSize',15)
legend( [p2],  ['Point (\nu,h,|\delta|) : ' '(',num2str(round(third_order_nu_best,3,'significant')),' , ',num2str(round(third_order_h_best,3,'significant')),' , ', num2str(round(third_order_delta_best,3,'significant')),')'],'Location','NorthEast')

xlabel('$\nu$ Poissons Ratio','Interpreter','LaTex','FontSize',15)
ylabel('h Step size','Interpreter','LaTex','FontSize',15)
zlabel('$|\delta|$ Error','Interpreter','LaTex','FontSize',15)


view(-58,30)
shading interp
camlight right

set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

subplot(2,2,4)
surf(fourth_order_nu,fourth_order_h,fourth_order_delta,log(fourth_order_delta));
hold on 
p2 = plot3(fourth_order_nu_best,fourth_order_h_best,fourth_order_delta_best,'r.', 'MarkerSize',40);
title({'Crossvalidation of parameters  on 4th order FDM: $\nu$, h. '; '$|\delta| = K_{approx} - K_{exact}$'},'Interpreter','LaTex','FontSize',15)
legend( [p2],  ['Point (\nu,h,|\delta|) : ' '(',num2str(round(fourth_order_nu_best,3,'significant')),' , ',num2str(round(fourth_order_h_best,3,'significant')),' , ', num2str(round(fourth_order_delta_best,3,'significant')),')'],'Location','NorthEast')

xlabel('$\nu$ Poissons Ratio','Interpreter','LaTex','FontSize',15)
ylabel('h Step size','Interpreter','LaTex','FontSize',15)
zlabel('$|\delta|$ Error','Interpreter','LaTex','FontSize',15)


view(-58,30)
shading interp
camlight right

set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

print -depsc approximations_nu
savefig('approximations_nu.fig')

hold off
