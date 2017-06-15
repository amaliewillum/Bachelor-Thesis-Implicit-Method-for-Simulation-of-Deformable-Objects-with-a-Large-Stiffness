close all

dut = load('CV_E_1storder');

first_order_CV_E1       = dut.CV_E1;
first_order_E           = first_order_CV_E1.E;
first_order_E_best      = first_order_CV_E1.E_best;
first_order_E_best1     = first_order_CV_E1.E_best1;
first_order_h           = first_order_CV_E1.h;
first_order_h_best      = first_order_CV_E1.h_best;
first_order_h_best1     = first_order_CV_E1.h_best1;
first_order_delta       = first_order_CV_E1.delta;
first_order_delta_best  = first_order_CV_E1.delta_best;
first_order_delta1      = first_order_CV_E1.delta1;
first_order_delta_best1 = first_order_CV_E1.delta_best1;




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

%%

dut = load('CV_E_2ndorder');

second_order_CV_E2       = dut.CV_E2;
second_order_E           = second_order_CV_E2.E;
second_order_E_best      = second_order_CV_E2.E_best;
second_order_E_best1     = second_order_CV_E2.E_best1;
second_order_h           = second_order_CV_E2.h;
second_order_h_best      = second_order_CV_E2.h_best;
second_order_h_best1     = second_order_CV_E2.h_best1;
second_order_delta       = second_order_CV_E2.delta;
second_order_delta_best  = second_order_CV_E2.delta_best;
second_order_delta1      = second_order_CV_E2.delta1;
second_order_delta_best1 = second_order_CV_E2.delta_best1;

dut = load('CV_E_3rdorder');

third_order_CV_E3       = dut.CV_E3;
third_order_E           = third_order_CV_E3.E;
third_order_E_best      = third_order_CV_E3.E_best;
third_order_E_best1     = third_order_CV_E3.E_best1;
third_order_h           = third_order_CV_E3.h;
third_order_h_best      = third_order_CV_E3.h_best;
third_order_h_best1     = third_order_CV_E3.h_best1;
third_order_delta       = third_order_CV_E3.delta;
third_order_delta_best  = third_order_CV_E3.delta_best;
third_order_delta1      = third_order_CV_E3.delta1;
third_order_delta_best1 = third_order_CV_E3.delta_best1;

dut = load('CV_E_4thorder');

fourth_order_CV_E4       = dut.CV_E4;
fourth_order_E           = fourth_order_CV_E4.E;
fourth_order_E_best      = fourth_order_CV_E4.E_best;
fourth_order_E_best1     = fourth_order_CV_E4.E_best1;
fourth_order_h           = fourth_order_CV_E4.h;
fourth_order_h_best      = fourth_order_CV_E4.h_best;
fourth_order_h_best1     = fourth_order_CV_E4.h_best1;
fourth_order_delta       = fourth_order_CV_E4.delta;
fourth_order_delta_best  = fourth_order_CV_E4.delta_best;
fourth_order_delta1      = fourth_order_CV_E4.delta1;
fourth_order_delta_best1 = fourth_order_CV_E4.delta_best1;



%%   HISTOGRAM plots for E

close all
clc
% val = delta1(:);

figure(1)
one = semilogxhist(first_order_delta1(:),'b');
hold on 
two = semilogxhist(second_order_delta1(:),'m');
three = semilogxhist(third_order_delta1(:),'g');
four = semilogxhist(fourth_order_delta1(:),'y');
legend([one,two,three,four],'1st order','2nd order','3rd order','4th order','Location','NorthWest')
title({'Error of FDMs of different orders from crossvalidation of $h$ and $E$'},'Interpreter','LaTex','FontSize',15)
xlabel('$|\delta| \frac{K_{approx}}{K_{exact}}$','FontSize',20,'Interpreter','LaTex')
ylabel('Number of occurences','Interpreter','LaTex')

print -depsc hist_approximations
savefig('hist_approximations.fig')
hold off

%%   SURFACE plots for E 

close all
clc

figh = figure(1);
pos = get(figh,'position');
set(figh,'position',[pos(1:2)/4 pos(3:4)*2])

subplot(2,2,1)
p1 = surf(first_order_E,first_order_h,first_order_delta,log(first_order_delta));
hold on 
p2 = plot3(first_order_E_best,first_order_h_best,first_order_delta_best,'r.', 'MarkerSize',40);
% text(E_best,h_best,delta_best,['(',num2str(E_best),',',num2str(h_best),',',num2str(delta_best),')'],'Color','red','FontSize',15)

title({'Crossvalidation of parameters  on 1st order FDM: E, h. '; '$|\delta| = K_{approx} - K_{exact}$'},'Interpreter','LaTex','FontSize',15)

% title('Crossvalidation of parameters  on 1st order FDM: E, h. Point (E,h,$|\delta|$)','Interpreter','LaTex','FontSize',15)
legend( [p2],  ['Point (E,h,|\delta|) : ' '(',num2str(round(first_order_E_best,3,'significant')),' , ',num2str(round(first_order_h_best,3,'significant')),' , ', num2str(round(first_order_delta_best,3,'significant')),')'],'Location','NorthEast')


xlabel('E Youngs Modulus','Interpreter','LaTex','FontSize',15)
ylabel('h Step size','Interpreter','LaTex','FontSize',15)
zlabel('$|\delta|$ Error','Interpreter','LaTex','FontSize',15)


view(-75,20)
colormap spring 
shading interp
camlight left
% camlight

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')



subplot(2,2,2)
p1 = surf(second_order_E,second_order_h,second_order_delta,log(second_order_delta));
hold on 
p2 = plot3(second_order_E_best,second_order_h_best,second_order_delta_best,'r.', 'MarkerSize',40);
title({'Crossvalidation of parameters  on 2nd order FDM: E, h. '; '$|\delta| = K_{approx} - K_{exact}$'},'Interpreter','LaTex','FontSize',15)
legend( [p2],  ['Point (E,h,|\delta|) : ' '(',num2str(round(second_order_E_best,3,'significant')),' , ',num2str(round(second_order_h_best,3,'significant')),' , ', num2str(round(second_order_delta_best,3,'significant')),')'],'Location','NorthEast')


xlabel('E Youngs Modulus','Interpreter','LaTex','FontSize',15)
ylabel('h Step size','Interpreter','LaTex','FontSize',15)
zlabel('$|\delta|$ Error','Interpreter','LaTex','FontSize',15)


view(-75,20)
colormap spring 
shading interp
camlight left
% camlight


set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')


subplot(2,2,3)
p1 = surf(third_order_E,third_order_h,third_order_delta,log(third_order_delta));
hold on 
p2 = plot3(third_order_E_best,third_order_h_best,third_order_delta_best,'r.', 'MarkerSize',40);
title({'Crossvalidation of parameters  on 3rd order FDM: E, h. '; '$|\delta| = K_{approx} - K_{exact}$'},'Interpreter','LaTex','FontSize',15)
legend( [p2],  ['Point (E,h,|\delta|) : ' '(',num2str(round(third_order_E_best,3,'significant')),' , ',num2str(round(third_order_h_best,3,'significant')),' , ', num2str(round(third_order_delta_best,3,'significant')),')'],'Location','NorthEast')

xlabel('E Youngs Modulus','Interpreter','LaTex','FontSize',15)
ylabel('h Step size','Interpreter','LaTex','FontSize',15)
zlabel('$|\delta|$ Error','Interpreter','LaTex','FontSize',15)


view(-75,20)
colormap spring 
shading interp
camlight left
% camlight


set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

subplot(2,2,4)
p1 = surf(fourth_order_E,fourth_order_h,fourth_order_delta,log(fourth_order_delta));
hold on 
p2 = plot3(fourth_order_E_best,fourth_order_h_best,fourth_order_delta_best,'r.', 'MarkerSize',40);
title({'Crossvalidation of parameters  on 4th order FDM: E, h. '; '$|\delta| = K_{approx} - K_{exact}$'},'Interpreter','LaTex','FontSize',15)
legend( [p2],  ['Point (E,h,|\delta|) : ' '(',num2str(round(fourth_order_E_best,3,'significant')),' , ',num2str(round(fourth_order_h_best,3,'significant')),' , ', num2str(round(fourth_order_delta_best,3,'significant')),')'],'Location','NorthEast')

xlabel('E Youngs Modulus','Interpreter','LaTex','FontSize',15)
ylabel('h Step size','Interpreter','LaTex','FontSize',15)
zlabel('$|\delta|$ Error','Interpreter','LaTex','FontSize',15)


view(-75,20)
colormap spring 
shading interp
camlight left
% camlight

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

print -depsc approximations
savefig('approximations.fig')

hold off