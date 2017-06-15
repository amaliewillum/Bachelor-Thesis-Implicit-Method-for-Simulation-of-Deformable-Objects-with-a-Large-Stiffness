clear all 
close all 

% Fast E = 1000
dut1 = load('CV_E_lin1_e.mat'); %elapsed time: 208.105272

lin1    = dut1.CV_E_lin1_e;
h1      = lin1.h;
delta1  = lin1.delta;


dut2 = load('CV_E_lin1_e_2nd.mat'); %elapsed time: 411.660101

lin2    = dut2.CV_E_lin1_e;
h2      = lin2.h;
delta2  = lin2.delta;

% Fast E = 10.000
dut1_10000 = load('CV_E_lin1_e5000.mat'); %elapsed time: 203.586806

lin1_10000    = dut1_10000.CV_E_lin1_e5000;
h1_10000      = lin1_10000.h;
delta1_10000  = lin1_10000.delta;

dut2_10000 = load('CV_E_lin1_e_2nd5000.mat'); %elapsed time: 446.628374

lin2_10000    = dut2_10000.CV_E_lin1_e5000;
h2_10000      = lin2_10000.h;
delta2_10000  = lin2_10000.delta;

close all

figure(272)
loglog((h1),delta1,'-b','LineWidth',10)
hold on 
loglog((h2),delta2,'-r','LineWidth',5)

loglog((h1_10000),delta1_10000,'-g','LineWidth',10)
loglog((h2_10000),delta2_10000,'-y','LineWidth',5)

legend('1st order: E = 1000','2nd order: E = 1000','1st order: E = 10000','2nd order: E = 10000')
title('Fast E = 1000. Første og anden ordens approximation. logaritmiske akser')
xlabel('h')
ylabel('K_{exact} - K_{approx}')
set(gca, 'xdir','reverse')

print -depsc approx_test
savefig('approx_test.fig')

%% E = 1000 og E = 10.000

dut1 = load('CV_E_lin1_rand1100_h1e50.mat');

lin1_e       = dut1.CV_E_lin1_e5000_h1e50;
h1_rand1        = lin1_e.h;
delta1_rand1    = lin1_e.delta;

dut2 = load('CV_E_lin1_rand_2nd1100_h1e50.mat');

lin1_rand2       = dut2.CV_E_lin1_e5000_h1e50;
h1_rand2        = lin1_e.h;
delta1_rand2    = lin1_e.delta;

close all


figure(233)
loglog((h1_rand1),delta1_rand1(:,1),'-b','LineWidth',10)
hold on 
loglog((h1_rand2),delta1_rand2(:,1),'-r','LineWidth',5)
loglog((h1_rand1),delta1_rand1(:,2),'-b','LineWidth',10)
loglog((h1_rand2),delta1_rand2(:,2),'-r','LineWidth',5)

set(gca, 'xdir','reverse')

%% Lin  tau = 1.  100 E 100 h with random vectors
close all

dut1 = load('CV_E_lin1_rand.mat'); %first order. lin material. E: elapsed time 5740.131956

lin1_e       = dut1.CV_E_lin1_rand;
h_e_lin1        = lin1_e.h;
delta_e_lin1    = lin1_e.delta;
E_e_lin1         = lin1_e.E;

dut2 = load('CV_E_lin1_rand_2nd.mat'); %second order. lin material. E: elapsed time 8934.406306

lin1_e_lin2      = dut2.CV_E_lin1_rand;
h_rand_lin2        = lin1_e.h;
delta_e_lin2    = lin1_e.delta;
E_e_lin2         = lin1_e.E;




figure(233) % first order. lin material. E 
surf(E_e_lin2,h_e_lin1,delta_e_lin1,log(delta_e_lin1))
hold on 

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')
view(-60,25)

colormap spring 
shading interp
camlight headlight

hold off

figure(76) % second order. lin material. E 
surf(E_e_lin2,h_rand_lin2,delta_e_lin2,log(delta_e_lin2))
hold on 

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')
view(-60,25)

colormap winter
shading interp
camlight headlight

hold off


%% Svk tau = 1. 100 E 100 h with random vectors
close all

dut1 = load('CV_E_svk1_rand_1st.mat'); %first order. lin material. E: elapsed time 6149.574477

svk1_rand       = dut1.CV_E_svk1_rand;
h_rand_svk1        = svk1_rand.h;
delta_rand_svk1    = svk1_rand.delta;
E_rand_svk1         = svk1_rand.E;

dut2 = load('CV_E_svk1_rand_2nd.mat'); %second order. lin material. E: elapsed time 9434.781194

svk1_e2      = dut2.CV_E_svk1_rand;
h_e_svk2        = svk1_rand.h;
delta_e_svk2    = svk1_rand.delta;
E_rand_svk2         = svk1_rand.E;




figure(23333) % first order. lin material. E 
surf(E_rand_svk2,h_rand_svk1,delta_rand_svk1,log(delta_rand_svk1))
hold on 

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')
view(-60,25)

colormap spring 
shading interp
camlight headlight

hold off

figure(7226) % second order. lin material. E 
surf(E_rand_svk2,h_e_svk2,delta_e_svk2,log(delta_e_svk2))
hold on 

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')
view(-60,25)

colormap winter
shading interp
camlight headlight

hold off

%% This plot shows that Lin is a better result with random vectors
figure(098765432)
loglog(h_e_lin1,delta_e_lin1(:,1),'b','LineWidth',5)
hold on 
loglog(h_rand_lin2,delta_e_lin2(:,1),'g','LineWidth',3)
loglog(h_rand_svk1,delta_rand_svk1(:,1),'m','LineWidth',5)
loglog(h_e_svk2,delta_e_svk2(:,1),'c','LineWidth',3)
set(gca, 'xdir','reverse')
hold off

%% Lin  tau = 1.  100 E 100 h
close all

dut1 = load('CV_E_lin1_e_1st.mat'); %first order. lin material. E: elapsed time 5740.131956

lin1_e       = dut1.CV_E_lin1_e;
h_e_lin1        = lin1_e.h;
delta_e_lin1    = lin1_e.delta;
E_e_lin1         = lin1_e.E;

dut2 = load('CV_E_lin1_e_2nd.mat'); %second order. lin material. E: elapsed time 8934.406306

lin1_e_lin2      = dut2.CV_E_lin1_e;
h_e_lin2        = lin1_e.h;
delta_e_lin2    = lin1_e.delta;
E_e_lin2         = lin1_e.E;




figure(233) % first order. lin material. E 
surf(E_e_lin2,h_e_lin1,delta_e_lin1,log(delta_e_lin1))
hold on 

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')
view(-60,25)

colormap spring 
shading interp
camlight headlight

hold off

figure(76) % second order. lin material. E 
surf(E_e_lin2,h_e_lin2,delta_e_lin2,log(delta_e_lin2))
hold on 

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')
view(-60,25)

colormap winter
shading interp
camlight headlight

hold off


%% Svk tau = 1. 100 E 100 h with e-vectors
close all

dut1 = load('CV_E_svk1_e_1st.mat'); %first order. lin material. E: elapsed time 6149.574477

svk1_e       = dut1.CV_E_svk1_e;
h_e_svk1        = svk1_e.h;
delta_e_svk1    = svk1_e.delta;
E_e_svk1         = svk1_e.E;

dut2 = load('CV_E_svk1_e_2nd.mat'); %second order. lin material. E: elapsed time 9434.781194

svk1_e2      = dut2.CV_E_svk1_e;
h_e_svk2        = svk1_e.h;
delta_e_svk2    = svk1_e.delta;
E_e_svk2         = svk1_e.E;




figure(23333) % first order. lin material. E 
surf(E_e_svk2,h_e_svk1,delta_e_svk1,log(delta_e_svk1))
hold on 

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')
view(-60,25)

colormap spring 
shading interp
camlight headlight

hold off

figure(7226) % second order. lin material. E 
surf(E_e_svk2,h_e_svk2,delta_e_svk2,log(delta_e_svk2))
hold on 

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')
view(-60,25)

colormap winter
shading interp
camlight headlight

hold off

%% This plot shows that Lin is a better result with e-vectors
figure(098765432)
loglog(h_e_lin1,delta_e_lin1(:,1),'b','LineWidth',10)
hold on 
loglog(h_e_lin2,delta_e_lin2(:,1),'g','LineWidth',8)
loglog(h_e_svk1,delta_e_svk1(:,1),'m','LineWidth',5)
loglog(h_e_svk2,delta_e_svk2(:,1),'c','LineWidth',3)
set(gca, 'xdir','reverse')
hold off

%%
close all
derp = load('kappa_approx_svk_rand.mat');

kappa = derp.kappa;
kappa2 = derp.kappa2;
tau = derp.tau;

figure(8765)
plot(tau,kappa,'-b','linewidth',3)
hold on 
plot(tau,kappa2,'-c','linewidth',3)
title({'$ \kappa(\tau) = \overline{k}\left(\overline{x} + \tau \overline{u} \right)$';'Without logarithmic scaling.'},'Interpreter','Latex','FontSize',20)
xlabel('$\tau$','FontSize',15,'Interpreter','Latex')
ylabel('$\kappa (\tau)$','FontSize',15,'Interpreter','Latex')

grid minor
set(gca, 'xdir','reverse')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

% print -depsc kappa
% savefig('kappa.fig')
hold off



%% E = 1000 og E = 10.000
close all



dut1 = load('CV_E_lin1_rand1100.mat'); %elapsed time: 512.300011

lin1_e     = dut1.CV_E_lin1_e5000;
h1_rand      = lin1_e.h;
delta1_rand  = lin1_e.delta;

dut2 = load('CV_E_lin1_rand_2nd1100.mat'); %elapsed time: 1025.154811

lin2_rand     = dut2.CV_E_lin1_e5000;
h2_rand      = lin1_e.h;
delta2_rand  = lin1_e.delta;




figure(2)
loglog((h1_rand),delta1_rand(:,1),'-b','LineWidth',10)
hold on 
loglog((h2_rand),delta2_rand(:,1),'-r','LineWidth',5)
loglog((h1_rand),delta1_rand(:,2),'-b','LineWidth',10)
loglog((h2_rand),delta2_rand(:,2),'-r','LineWidth',5)

set(gca, 'xdir','reverse')

