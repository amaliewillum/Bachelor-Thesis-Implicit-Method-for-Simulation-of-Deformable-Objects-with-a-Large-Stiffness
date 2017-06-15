%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   SYMMETRY OF J AND A
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Symmetri J vs Approx svk 1st
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Test for approximation of J - symmetrisk med approx 
% St. Venant Kirchhoff 1st order, 2500 testvectors.
% Error percentage 5.48 % 
% 3240000 test vectors (indentity vectors) Error percentage: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear all 
close all 

addpath('../fem');
tic 
Jacobian    = load('Jacobian_svk.mat');
J           = Jacobian.J;
dt          = Jacobian.dt;
mesh        = Jacobian.mesh;
state       = Jacobian.state;
params      = Jacobian.params;
M           = Jacobian.M;
C           = Jacobian.C;
free        = Jacobian.free;
k           = Jacobian.k;

K_file = load('stiffness_matrix_svk.mat');
K = K_file.K;

material    = 'svk';

order       = '1st';
tau         = 1;

dt = 1e-5;

A   = @(xx) afunAx(xx, order, dt, state, mesh, params, material, tau);


dont = 0;
lal = 0;

J = (M + dt*C + dt^2*K);

dut = load('twist.mat');
vectors2 = dut.test_vectors;
dut1 = load('squeeze1.mat');
vectors1 = dut1.test_vectors;

l = size(J,1)/50;
for i = 1:size(J,1)
    if mod(i,l) == 0
%         ei = zeros(size(J,1),1);
%         ei(i) = 1;
        ei = mean2(vectors1).*randn(size(J,1),1) + std2(vectors1);
%         ei = randn(size(J,1),1);
        parfor j = 1:size(J,1)
            if mod(j,l) == 0
                lal = lal +1;
%                 ej = zeros(size(J,1),1);
%                 ej(j) = 1;
                ej = mean2(vectors1).*randn(size(J,1),1) + std2(vectors1);
%                 ej = randn(size(J,1),1);

                Jij = ei'*J*ej;
                Aij = ei'*A(ej);

                if abs(Aij - Jij) > 1e-5
                    dont = dont + 1;
%                     display(['Not a good approximation for J. ', num2str(i), ' ' ,num2str(j)])
                end
            end
        end
    end
end


err_svk_1st = dont/lal * 100;
disp(['Error Percentage svk 1st: ', num2str(err_svk_1st)])
toc
%% Symmetri J vs Approx svk 2nd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Test for approximation of J - symmetrisk med approx 
% St. Venant Kirchhoff 2nd order, 2500 testvectors.
% Error percentage  5.48 % 
% 3240000 test vectors (indentity vectors) Error percentage: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear all 
close all 

addpath('../fem');

Jacobian    = load('Jacobian_svk.mat');
J           = Jacobian.J;
dt          = Jacobian.dt;
mesh        = Jacobian.mesh;
state       = Jacobian.state;
params      = Jacobian.params;
M           = Jacobian.M;
C           = Jacobian.C;
free        = Jacobian.free;
k           = Jacobian.k;
% e         = eye(size(J));

material    = 'svk';

order       = '2nd';
tau         = 1;
final_new   = @(xx) afunAx(xx, order, dt, state, free, mesh, params, material, k, tau, [], M, C);

J = (M + dt*C + dt^2*K);

dont = 0;
lal = 0;

l = size(J,1)/10;
% l = size(J,1)/1800;
Aijj = [];
for i = 1:size(J,1)
    if mod(i,l) == 0
%         ei = zeros(size(J,1),1);
%         ei(i) = 1;
        ei = randn(size(J,1),1);
        for j = 1:size(J,1)
            if mod(j,l) == 0
                lal = lal +1;
%                 ej = zeros(size(J,1),1);
%                 ej(j) = 1;
                ej = randn(size(J,1),1);
        
                Jij = ei'*J*ej;
                Aijj = [Aijj; ei'*final_new(ej)];
        
                if abs(Aijj - Jij) > 1e-5
                    dont = dont + 1;
%                     display(['Not a good approximation for J. ', num2str(i), ' ' ,num2str(j)])
                end
            end
        end
    end
end


err_svk_2nd = dont/lal * 100;
disp(['Error Percentage svk 2nd: ', num2str(err_svk_2nd)])
%% Symmetri J vs Approx Lin 1st
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Test for approximation of J - symmetrisk med approx 
% Linear Isotropic Material 1st order, 2500 testvectors.
% Error percentage  5.48 % 
% 3240000 test vectors (indentity vectors) Error percentage: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all 
close all 

addpath('../fem');

Jacobian    = load('Jacobian_lin.mat');
J           = Jacobian.J;
dt          = Jacobian.dt;
mesh        = Jacobian.mesh;
state       = Jacobian.state;
params      = Jacobian.params;
M           = Jacobian.M;
C           = Jacobian.C;
% e         = eye(size(J));

K_file      = load('stiffness_matrix_svk.mat');
K           = K_file.K;

material    = 'lin';

order       = '1st';
tau         = 1;
A           = @(xx) afunAx(xx, order, dt, state, mesh, params, material, tau);


dont = 0;
lal = 0;

% dt = 1e-3;
% dt = 1e-5;

J = (M + dt*C + dt^2*K);

dut = load('twist.mat');
vectors2 = dut.test_vectors;
dut1 = load('squeeze1.mat');
vectors1 = dut1.test_vectors;

l = size(J,1)/50;
% l = size(J,1)/1800;

for i = 1:size(J,1)
    if mod(i,l) == 0
%         ei = zeros(size(J,1),1);
%         ei(i) = 1;
        ei = mean2(vectors1).*randn(size(J,1),1) + 3*std2(vectors1);
%         ei = randn(size(J,1),1);
        parfor j = 1:size(J,1)
            if mod(j,l) == 0
                lal = lal +1;
%                 ej = zeros(size(J,1),1);
%                 ej(i) = 1;
                ej = mean2(vectors1).*randn(size(J,1),1) + 3*std2(vectors1);
%                 ej = randn(size(J,1),1);
        
                Jij = ei'*J*ej;
                Aij = ei'*A(ej);
        
                if abs(Aij - Jij) > 1e-5
                    dont = dont + 1;
%                     display(['Not a good approximation for J. ', num2str(i), ' ' ,num2str(j)])
                end
            end
        end
    end
end


err_lin_1st = dont/lal * 100;
disp(['Error Percentage lin 1st: ', num2str(err_lin_1st)])
%% Symmetri J vs Approx Lin 2nd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Test for approximation of J - symmetrisk med approx 
% Linear Isotropic Material 2nd order, 2500 testvectors.
% Error percentage   5.48% 
% 3240000 test vectors (indentity vectors) Error percentage:  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all 
close all 

addpath('../fem');

Jacobian    = load('Jacobian_lin.mat');
J           = Jacobian.J;
dt          = Jacobian.dt;
mesh        = Jacobian.mesh;
state       = Jacobian.state;
params      = Jacobian.params;
M           = Jacobian.M;
C           = Jacobian.C;
free        = Jacobian.free;
k           = Jacobian.k;
% e         = eye(size(J));

material    = 'lin';

order       = '2nd';
tau         = 1;
final_new   = @(xx) afunAx(xx, order, dt, state, free, mesh, params, material, k, tau, [], M, C);


dont = 0;
lal = 0;

l = size(J,1)/50;
% l = size(J,1)/1800;

for i = 1:size(J,1)
    if mod(i,l) == 0
        ei = zeros(size(J,1),1);
        ei(i) = 1;
        for j = 1:size(J,1)
            if mod(j,l) == 0
                lal = lal +1;
                ej = zeros(size(J,1),1);
                ej(j) = 1;
        
                Jij = ei'*J*ej;
                Aij = ei'*final_new(ej);
        
                if abs(Aij - Jij) > 1e-5
                    dont = dont + 1;
%                     display(['Not a good approximation for J. ', num2str(i), ' ' ,num2str(j)])
                end
            end
        end
    end
end


err_lin_2nd = dont/lal * 100;
disp(['Error Percentage lin 2nd: ', num2str(err_lin_2nd)])
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   SYMMETRY OF K
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% K symmetrisk  YES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check if K (stiffness matrix) is symmetric  - K er symmetrisk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all 
close all 

addpath('../fem');

K_file = load('stiffness_matrix_svk.mat');
K_svk = K_file.K;

K_file = load('stiffness_matrix_lin.mat');
K_lin = K_file.K;

z = rand(length(K_svk),1);
y = rand(length(K_svk),1);

svk1 = z'*K_svk*y;
svk2 = y'*K_svk*z;

lin1 = z'*K_lin*y;
lin2 = y'*K_lin*z;

diff_svk = abs(svk1 - svk2);
diff_lin = abs(lin1 - lin2);

disp(['Error svk : ', num2str(diff_svk)])
disp(['Error lin : ', num2str(diff_lin)])
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   SYMMETRI OF K EXACT AND K APPROX
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Symmetri K vs Approx svk 1st
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Test for approximation of K - symmetrisk med approx 
% St. Venant Kirchhoff 1st order, 2500 testvectors .
% Error percentage   10.72% 
% 10.000 testvectors (identity vectors). Error percentage 5.62 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all 
close all 
tic
addpath('../fem');

K_file = load('stiffness_matrix_svk.mat');
K_svk = K_file.K;


material    = 'svk';
order       = '1st';

% l = size(K_svk,1)/100;
l = size(K_svk,1)/50;

hallo = @(ee) K_k_loop(ee,material,order,l);

dont = 0;
lal = 0;

params = create_params([],'implicit_amalie');
params.E = 1e+10;
params.h__max = 1e-5;

parfor i = 1:size(K_svk,1)
    if mod(i,l) == 0
%         ei = zeros(size(K_svk,1),1);
%         ei(i) = 1;
        ei = randn(size(K_svk,1),1);
        
        
        [lal1,dont1] = hallo(ei);

        
        lal = lal + lal1;
        dont = dont + dont1;
        
    end
end


err_K_svk_1st = dont/lal * 100;
disp(['Error Percentage K svk 1st: ', num2str(err_K_svk_1st)])
toc
%% Symmetri K vs Approx svk 2nd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Test for approximation of K - symmetrisk med approx 
% St. Venant Kirchhoff 2nd order, 2500 testvectors.
% Error percentage   10.72% 
% 10.000 testvectors (identity vectors). Error percentage 5.62 %
% random vectors. error percentage : 100 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all 
close all 

addpath('../fem');

K_file = load('stiffness_matrix_svk.mat');
K_svk = K_file.K;


material    = 'svk';
order       = '2nd';

l = size(K_svk,1)/10;

hallo = @(ee) K_k_loop(ee,material,order,l);

dont = 0;
lal = 0;


parfor i = 1:size(K_svk,1)
    if mod(i,l) == 0
%         ei = zeros(size(K_svk,1),1);
%         ei(i) = 1;
        ei = randn(size(K_svk,1),1);
        
        
        [lal1,dont1] = hallo(ei);

        
        lal = lal + lal1;
        dont = dont + dont1;
        
    end
end


err_K_svk_2nd = dont/lal * 100;
disp(['Error Percentage K svk 2nd: ', num2str(err_K_svk_2nd)])
%% Symmetri K vs Approx lin 1st
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Test for approximation of K if symmetric with approx 
% Linear Isotropic Material 1st order, 2500 testvectors.
% Error percentage   10.72% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all 
close all 

addpath('../fem');

K_file = load('stiffness_matrix_lin.mat');
K_lin = K_file.K;

material    = 'lin';
order       = '1st';

l = size(K_lin,1)/10;

hallo = @(ee) K_k_loop(ee,material,order,l);

dont = 0;
lal = 0;
params = create_params([],'implicit_amalie');
params.E = 1e+4;
params.h_max = 1e-10;

parfor i = 1:size(K_lin,1)
    if mod(i,l) == 0
%         ei = zeros(size(K_lin,1),1);
%         ei(i) = 1;
        ei = randn(size(K_lin,1),1);
        
        [lal1,dont1] = hallo(ei);

        
        lal = lal + lal1;
        dont = dont + dont1;
        
    end
end


err_K_lin_1st = dont/lal * 100;
disp(['Error Percentage K lin 1st: ', num2str(err_K_lin_1st)])
%% Symmetri K vs Approx lin 2nd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Test for approximation of K if symmetric with approx 
% Linear Isotropic Material 2nd order, 2500 testvectors.
% Error percentage   10.72% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all 
close all 

tic
addpath('../fem');

K_file = load('stiffness_matrix_lin.mat');
K_lin = K_file.K;

material    = 'lin';
order       = '2nd';

l = size(K_lin,1)/50;

hallo = @(ee) K_k_loop(ee,material,order,l);

dont = 0;
lal = 0;


parfor i = 1:size(K_lin,1)
    if mod(i,l) == 0
%         ei = zeros(size(K_lin,1),1);
%         ei(i) = 1;
        ei = randn(size(K_lin,1),1);
        
        [lal1,dont1] = hallo(ei);

        
        lal = lal + lal1;
        dont = dont + dont1;
        
    end
end

err_K_lin_2nd = dont/lal * 100;
disp(['Error Percentage K lin 2nd: ', num2str(err_K_lin_2nd)])
toc
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   IS KAPPA LINEAR??? 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test Kappa Svk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test if there is a reason to believe that the material is linear, and
% therefore no reason to do a 2nd order FDM - by step-reduction.
%If the 1st order FDM is linear, then it is not necessary to do a 2nd order
% St. Venant-Kirchhoff Material Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all 

addpath('../fem');

Jacobian    = load('Jacobian_svk.mat');
mesh        = Jacobian.mesh;
state       = Jacobian.state;
params      = Jacobian.params;
free        = Jacobian.free;
k           = Jacobian.k;

material    = 'svk';


tau = linspace(0,10,1000);
% tau = logspace(-10,0,100);

kappa_2_svk     = zeros(size(tau));
kappa_inf_svk   = zeros(size(tau));


dut = load('squeeze1.mat');
test_vectors = dut.test_vectors;

mu = mean2(test_vectors);
sigma = std2(test_vectors);

vector = mu.*randn(900,1) + sigma;


params.E = 500;
 

parfor i = 1:length(tau)
    tmp = state;
    p = [tmp.x;tmp.y;tmp.z];
    

    p1              = p + tau(i) * vector;
    tmp.x           = p1(1:length(p1)/3);
    tmp.y           = p1(length(p1)/3 +1:length(p1)/3 * 2);
    tmp.z           = p1(length(p1)/3 * 2 +1:end);

    ke              = fem_compute_elastic_force_elements(mesh, tmp, params, material,[]);
    k1              = fem_assemble_global_vector(mesh, ke);    
        
    kappa_2_svk(i)      = norm(k1,2);
    kappa_inf_svk(i)    = norm(k1,inf);
    
    
end
%%
close all

figure('units','normalized','position',[.25,.35,.45,.55])

p1=plot(tau(1:4:end),kappa_2(1:4:end),'b.-','MarkerSize',20);
hold on
p2 = plot(tau(1:4:end),kappa_inf(1:4:end),'r.-','MarkerSize',20);

L = legend([p1,p2],'$\Vert \kappa \Vert_2$', '$\Vert \kappa \Vert_\infty$');
set(L,'Interpreter','latex','FontSize',15)

% th = title('Elastic Forces Fucntion $\kappa(\tau) =   \vec{k}(\vec{u} + \tau \Delta \vec x )$ for St. Venant-Kirchhoff Material','Interpreter','Latex','FontSize',18);
% titlePos = get( th , 'position');
% titlePos(2) = 47000;
% set( th , 'position' , titlePos);

ylabel('$\Vert \kappa \Vert$','Interpreter','Latex','FontSize',15)
xlabel('$-\tau$','Interpreter','Latex','FontSize',15)


set(gca, 'xdir','reverse')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

% print( '-deps' ,'ku_svk')

print -depsc ku_svkloglog
savefig('ku_svkloglog.fig');
hold off         

%% Test Kappa Lin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test if there is a reason to believe that the material is linear, and
% therefore no reason to do a 2nd order FDM - by step-reduction.
%If the 1st order FDM is linear, then it is not necessary to do a 2nd order
% Linear Isotropic Material Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

addpath('../fem');

Jacobian    = load('Jacobian_lin.mat');
mesh        = Jacobian.mesh;
state       = Jacobian.state;
params      = Jacobian.params;
free        = Jacobian.free;
k           = Jacobian.k;
dt          = Jacobian.dt;

material    = 'lin';

tau = linspace(0,10,1000);
% tau = logspace(-10,0,100);

kappa_2     = zeros(size(tau));
kappa_inf   = zeros(size(tau));


dut = load('squeeze1.mat');
test_vectors = dut.test_vectors;

mu = mean2(test_vectors);
sigma = std2(test_vectors);

vector = mu.*randn(900,1) + sigma;


params.E = 500;
 

parfor i = 1:length(tau)
    tmp = state;
    p = [tmp.x;tmp.y;tmp.z];
    

    p1              = p + tau(i) * vector;
    tmp.x           = p1(1:length(p1)/3);
    tmp.y           = p1(length(p1)/3 +1:length(p1)/3 * 2);
    tmp.z           = p1(length(p1)/3 * 2 +1:end);

    ke              = fem_compute_elastic_force_elements(mesh, tmp, params, material,[]);
    k1              = fem_assemble_global_vector(mesh, ke);    
        
    kappa_2(i)      = norm(k1,2);
    kappa_inf(i)    = norm(k1,inf);
    
    
end

%%

A1 = (sum(tau.^2)*sum(kappa_inf_svk) - sum(tau)*sum(tau.*kappa_inf_svk))/(1000*sum(tau.^2) - sum(tau)^2);
A2 = (sum(tau.^2)*sum(kappa_inf) - sum(tau)*sum(tau.*kappa_inf))/(1000*sum(tau.^2) - sum(tau)^2);

B1 = (1000* sum(tau.*kappa_inf_svk) - sum(tau)*sum(kappa_inf_svk))/(1000*sum(tau.^2) - sum(tau)^2);
B2 = (1000* sum(tau.*kappa_inf) - sum(tau)*sum(kappa_inf))/(1000*sum(tau.^2) - sum(tau)^2);

%%
sigma1 = sqrt(1/1000 * sum( kappa_inf_svk - A1 - B1*tau).^2)
sigma2 = sqrt(1/1000 * sum( kappa_inf - A2 - B2*tau).^2)

%%
close all
% mdl = fitlm( tau, kappa_2,'linear')
mdl1 = fitlm( tau, kappa_inf,'linear')
% mdl = fitlm( tau, kappa_2_svk,'linear')
mdl2 = fitlm( tau, kappa_inf_svk,'linear')
%%
figure; 
plotResiduals(mdl,'probability');

%%
close all

colors = fliplr(hsv(8));

figure('units','normalized','position',[.25,.35,.45,.55])

% p1 = plot(tau,kappa_2,'Color',colors(8,:),'linewidth',5);
hold on
p3 = plot(tau,kappa_inf,'Color',colors(7,:),'linewidth',6);
% p4 = plot(tau, 3.4057* tau + (-1.0028e-13) ,'--','color',colors(3,:),'linewidth',4);
p4 = plot(tau, B2*tau + A2 ,'--','color',colors(3,:),'linewidth',4);
% p2 = plot(tau, 24.518* tau + (-1.44e-13) ,'--','color',colors(4,:),'linewidth',2);

% p5 = plot(tau,kappa_2_svk,'Color',colors(7,:),'linewidth',5);
hold on
p7 = plot(tau,kappa_inf_svk,'Color',colors(2,:),'linewidth',6);
% p6 = plot(tau, 25.767* tau + (-1.2281) ,'--','color',colors(5,:),'linewidth',2);
% p8 = plot(tau, 4.5269* tau + (0.084999) ,'--','color',colors(6,:),'linewidth',4);
p8 = plot(tau, B1* tau + A1 ,'--','color',colors(6,:),'linewidth',4);

title('The elastic forces as a funtion of $\tau$ with linear least squares (LLS) fit $\kappa = \alpha + \beta \tau$ ','interpreter','latex','fontsize',15)

axis([0,10,0,50])
ylabel('$\Vert \kappa \Vert$','Interpreter','Latex','FontSize',20)
xlabel('$-\tau$','Interpreter','Latex','FontSize',15)
box on
grid on
l1 = legend([p3,p4,p7,p8],...
    '$\Vert \kappa_{lin} \Vert_\infty$','LLS fit ',...
    '$\Vert \kappa_{svk} \Vert_\infty$','LLS fit'...
    );
set(l1, 'interpreter','latex','location','northeast','fontsize',15)
set(gca, 'xdir','reverse')

ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

fname= '/users/amaliewillum/dropbox/2016_WILLUM/rapport';
filename = 'kappa_u_linear';


print -depsc kappa_u_linear
print( '-depsc' ,fullfile(fname, filename))
hold off 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Den relative fejl
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%














%% E x h: CV_E_nh_real_squeeze_1st_h301
clear all 
close all

addpath('../fem');

load('CV_E_nh_real_squeeze_1st_h301.mat')

% Get the field names of the structure.
fields = fieldnames(CV_E_nh_real, '-full');
% Find out how many there are - for our loop.
numberOfFields = length(fields);
for f = 1 : numberOfFields
    thisField = fields{f};
    commandLine = sprintf('%s = CV_E_nh_real.%s', thisField, thisField);
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
%% Plot CV_E_nh_real_squeeze_1st_h301_delta.png
close all

figure('units','normalized','position',[.22,.3,.42,.5])
p1 = surf(E,h,delta,log(delta));
hold on 
p2 = plot3([E(1) E(end)],[h_thing(1) h_thing(end)],[delta_thing(1) delta_thing(end)],'k-','LineWidth',3);
% p2 = plot3(E,h_thing,delta_thing,'k-','LineWidth',3);
p3 = plot3(E(1),h_thing(1),delta_thing(1),'r.','MarkerSize',40);
p4 = plot3(E(end),h_thing(end),delta_thing(end),'b.','MarkerSize',40);

L = legend([p3,p4],['( h , E ) = ( ' num2str(h_thing(1),'%.2e') ' , ' num2str(E(1),'%.1e') ' )'   ],...
    ['( h , E ) = ( ' num2str(h_thing(end),'%.2e') ' , ' num2str(E(end),'%.1e') ' )' ],...
    'Location','NorthEast');
set(L,'Interpreter','latex','FontSize',12)

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

title('Absolute Error','interpreter','latex','fontsize',18)
xlabel('Youngs Modulus $E$','interpreter','latex','fontsize',15)
ylabel('Stepsize $h$','interpreter','latex','fontsize',15)
zlabel('Error $\delta$','interpreter','latex','fontsize',15)

% axis([E(1),E(end),h(34),h(1),0,10^8])


colormap spring
% colormap cool 
% colorbar eastoutside
% shading flat

box on

shading interp 
lighting gouraud
camlight infinite
material dull

view(-50,30)

print( '-dpng' ,'CV_E_nh_real_squeeze_1st_h301_delta')
hold off
%% Plot CV_E_nh_real_squeeze_1st_h301_delta1.png
close all


% gennemsnit af relative error 538.1714
nr = 4;
E_tmp = E(1:nr:end);
h_tmp = h(1:nr:end);
delta1_tmp = delta1(1:nr:end,1:nr:end);

figure('units','normalized','position',[.22,.3,.42,.5])
p1 = surf(E_tmp,h_tmp,delta1_tmp,log(delta1_tmp));
hold on 
% p2 = plot3([E(1) E(end)],[h_thing(1) h_thing(end)],[delta_thing(1) delta_thing(end)],'k-','LineWidth',3);
% % p2 = plot3(E,h_thing,delta_thing,'k-','LineWidth',3);
% p3 = plot3(E(1),h_thing(1),delta_thing(1),'r.','MarkerSize',40);
% p4 = plot3(E(end),h_thing(end),delta_thing(end),'b.','MarkerSize',40);
% 
% L = legend([p3,p4],['( h , E ) = ( ' num2str(h_thing(1),'%.2e') ' , ' num2str(E(1),'%.1e') ' )'   ],...
%     ['( h , E ) = ( ' num2str(h_thing(end),'%.2e') ' , ' num2str(E(end),'%.1e') ' )' ],...
%     'Location','NorthEast');
% set(L,'Interpreter','latex','FontSize',12)

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

title('Relative Error','interpreter','latex','fontsize',18)
xlabel('Youngs Modulus $E$','interpreter','latex','fontsize',15)
ylabel('Stepsize $h$','interpreter','latex','fontsize',15)
zlabel('Error $\delta$','interpreter','latex','fontsize',15)

% axis([E(1),E(end),h(34),h(1),0,10^8])


colormap spring
% colormap cool 
% colorbar eastoutside
% shading flat

box on

shading interp 
lighting gouraud
camlight('headlight','infinite')
material dull

view(-50,30)

print( '-dpng' ,'CV_E_nh_real_squeeze_1st_h301_delta1')
hold off


%% E x h: CV_E_nh_e_squeeze_1st_h301
clear all 
close all

addpath('../fem');

load('CV_E_nh_e_squeeze_1st_h301.mat')

% Get the field names of the structure.
fields = fieldnames(CV_E_nh_real, '-full');
% Find out how many there are - for our loop.
numberOfFields = length(fields);
for f = 1 : numberOfFields
    thisField = fields{f};
    commandLine = sprintf('%s = CV_E_nh_real.%s', thisField, thisField);
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
%% Plot CV_E_nh_e_squeeze_1st_h301_delta.png
close all

figure('units','normalized','position',[.22,.3,.42,.5])
p1 = surf(E,h,delta,log(delta));
hold on 
p2 = plot3([E(1) E(end)],[h_thing(1) h_thing(end)],[delta_thing(1) delta_thing(end)],'k-','LineWidth',3);
% p2 = plot3(E,h_thing,delta_thing,'k-','LineWidth',3);
p3 = plot3(E(1),h_thing(1),delta_thing(1),'r.','MarkerSize',40);
p4 = plot3(E(end),h_thing(end),delta_thing(end),'b.','MarkerSize',40);

L = legend([p3,p4],['( h , E ) = ( ' num2str(h_thing(1),'%.2e') ' , ' num2str(E(1),'%.1e') ' )'   ],...
    ['( h , E ) = ( ' num2str(h_thing(end),'%.2e') ' , ' num2str(E(end),'%.1e') ' )' ],...
    'Location','NorthEast');
set(L,'Interpreter','latex','FontSize',12)

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

title('Absolute Error','interpreter','latex','fontsize',18)
xlabel('Youngs Modulus $E$','interpreter','latex','fontsize',15)
ylabel('Stepsize $h$','interpreter','latex','fontsize',15)
zlabel('Error $\delta$','interpreter','latex','fontsize',15)

% axis([E(1),E(end),h(34),h(1),0,10^8])


colormap spring
% colormap cool 
% colorbar eastoutside
% shading flat

box on

shading interp 
lighting gouraud
camlight headlight infinite
material dull

view(-50,30)

print( '-dpng' ,'CV_E_nh_e_squeeze_1st_h301_delta')
hold off
%% Plot CV_E_nh_e_squeeze_1st_h301_delta1.png
close all

figure('units','normalized','position',[.22,.3,.42,.5])
p1 = surf(E,h,delta1,log(delta1));
hold on 
% p2 = plot3([E(1) E(end)],[h_thing(1) h_thing(end)],[delta_thing(1) delta_thing(end)],'k-','LineWidth',3);
% % p2 = plot3(E,h_thing,delta_thing,'k-','LineWidth',3);
% p3 = plot3(E(1),h_thing(1),delta_thing(1),'r.','MarkerSize',40);
% p4 = plot3(E(end),h_thing(end),delta_thing(end),'b.','MarkerSize',40);
% 
% L = legend([p3,p4],['( h , E ) = ( ' num2str(h_thing(1),'%.2e') ' , ' num2str(E(1),'%.1e') ' )'   ],...
%     ['( h , E ) = ( ' num2str(h_thing(end),'%.2e') ' , ' num2str(E(end),'%.1e') ' )' ],...
%     'Location','NorthEast');
% set(L,'Interpreter','latex','FontSize',12)

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

title('Absolute Error','interpreter','latex','fontsize',18)
xlabel('Youngs Modulus $E$','interpreter','latex','fontsize',15)
ylabel('Stepsize $h$','interpreter','latex','fontsize',15)
zlabel('Error $\delta$','interpreter','latex','fontsize',15)

% axis([E(1),E(end),h(34),h(1),0,10^8])


colormap spring
% colormap cool 
% colorbar eastoutside
% shading flat

box on

shading interp 
lighting gouraud
camlight infinite
material shiny

view(-50,30)

print( '-dpng' ,'CV_E_nh_e_squeeze_1st_h301_delta1')
hold off


%% E x h: CV_E_nh_rand_squeeze_1st_h301
clear all 
close all

addpath('../fem');

load('CV_E_nh_rand_squeeze_1st_h301.mat')

% Get the field names of the structure.
fields = fieldnames(CV_E_nh_real, '-full');
% Find out how many there are - for our loop.
numberOfFields = length(fields);
for f = 1 : numberOfFields
    thisField = fields{f};
    commandLine = sprintf('%s = CV_E_nh_real.%s', thisField, thisField);
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
%% Plot CV_E_nh_rand_squeeze_1st_h301_delta.png
close all

figure('units','normalized','position',[.22,.3,.42,.5])
p1 = surf(E,h,delta,log(delta));
hold on 
p2 = plot3([E(1) E(end)],[h_thing(1) h_thing(end)],[delta_thing(1) delta_thing(end)],'k-','LineWidth',3);
% p2 = plot3(E,h_thing,delta_thing,'k-','LineWidth',3);
p3 = plot3(E(1),h_thing(1),delta_thing(1),'r.','MarkerSize',40);
p4 = plot3(E(end),h_thing(end),delta_thing(end),'b.','MarkerSize',40);

L = legend([p3,p4],['( h , E ) = ( ' num2str(h_thing(1),'%.2e') ' , ' num2str(E(1),'%.1e') ' )'   ],...
    ['( h , E ) = ( ' num2str(h_thing(end),'%.2e') ' , ' num2str(E(end),'%.1e') ' )' ],...
    'Location','NorthEast');
set(L,'Interpreter','latex','FontSize',12)

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

title('Absolute Error','interpreter','latex','fontsize',18)
xlabel('Youngs Modulus $E$','interpreter','latex','fontsize',15)
ylabel('Stepsize $h$','interpreter','latex','fontsize',15)
zlabel('Error $\delta$','interpreter','latex','fontsize',15)

% axis([E(1),E(end),h(34),h(1),0,10^8])


colormap spring
% colormap cool 
% colorbar eastoutside
% shading flat

box on

shading interp 
lighting gouraud
camlight headlight infinite
material dull

view(-50,30)

print( '-dpng' ,'CV_E_nh_rand_squeeze_1st_h301_delta')
hold off
%% Plot CV_E_nh_rand_squeeze_1st_h301_delta1.png
close all

% gennemsnit delta1 = 1.8904e + 6

nr = 4;
E_tmp = E(1:nr:end);
h_tmp = h(1:nr:end);
delta1_tmp = delta1(1:nr:end,1:nr:end);


figure('units','normalized','position',[.22,.3,.42,.5])
p1 = surf(E_tmp,h_tmp,delta1_tmp,log(delta1_tmp));
hold on 
% p2 = plot3([E(1) E(end)],[h_thing(1) h_thing(end)],[delta_thing(1) delta_thing(end)],'k-','LineWidth',3);
% % p2 = plot3(E,h_thing,delta_thing,'k-','LineWidth',3);
% p3 = plot3(E(1),h_thing(1),delta_thing(1),'r.','MarkerSize',40);
% p4 = plot3(E(end),h_thing(end),delta_thing(end),'b.','MarkerSize',40);
% 
% L = legend([p3,p4],['( h , E ) = ( ' num2str(h_thing(1),'%.2e') ' , ' num2str(E(1),'%.1e') ' )'   ],...
%     ['( h , E ) = ( ' num2str(h_thing(end),'%.2e') ' , ' num2str(E(end),'%.1e') ' )' ],...
%     'Location','NorthEast');
% set(L,'Interpreter','latex','FontSize',12)

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')

title('Relative Error','interpreter','latex','fontsize',18)
xlabel('Youngs Modulus $E$','interpreter','latex','fontsize',15)
ylabel('Stepsize $h$','interpreter','latex','fontsize',15)
zlabel('Error $\delta$','interpreter','latex','fontsize',15)

% axis([E(1),E(end),h(34),h(1),0,10^8])


colormap spring
% colormap cool 
% colorbar eastoutside
% shading flat

box on

shading interp 
lighting gouraud
camlight headlight infinite
material dull

view(-50,30)

print( '-dpng' ,'CV_E_nh_rand_squeeze_1st_h301_delta1')
hold off






