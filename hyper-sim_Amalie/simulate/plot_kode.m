%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1st order FDM E %%%%%%%%%%%%%%%%%%%%%%%%%




% h = CV_E.h;
% h_best = CV_E.h_best;
% h_best1 = CV_E.h_best1;
% 
% E = CV_E.E;
% E_best = CV_E.E_best;
% E_best1 = CV_E.E_best1;
% 
% delta = CV_E.delta;
% delta_best = CV_E.delta_best;
% 
% delta1 = CV_E.delta1;
% delta_best1 =  CV_E.delta_best1;


% assignin('base', 'delta', delta)
% assignin('base', 'delta_best', delta_best)
% assignin('base', 'delta1', delta1)
% assignin('base', 'delta_best1', delta_best1)
% 
% assignin('base', 'E', E)
% assignin('base', 'E_best', E_best)
% assignin('base', 'E_best1', E_best1)
% 
% assignin('base', 'h', h)
% assignin('base', 'h_best', h_best)
% assignin('base', 'h_best1', h_best1)


% assignin('base', 'nu_best', nu_best)
% assignin('base', 'nu', nu)

% disp(['h = ',num2str(h_best),'. E = ',num2str(E_best),'. delta = ',num2str(delta_best)])
% 
% figure(1)
% surf(E,h,delta)
% hold on 
% plot3(E_best,h_best,delta_best,'r.', 'MarkerSize',40)
% text(E_best,h_best,delta_best,['(',num2str(E_best),',',num2str(h_best),',',num2str(delta_best),')'],'Color','red','FontSize',15)
% 
% title('Crossvalidation of parameters  on 1st order FDM: E, h. Point (E,h,$|\delta|$)','Interpreter','LaTex','FontSize',15)
% 
% xlabel('E Youngs Modulus','Interpreter','LaTex','FontSize',15)
% ylabel('h Step size','Interpreter','LaTex','FontSize',15)
% zlabel('$|\delta|$ Error','Interpreter','LaTex','FontSize',15)
% 
% view(233,26)
% 
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
% set(gca, 'ZScale', 'log')
% set(gca,'Ydir','reverse')
% 
% print -depsc crossvalE_1stFDM_neo
% savefig('crossvalE_1stFDM_neo.fig');
% hold off
% 
% kjh=lkjh;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2nd order FDM E %%%%%%%%%%%%%%%%%%%%%%%%%



% h = CV_E.h;
% h_best = CV_E.h_best;
% h_best1 = CV_E.h_best1;
% 
% E = CV_E.E;
% E_best = CV_E.E_best;
% E_best1 = CV_E.E_best1;
% 
% delta = CV_E.delta;
% delta_best = CV_E.delta_best;
% 
% delta1 = CV_E.delta1;
% delta_best1 =  CV_E.delta_best1;
% 
% 
% assignin('base', 'delta', delta)
% assignin('base', 'delta_best', delta_best)
% assignin('base', 'delta1', delta1)
% assignin('base', 'delta_best1', delta_best1)
% 
% assignin('base', 'E', E)
% assignin('base', 'E_best', E_best)
% assignin('base', 'E_best1', E_best1)
% 
% assignin('base', 'h', h)
% assignin('base', 'h_best', h_best)
% assignin('base', 'h_best1', h_best1)
% 
% 
% disp(['h = ',num2str(h_best),'. E = ',num2str(E_best),'. delta = ',num2str(delta_best)])
% 
% figure(2)
% surf(E,h,delta)
% hold on 
% plot3(E_best,h_best,delta_best,'r.', 'MarkerSize',40)
% text(E_best,h_best,delta_best,['(',num2str(E_best),',',num2str(h_best),',',num2str(delta_best),')'],'Color','red','FontSize',15)
% 
% title('Crossvalidation of parameters  on 2nd order FDM: E, h. Point (E,h,$|\delta|$)','Interpreter','LaTex','FontSize',15)
% 
% xlabel('E Youngs Modulus','Interpreter','LaTex','FontSize',15)
% ylabel('h Step size','Interpreter','LaTex','FontSize',15)
% zlabel('$|\delta|$ Error','Interpreter','LaTex','FontSize',15)
% 
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
% set(gca, 'ZScale', 'log')
% set(gca,'Ydir','reverse')
% 
% print -depsc crossvalE_2ndFDM_neo
% savefig('crossvalE_2ndFDM_neo.fig');
% hold off




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3rd order FDM E %%%%%%%%%%%%%%%%%%%%%%%%%


% h = CV_E.h;
% E = CV_E.E;
% delta = CV_E.delta;
% 
% E_best = CV_E.E_best;
% h_best = CV_E.h_best;
% delta_best = CV_E.delta_best;
% 
% 
% 
% disp(['h = ',num2str(h_best),'. E = ',num2str(E_best),'. delta = ',num2str(delta_best)])
% 
% figure(3)
% surf(E,h,delta)
% hold on 
% plot3(E_best,h_best,delta_best,'r.', 'MarkerSize',40)
% text(E_best,h_best,delta_best,['(',num2str(E_best),',',num2str(h_best),',',num2str(delta_best),')'],'Color','red','FontSize',15)
% 
% title('Crossvalidation of parameters  on 3rd order FDM: E, h. Point (E,h,$|\delta|$)','Interpreter','LaTex','FontSize',15)
% 
% xlabel('E Youngs Modulus','Interpreter','LaTex','FontSize',15)
% ylabel('h Step size','Interpreter','LaTex','FontSize',15)
% zlabel('$|\delta|$ Error','Interpreter','LaTex','FontSize',15)
% 
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
% set(gca, 'ZScale', 'log')
% set(gca,'Ydir','reverse')
% 
% print -depsc crossvalE_3rdFDM_neo
% savefig('crossvalE_3rdFDM_neo.fig');
% hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 4th order FDM E %%%%%%%%%%%%%%%%%%%%%%%%%

% h = CV_E.h;
% E = CV_E.E;
% delta = CV_E.delta;
% 
% E_best = CV_E.E_best;
% h_best = CV_E.h_best;
% delta_best = CV_E.delta_best;
% 
% disp(['h = ',num2str(h_best),'. E = ',num2str(E_best),'. delta = ',num2str(delta_best)])
% 
% figure(4)
% surf(E,h,delta)
% hold on 
% plot3(E_best,h_best,delta_best,'r.', 'MarkerSize',40)
% text(E_best,h_best,delta_best,['(',num2str(E_best),',',num2str(h_best),',',num2str(delta_best),')'],'Color','red','FontSize',15)
% 
% title('Crossvalidation of parameters  on 4th order FDM: E, h. Point (E,h,$|\delta|$)','Interpreter','LaTex','FontSize',15)
% 
% xlabel('E Youngs Modulus','Interpreter','LaTex','FontSize',15)
% ylabel('h Step size','Interpreter','LaTex','FontSize',15)
% zlabel('$|\delta| Error$','Interpreter','LaTex','FontSize',15)
% 
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
% set(gca, 'ZScale', 'log')
% set(gca,'Ydir','reverse')
% 
% print -depsc crossvalE_4thFDM_neo
% savefig('crossvalE_4thFDM_neo.fig');
% hold off



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1st order FDM nu %%%%%%%%%%%%%%%%%%%%%%%%


% h = CV_nu.h;
% nu = CV_nu.nu;
% delta = CV_nu.delta;
% 
% nu_best = CV_nu.nu_best;
% h_best = CV_nu.h_best;
% delta_best = CV_nu.delta_best;
% 
% disp(['h = ',num2str(h_best),'. nu = ',num2str(nu_best),'. delta = ',num2str(delta_best)])
% 
% figure(21)
% surf(nu,h,delta)
% hold on 
% plot3(nu_best,h_best,delta_best,'r.', 'MarkerSize',40)
% text(nu_best,h_best,delta_best,['(',num2str(nu_best),',',num2str(h_best),',',num2str(delta_best),')'],'Color','red','FontSize',20)
% 
% title('Crossvalidation of parameters on 1st order FDM: $\nu$, h. Point ($\nu$,h,$|\delta|$)','Interpreter','LaTex','FontSize',15)
% 
% xlabel('$\nu$ Poissons Ratio','Interpreter','LaTex','FontSize',15)
% ylabel('h Step size','Interpreter','LaTex','FontSize',15)
% zlabel('$|\delta|$ Error','Interpreter','LaTex','FontSize',15)
% 
% set(gca, 'ZScale', 'log')
% set(gca, 'YScale', 'log')
% set(gca,'Ydir','reverse')
% 
% print -depsc crossvalnu_1stFDM_neo
% savefig('crossvalnu_1stFDM_neo.fig');
% hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2nd order FDM nu %%%%%%%%%%%%%%%%%%%%%%%%




% h = CV_nu.h;
% nu = CV_nu.nu;
% delta = CV_nu.delta;
% 
% nu_best = CV_nu.nu_best;
% h_best = CV_nu.h_best;
% delta_best = CV_nu.delta_best;
% 
% disp(['h = ',num2str(h_best),'. nu = ',num2str(nu_best),'. delta = ',num2str(delta_best)])
% 
% figure(22)
% surf(nu,h,delta)
% hold on 
% plot3(nu_best,h_best,delta_best,'r.', 'MarkerSize',40)
% text(nu_best,h_best,delta_best,['(',num2str(nu_best),',',num2str(h_best),',',num2str(delta_best),')'],'Color','red','FontSize',20)
% 
% title('Crossvalidation of parameters on 2nd order FDM: $\nu$, h. Point ($\nu$,h,$|\delta|$)','Interpreter','LaTex','FontSize',15)
% 
% xlabel('$\nu$ Poissons Ratio','Interpreter','LaTex','FontSize',15)
% ylabel('h Step size','Interpreter','LaTex','FontSize',15)
% zlabel('$|\delta|$ Error','Interpreter','LaTex','FontSize',15)
% 
% set(gca, 'ZScale', 'log')
% set(gca, 'YScale', 'log')
% set(gca,'Ydir','reverse')
% 
% print -depsc crossvalnu_2ndFDM_neo
% savefig('crossvalnu_2ndFDM_neo.fig');
% hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3rd order FDM nu %%%%%%%%%%%%%%%%%%%%%%%%



% 
% h = CV_nu.h;
% nu = CV_nu.nu;
% delta = CV_nu.delta;
% 
% nu_best = CV_nu.nu_best;
% h_best = CV_nu.h_best;
% delta_best = CV_nu.delta_best;
% 
% disp(['h = ',num2str(h_best),'. nu = ',num2str(nu_best),'. delta = ',num2str(delta_best)])
% 
% figure(23)
% surf(nu,h,delta)
% hold on 
% plot3(nu_best,h_best,delta_best,'r.', 'MarkerSize',40)
% text(nu_best,h_best,delta_best,['(',num2str(nu_best),',',num2str(h_best),',',num2str(delta_best),')'],'Color','red','FontSize',20)
% 
% title('Crossvalidation of parameters on 3rd order FDM: $\nu$, h. Point ($\nu$,h,$|\delta|$)','Interpreter','LaTex','FontSize',15)
% 
% xlabel('$\nu$ Poissons Ratio','Interpreter','LaTex','FontSize',15)
% ylabel('h Step size','Interpreter','LaTex','FontSize',15)
% zlabel('$|\delta|$ Error','Interpreter','LaTex','FontSize',15)
% 
% set(gca, 'ZScale', 'log')
% set(gca, 'YScale', 'log')
% set(gca,'Ydir','reverse')
% 
% print -depsc crossvalnu_3rdFDM_neo
% savefig('crossvalnu_3rdFDM_neo.fig');
% hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 4th order FDM nu %%%%%%%%%%%%%%%%%%%%%%%%



% h = CV_nu.h;
% nu = CV_nu.nu;
% delta = CV_nu.delta;
% 
% nu_best = CV_nu.nu_best;
% h_best = CV_nu.h_best;
% delta_best = CV_nu.delta_best;
% 
% disp(['h = ',num2str(h_best),'. nu = ',num2str(nu_best),'. delta = ',num2str(delta_best)])
% 
% figure(24)
% surf(nu,h,delta)
% hold on 
% plot3(nu_best,h_best,delta_best,'r.', 'MarkerSize',40)
% text(nu_best,h_best,delta_best,['(',num2str(nu_best),',',num2str(h_best),',',num2str(delta_best),')'],'Color','red','FontSize',20)
% 
% title('Crossvalidation of parameters on 4th order FDM: $\nu$, h. Point ($\nu$,h,$|\delta|$)','Interpreter','LaTex','FontSize',15)
% 
% xlabel('$\nu$ Poissons Ratio','Interpreter','LaTex','FontSize',15)
% ylabel('h Step size','Interpreter','LaTex','FontSize',15)
% zlabel('$|\delta|$ Error','Interpreter','LaTex','FontSize',15)
% 
% set(gca, 'ZScale', 'log')
% set(gca, 'YScale', 'log')
% set(gca,'Ydir','reverse')
% 
% print -depsc crossvalnu_4thFDM_neo
% savefig('crossvalnu_4thFDM_neo.fig');