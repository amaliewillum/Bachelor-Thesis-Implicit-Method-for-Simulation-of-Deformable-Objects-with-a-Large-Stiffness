clear all 
close all 
clc

res         = load('res_vectors.mat');
res_vecs    = res.res_vectors;
gmres_LU    = res_vecs.resvec_gmres_LU;
gmres       = res_vecs.resvec_gmres;
pcg_LU      = res_vecs.resvec_pcg_LU;
pcg         = res_vecs.resvec_pcg;

norm_phi = 0.1241;
%%
close all 

    % divideret med normen af phi
    figure(200)
    semilogy(0:length(gmres_LU)-1,gmres_LU/norm_phi,'-o');
    hold on
    semilogy(0:length(gmres)-1,gmres/norm_phi,'-o');
	semilogy(0:length(pcg_LU)-1,pcg_LU/norm_phi,'-o');
  	semilogy(0:length(pcg)-1,pcg/norm_phi,'-o');
    legend('gmres LU','gmres','pcg LU','pcg')
    grid on
    xlabel('Iteration number');
    ylabel('Relative residual');
    print -depsc dutte
    savefig('dutte.fig')
    
    hold off
    
    
    % ikke divideret med normen af phi
    figure(300)
    semilogy(0:length(gmres_LU)-1,gmres_LU,'-o');
    hold on
    semilogy(0:length(gmres)-1,gmres,'-o');
	semilogy(0:length(pcg_LU)-1,pcg_LU,'-o');
  	semilogy(0:length(pcg)-1,pcg,'-o');
    legend('gmres LU','gmres','pcg LU','pcg')
    grid on
    xlabel('Iteration number');
    ylabel('Relative residual');
    print -depsc duttel
    savefig('duttel.fig')
    
    hold off
   