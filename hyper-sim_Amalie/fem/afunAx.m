function Ax = afunAx(xx, order, dt, state, mesh, params, material, tau)
%         epsilon = 1e-6;
         


        M = state.M;
        C = state.C;
        tmp = state;
        ke = fem_compute_elastic_force_elements(mesh, tmp, params, material, []);
        k  = fem_assemble_global_vector(mesh, ke);

% %         V = length(state.x);
% %         dx = xx(1:3*V);
% %         dv = xx(3*V+1:end);

        dx = xx;
                  
        % k(u+delta(u))
        approx = @(dx,order,tau) approximation_test(dx, order, state, mesh, params, material, k, tau, []); 

        
% %         tau = adaptive_step(dx, order, state, mesh, free, params, material, k, 10, []);
        
        approx_k = approx(dx,order,tau);

% %         Iff = eye(size(M(free,free)));
% % 
% %         Ax1 = eye(size(M));
% % %         Ax2 = zeros(length(M),1);
% %         Ax3 = zeros(size(M));
% %         Ax4 = eye(size(M));
% % 
% %         Ax1(free,free) = Iff;
% %         A1 = Ax1 * dx;
% % 
% %         Ax2 = approx_k*dt;
% % 
% %         Ax3(free,free) = -dt * Iff ;
% %         A3 = Ax3 * dv;
% % 
% %         Ax4(free,free) = (M(free,free) + dt * Cff); 
% %         A4 = Ax4 * dv;
% % 
% %         Ax13 = A1 + A3;
% %         Ax24 = Ax2 + A4;
% % 
% % %         Ax = [Ax13 ; Ax24];
%         Ke  = fem_compute_stiffness_elements(mesh, tmp, params, material, []);
%         K   = fem_assemble_global_matrix( mesh, Ke );
        
        
%         Ax = M*dx + dt*C*dx + dt^2*approx_k;
        Ax = (M + dt*C)*dx + dt^2*approx_k;
end
    