function [approx] = approximation_test(xx, order, state, mesh, params, material, k, tau, h, E, nu, rho)
        epsilon = 0;
       
    if h
        h = params.h_max;
        E = params.E;
        nu = params.nu;
        rho = params.rho;
    end  
        
        V = length(state.x);
        dx = xx;
        p = [state.x ; state.y ; state.z];
        
        switch order
            case '1st'
                %%% first derivative - 1st order FDM... O(h)
                % k(u+delta(u))
                p1 = p + tau*dx;
                
                tmp = state;
                tmp.x = p1(1:V);
                tmp.y = p1(V+1:2*V);
                tmp.z = p1((2*V)+1:3*V);
                
                if h
                    h = params.h_max;
                    E = params.E;
                    nu = params.nu;
                    rho = params.rho;
                    
                    ke1 = fem_compute_elastic_force_elements(mesh, tmp, params, material, h, E, nu, rho);
                    k1 = fem_assemble_global_vector(mesh, ke1);
                else
                    ke1 = fem_compute_elastic_force_elements(mesh, tmp, params, material, []);
                    k1 = fem_assemble_global_vector(mesh, ke1);
                end
                
                

%                 approx = (norm(dx,2)./tau)*((k1 - k)./(tau*norm(dx,2)));
%                 approx = (k1-k)./(tau^2);
                approx = (k1-k)./(norm(dx,2));
%                 approx = ((k1 - k)./(tau*norm(dx,2)));
                

            case '2nd'
                %%% first derivative - 2nd order FDM... O(h^2)
                % k(u+delta(u))
                p1 = p + tau*dx;
                tmp = state;
                tmp.x = p1(1:V);
                tmp.y = p1(V+1:2*V);
                tmp.z = p1((2*V)+1:3*V);
                
                if h
                    h = params.h_max;
                    E = params.E;
                    nu = params.nu;
                    rho = params.rho;
                    
                    ke1 = fem_compute_elastic_force_elements(mesh, tmp, params, material, h, E, nu, rho);
                    k1 = fem_assemble_global_vector(mesh, ke1);
                else
                    ke1 = fem_compute_elastic_force_elements(mesh, tmp, params, material, []);
                    k1 = fem_assemble_global_vector(mesh, ke1);
                end
                    
                
                % k(u+ 2*delta(u))
                p2 = p + tau*2*dx;
                tmp = state;
                tmp.x = p2(1:V);
                tmp.y = p2(V+1:2*V);
                tmp.z = p2((2*V)+1:3*V);
                
                if h
                    h = params.h_max;
                    E = params.E;
                    nu = params.nu;
                    rho = params.rho;
                    
                    ke2 = fem_compute_elastic_force_elements(mesh, tmp, params, material, h, E, nu, rho);
                    k2 = fem_assemble_global_vector(mesh, ke2);
                else
                    ke2 = fem_compute_elastic_force_elements(mesh, tmp, params, material, []);
                    k2 = fem_assemble_global_vector(mesh, ke2);
                end
                
%                 approx = (norm(dx,2)./tau) .* ((-k2 + 4*k1 -3*k)./(2*tau*norm(dx,2)));
                approx = (-k2 + 4*k1 -3*k)./(2*tau^2);
%                 approx = ((-k2 + 4*k1 -3*k)./(2*tau*norm(dx,2)));
                
            case '3rd'
                %%% first derivative - 3rd order FDM... O(h^3)
                % k(u+delta(u))
                p = p + tau*dx;
                tmp.x = p(1:V);
                tmp.y = p(V+1:2*V);
                tmp.z = p((2*V)+1:3*V);
                tmp = tmp;
                ke1 = fem_compute_elastic_force_elements(mesh, tmp, params, h, E, nu);
                k1 = fem_assemble_global_vector(mesh, ke1);

                % k(u-delta(u))
%                 p2 = p - tau*dx;
%                 state.x = p2(1:V);
%                 state.y = p2(V+1:2*V);
%                 state.z = p2((2*V)+1:3*V);
%                 tmp = state;
%                 ke2 = fem_compute_elastic_force_elements(mesh, tmp, params, h, E, nu);
%                 k2 = fem_assemble_global_vector(mesh, ke2);

                % k(u+2delta(u))
                p3 = p + 2*tau*dx;
                tmp.x = p3(1:V);
                tmp.y = p3(V+1:2*V);
                tmp.z = p3((2*V)+1:3*V);
                tmp = tmp;
                ke3 = fem_compute_elastic_force_elements(mesh, tmp, params, h, E, nu);
                k3 = fem_assemble_global_vector(mesh, ke3);

                % k(u+3delta(u))
                p4 = p + 3*tau*dx;
                tmp.x = p4(1:V);
                tmp.y = p4(V+1:2*V);
                tmp.z = p4((2*V)+1:3*V);
                tmp = tmp;
                ke4 = fem_compute_elastic_force_elements(mesh, tmp, params, h, E, nu);
                k4 = fem_assemble_global_vector(mesh, ke4);
                
%                 approx = (k3 - 2*k1 + 2*k2 - k4) / (2*tau*norm(dx,2)^3+ epsilon);
                approx = (-11/6*k + 3*k1 - 3/2*k3 + 1/3*k4) / (tau * norm(dx,2) + epsilon);
            case '4th'
                %%% first derivative - 4th order FDM... O(h^4)
                % k(u+delta(u))
                p = p + tau*dx;
                tmp.x = p(1:V);
                tmp.y = p(V+1:2*V);
                tmp.z = p((2*V)+1:3*V);
                tmp = tmp;
                ke1 = fem_compute_elastic_force_elements(mesh, tmp, params, h, E, nu);
                k1 = fem_assemble_global_vector(mesh, ke1);

                % k(u-delta(u))
                p2 = p - tau*dx;
                tmp.x = p2(1:V);
                tmp.y = p2(V+1:2*V);
                tmp.z = p2((2*V)+1:3*V);
                tmp = tmp;
                ke2 = fem_compute_elastic_force_elements(mesh, tmp, params, h, E, nu);
                k2 = fem_assemble_global_vector(mesh, ke2);

                % k(u+2delta(u))
                p3 = p + 2*tau*dx;
                tmp.x = p3(1:V);
                tmp.y = p3(V+1:2*V);
                tmp.z = p3((2*V)+1:3*V);
                tmp = tmp;
                ke3 = fem_compute_elastic_force_elements(mesh, tmp, params, h, E, nu);
                k3 = fem_assemble_global_vector(mesh, ke3);

                % k(u-2delta(u))
                p4 = p - 2*tau*dx;
                tmp.x = p4(1:V);
                tmp.y = p4(V+1:2*V);
                tmp.z = p4((2*V)+1:3*V);
                tmp = tmp;
                ke4 = fem_compute_elastic_force_elements(mesh, tmp, params, h, E, nu);
                k4 = fem_assemble_global_vector(mesh, ke4);
                
                
                approx = (-k3 + 8*k1 - 8*k2 + k4) / (tau * 12*norm(dx,2) + epsilon);
        end

end 


