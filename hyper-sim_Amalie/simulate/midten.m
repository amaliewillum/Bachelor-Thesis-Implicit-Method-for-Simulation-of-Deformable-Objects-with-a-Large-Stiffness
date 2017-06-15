function [eps_inf1, eps_inf2] = midten(string, mesh, mesh_no, order, state, params, material, nu, h, E, rho, l, method, vector)
    addpath('../fem');    

    switch string
        case 'E'
            nu = 0.3; rho = 1000;
            eps_inf1 = zeros(size(E));
            eps_inf2 = zeros(size(E));
            parfor j = 1:length(E)
                ke = fem_compute_elastic_force_elements(mesh, state, params, material, h, E(j), nu, rho);
                k  = fem_assemble_global_vector(mesh, ke);

                Ke = fem_compute_stiffness_elements(mesh, state, params, material, h, E(j), nu, rho);
                K  = fem_assemble_global_matrix( mesh, Ke );


                [epso1, epso2] = indre(order, state, mesh, params, material, k, l, K, h, E(j), nu, rho, vector);

                eps_inf1(j) = norm(epso1,length(k));
                eps_inf2(j) = norm(epso2,length(k));
            end
        case 'nu'
            E = 1;
            rho = 1;
            eps_inf1 = zeros(size(nu));
            eps_inf2 = zeros(size(nu));
            parfor j = 1:length(nu)

                ke = fem_compute_elastic_force_elements(mesh, state, params, material, h, E, nu(j),rho);
                k  = fem_assemble_global_vector(mesh, ke);

                Ke = fem_compute_stiffness_elements(mesh, state, params, material, h, E, nu(j), rho);
                K  = fem_assemble_global_matrix( mesh, Ke );


                [epso1, epso2] = indre(order, state, mesh, params, material, k, l, K, h, E, nu(j), rho, vector);

                eps_inf1(j) = norm(epso1,length(k));
                eps_inf2(j) = norm(epso2,length(k));
            end
    
        case 'L'
            nu = 0.3;
            E = 1;
            rho = 1;

        %     min_edge_length = zeros(size(mesh_no));

            eps_inf1 = zeros(size(mesh_no));
            eps_inf2 = zeros(size(mesh_no));


            parfor j = 1:length(mesh_no)
                scene   = twist_create_scene(mesh_no(j),[]);
                dut = load(scene.meshfile);
                T = dut.T; X = dut.X; Y = dut.Y; Z = dut.Z;
                mesh    = method.create_mesh(T,X,Y,Z);
        %         min_edge_length(j) = dut.min_edge_length;
                state   = method.create_state(mesh, params);
                l       = (length(state.x)*3)/10;

                ke = fem_compute_elastic_force_elements(mesh, state, params, material, h, E, nu,rho);
                k  = fem_assemble_global_vector(mesh, ke);

                Ke = fem_compute_stiffness_elements(mesh, state, params, material, h, E, nu, rho);
                K  = fem_assemble_global_matrix( mesh, Ke );


                [epso1, epso2] = indre(order, state, mesh, params, material, k, l, K, h, E, nu, rho, vector);

                eps_inf1(j) = norm(epso1,inf);
                eps_inf2(j) = norm(epso2,inf);

            end
        case 'rho'
            E = 10e5;
            nu = 0.3;
            eps_inf1 = zeros(size(rho));
            eps_inf2 = zeros(size(rho));

            parfor j = 1:length(rho)

                ke = fem_compute_elastic_force_elements(mesh, state, params, material, h, E, nu, rho(j));
                k  = fem_assemble_global_vector(mesh, ke);

                Ke = fem_compute_stiffness_elements(mesh, state, params, material, h, E, nu, rho(j));
                K  = fem_assemble_global_matrix( mesh, Ke );


                [epso1, epso2] = indre(order, state, mesh, params, material, k, l, K, h, E, nu, rho(j), vector);

                eps_inf1(j) = norm(epso1,inf);
                eps_inf2(j) = norm(epso2,inf);
            end
    end
end


        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        