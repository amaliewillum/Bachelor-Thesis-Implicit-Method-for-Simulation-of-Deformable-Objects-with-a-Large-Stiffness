function [lal,dont] = K_k_loop(ei,material,order,l)

    addpath('../fem');


    switch material 
        case 'lin';
            Jacobian    = load('Jacobian_lin.mat');
%             K_file = load('stiffness_matrix_lin.mat');
%             K_lin = K_file.K;
            mesh        = Jacobian.mesh;
            state       = Jacobian.state;
            params      = Jacobian.params;
            k           = Jacobian.k;
        case 'svk';
            Jacobian    = load('Jacobian_svk.mat');
%             K_file = load('stiffness_matrix_svk.mat');
%             K_lin = K_file.K;
            mesh        = Jacobian.mesh;
            state       = Jacobian.state;
            params      = Jacobian.params;
            k           = Jacobian.k;
    end
    Jacobian    = load('Jacobian_svk.mat');
    mesh        = Jacobian.mesh;
    state       = Jacobian.state;
    params      = Jacobian.params;
    k           = Jacobian.k;
%     method = fem_method();

%     mesh_no = 1;
%     scene  = twist_create_scene(mesh_no);
    
    params = create_params([],'implicit_amalie'); 
    params.E = 1e+3;
    params.h_max = 1e-11;
    
    tau = 1;
    approx_K = @(xx) approximation_test(xx, order, state, mesh, params, material, k, tau, []);
    
    dont = 0;
    lal = 0;
    Ke     = fem_compute_stiffness_elements(mesh, state, params, material, []);
    K_lin  = fem_assemble_global_matrix( mesh, Ke );

    
    
    dt = params.h_max;
    
    parfor j = 1:size(K_lin,1)
        if mod(j,l) == 0
            lal = lal + 1;
%             ej = zeros(size(K_lin,1),1);
%             ej(j) = 1;
            ej = randn(size(K_lin,1),1);

%             tau = adaptive_step(ej, order, state, mesh, free, params, material, k, 10, []);

            Kij =  dt*(ei'*K_lin*ej);
            kij =  dt*(ei'*approx_K(ej));
            
            disp(abs(Kij - kij))
            
            if abs(Kij - kij) > 1e-5
                dont = dont + 1;
                
                
    %                     display(['Not a good approximation for J. ', num2str(i), ' ' ,num2str(j)])
            end
        end
    end
end