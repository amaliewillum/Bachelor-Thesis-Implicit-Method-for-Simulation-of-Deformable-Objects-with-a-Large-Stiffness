function [crossval] = CrossvalValidation(string, order, mesh, params, material, state, vector)
    addpath('../fem');
    
    h                   = logspace(-1,-30,100);
    l                   = (length(state.x)*3)/10;
    E                   = 0;
    E_best              = 0;
    E_best1             = 0;
    nu                  = 0;
    nu_best             = 0;
    nu_best1            = 0;
    rho                 = 0;
    rho_best            = 0;
    rho_best1           = 0;
    min_edge_length     = 0;
    delta               = 0;
    delta1              = 0;
    delta_best          = 0;
    delta_best1         = 0;

    method = fem_method();
    
    
    switch string
        case 'E'
            E = logspace(3,15,10);
%             E = 1000;
            delta = zeros(length(h),length(E));
            delta1 = zeros(length(h),length(E));
            mesh_no = 0;
            parfor i = 1:length(h)
                [eps_inf1, eps_inf2] = midten(string, mesh, mesh_no, order, state, params, material, nu, h(i), E, rho, l, method, vector);

                delta(i,:) = eps_inf1;
                delta1(i,:) = eps_inf2;
            end
            
            [~, index] = min(delta(:));
            [row, col] = ind2sub(size(delta),index);
            
            E_best      = E(col);
            h_best      = h(row);
            delta_best  = delta(row,col);
            
            [~, index] = min(delta1(:));
            [row, col] = ind2sub(size(delta1),index);
            
            E_best1     = E(col);
            h_best1      = h(row);
            delta_best1  = delta1(row,col);
            
        case 'nu'
            nu = linspace(0,0.5,100);
            delta = zeros(length(h),length(nu));
            delta1 = zeros(length(h),length(nu));
            mesh_no = 0;

            parfor i = 1:length(h)
                [eps_inf1, eps_inf2] = midten(string, mesh, mesh_no, order, state, params, material, nu, h(i), E, rho, l, method, vector);

                delta(i,:) = eps_inf1;
                delta1(i,:) = eps_inf2;
            end
            
            [~, index] = min(delta(:));
            [row, col] = ind2sub(size(delta),index);
            
            nu_best     = nu(col);
            h_best      = h(row);
            delta_best  = delta(row,col);
            
            [~, index] = min(delta1(:));
            [row, col] = ind2sub(size(delta1),index);
            
            nu_best1     = nu(col);
            h_best1      = h(row);
            delta_best1  = delta1(row,col);
            
        case 'L'
            mesh_no = linspace(1,9,9);
            delta = zeros(length(h),length(mesh_no));
            delta1 = zeros(length(h),length(mesh_no));

            parfor i = 1:length(h)

                [eps_inf1, eps_inf2] = midten(string,mesh, mesh_no, order, state, params, material, nu, h(i), E, rho, l, method, vector);
 
                delta(i,:) = eps_inf1;
                delta1(i,:) = eps_inf2;

            end
            [~, index] = min(delta(:));
            [row, col] = ind2sub(size(delta),index);
            
            h_best      = h(row);
            delta_best  = delta(row,col);
            
            [~, index] = min(delta1(:));
            [row, col] = ind2sub(size(delta1),index);
            
            h_best1      = h(row);
            delta_best1  = delta1(row,col);
        case 'rho'
            rho = linspace(100,5000,100);
            delta = zeros(length(h),length(rho));
            delta1 = zeros(length(h),length(rho));
            mesh_no = 0;

            for i = 1:length(h);
                [eps_inf1, eps_inf2] = midten(string,mesh, mesh_no, order, state, params, material, 0, h(i), 0, rho, l, method, vector);
                
                delta(i,:) = eps_inf1;
                delta1(i,:) = eps_inf2;
            end
            [~, index] = min(delta(:));
            [row, col] = ind2sub(size(delta),index);
            
            rho_best     = rho(col);
            h_best      = h(row);
            delta_best  = delta(row,col);
            
            [~, index] = min(delta1(:));
            [row, col] = ind2sub(size(delta1),index);
            
            rho_best1     = rho(col);
            h_best1      = h(row);
            delta_best1  = delta1(row,col);
        
    end



    crossval = struct( 'E', E,...
                'E_best', E_best,...
                'nu', nu,...
                'nu_best', nu_best,...
                'h', h,...
                'h_best', h_best,...
                'delta', delta,...
                'delta_best', delta_best,...
                'nu_best1', nu_best1,...
                'E_best1', E_best1,...
                'h_best1', h_best1,...
                'rho', rho,...
                'rho_best', rho_best,...
                'rho_best1', rho_best1,...
                'delta1', delta1,...
                'min_edge_length',min_edge_length,...
                'delta_best1', delta_best1...
                );
end