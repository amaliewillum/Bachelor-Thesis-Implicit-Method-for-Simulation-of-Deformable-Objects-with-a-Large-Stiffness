function [epso1, epso2] = indre(order, state, mesh, params, material, k, l, K, h, E, nu, rho, vector)
    addpath('../fem');

    params.h = h;
    params.E = E;
    params.nu = nu;
    params.rho = rho;
        
    epso1 = zeros(size(k));
    tau = 1;
    epso2 = zeros(size(k));
    
    derp = load('squeeze1.mat');
    test_vectors = derp.test_vectors;
    mu = mean2(test_vectors);
    sigma = std2(test_vectors);
    
    parfor q = 1:length(k)
        if mod(q,l) == 0
            
            switch vector
                case 'e'
                    e = zeros(size(k));
                    e(q) = 1;
                case 'rand'
                    e = randn(size(k));
                case 'real'
                    e = mu.*randn(size(k)) + sigma;
            end
                        
            approx = @(xx,order,tau) approximation_test(xx, order, state, mesh, params, material, k, tau, h, E, nu, rho);
%             maxTry = 100;
%             tau = adaptive_step(e*h, order, state, mesh, free, params, material, k, maxTry, h, E, nu, rho);
%             tau = 1;  
            K_approx = approx(e*h, order, tau);
            K_exact = K*h*e;

            dut = abs(K_exact - K_approx);
            dut1 = abs(1-( K_approx./K_exact));

            epso1(q) = norm(dut,inf);
            epso2(q) = norm(dut1,inf);
        end
    end
end