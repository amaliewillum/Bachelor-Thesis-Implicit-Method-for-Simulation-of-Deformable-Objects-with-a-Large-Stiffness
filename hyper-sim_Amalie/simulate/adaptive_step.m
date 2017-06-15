function [tau] = adaptive_step(dx, order, state, mesh, free, params, material, k, maxTry, h, E, nu, rho)

        safe1 = 0.9;
        safe2 = 4.;
        err = 1e-5;
        tau = 1;
        
        for i = 1:maxTry
            approx = @(dx,order,tau) approximation_test(dx, order, state, mesh, free, params, material, k, tau, []); 

            half_tau = 0.5 * tau;

            xSmall = approx(dx,order,half_tau);
            xBig = approx(dx,order,tau);

            %* Compute the estimated truncation error
            scale = err * (abs(xSmall) + abs(xBig))/2.;
            xDiff = xSmall - xBig;
                        
            errorRatio = max( abs(xDiff)./(scale + eps) );

            %* Estimate new tau value (including safety factors)
            tau_old = tau;
            tau = safe1*tau_old*errorRatio^(-0.20);
            tau = max(tau,tau_old/safe2);
            tau = min(tau,safe2*tau_old);
            
            %* If error is acceptable, return computed values
            if (errorRatio < 1)  
%                 disp([errorRatio errorRatio^(-0.20) tau])
                return  
            end 
        end
        
        %         disp(tau)
        
        %* Issue error message if error bound never satisfied
%         error('ERROR: Adaptive Runge-Kutta routine failed');
end