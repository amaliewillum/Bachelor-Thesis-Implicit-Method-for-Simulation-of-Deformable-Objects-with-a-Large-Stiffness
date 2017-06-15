function [lal,dont] = duttelut(order,dt,state,mesh,params,material,tau,J,l,ei,vectors1)

        A  = @(xx) afunAx(xx, order, dt, state, mesh, params, material, tau);

        
        lal = 0;
        dont = 0;
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