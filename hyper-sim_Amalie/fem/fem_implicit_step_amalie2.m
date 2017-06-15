function [ state, conv, flag, iter, message ] = fem_implicit_step_amalie(dt, mesh, state, params, bcon, profile, iter_max, abs_tol, rel_tol, dir_tol, stag_tol)
% Copyright 2011, Kenny Erleben

%--- Check if arguments are valid or assing defaults ----------------------
if (nargin < 6)
  error('Not enough input arguments');
end

if(nargin < 7)
  iter_max = 30;
end
if (iter_max <= 0)
  iter_max = 30;
  display( strcat( 'Forcing iter max =', num2str(iter_max) ) );
end

if(nargin < 8)
  abs_tol = 1e-2;
end
if (abs_tol <= 0)
  abs_tol = 1e-2;
  display( strcat( 'Forcing abs. tol =', num2str(abs_tol) ) );
end

if(nargin < 9)
  rel_tol = 1e-5;
end
if (rel_tol <= 0)
  rel_tol = 1e-5;
  display( strcat( 'Forcing rel. tol =', num2str(rel_tol) ) );
end

if(nargin < 10)
  dir_tol = 1e-7;
end
if (dir_tol <= 0)
  dir_tol = 1e-7;
  display( strcat( 'Forcing dir. tol =', num2str(dir_tol) ) );
end

if(nargin < 11)
  stag_tol = eps*100;
end
if (stag_tol <= 0)
  stag_tol = eps*100;
  display( strcat( 'Forcing stag. tol =', num2str(stag_tol) ) );
end

conv            = [];

ITERATING       = 0;
ABSOLUTE        = 1;
RELATIVE        = 2;
STAGNATION      = 3;
NONDESCENT      = 4;
SMALL_DIRECTION = 5;
MAXITER         = 6;
LOCALMIN        = 7;

%--- Assemble all vector and matrices needed ------------------------------


fx    = state.fx;
fy    = state.fy;
fz    = state.fz;

p0    = [ state.x;  state.y;  state.z];   % Current spatial position
v0    = [ state.vx; state.vy; state.vz];  % Current spatial velocity

x0    = [p0;v0];
x     = x0;

f       = [ fx; fy; fz];                    % External forces
M       = state.M;                          % Mass matrix
C       = state.C;                          % Lumped matrix

%--- Get information about boundary conditions ----------------------------

V      = length( state.x );      % Number of vertices in mesh
idx    = bcon.idx;               % Get indices of boundary conditions
values = bcon.values;            % Get values of the boundary conditions              
free   = setdiff( 1:3*V, idx );  % Get indices of non-boundary conditions

tmp = state;
ke = fem_compute_elastic_force_elements(mesh, tmp, params);
k  = fem_assemble_global_vector(mesh, ke);

Ke     = fem_compute_stiffness_elements(mesh, tmp, params);
K      = fem_assemble_global_matrix( mesh, Ke );


p  = x(1:length(p0));
v  = x(length(p0)+1:end);

%--- Apply boundary conditions --------------------------------------------

v(idx) = 0;
p(idx) = values;

%--- Partioning of system -------------------------------------------------
Mff     = M(free,free);
Cff     = C(free,free);
Kff    = K(free,free);

vf      = v(free);
pf      = p(free);

v0f     = v0(free);
p0f     = p0(free);

Ff      = f(free);
kf      = k(free);

Iff     = eye(size(Mff));

%--- Assemble Jacobian ----------------------------------------------------

    function J = jacobian(dt,M,Mff,Cff,Iff,Kff,free)

        Jpp = eye(size(M));
        Jpv = zeros(size(M));
        Jvp = zeros(size(M));
        Jvv = eye(size(M));

        Jpp(free,free) = Iff;        Jpv(free,free) = -dt*Iff;
        Jvp(free,free) = dt*Kff;     Jvv(free,free) = (Mff + dt*Cff);

        J = [Jpp, Jpv; Jvp, Jvv];
    end

J = jacobian( dt, M, Mff, Cff, Iff, Kff, free);
J0 = jacobian( 0, M, Mff, Cff, Iff, Kff, free);

%--- Second Part ---------------------------------------------------------- 

J1 = (Mff + dt * Cff + dt^2 * Kff); 
J2 = (1/dt^2 * M + 1/dt * C + K);

%--- Assemble Phi ---------------------------------------------------------

phi_p       = zeros(size(p));
phi_p(free) = pf - p0f - dt*vf;
phi_v       = zeros(size(v));
phi_v(free) = (Mff + dt*Cff)*vf - Mff*v0f - dt*Ff + dt*kf;


y = 1/dt^2*[phi_p;phi_v];
% y = [phi_p;phi_v];


%--- Second Part ---------------------------------------------------------- 
y1 = dt * Kff * (pf - p0f - dt * vf) - ((Mff + dt * Cff) * vf - Mff * vf - dt * Ff + dt * kf);
y2 = -1/dt * ((Mff + dt * Cff) * vf - Mff * vf - dt * Ff + dt * kf);

y1 = (dt*K*phi_p) - phi_v;
y2 = -1/dt*phi_v;

%--- Matrix free ---------------------------------------------------------- 
% 
%     function Ax = afun(x,Iff,V,Ff,kf,Mff,Cff,dt,free)
%         
%         dx = x(1:3*V);
%         dv = x(3*V+1:end);
%         Ax13 = zeros(length(x(1:3*V)),1);
%         Ax24 = zeros(length(x(1:3*V)),1);
%         
% %         Ax1 = Iff * x(1:3*V);
%         Ax1 = Iff * dx(free);
%         Ax2 = dt*(Ff - kf);
% %         Ax3 = -dt * Iff * x(3*V+1:end);
%         Ax3 = -dt * Iff * dv(free);
% %         Ax4 = (Mff + dt * Cff) * x(3*V+1:end);
%         Ax4 = (Mff + dt * Cff) * dv(free);
%         
%         Ax13(free) = Ax1 + Ax3;
%         Ax24(free) = Ax2 + Ax4;
%         
%         Ax = [Ax13 ; Ax24];
%         Ax = Ax * 1/dt^2;
%         
%     end 
%--- Second Part ---------------------------------------------------------- 




% final_new = @(x) afun(x,Iff,V,Ff,kf,Mff,Cff,dt,free);

% D =  sparse( diag(diag(J)) );

% dut = pcg(final_new,-y,[],1000,D); %check-up - not working


% dut1 = final_new(x);

% dutdut = J\(-y);

% disp([dut(1000:1400) dutdut(1000:1400) dut1(1000:1400) -y(1000:1400)])

% ljsdlks=dsaljdk;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--- Old function
%--- Matrix_vector_product
%     Ax = c1 * Mff * x;
%     Ax = Ax + c2 * Cff * x; 
%     Ax = Ax + c3 * (Ff - kf);


c1 = 1; c2 = dt; c3 = dt^2;

D =  sparse( diag(diag(J2)) );

disp(size(y2))

final = @(x) matrix_vector_mul(x,dt,c1,c2,c3,M,C,f,k);
[delta2,flag,~,iter] = pcg(final,y1,[],2000,D);
derp=  pcg(final,y1,[],2000,D);
disp([flag iter])
% 
ljsdlks=dsaljdk;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



outer_iter     = 1;
max_outer_iter = 1;



%--- Incremental loading --------------------------------------------------
while 1,
  
  if( outer_iter > max_outer_iter)
    break;
  end
  
  load = outer_iter / max_outer_iter;
  
  state.fx = load*fx;
  state.fy = load*fy;
  state.fz = load*fz;
  
  outer_iter = outer_iter + 1;
  
  %--- Make a good prediction for starting iterate ------------------------
%   display('Gradient descent method...');
%   for i=1:0,
%       
%    b = fem_phi( x,  dt, mesh, state, params, bcon);
%    y2 = -1/dt * ((Mff + dt * Cff) * vf - Mff * vf - dt * Ff + dt * kf);
% 
% %    J = fem_nabla_phi( x,  0, mesh, state, params, bcon);
% %     b = y;
%     nabla_psi = J2' * y2;
%    
%    % Do a line search along the gradient descent direction
%    % in order to not over step
%    psi     = y2'*y2;
%    tau     = 1;
%    p0      = x(1:3*V);
%    v0      = x(3*V+1:end);
%    
%      pf = p0(free) - tau*nabla_psi  
%      x = x0 - tau*nabla_psi;
%      disp(size(J2))
%      disp(size(y2))
%      disp(size(nabla_psi))
%      disp(size(pf))
%      disp(size(x))
%      
%      
%      ljsdlks=dsaljdk;   
%    while 1,
%      pf = p0(free) - tau*nabla_psi  
%      x = x0 - tau*nabla_psi;
%      disp(size(J2))
%      disp(size(y2))
%      disp(size(nabla_psi))
%      disp(size(pf))
%      disp(size(x))
%      
%      
%      ljsdlks=dsaljdk;
%      b = fem_phi( x,  dt, mesh, state, params, bcon);
%      psi_tau = b'*b;
%      if(psi_tau < psi)
%        %display( strcat('  step length =', num2str(tau) ) );
%        %display( strcat('  psi =', num2str(psi_tau) ) );
%        break;
%      end
%      tau = tau*0.5;
%    end
%    
%    
%   end
  
  %--- Call the Newton method ---------------------------------------------
  display('Newton method...');
  flag  = ITERATING;
  iter  = 1;           % Iteration counter
  psi   = inf;
  while 1,
    
    %--- Test if we exceeded maximum iteration count ----------------------
    if(iter > iter_max)
      message = 'Exceeded max iteration limit';
      iter    = iter_max;
      flag    = MAXITER;
      break
    end
    
    %--- Get value of current iterative -----------------------------------
%     phi = fem_phi( x,  dt, mesh, state, params, bcon);
    phi = -1/dt * ((Mff + dt * Cff) * vf - Mff * vf - dt * Ff + dt * kf);
    
    %--- Test for absolute convergence ------------------------------------
    psi_old  = psi;
    psi      = norm(phi,2);
    if( profile )
      conv = [conv; psi];
    end
    if( psi < abs_tol )
      message = 'Absolute convergence';
      flag    = ABSOLUTE;
      break
    end
    
    %--- Test for relative convergence ------------------------------------
    rel_tst = abs(psi_old - psi)/ abs(psi);
    if( rel_tst < rel_tol )
%       display( strcat('  Relative convergence =', num2str( rel_tst )) );
      message = 'Relative convergence';
      flag    = RELATIVE;
      break
    end
        
    %--- Compute Jacobian -------------------------------------------------
%     nabla_phi = fem_nabla_phi( x,  dt, mesh, state, params, bcon);
    nabla_phi = J2; %just renaming
    
    %--- Compute Newton Direction -----------------------------------------
%     [delta ~] = gmres( nabla_phi, -phi);
    [delta,flag,~,iter] = pcg(final,y2,[],1000,D);
    disp([flag iter])
    
%     disp(size(delta))
%     disp(size(nabla_phi))
    
    
    %--- Test if we have a sufficient large Newton direction
    delta_norm = norm(delta, inf );
    if( delta_norm < dir_tol )  
%         disp([delta_norm dir_tol])
      message = 'Newton direction too small -- Im giving up';
      flag    = SMALL_DIRECTION;
      break
    end
        
    %--- Test for local minimum -------------------------------------------
    nabla_psi = 2 * delta' * nabla_phi;
    if( norm(nabla_psi,2) < rel_tol )
      message = 'Local minimum';
      flag    = LOCALMIN;
      break
    end
    
    %--- Test for descent direction ---------------------------------------
    dir_deriv = phi' * nabla_phi * delta;
    if (dir_deriv > 0)
      message = 'Non-descent direction -- Im giving up';
      flag    = NONDESCENT;
      break
    end
    
    %--- Do a line search -------------------------------------------------
    tau     = 1;        % Initial step length value
    beta    = 0.001;    % Sufficient decrease parameter
    alpha   = 0.5;      % Step length reduction parameter
    tau_min = eps*100;  % Minimum allowed step length
    dpsi    = dir_deriv*beta;
    %vol     = fem_compute_volume(x, mesh);
    while 1,
      
%       x_tau   = x + delta*tau;
    
      pf_tau = pf + tau * delta;
      vf_tau = vf + tau * (delta/dt);
      
      phi_tau = -1/dt * ((Mff + dt * Cff) * vf - Mff * vf - dt * Ff + dt * kf);
%       phi_tau = fem_phi( x_tau,  dt, mesh, state, params, bcon);
      psi_tau = norm(phi_tau,2);
      x_tau = [pf_tau;vf_tau];
      vol_tau = fem_compute_volume(x_tau, mesh);
      min_vol = min(vol_tau(:))>0;
      % ratio   = abs(vol_tau - vol) ./ vol;
      
      if (psi_tau < (psi + dpsi*tau)) && (min_vol > 0), % && (ratio < 0.1),
        display( strcat('  line-search: step length=', num2str(tau)) );
        break;
      end
      
      tau = tau*alpha;
      
      if tau < tau_min,
        display('  line-search: too small step length');
        break;
      end
      
    end
    
    %--- Do a Newton update -----------------------------------------------
%     x_old = x;
    pf_old = pf;
    vf_old = vf;

%     x     = x_tau;
    pf = pf_tau;
    vf = vf_tau;
    
    %--- Test if we have stagnation ---------------------------------------
    diff_p = max( norm( pf - pf_old, inf ) );
    diff_v = max( norm( vf - vf_old, inf ) );
    if( diff_p < stag_tol && diff_v < stag_tol )
      message = 'Stagnation of updates -- Im giving up';
      flag    = STAGNATION;
      break
    end
    
    iter = iter + 1;
  end
  
  %--- Debug output -------------------------------------------------------
  if( profile )
    display(message)
  end
  
end

%--- Put back the solution into state -------------------------------------
p(free)  = pf;
v(free)  = vf;

cntV     = length(mesh.x0);

state.x  = p(1:cntV);
state.y  = p(cntV+1:2*cntV);
state.z  = p(2*cntV+1:end);

state.vx = v(1:cntV);
state.vy = v(cntV+1:2*cntV);
state.vz = v(2*cntV+1:end);

end



