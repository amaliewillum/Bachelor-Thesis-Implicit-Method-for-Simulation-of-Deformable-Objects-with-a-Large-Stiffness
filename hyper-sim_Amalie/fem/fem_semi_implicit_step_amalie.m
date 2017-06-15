function [ state, conv ] = fem_semi_implicit_step_amalie(dt, state, bcon, profile)
% Copyright 2011, Kenny Erleben

%--- Assemble all vector and matrices needed ------------------------------
k    =   state.k;                        % Elastic forces
f    = [ state.fx; state.fy; state.fz];  % External forces
p    = [ state.x; state.y; state.z];     % Current spatial position
v    = [ state.vx; state.vy; state.vz];  % Current spatial velocity
M    = state.M;
C    = state.C;



%--- Get information about boundary conditions ----------------------------
V      = length( state.x );      % Number of vertices in mesh
idx    = bcon.idx;               % Get indices of boundary conditions
values = bcon.values;            % Get values of the boundary conditions              
free   = setdiff( 1:3*V, idx );  % Get indices of non-boundary conditions


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--- Initialize 
c1      = 1;
c2      = dt;
c3      = dt^2;

%--- Apply boundary conditions --------------------------------------------
v(idx) = 0;
p(idx) = values;

%--- Partioning of system -------------------------------------------------
Mff     = M(free,free);
Cff     = C(free,free);

kf      = k(free);
Ff      = f(free);

Atmp    = M;
btmp    = M*v + dt*(f - k - C*v );

b       =  btmp(free)  - Atmp(free,idx)*v(idx);
% b       =  btmp  - Atmp*v;
A       =  sparse( Atmp(free,free) );
% A       =  sparse( Atmp );

% J = 1/(dt^2*beta) * M + gamma/(dt*beta) * C + K;
% Jff = 1/(dt^2*beta) * Mff + gamma/(dt*beta) * Cff + Kff;

% y = dt*kf.*phi_p(free) - phi_v(free);

%--- Matrix_vector_product
%     Ax = c1 * Mff * x;
%     Ax = Ax + c2 * Cff * x; 
%     Ax = Ax + c3 * (Ff - kf);


final = @(x) matrix_vector_mul(x,c1,c2,c3,Mff,Cff,Ff,kf);

% final = @(x) matrix_vector_mul(x,c1,c2,c3,M,C,f,k);


%--- Do velocity update ---------------------------------------------------
% v(free)    = A \ b;       % Direct method

% Use plain old conjugate gradient method
D =  sparse( diag(diag(A)) );

% [R,p] = chol(A); % p = 0 => A is positive definite
% all the eigenvalues of A is greater than zero.
% sym = issymmetric(A); % true => A is symmetric

% thing = v(free);

if profile
  % PCG return values:  X,FLAG,RELRES,ITER,RESVEC
    [v(free), FLAG ,~,iter,conv] = pcg(final, b, [], [], D);
%     disp([FLAG iter])
     
else
    [v(free), FLAG ,~,iter,~] = pcg(final, b, [], [], D); 
%     disp([FLAG iter])
    conv    = [];
end

% thing1 = v(free);
% disp([thing(1:10) thing1(1:10)])


%--- Do position update ---------------------------------------------------
% dut = p(free);

p(free)    = p(free) + dt*v(free);

% dut1 = p(free);
% disp([dut1(1:10) dut(1:10)])

%--- Store the updated values in state structure --------------------------
state.vx = v(1:V);
state.vy = v(V+1:2*V);
state.vz = v(2*V+1:end);
state.x  = p(1:V);
state.y  = p(V+1:2*V);
state.z  = p(2*V+1:end);

end

