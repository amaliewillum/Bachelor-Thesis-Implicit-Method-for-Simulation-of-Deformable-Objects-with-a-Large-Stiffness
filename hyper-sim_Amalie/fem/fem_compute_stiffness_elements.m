function [ Ke ] = fem_compute_stiffness_elements(mesh, state, params, material, h, E, nu, rho)
% Copyright 2011, Kenny Erleben
if h
    params.h_max = h;
    params.E = E;
    params.nu = nu;
    params.rho = rho;
end

%--- Get number of elements
T = length(mesh.T(:,1));

%--- Allocate space for local elastic force vector
Ke = repmat(zeros(3,3),4,4*T);

%--- Get Lamé coefficients
lambda = (params.E*params.nu) / ( (1 + params.nu)*(1 - 2*params.nu) );
mu     =  params.E            / (        2*(1+params.nu)     );

%--- Now compute local elastic forces
for e=1:T
  
  %--- Pre computation of block indices -- these will be used again and again
  from = (e-1)*3 + 1;
  to   = e*3;
  
  i  = mesh.T(e,1);                              % Tetrahedron indices
  j  = mesh.T(e,2);
  k  = mesh.T(e,3);
  m  = mesh.T(e,4);
  
%   vol = mesh.V(e);                                % Volume
  
  Pi = [ state.x(i);  state.y(i);  state.z(i) ]; % Spatial vertex coords
  Pj = [ state.x(j);  state.y(j);  state.z(j) ];
  Pk = [ state.x(k);  state.y(k);  state.z(k) ];
  Pm = [ state.x(m);  state.y(m);  state.z(m) ];
  
  %--- Define current spatial edge matrix
  D  = [  Pj-Pi, Pk-Pi, Pm-Pi ];
  
  %--- Compute deformation gradient
  %
  %    D = F * D0
  %
  % note inv(D0) is precomputed
  invD0 = mesh.invE0(:,from:to);
  
  F = D*invD0;                     % The Deformation gradient
  
  
  % By definition A = d^2 psi/dF dF, we used symbolic differentiation and ccode to find a closed form formula for A
  F11 = F(1,1); F12 = F(1,2); F13 = F(1,3);
  F21 = F(2,1); F22 = F(2,2); F23 = F(2,3);
  F31 = F(3,1); F32 = F(3,2); F33 = F(3,3);
  %A = fem_compute_A(F11,F12,F13,F21,F22,F23,F31,F32,F33,lambda,mu);
  
  % Given  Ga, Gb, V, F, S, lambda, and nu then we need to compute....
  %
  % C(R,M,V,W)  = lambda*kron(R,M)*kron(V,W)
  %                  + mu*( kron(R,V)*kron(M,W) + kron(R,V)*kron(M,W) )
  %
  % J(V,W,I,J)  = kron(V,J)*F(W,I) + kron(W,J)*F(V,I)
  %
  % A(I,J,K,M)  = kron(I,K)*S(J,M) + (1/2) * sum_R  F(K,R)*(sum_V sum_W C(R,M,V,W)*J(V,W,I,J))
  %
  % Kab(I,K)    = sum_J sum_M  A(I,J,K,M) * Gb(M) * Ga(J) * V
  K = repmat(zeros(3,3),4,4);
  for i=1:4
    for j=1:4
      %      a  = mesh.T(e,i);
      %      b  = mesh.T(e,j);
      Ga =  mesh.nabla_Ne(from:to,i);
      Gb =  mesh.nabla_Ne(from:to,j);
      

      Ga1 = Ga(1); Ga2 = Ga(2); Ga3 = Ga(3);
      Gb1 = Gb(1); Gb2 = Gb(2); Gb3 = Gb(3);
      switch material 
          case 'svk'
            Kab = fem_compute_K_svk(F11,F12,F13,F21,F22,F23,F31,F32,F33,Ga1,Ga2,Ga3,Gb1,Gb2,Gb3,lambda,mu);
          case 'lin'
            Kab = fem_compute_K_lin(Ga1,Ga2,Ga3,Gb1,Gb2,Gb3,lambda,mu);
          case 'nh'
            Kab = fem_compute_K_nh(F11,F12,F13,F21,F22,F23,F31,F32,F33,Ga1,Ga2,Ga3,Gb1,Gb2,Gb3,lambda,mu);
      end

      i_from = (i-1)*3 + 1;
      j_from = (j-1)*3 + 1;
      i_to   = i*3;
      j_to   = j*3;
      
      K(i_from:i_to, j_from:j_to) =  Kab(:,:);
      
    end % Next j
  end % Next i
  
  % Store tangent stiffness matrices in element array
  from = (e-1)*12 + 1;
  to   = e*12;
  
  Ke(:, from:to) =  K;
end

end





