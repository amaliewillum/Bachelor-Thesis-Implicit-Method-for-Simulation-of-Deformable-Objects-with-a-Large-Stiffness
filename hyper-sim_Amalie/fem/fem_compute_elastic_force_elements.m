function [ ke ] = fem_compute_elastic_force_elements(mesh, state, params, material, h, E, nu,rho)
% Copyright 2011, Kenny Erleben

if h
    params.h_max = h;
    params.E = E;
    params.nu = nu;
    params.rho = rho;
end

%--- Get number of elements
cntT = length(mesh.T(:,1));

%--- Allocate space for local elastic force vector
ke = repmat(zeros(3,1),4,cntT);

%--- Get elastic material parameters
E  = params.E;    % Young modulus
nu = params.nu;   % Poisson ratio
%--- Convert E and nu into Lamé coefficients
lambda = (E*nu) / ( (1 + nu)*(1 - 2*nu) );
mu     =  E     / (        2*(1+nu)     );

%--- Now compute local elastic forces
for e=1:cntT

  %--- Pre computation of block indices -- these will be used again and again
  from = (e-1)*3 + 1;
  to   = e*3;
  
  %--- Get tetrahedron indices
  i = mesh.T(e,1);
  j = mesh.T(e,2);
  k = mesh.T(e,3);
  m = mesh.T(e,4);
  
  %--- Get spatial vertex coordinates
  Pei = [ state.x(i);  state.y(i);  state.z(i) ];
  Pej = [ state.x(j);  state.y(j);  state.z(j) ];
  Pek = [ state.x(k);  state.y(k);  state.z(k) ];
  Pem = [ state.x(m);  state.y(m);  state.z(m) ];
  
  %--- Define current spatial edge matrix
  E  = [  Pej-Pei, Pek-Pei, Pem-Pei ];
  
  %--- Compute deformation gradient
  %
  %    E = F * E0
  %
  % note inv(E0) is precomputed
  invE0 = mesh.invE0(:,from:to);
  
  Fe    = E*invE0;                     % The Deformation gradient
  
  F11 = Fe(1,1);    F12 = Fe(1,2);  F13 = Fe(1,3);
  F21 = Fe(2,1);    F22 = Fe(2,2);  F23 = Fe(2,3);
  F31 = Fe(3,1);    F32 = Fe(3,2);  F33 = Fe(3,3);
  
  switch material
      case 'svk'    % St. Venant Kirchhoff
          Pe = fem_compute_P_svk(F11,F12,F13,F21,F22,F23,F31,F32,F33,lambda,mu);
      case 'lin'    % Linear Isotropic 
          Pe = fem_compute_P_lin(F11,F12,F13,F21,F22,F23,F31,F32,F33,lambda,mu);
      case 'nh'     % Neo-Hookean
          Pe = fem_compute_P_nh(F11,F12,F13,F21,F22,F23,F31,F32,F33,lambda,mu);
  end
  
  kei = mesh.V(e) .* Pe * mesh.nabla_Ne(from:to,1);
  kej = mesh.V(e) .* Pe * mesh.nabla_Ne(from:to,2);
  kek = mesh.V(e) .* Pe * mesh.nabla_Ne(from:to,3);
  kem = mesh.V(e) .* Pe * mesh.nabla_Ne(from:to,4);
          
  ke(:, e) =  [kei; kej; kek; kem ];
end

end
