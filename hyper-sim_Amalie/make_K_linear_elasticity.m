syms lambda mu real;
syms F11 F12 F13...
    F21 F22 F23...
    F31 F32 F33 real;
F = [F11 F12 F13;...
    F21 F22 F23;...
    F31 F32 F33...
    ];
C = (F')*F;
E = (1/2) * (C-eye(3,3));  % Green Strain Tensor
epsilon =  (1/2) * (F' + F ) - eye(3,3); % Cauchy Strain Tensor

% Saint Venant Kirchoff Strain energy function
psi_svk = (1/2) * lambda* trace(E)^2 +...
    mu * trace(E*E);

% Compressible Linear Elasticity strain energy funciton
psi_lin = (1/2) * lambda* trace(epsilon)^2 +...
    mu * trace(epsilon*epsilon);

% First Piola--Kirchhoff stress tensor
P = sym('TEMP')*ones(3,3);
for i=1:3
    for j=1:3
        P(i,j) = diff(psi_svk, F(i,j));
    end
end
P = simplify(P)
matlabFunction(P,'file','compute_P_lin.m')

% Tangent stiffness element
syms Ga1 Ga2 Ga3 real;
syms Gb1 Gb2 Gb3 real;
Ga = [Ga1; Ga2; Ga3];
Gb = [Gb1; Gb2; Gb3];

Kab = sym('TEMP')*ones(3,3);
Kab(:,:) = 0;
for i=1:3
    for k=1:3
        for j=1:3
            for m=1:3
                r = (i-1)*3 + k;
                c = (j-1)*3 + m;
                Kab( i, k ) = Kab( i, k ) +...
                    A(r,c)*Gb(m)*Ga(j);
            end
        end
    end
end
Kab = simplify(Kab)
matlabFunction(Kab,'file','compute_K_svk.m');
