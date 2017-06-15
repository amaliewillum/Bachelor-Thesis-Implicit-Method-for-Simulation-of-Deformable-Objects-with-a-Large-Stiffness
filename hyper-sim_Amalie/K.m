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
matlabFunction(Kab,'file','compute_K.m');