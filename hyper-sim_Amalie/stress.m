% First Piola--Kirchhoff stress tensor
P = sym('TEMP')*ones(3,3);
for i=1:3
    for j=1:3
        P(i,j) = diff(psi F(i,j));
    end
end
P = simplify(P)
matlabFunction(P,'file','compute_P.m')