% Material Elasticity Tensor
A = sym('TEMP')*ones(9,9);
for i=1:3
    for j=1:3
        for k=1:3
            for m=1:3
                r = (i-1)*3 + k;
                c = (j-1)*3 + m;
                A( r, c ) = diff( P(i,j) , F(k,m));
            end
        end
    end
end
A = simplify(A)
matlabFunction(A,'file','compute_A.m');