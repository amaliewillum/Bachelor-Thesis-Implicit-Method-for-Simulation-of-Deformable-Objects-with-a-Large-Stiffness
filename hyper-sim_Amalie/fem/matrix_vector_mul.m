function [Ax] = matrix_vector_mul(x,dt,c1,c2,c3,M,C,f,k)
    Ax = c1 * M * x;
    Ax = Ax + c2 * C * x; 
    Ax = Ax + c3 * (f - k);
    Ax = Ax * 1/dt^2;
end
