function A = fem_compute_A_lin(lambda,mu)
%COMPUTE_A_LIN
%    A = COMPUTE_A_LIN(LAMBDA,MU)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    06-Jul-2016 17:00:24

t2 = mu.*2.0;
t3 = lambda+t2;
A = reshape([t3,0.0,0.0,0.0,mu,0.0,0.0,0.0,mu,0.0,lambda,0.0,mu,0.0,0.0,0.0,0.0,0.0,0.0,0.0,lambda,0.0,0.0,0.0,mu,0.0,0.0,0.0,mu,0.0,lambda,0.0,0.0,0.0,0.0,0.0,mu,0.0,0.0,0.0,t3,0.0,0.0,0.0,mu,0.0,0.0,0.0,0.0,0.0,lambda,0.0,mu,0.0,0.0,0.0,mu,0.0,0.0,0.0,lambda,0.0,0.0,0.0,0.0,0.0,0.0,0.0,mu,0.0,lambda,0.0,mu,0.0,0.0,0.0,mu,0.0,0.0,0.0,t3],[9,9]);