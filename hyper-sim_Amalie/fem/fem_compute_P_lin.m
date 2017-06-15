function P = fem_compute_P_lin(F11,F12,F13,F21,F22,F23,F31,F32,F33,lambda,mu)
%COMPUTE_P_LIN
%    P = COMPUTE_P_LIN(F11,F12,F13,F21,F22,F23,F31,F32,F33,LAMBDA,MU)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    06-Jul-2016 17:00:21

t2 = F11.*2.0;
t3 = F12+F21;
t4 = mu.*t3;
t5 = F22.*2.0;
t6 = F33.*2.0;
t7 = t2+t5+t6-6.0;
t8 = lambda.*t7.*(1.0./2.0);
t9 = F13+F31;
t10 = mu.*t9;
t11 = F23+F32;
t12 = mu.*t11;
P = reshape([t8+mu.*(t2-2.0),t4,t10,t4,t8+mu.*(t5-2.0),t12,t10,t12,t8+mu.*(t6-2.0)],[3,3]);
