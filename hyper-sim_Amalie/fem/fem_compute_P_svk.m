function P = fem_compute_P_svk(F11,F12,F13,F21,F22,F23,F31,F32,F33,lambda,mu)
%COMPUTE_P_SVK
%    P = COMPUTE_P_SVK(F11,F12,F13,F21,F22,F23,F31,F32,F33,LAMBDA,MU)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    06-Jul-2016 16:45:17

t2 = F11.^2;
t3 = t2.*(1.0./2.0);
t4 = F21.^2;
t5 = t4.*(1.0./2.0);
t6 = F31.^2;
t7 = t6.*(1.0./2.0);
t8 = F12.^2;
t9 = t8.*(1.0./2.0);
t10 = F22.^2;
t11 = t10.*(1.0./2.0);
t12 = F32.^2;
t13 = t12.*(1.0./2.0);
t14 = F11.*F12.*(1.0./2.0);
t15 = F21.*F22.*(1.0./2.0);
t16 = F31.*F32.*(1.0./2.0);
t17 = t14+t15+t16;
t18 = F13.^2;
t19 = t18.*(1.0./2.0);
t20 = F23.^2;
t21 = t20.*(1.0./2.0);
t22 = F33.^2;
t23 = t22.*(1.0./2.0);
t24 = t3+t5+t7+t9+t11+t13+t19+t21+t23-3.0./2.0;
t25 = F11.*F13.*(1.0./2.0);
t26 = F21.*F23.*(1.0./2.0);
t27 = F31.*F33.*(1.0./2.0);
t28 = t25+t26+t27;
t29 = F12.*F13.*(1.0./2.0);
t30 = F22.*F23.*(1.0./2.0);
t31 = F32.*F33.*(1.0./2.0);
t32 = t29+t30+t31;
t33 = t3+t5+t7-1.0./2.0;
t34 = t9+t11+t13-1.0./2.0;
t35 = t19+t21+t23-1.0./2.0;
P = reshape([mu.*(F12.*t17.*2.0+F13.*t28.*2.0+F11.*t33.*2.0)+F11.*lambda.*t24,mu.*(F22.*t17.*2.0+F23.*t28.*2.0+F21.*t33.*2.0)+F21.*lambda.*t24,mu.*(F32.*t17.*2.0+F33.*t28.*2.0+F31.*t33.*2.0)+F31.*lambda.*t24,mu.*(F11.*t17.*2.0+F13.*t32.*2.0+F12.*t34.*2.0)+F12.*lambda.*t24,mu.*(F21.*t17.*2.0+F23.*t32.*2.0+F22.*t34.*2.0)+F22.*lambda.*t24,mu.*(F31.*t17.*2.0+F33.*t32.*2.0+F32.*t34.*2.0)+F32.*lambda.*t24,mu.*(F11.*t28.*2.0+F12.*t32.*2.0+F13.*t35.*2.0)+F13.*lambda.*t24,mu.*(F21.*t28.*2.0+F22.*t32.*2.0+F23.*t35.*2.0)+F23.*lambda.*t24,mu.*(F31.*t28.*2.0+F32.*t32.*2.0+F33.*t35.*2.0)+F33.*lambda.*t24],[3,3]);