clear all; close all; clc
%% Following works well for obtaining symbolic TF
% Define symbols for components in the suspension system
tic
syms ms mb mw ks kb kw cw z s cb cs x y t;
syms x1  x2 real
x = [x1;x2];
% Parameters for K1 and K2
K1 = x1;
K2 = x2;
% Suspension admittances
Q1 = (ks/s) + K1;
Q2 = (kb/s) + K2;
% State space representation of the system
A2 = [ms*s^2+Q1*s, -s*Q1, 0; -s*Q1, s^2*mb+s*(Q1+Q2), -s*Q2; 0, -s*Q2, mw*s^2+Q2*s+cw*s+kw];
B2 = [0; 0; cw*s+kw];
% Extract transfer function
TFT = [s, 0, 0]*inv(A2)*B2;
TFTF = TFT(1,1);
% toc
[n, d] = numden(TFTF);

%% Evaluate symbolic TF so they're functions of x and y
symbolic = 1;
if symbolic == 1
    ms = 3500; mb = 250; mw = 350; ks = 141E3; kb = 1260E3; cb = 7100; kw = 8E9; cw = 670E3; z = 200; cs = 8870;
    cd = coeffs(subs(d),s, 'All');
    cn = coeffs(subs(n),s,'All');
    cn = [zeros(1, length(cd) - length(cn) - 1) cn];
else
    ms = 3500; mb = 250; mw = 350; ks = 141E3; kb = 1260E3; cb = 7100; kw = 8E9; cw = 670E3; z = 200; cs = 8870;
    x1 = 7616;
    x2 = 11844;
    cd = sym2poly(subs(d));
    cn = sym2poly(subs(n));
    cn = [zeros(1, length(cd) - length(cn) - 1) cn];
end

%%
% for x1 = 1000:10:10000
%     for x2 = 10000:10:20000
%         ms = 3500; mb = 250; mw = 350; ks = 141E3; kb = 1260E3; cb = 7100; kw = 8E9; cw = 670E3; z = 200; cs = 8870; %x1 = 7616; x2 = 11844;
%         cd = sym2poly(subs(d));
%         cn = sym2poly(subs(n));
%         H2n2(x1-999, x2-9999) = norm(tf(cn,cd),2)^2;
%     end
% end
% toc

%%
time = [];
for i = 1:100
    tic
    x0 = [7000,11000];
    fun = @inbuilt_optim_H2;
    options = optimoptions('fminunc','SpecifyObjectiveGradient',false);
    [x,fval] = fminunc(fun,x0,options);
    time(i) = toc;
end