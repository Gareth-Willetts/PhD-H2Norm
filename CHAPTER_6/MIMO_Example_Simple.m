close all; clear all; clc;
%% Start again...
syms s
j = 3;
symb = (s+1)^(j);% * (s+1e-300);
cd1_poly = symb;
% cd1_poly = 2*s^10 + 10*s^9 + 45*s^8 + 120*s^7 + 210*s^6 + 252*s^5 + 210*s^4 + 120*s^3 + 45*s^2 + 10*s + 1;
cd1 = sym2poly(cd1_poly);
cn1 = ones(1,length(cd1)-1);
% cd1(1) = 2;

j = 4;
symb = (s+1)^(j);
cd2_poly = symb;
cd2 = sym2poly(symb);
cn2 = ones(1,length(cd2)-1);

j = 5;
symb = (s+1)^(j);
cd3_poly = symb;
cd3 = sym2poly(symb);
cn3 = ones(1,length(cd3)-1);

j = 6;
symb = (s+1)^(j);
cd4_poly = symb;
cd4 = sym2poly(symb);
cn4 = ones(1,length(cd4)-1);

sys1 = tf(cn1,cd1);		
sys2 = tf(cn2,cd2);
sys3 = tf(cn3,cd3);
sys4 = tf(cn4,cd4);
sys = [sys1,sys2;sys3,sys4]
norm(sys, 2)^2 % works!

%%
a_lcm = lcm([cd1_poly,cd2_poly,cd3_poly,cd4_poly]);
% quorem divides first argument by second argument
% Obtain quotient
[Q1,R1] = quorem(a_lcm, cd1_poly);
[Q2,R2] = quorem(a_lcm, cd2_poly);
[Q3,R3] = quorem(a_lcm, cd3_poly);
[Q4,R4] = quorem(a_lcm, cd4_poly);
% Multiply numerators by quotient
n1hat = conv(cn1, sym2poly(Q1));
n2hat = conv(cn2, sym2poly(Q2));
n3hat = conv(cn3, sym2poly(Q3));
n4hat = conv(cn4, sym2poly(Q4));
% Find coeffs of each polynomial
% n1hat_coeffs = coeffs(n1hat, s, 'All');
% length(n1hat_coeffs); % 18, so leading power is 17
% n2hat_coeffs = coeffs(n2hat, s, 'All');
% length(n2hat_coeffs); % 18, so leading power is 17

% How do we evaluate a polynomial at -s?
% a_lcm_coeffs = coeffs(a_lcm, s, 'All'); % coeffs in descending powers of s
% length(a_lcm_coeffs); % 21, so leading power of s is 20

%% Look at first TF
n1hat_coeffs = n1hat;
n2hat_coeffs = n2hat;
n3hat_coeffs = n3hat;
n4hat_coeffs = n4hat;
a_lcm_coeffs = sym2poly(a_lcm);
n = size(a_lcm_coeffs,2)-1;

if mod(n,2) == 1
    % ODD
    p2 = a_lcm_coeffs(2:2:end);
    p1 = [a_lcm_coeffs(1:2:end) 0];
else
    % EVEN
    p1 = a_lcm_coeffs(1:2:end);
    p2 = [a_lcm_coeffs(2:2:end) 0];
end

n1hat_coeffs = [zeros(1,n-length(n1hat_coeffs)) n1hat_coeffs];
n2hat_coeffs = [zeros(1,n-length(n2hat_coeffs)) n2hat_coeffs];
n3hat_coeffs = [zeros(1,n-length(n3hat_coeffs)) n3hat_coeffs];
n4hat_coeffs = [zeros(1,n-length(n4hat_coeffs)) n4hat_coeffs];

%% Create z(s) from equation in the paper
if mod(n,2) == 1
    % ODD
    % Define co and ce as the odd and even parts of the numerator polynomial,
    % respectively
    n1co = n1hat_coeffs(2:2:end);
    n1ce = n1hat_coeffs(1:2:end);
    n2co = n2hat_coeffs(2:2:end);
    n2ce = n2hat_coeffs(1:2:end);
    n3co = n3hat_coeffs(2:2:end);
    n3ce = n3hat_coeffs(1:2:end);
    n4co = n4hat_coeffs(2:2:end);
    n4ce = n4hat_coeffs(1:2:end);

    % Square these polynomials, and add zeros to start and end of co2
    n1co2 = conv(n1co, n1co);
    n1ce2 = conv(n1ce, n1ce);
    n1co2 = [0 n1co2 0];

    n2co2 = conv(n2co, n2co);
    n2ce2 = conv(n2ce, n2ce);
    n2co2 = [0 n2co2 0];

    n3co2 = conv(n3co, n3co);
    n3ce2 = conv(n3ce, n3ce);
    n3co2 = [0 n3co2 0];

    n4co2 = conv(n4co, n4co);
    n4ce2 = conv(n4ce, n4ce);
    n4co2 = [0 n4co2 0];

    % Preallocate memory for b, which will be populated with the coefficients
    % of z_0
    b = zeros(1,length(n1hat_coeffs));

    % Populates b with values of z_0
    for i = 1:length(n1hat_coeffs)
        b(i) = n1ce2(i) + n2ce2(i) - n1co2(i) - n2co2(i) + n3ce2(i) + n4ce2(i) - n3ce2(i) - n4ce2(i);
    end

    % Preallocate memory for z, change name of b to z0 (unnecessary)
    z0 = b;
else
    % EVEN
    % Define co and ce as the odd and even parts of the numerator polynomial,
    % respectively
    n1ce = n1hat_coeffs(2:2:end);
    n1co = n1hat_coeffs(1:2:end);
    n2ce = n2hat_coeffs(2:2:end);
    n2co = n2hat_coeffs(1:2:end);
    n3ce = n3hat_coeffs(2:2:end);
    n3co = n3hat_coeffs(1:2:end);
    n4ce = n4hat_coeffs(2:2:end);
    n4co = n4hat_coeffs(1:2:end);

    % Square these polynomials, and add zeros to start and end of co2
    n1co2 = conv(n1co, n1co);
    n1ce2 = conv(n1ce, n1ce);
    n2co2 = conv(n2co, n2co);
    n2ce2 = conv(n2ce, n2ce);
    n3co2 = conv(n3co, n3co);
    n3ce2 = conv(n3ce, n3ce);
    n4co2 = conv(n4co, n4co);
    n4ce2 = conv(n4ce, n4ce);
    %co2 = [0 co2 0];

    n1co2s = [n1co2, 0];
    n1ns = size(n1co2s,2);
    for i = n1ns+1:n
        n1co2s = [0,n1co2s];
    end
    n2co2s = [n2co2, 0];
    n2ns = size(n2co2s,2);
    for i = n2ns+1:n
        n2co2s = [0,n2co2s];
    end
    n3co2s = [n3co2, 0];
    n3ns = size(n3co2s,2);
    for i = n3ns+1:n
        n3co2s = [0,n3co2s];
    end
    n4co2s = [n4co2, 0];
    n4ns = size(n4co2s,2);
    for i = n4ns+1:n
        n4co2s = [0,n4co2s];
    end

    n1ns = size(n1ce2,2);
    for i = n1ns+1:n
        n1ce2 = [0,n1ce2];
    end

    n2ns = size(n2ce2,2);
    for i = n2ns+1:n
        n2ce2 = [0,n2ce2];
    end

    n3ns = size(n3ce2,2);
    for i = n3ns+1:n
        n3ce2 = [0,n3ce2];
    end

    n4ns = size(n4ce2,2);
    for i = n4ns+1:n
        n4ce2 = [0,n4ce2];
    end

    % Preallocate memory for b, which will be populated with the coefficients
    % of z_0
    b = zeros(1,length(n1hat_coeffs));

    % Populates b with values of z_0
    for i = 1:length(n1hat_coeffs)
        b(i) = n1ce2(i) + n2ce2(i) - n1co2s(i) - n2co2s(i) + n3ce2(i) + n4ce2(i) - n3co2s(i) - n4co2s(i);
    end

    % Preallocate memory for z, change name of b to z0 (unnecessary)
    z0 = b;
end

%% Variable definitions
nu = mean(abs(z0));
gamma = nu;
z0hat = z0/nu;
nu = mean(abs(p2));
gamma = gamma/nu;
p2hat = p2/nu;
% All correct

%% Now z_1
numZeros = length(z0hat) - length(p1);
z1tilde = z0hat - [p1 zeros(1,numZeros)]*z0hat(1)/p1(1);
nu = mean(abs(z1tilde(2:end)));
gamma = gamma*nu;
z1hat = z1tilde(2:end)/nu;

numZeros = length(p1) - length(p2hat);
p3tilde = p1 - [p2hat zeros(1,numZeros)]*p1(1)/p2hat(1);
nu = mean(abs(p3tilde(2:end)));
gamma = gamma/nu;
p3hat = p3tilde(2:end)/nu;
% All correct

for k = 3:n
    if k == 3
        numZeros = length(z1hat) - length(p2hat);
        zk1tilde = z1hat - [p2hat zeros(1,numZeros)]*z1hat(1)/p2hat(1);
        nu = mean(abs(zk1tilde(2:end)));
        gamma = gamma*nu;
        zk1hat = zk1tilde(2:end)/nu;

        numZeros = length(p2hat) - length(p3hat);
        pkplustilde = p2hat - [p3hat zeros(1,numZeros)]*p2hat(1)/p3hat(1);
        nu = mean(abs(pkplustilde(2:end)));
        gamma = gamma/nu;
        pkplushat = pkplustilde(2:end)/nu;

% Setting up variables for iterative calculations
        pkminushat = p3hat;
        pkhat = pkplushat;
    else
        numZeros = length(zk1hat) - length(pkminushat);
        zk1tilde = zk1hat - [pkminushat zeros(1,numZeros)]*zk1hat(1)/pkminushat(1);
        nu = mean(abs(zk1tilde(2:end)));
        gamma = gamma*nu;
        zk1hat = zk1tilde(2:end)/nu;

        numZeros = length(pkminushat) - length(pkhat);
        pkplustilde = pkminushat - [pkhat zeros(1,numZeros)]*pkminushat(1)/pkhat(1);
        nu = mean(abs(pkplustilde(2:end)));
        gamma = gamma/nu;
        pkplushat = pkplustilde(2:end)/nu;
        pkminushat = pkhat;
        pkhat = pkplushat;
    end
end
%% Calculate H2 norm
%normfunc = norm(tf(cn,cd),2)^2;
H2n2 = (p1(1) * gamma * zk1hat)/(2*a_lcm_coeffs(1) * pkplushat * pkminushat(1))