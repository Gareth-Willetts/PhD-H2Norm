%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc;                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define symbols for mechanical system
syms d1 d2 d3 d4 % Deltas -> these are our inputs
% syms zv theta zb1 zb2 % polynomial coefficients
syms s %k b c % mv Ivz lvx mb k bs c % G1 G2 % mechanical variables
% syms Kpz Cpz Ksz Krz Crz Kaz
% Define values for variables - symbolic version takes too long
k = 5*10^6; b = 92; c = 11.049;
mv = 38000; mb = 2500;
lvx = 9.5; Ivz = 2.31*10^6;
Kpz = 4.935*10^6; Ksz = 1.016*10^6; Krz = 0.508*10^6; Kaz = 0;
Cpz = 50.74*10^3; Crz = 64.11*10^3;
V = 55; Av = 2.5*10^(-7); lbx = 1.25;

G1 = 2*Kpz/s + 2*Cpz;
G0 = 1/((s/k) + 1/(b*s+c)); % ??? % 1/((s/k) + 1/(bs + c)) Layout S4

gamma_k = (Krz*(Ksz + Kaz) + Ksz*Kaz)/Crz;
gamma_c = Ksz + Kaz;
alpha = (Ksz + Krz)/Crz;

Gairspring = (gamma_k/s + gamma_c)/(s + alpha);
G2 = Gairspring + G0;
% Rearrange mechanical equations and put into matrix form
F = [(mv*s^2 + 2*G2*s), 0, -G2*s, -G2*s;
    0, (Ivz*s^2 + 2*lvx^2*G2*s), -lvx*G2*s, lvx*G2*s;
    -G2*s, (-lvx*G2*s), (mb*s^2 + G2*s + 2*G1*s), 0;
    -G2*s, (lvx*G2*s), 0, (mb*s^2 + 2*G1*s + G2*s)];

G = [0, 0, 0, 0;
    0, 0, 0, 0;
    G1*s, G1*s, 0, 0;
    0, 0, G1*s, G1*s];
% Obtain TF
% Plot bode diagram of exponential approximations against e^(-s*tau4) in
% frequency range 0-40*pi
TFs = [1 0; 1 lbx+lvx; 1 -lbx-lvx] * [1 0 0 0; 0 1 0 0] * (F\G) * [1; 1/(1+(s*2*lbx/V)); 1/(1+(s*2*lvx/V)); 1/(1+(s*(2*(lbx+lvx))/V))] * 1/(1 + (s/40*pi)) * pi*sqrt(2*Av*V); %*[d1; d2; d3; d4];
%%
% Split outputs
TF1 = TFs(1);
TF2 = TFs(2);
% Find numerator and denominator
[n1, d1] = numden(TF1);
[n2, d2] = numden(TF2);
% Find LCM of denominators
a_lcm = lcm(d1,d2);
% quorem divides first argument by second argument
% Obtain quotient - in this case, it's 1
[Q1,R1] = quorem(a_lcm, d1);
[Q2,R2] = quorem(a_lcm, d2);
% Multiply numerators by quotient
n1hat = n1*Q1;
n2hat = n2*Q2;
% Find coeffs of each polynomial
% n1hat_coeffs = coeffs(n1hat, s, 'All');
% length(n1hat_coeffs); % 18, so leading power is 17
% n2hat_coeffs = coeffs(n2hat, s, 'All');
% length(n2hat_coeffs); % 18, so leading power is 17

% How do we evaluate a polynomial at -s?
% a_lcm_coeffs = coeffs(a_lcm, s, 'All'); % coeffs in descending powers of s
% length(a_lcm_coeffs); % 21, so leading power of s is 20

%% Look at first TF
n1hat_coeffs = sym2poly(n1hat);
n2hat_coeffs = sym2poly(n2hat);
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

%% Create z(s) from equation in the paper
if mod(n,2) == 1
    % ODD
    % Define co and ce as the odd and even parts of the numerator polynomial,
    % respectively
    n1co = n1hat_coeffs(2:2:end);
    n1ce = n1hat_coeffs(1:2:end);
    n2co = n2hat_coeffs(2:2:end);
    n2ce = n2hat_coeffs(1:2:end);

    % Square these polynomials, and add zeros to start and end of co2
    n1co2 = conv(n1co, n1co);
    n1ce2 = conv(n1ce, n1ce);
    n1co2 = [0 n1co2 0];

    n2co2 = conv(n2co, n2co);
    n2ce2 = conv(n2ce, n2ce);
    n2co2 = [0 n2co2 0];

    % Preallocate memory for b, which will be populated with the coefficients
    % of z_0
    b = zeros(1,length(cn));

    % Populates b with values of z_0
    for i = 1:length(cn)
        b(i) = n1ce2(i) + n2ce2(i) - n1co2(i) - n2co2(i);
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

    % Square these polynomials, and add zeros to start and end of co2
    n1co2 = conv(n1co, n1co);
    n1ce2 = conv(n1ce, n1ce);
    n2co2 = conv(n2co, n2co);
    n2ce2 = conv(n2ce, n2ce);
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

    n1ns = size(n1ce2,2);
    for i = n1ns+1:n
        n1ce2 = [0,n1ce2];
    end

    n2ns = size(n2ce2,2);
    for i = n2ns+1:n
        n2ce2 = [0,n2ce2];
    end

    % Preallocate memory for b, which will be populated with the coefficients
    % of z_0
    b = zeros(1,length(n1hat_coeffs));

    % Populates b with values of z_0
    for i = 1:length(n1hat_coeffs)
        b(i) = n1ce2(i) + n2ce2(i) - n1co2s(i) - n2co2s(i);
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

%% MIMO Test
% https://uk.mathworks.com/help/control/ug/mimo-transfer-function-models.html
% https://uk.mathworks.com/help/ident/ref/dynamicsystem.norm.html
% https://uk.mathworks.com/help/control/ref/tf.html
%openExample('control/ConcatenateSISOTFIntoMIMOTransferFunctionModelExample')
%In this example, you create a MIMO transfer function model by concatenating SISO transfer function models. Consider the following single-input, two-output transfer function:
% sys1 = tf([1 -1],[1 1]);		
% sys2 = tf([1 2],[1 4 5]);
% sys = [sys1;sys2]
% norm(sys, 2)

sys1 = tf(n1hat_coeffs, a_lcm_coeffs);
sys2 = tf(n2hat_coeffs, a_lcm_coeffs);
sys = [sys1; sys2];
norm(sys, 2)

%% Old stuff!
% Multiply odd powers of s by -1 to simulate evaluating the polynomial at
% -s
% for i = 1:1:length(a_lcm_coeffs)
%     if i == 1
%         a_lcm_neg = a_lcm_coeffs(i);
%     else
%         if mod(i,2) == 0% && i < 20
%             a_lcm_neg = cat(3, a_lcm_neg, a_lcm_coeffs(i)*-1);
%         else
%             a_lcm_neg = cat(3, a_lcm_neg, a_lcm_coeffs(i));
%         end
%     end
% end
% % Repeat for numerators
% for i = 1:1:length(n1hat_coeffs)
%     if i == 1
%         n1hat_neg = n1hat_coeffs(i)*-1;
%     else
%         if mod(i,2) == 0
%             n1hat_neg = cat(3, n1hat_neg, n1hat_coeffs(i));
%         else
%             n1hat_neg = cat(3, n1hat_neg, n1hat_coeffs(i)*-1);
%         end
%     end
% end
% 
% for i = 1:1:length(n2hat_coeffs)
%     if i == 1
%         n2hat_neg = n2hat_coeffs(i)*-1;
%     else
%         if mod(i,2) == 0
%             n2hat_neg = cat(3, n2hat_neg, n2hat_coeffs(i));
%         else
%             n2hat_neg = cat(3, n2hat_neg, n2hat_coeffs(i)*-1);
%         end
%     end
% end

%%
% Reshape into vectors so that they can be multiplied together
% n1hat_neg = reshape(n1hat_neg, 1, 18);
% n2hat_neg = reshape(n2hat_neg, 1, 18);
% a_lcm_neg = reshape(a_lcm_neg, 1, 21);
% Multiply together to obtain TF - check with Tim!
% finaln = (n1hat_neg.*n1hat_coeffs + n2hat_neg.*n2hat_coeffs);
% finald = (a_lcm_coeffs.*a_lcm_neg);

% subs_syms = [mv, mb, lvx, Kpz, Cpz, Ksz, Krz, Crz, Kaz, Ivz];
% subs_vals = [38000, 2500, 9.5, 4.935*10^6, 50.74*10^3, 1.016*10^6, 0.508*10^6, 64.11*10^3, 0, 2.31*10^6];
% mv = 38000;
% mb = 2500;
% Ivz = 2.31*10^6 (use value for Ivy)
% lvx = 9.5;
% Kpz = 4.935 * 10^6;
% Cpz = 50.74 * 10^3;
% Ksz = 1.016 * 10^6;
% Krz = 0.508 * 10^6;
% Crz = 64.11 * 10^3;
% Kaz = 0;

% cn1 = coeffs(subs(n1, subs_syms, subs_vals),s, 'All');
% cd1 = coeffs(subs(d1, subs_syms, subs_vals),s, 'All');
% cn1 = [zeros(1, length(cd1) - length(cn1) - 1) cn1];
% cn2 = coeffs(subs(n2, subs_syms, subs_vals),s, 'All');
% cd2 = coeffs(subs(d2, subs_syms, subs_vals),s, 'All');
% cn2 = [zeros(1, length(cd2) - length(cn2) - 1) cn2];


%[zhatv; thetahatv] = [1 0 0 0; 0 1 0 0] * F^(-1)*G*[d1; d2; d3; d4];
