%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc;                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define symbols for mechanical system
syms d1 d2 d3 d4 % Deltas -> these are our inputs
% syms zv theta zb1 zb2 % polynomial coefficients
syms s %k b c % mv Ivz lvx mb k bs c % G1 G2 % mechanical variables
% syms Kpz Cpz Ksz Krz Crz Kaz
% Define values for variables - symbolic version takes too long
% s = tf('s');
% Values taken from paper
k = 5*10^6; b = 92; c = 11.049;
mv = 38000; mb = 2500;
lvx = 9.5; Ivz = 2.31*10^6;
Kpz = 4.935*10^6; Ksz = 1.016*10^6; Krz = 0.508*10^6; Kaz = 0;
Cpz = 50.74*10^3; Crz = 64.11*10^3;
V = 55; Av = 2.5*10^(-7); lbx = 1.25;
% Suspension layouts
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
    -G2*s, (lvx*G2*s), 0, (mb*s^2 + 2*G1*s + 2*G2*s)];

K = [0, 0, 0, 0;
    0, 0, 0, 0;
    G1*s, G1*s, 0, 0;
    0, 0, G1*s, G1*s];
% Define time delays
tau2 = 2*lbx/V;
tau3 = 2*lvx/V;
tau4 = 2*(lbx+lvx)/V;
% Pade approximations
% (1680 - 840*(s*tau2) + 180*(s*tau2)^2 - 20*(s*tau2)^3 + (s*tau2)^4)/(1680 + 840*(s*tau2) + 180*(s*tau2)^2 + 20*(s*tau2)^3 + (s*tau2)^4)
% TD2 = (1680 - 840*s*tau2 + 180*(s*tau2)^2 - 20*(s*tau2)^3 + (s*tau2)^4)/(1680 + 840*s*tau2 + 180*(s*tau2)^2 + 20*(s*tau2)^3 + (s*tau2)^4);
% TD3 = (1680 - 840*s*tau3 + 180*(s*tau3)^2 - 20*(s*tau3)^3 + (s*tau3)^4)/(1680 + 840*s*tau3 + 180*(s*tau3)^2 + 20*(s*tau3)^3 + (s*tau3)^4);
% TD4 = (1680 - 840*s*tau4 + 180*(s*tau4)^2 - 20*(s*tau4)^3 + (s*tau4)^4)/(1680 + 840*s*tau4 + 180*(s*tau4)^2 + 20*(s*tau4)^3 + (s*tau4)^4);
TD2 = (30240 - 15120*s*tau2 + 3360*(s*tau2)^2 - 420*(s*tau2)^3 + 30*(s*tau2)^4 - (s*tau2)^5)/(30240 + 15120*s*tau2 + 3360*(s*tau2)^2 + 420*(s*tau2)^3 + 30*(s*tau2)^4 + (s*tau2)^5);
TD3 = (30240 - 15120*s*tau3 + 3360*(s*tau3)^2 - 420*(s*tau3)^3 + 30*(s*tau3)^4 - (s*tau3)^5)/(30240 + 15120*s*tau3 + 3360*(s*tau3)^2 + 420*(s*tau3)^3 + 30*(s*tau3)^4 + (s*tau3)^5);
TD4 = (30240 - 15120*s*tau4 + 3360*(s*tau4)^2 - 420*(s*tau4)^3 + 30*(s*tau4)^4 - (s*tau4)^5)/(30240 + 15120*s*tau4 + 3360*(s*tau4)^2 + 420*(s*tau4)^3 + 30*(s*tau4)^4 + (s*tau4)^5);
% TD2 = (120 - 60*s*tau2 + 12*(s*tau2)^2 - (s*tau2)^3)/(120 + 60*s*tau2 + 12*(s*tau2)^2 + (s*tau2)^3);
% TD3 = (120 - 60*s*tau3 + 12*(s*tau3)^2 - (s*tau3)^3)/(120 + 60*s*tau3 + 12*(s*tau3)^2 + (s*tau3)^3);
% TD4 = (120 - 60*s*tau4 + 12*(s*tau4)^2 - (s*tau4)^3)/(120 + 60*s*tau4 + 12*(s*tau4)^2 + (s*tau4)^3);
P = [1; TD2; TD3; TD4];
% Moments
Q = [1 0 0 0; 1 -(lbx+lvx) 0 0; 1 (lbx+lvx) 0 0];
% Power spectral density of delta_1
% delta1dot = sqrt(2*Av*V);
% Obtain TF
% TFs = [1 0; 1 lbx+lvx; 1 -lbx-lvx] * [1 0 0 0; 0 1 0 0] * (F\K) * P * 1/(1 + (s/40*pi)) * delta1dot; %*[d1; d2; d3; d4];
TFs = s * Q * (F\K) * P * 1/(1 + (s/40*pi));
% Final TF should be
% TFs = s * Q * (F\G) * P * delta1hat
%%
% Split outputs
TF1 = TFs(1);
TF2 = TFs(2);
TF3 = TFs(3);
% Find numerator and denominator
[n1, d1] = numden(TF1);
[n2, d2] = numden(TF2);
[n3, d3] = numden(TF3);
% Find LCM of denominators
a_lcm = lcm(lcm(d1,d2), d3);
% quorem divides first argument by second argument
% Obtain quotient
[Q1,R1] = quorem(a_lcm, d1);
[Q2,R2] = quorem(a_lcm, d2);
[Q3,R3] = quorem(a_lcm, d3);
% Multiply numerators by quotient
n1hat = n1*Q1;
n2hat = n2*Q2;
n3hat = n3*Q3;

%% Look at first TF
n1hat_coeffs = sym2poly(n1hat);
n2hat_coeffs = sym2poly(n2hat);
n3hat_coeffs = sym2poly(n3hat);
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

    % Preallocate memory for b, which will be populated with the coefficients
    % of z_0
    b = zeros(1,length(n1hat_coeffs));

    % Populates b with values of z_0
    for i = 1:length(n1hat_coeffs)
        b(i) = n1ce2(i) + n2ce2(i) + n3ce2(i) - n1co2(i) - n2co2(i) - n3co2(i);
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

    % Square these polynomials, and add zeros to start and end of co2
    n1co2 = conv(n1co, n1co);
    n1ce2 = conv(n1ce, n1ce);
    n2co2 = conv(n2co, n2co);
    n2ce2 = conv(n2ce, n2ce);
    n3co2 = conv(n3co, n3co);
    n3ce2 = conv(n3ce, n3ce);
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

    % Preallocate memory for b, which will be populated with the coefficients
    % of z_0
    b = zeros(1,length(n1hat_coeffs));

    % Populates b with values of z_0
    for i = 1:length(n1hat_coeffs)
        b(i) = n1ce2(i) + n2ce2(i) + n3ce2(i) - n1co2s(i) - n2co2s(i) - n3co2s(i);
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
sys3 = tf(n3hat_coeffs, a_lcm_coeffs);
sys = [sys1; sys2; sys3];
norm(sys, 2)

%%
% subs_syms = [mv, mb, lvx, Kpz, Cpz, Ksz, Krz, Crz, Kaz, Ivz];
% subs_vals = [38000, 2500, 9.5, 4.935*10^6, 50.74*10^3, 1.016*10^6, 0.508*10^6, 64.11*10^3, 0, 2.31*10^6];

% cn1 = coeffs(subs(n1, subs_syms, subs_vals),s, 'All');
% cd1 = coeffs(subs(d1, subs_syms, subs_vals),s, 'All');
% cn1 = [zeros(1, length(cd1) - length(cn1) - 1) cn1];
% cn2 = coeffs(subs(n2, subs_syms, subs_vals),s, 'All');
% cd2 = coeffs(subs(d2, subs_syms, subs_vals),s, 'All');
% cn2 = [zeros(1, length(cd2) - length(cn2) - 1) cn2];

%% Bode plots for Taylor series approximations
% s = tf('s');
% f = [1/(1+(s*tau2)); exp(-s*tau2)];
% [mag,phase,wout] = bode(f, {1E-1 1E3}); % expects frequency in radians so 40*pi/pi
% figure;
% subplot(2,3,1)
% loglog(squeeze(wout), squeeze(mag), 'LineWidth',2); ylabel("Magnitude"); title("$\tau_2$", "Interpreter","latex");
% legend("Taylor Series", "Exponential", 'Location','southwest');
% set(gca, 'YLim',[1E-3  40], 'FontSize', 14)
% grid
% subplot(2,3,4)
% semilogx(squeeze(wout), squeeze(phase), 'LineWidth',2); xlabel("Frequency"); ylabel("Phase"); legend("Taylor Series", "Exponential", 'Location','southwest');
% set(gca, 'FontSize', 14);
% 
% s = tf('s');
% f = [1/(1+(s*tau3)); exp(-s*tau3)];
% [mag,phase,wout] = bode(f, {1E-1 1E3}); % expects frequency in radians so 40*pi/pi
% subplot(2,3,2)
% loglog(squeeze(wout), squeeze(mag), 'LineWidth',2); ylabel("Magnitude"); title("$\tau_3$", "Interpreter","latex");
% set(gca, 'YLim',[1E-3  40], 'FontSize', 14); legend("Taylor Series", "Exponential", 'Location','southwest');
% grid
% subplot(2,3,5)
% semilogx(squeeze(wout), squeeze(phase), 'LineWidth',2); xlabel("Frequency"); ylabel("Phase"); legend("Taylor Series", "Exponential", 'Location','southwest');
% set(gca, 'FontSize', 14);
% 
% s = tf('s');
% f = [1/(1+(s*tau4)); exp(-s*tau4)];
% [mag,phase,wout] = bode(f, {1E-1 1E3}); % expects frequency in radians so 40*pi/pi
% subplot(2,3,3)
% loglog(squeeze(wout), squeeze(mag), 'LineWidth',2); ylabel("Magnitude"); title("$\tau_4$", "Interpreter","latex"); legend("Taylor Series", "Exponential", 'Location','southwest');
% set(gca, 'YLim',[1E-3  40], 'FontSize', 14)
% grid
% subplot(2,3,6)
% semilogx(squeeze(wout), squeeze(phase), 'LineWidth',2); xlabel("Frequency"); ylabel("Phase"); legend("Taylor Series", "Exponential", 'Location','southwest');
% set(gca, 'FontSize', 14);

%% Bode plot for frequency filter
% f = 1/(1 + (s/40*pi));
% [mag,phase,wout] = bode(f, {1E-1 1E3}); % expects frequency in radians so 40*pi/pi
% figure(2)
% subplot(2,1,1)
% loglog(squeeze(wout), squeeze(mag)); ylabel("Magnitude");
% set(gca, 'YLim',[1E-3  40])
% grid
% subplot(2,1,2)
% semilogx(squeeze(wout), squeeze(phase)); xlabel("Frequency"); ylabel("Phase")

%% Bode plot for tau2
% s = tf('s');
% f = [(120 - 60*s*tau2 + 12*(s*tau2)^2 - (s*tau2)^3)/(120 + 60*s*tau2 + 12*(s*tau2)^2 + (s*tau2)^3); exp(-s*tau2)];
% % f1 = tf([-1 12 -60 120], [120 60 12 1]);
% [mag,phase,wout] = bode(f, {1E-1 40*pi}); % expects frequency in radians so 40*pi/pi
% figure
% subplot(2,3,1)
% loglog(squeeze(wout), squeeze(mag), 'LineWidth',2); ylabel("Magnitude"); title("$\tau_2$", "Interpreter","latex"); legend("Padé Approximant", "Exponential", 'Location','southwest');
% set(gca, 'YLim',[1E-3  40]); set(gca, 'FontSize', 14);
% grid
% subplot(2,3,4)
% phase = squeeze(phase);
% phase(1,:) = phase(1,:) - 720;
% semilogx(squeeze(wout), squeeze(phase), 'LineWidth',2); xlabel("Frequency"); ylabel("Phase"); legend("Padé Approximant", "Exponential", 'Location','southwest');
% set(gca, 'FontSize', 14);
% 
% % Bode plot for tau3
% s = tf('s');
% f = [(120 - 60*s*tau3 + 12*(s*tau3)^2 - (s*tau3)^3)/(120 + 60*s*tau3 + 12*(s*tau3)^2 + (s*tau3)^3); exp(-s*tau3)];
% [mag,phase,wout] = bode(f, {1E-1 40*pi}); % expects frequency in radians so 40*pi/pi
% subplot(2,3,2)
% loglog(squeeze(wout), squeeze(mag), 'LineWidth',2); ylabel("Magnitude"); title("$\tau_3$", "Interpreter","latex"); legend("Padé Approximant", "Exponential", 'Location','southwest');
% set(gca, 'YLim',[1E-3  40]); set(gca, 'FontSize', 14);
% grid
% subplot(2,3,5)
% phase = squeeze(phase);
% phase(1,:) = phase(1,:) - 720;
% semilogx(squeeze(wout), squeeze(phase), 'LineWidth',2); xlabel("Frequency"); ylabel("Phase"); legend("Padé Approximant", "Exponential", 'Location','southwest');
% set(gca, 'FontSize', 14);
% 
% % Bode plot for tau4
% s = tf('s');
% f = [(120 - 60*s*tau4 + 12*(s*tau4)^2 - (s*tau4)^3)/(120 + 60*s*tau4 + 12*(s*tau4)^2 + (s*tau4)^3); exp(-s*tau4)];
% [mag,phase,wout] = bode(f, {1E-1 40*pi}); % expects frequency in radians so 40*pi/pi
% subplot(2,3,3)
% loglog(squeeze(wout), squeeze(mag), 'LineWidth',2); ylabel("Magnitude"); title("$\tau_4$", "Interpreter","latex"); legend("Padé Approximant", "Exponential", 'Location','southwest');
% set(gca, 'YLim',[1E-3  40]); set(gca, 'FontSize', 14);
% grid
% subplot(2,3,6)
% phase = squeeze(phase);
% phase(1,:) = phase(1,:) - 720;
% semilogx(squeeze(wout), squeeze(phase), 'LineWidth',2); xlabel("Frequency"); ylabel("Phase"); legend("Padé Approximant", "Exponential", 'Location','southwest');
% set(gca, 'FontSize', 14);

%% Taylor Series Approximation
% %(1680 - 840*x(i)*T + 180*x(i)*T^2 - 20*x(i)*T^3 + x(i)*T^4)/(1680 + 840*x(i)*T + 180*x(i)*T^2 + 20*x(i)*T^3 + x(i)*T^4)
% %(30240 - 15120*x(i)*T + 3360*(x(i)*T)^2 - 420*(x(i)*T)^3 + 30*(x(i)*T)^4 - (x(i)*T)^5)/(30240 + 15120*x(i)*T + 3360*(x(i)*T)^2 + 420*(x(i)*T)^3 + 30*(x(i)*T)^4 + (x(i)*T)^5)
% x = 0:0.1:40*pi;
% y = [];
% y2 = [];
% for i = 1:length(x)
%     y = [y 1/(1+(x(i)*tau2))];
%     y2 = [y2 exp(-x(i)*tau2)];
% end
% figure;
% subplot(3,1,1)
% plot(x,y, 'LineWidth',2)
% hold on;
% plot(x,y2, 'LineWidth',2)
% legend("Taylor Series", "Exponential")
% ylabel("f(x)")
% title("$\tau_2$", "Interpreter", "latex")
% set(gca, 'FontSize', 14)
% 
% % Taylor Series Approximation - tau3
% x = 0:0.1:40*pi;
% y = [];
% y2 = [];
% for i = 1:length(x)
%     y = [y 1/(1+(x(i)*tau3))];
%     y2 = [y2 exp(-x(i)*tau3)];
% end
% subplot(3,1,2)
% plot(x,y, 'LineWidth',2)
% hold on;
% plot(x,y2, 'LineWidth',2)
% legend("Taylor Series", "Exponential")
% ylabel("f(x)")
% title("$\tau_3$", "Interpreter", "latex")
% set(gca, 'FontSize', 14)
% 
% % Taylor Series Approximation - tau4
% x = 0:0.1:40*pi;
% y = [];
% y2 = [];
% for i = 1:length(x)
%     y = [y 1/(1+(x(i)*tau4))];
%     y2 = [y2 exp(-x(i)*tau4)];
% end
% subplot(3,1,3)
% plot(x,y, 'LineWidth',2)
% hold on;
% plot(x,y2, 'LineWidth',2)
% legend("Taylor Series", "Exponential")
% title("$\tau_4$", "Interpreter", "latex")
% xlabel("x"); ylabel("f(x)"); set(gca, 'FontSize', 14);


%% 3rd order Pade Comparison - tau2
% %(1680 - 840*x(i)*T + 180*x(i)*T^2 - 20*x(i)*T^3 + x(i)*T^4)/(1680 + 840*x(i)*T + 180*x(i)*T^2 + 20*x(i)*T^3 + x(i)*T^4)
% %(30240 - 15120*x(i)*T + 3360*(x(i)*T)^2 - 420*(x(i)*T)^3 + 30*(x(i)*T)^4 - (x(i)*T)^5)/(30240 + 15120*x(i)*T + 3360*(x(i)*T)^2 + 420*(x(i)*T)^3 + 30*(x(i)*T)^4 + (x(i)*T)^5)
% x = 0:0.1:40*pi;
% y = [];
% y2 = [];
% y3 = [];
% for i = 1:length(x)
%     y = [y 1/(1+(x(i)*tau2))];
%     y2 = [y2 exp(-x(i)*tau2)];
%     y3 = [y3 (120 - 60*x(i)*tau2 + 12*(x(i)*tau2)^2 - (x(i)*tau2)^3)/(120 + 60*x(i)*tau2 + 12*(x(i)*tau2)^2 + (x(i)*tau2)^3);];
% end
% figure;
% subplot(3,1,1)
% plot(x,y, 'LineWidth',2)
% hold on;
% plot(x,y2, 'LineWidth',2)
% plot(x, y3, 'LineWidth',2)
% legend("Taylor Series", "Exponential", "3rd Order Padé Approximation")
% ylabel("f(x)")
% title("$\tau_2$", "Interpreter", "latex")
% set(gca, 'FontSize', 14)
% 
% % 3rd order Pade Comparison - tau3
% x = 0:0.1:40*pi;
% y = [];
% y2 = [];
% y3 = [];
% for i = 1:length(x)
%     y = [y 1/(1+(x(i)*tau3))];
%     y2 = [y2 exp(-x(i)*tau3)];
%     y3 = [y3 (120 - 60*x(i)*tau3 + 12*(x(i)*tau3)^2 - (x(i)*tau3)^3)/(120 + 60*x(i)*tau3 + 12*(x(i)*tau3)^2 + (x(i)*tau3)^3);];
% end
% subplot(3,1,2)
% plot(x,y, 'LineWidth',2)
% hold on;
% plot(x,y2, 'LineWidth',2)
% plot(x, y3, 'LineWidth',2)
% legend("Taylor Series", "Exponential", "3rd Order Padé Approximation")
% ylabel("f(x)")
% title("$\tau_3$", "Interpreter", "latex")
% set(gca, 'FontSize', 14)
% 
% % 3rd order Pade Comparison - tau4
% x = 0:0.1:40*pi;
% y = [];
% y2 = [];
% y3 = [];
% for i = 1:length(x)
%     y = [y 1/(1+(x(i)*tau4))];
%     y2 = [y2 exp(-x(i)*tau4)];
%     y3 = [y3 (120 - 60*x(i)*tau4 + 12*(x(i)*tau4)^2 - (x(i)*tau4)^3)/(120 + 60*x(i)*tau4 + 12*(x(i)*tau4)^2 + (x(i)*tau4)^3);];
% end
% subplot(3,1,3)
% plot(x,y, 'LineWidth',2)
% hold on;
% plot(x,y2, 'LineWidth',2)
% plot(x, y3, 'LineWidth',2)
% legend("Taylor Series", "Exponential", "3rd Order Padé Approximation")
% title("$\tau_4$", "Interpreter", "latex")
% xlabel("x"); ylabel("f(x)"); set(gca, 'FontSize', 14);

%%
% TF1 = TFs(1);
% TF2 = TFs(2);
% [n1, d1] = numden(TF1);
% [n2, d2] = numden(TF2);
% 
% TF1 = tf(sym2poly(n1), sym2poly(d1));
% TF2 = tf(sym2poly(n2), sym2poly(d2));
%% Bode plot for transfer functions - expand this to higher frequencies
% f = [TF1; TF2];
% [mag,phase,wout] = bode(f, {1E-1 40*pi}); % expects frequency in radians so 40*pi/pi
% figure(2)
% subplot(2,1,1)
% loglog(squeeze(wout), squeeze(mag)); ylabel("Magnitude");
% set(gca, 'YLim',[1E-3  40])
% grid
% subplot(2,1,2)
% phase = squeeze(phase);
% phase(1,:) = phase(1,:) - 720;
% semilogx(squeeze(wout), squeeze(phase)); xlabel("Frequency"); ylabel("Phase")