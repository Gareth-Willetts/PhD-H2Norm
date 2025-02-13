%% Create TF with exponential P matrix
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
% TD2 = (120 - 60*s*tau2 + 12*(s*tau2)^2 - (s*tau2)^3)/(120 + 60*s*tau2 + 12*(s*tau2)^2 + (s*tau2)^3);
% TD3 = (120 - 60*s*tau3 + 12*(s*tau3)^2 - (s*tau3)^3)/(120 + 60*s*tau3 + 12*(s*tau3)^2 + (s*tau3)^3);
% TD4 = (120 - 60*s*tau4 + 12*(s*tau4)^2 - (s*tau4)^3)/(120 + 60*s*tau4 + 12*(s*tau4)^2 + (s*tau4)^3);
% TD2 = (1680 - 840*(s*tau2) + 180*(s*tau2)^2 - 20*(s*tau2)^3 + (s*tau2)^4)/(1680 + 840*(s*tau2) + 180*(s*tau2)^2 + 20*(s*tau2)^3 + (s*tau2)^4);
% TD3 = (1680 - 840*(s*tau3) + 180*(s*tau3)^2 - 20*(s*tau3)^3 + (s*tau3)^4)/(1680 + 840*(s*tau3) + 180*(s*tau3)^2 + 20*(s*tau3)^3 + (s*tau3)^4);
% TD4 = (1680 - 840*(s*tau4) + 180*(s*tau4)^2 - 20*(s*tau4)^3 + (s*tau4)^4)/(1680 + 840*(s*tau4) + 180*(s*tau4)^2 + 20*(s*tau4)^3 + (s*tau4)^4);
TD2 = (30240 - 15120*s*tau2 + 3360*(s*tau2)^2 - 420*(s*tau2)^3 + 30*(s*tau2)^4 - (s*tau2)^5)/(30240 + 15120*s*tau2 + 3360*(s*tau2)^2 + 420*(s*tau2)^3 + 30*(s*tau2)^4 + (s*tau2)^5);
TD3 = (30240 - 15120*s*tau3 + 3360*(s*tau3)^2 - 420*(s*tau3)^3 + 30*(s*tau3)^4 - (s*tau3)^5)/(30240 + 15120*s*tau3 + 3360*(s*tau3)^2 + 420*(s*tau3)^3 + 30*(s*tau3)^4 + (s*tau3)^5);
TD4 = (30240 - 15120*s*tau4 + 3360*(s*tau4)^2 - 420*(s*tau4)^3 + 30*(s*tau4)^4 - (s*tau4)^5)/(30240 + 15120*s*tau4 + 3360*(s*tau4)^2 + 420*(s*tau4)^3 + 30*(s*tau4)^4 + (s*tau4)^5);
P = [1; TD2; TD3; TD4];
% Moments
Q = [1 0 0 0; 1 -(lbx+lvx) 0 0; 1 (lbx+lvx) 0 0];
% Power spectral density of delta_1
delta1dot = sqrt(2*Av*V);
% Obtain TF
TFs = s * Q * (F\K) * P * delta1dot * 1/(1 + (s/40*pi));
% TFs = [1 0; 1 lbx+lvx; 1 -lbx-lvx] * [1 0 0 0; 0 1 0 0] * (F\K) * P * 1/(1 + (s/40*pi)) * delta1dot; %*[d1; d2; d3; d4];

TF1 = TFs(1); TF2 = TFs(2); TF3 = TFs(3);
[n1, d1] = numden(TF1);
[n2, d2] = numden(TF2);
[n3, d3] = numden(TF3);

TF1 = tf(sym2poly(n1), sym2poly(d1));
TF2 = tf(sym2poly(n2), sym2poly(d2));
TF3 = tf(sym2poly(n3), sym2poly(d3));
%%
% Bode plot for transfer functions - expand this to higher frequencies
f = [TF1; TF2; TF3];
[mag,phase,wout] = bode(f, {15 1E3}); % expects frequency in radians so 40*pi/pi
figure
subplot(2,1,1)
loglog(squeeze(wout), squeeze(mag).^2, 'LineWidth',2); ylabel("Magnitude"); set(gca, 'FontSize', 14);
legend("Middle of train", "Leading edge of train", "Trailing edge of train", 'Location', 'southwest')
%set(gca, 'YLim',[1E-3  40])
grid
subplot(2,1,2)
phase = squeeze(phase);
phase(1,:) = phase(1,:) - 720;
semilogx(squeeze(wout), squeeze(phase), 'LineWidth',2); xlabel("Frequency"); ylabel("Phase"); grid; set(gca, 'FontSize', 14);
%legend("$J_{1Z(M)}^2$", "$J_{1Z(L)}^2$", "Interpreter", "latex");
legend("Middle of train", "Leading edge of train", "Trailing edge of train", 'Location', 'southwest')

% Trapz
test = squeeze(mag).^2;
trapz(squeeze(wout), test(3,:))/(2*pi)

%%
d1coeffs = sym2poly(d1);
n1coeffs = sym2poly(n1);

d2coeffs = sym2poly(d2);
n2coeffs = sym2poly(n2);

d3coeffs = sym2poly(d3);
n3coeffs = sym2poly(n3);
max([abs(n3coeffs(1)/d3coeffs(1)), abs(n2coeffs(1)/d2coeffs(1)), abs(n1coeffs(1)/d1coeffs(1))])
d1roots = roots(d1coeffs);
%%
s = tf('s');
f = [TF1; TF2; TF3; 10^(7/2)/s^3];
[mag,phase,wout] = bode(f, {1E3 1E9}); % expects frequency in radians so 40*pi/pi
figure
%subplot(2,1,1)
loglog(squeeze(wout), squeeze(mag).^2, 'LineWidth', 2); ylabel("Magnitude");
set(gca, 'FontSize', 14);
%set(gca, 'YLim',[1E-3  40])
% grid
% subplot(2,1,2)
% phase = squeeze(phase);
% phase(1,:) = phase(1,:) - 720;
% semilogx(squeeze(wout), squeeze(phase)); xlabel("Frequency"); ylabel("Phase"); grid
legend("Middle of train", "Leading edge of train", "Trailing edge of train", "Upper Bound", 'Location', 'southwest')

%% Bode plots for exponential comparison - 4TH ORDER
s = tf('s');
f = [(1680 - 840*(s*tau2) + 180*(s*tau2)^2 - 20*(s*tau2)^3 + (s*tau2)^4)/(1680 + 840*(s*tau2) + 180*(s*tau2)^2 + 20*(s*tau2)^3 + (s*tau2)^4); exp(-s*tau2)];
% f1 = tf([-1 12 -60 120], [120 60 12 1]);
[mag,phase,wout] = bode(f, {1E-1 40*pi}); % expects frequency in radians so 40*pi/pi %1E3
figure
subplot(2,3,1)
loglog(squeeze(wout), squeeze(mag), 'LineWidth',2); ylabel("Magnitude"); title("$\tau_2$", "Interpreter","latex"); legend("Padé Approximant", "Exponential", 'Location','southwest');
set(gca, 'YLim',[1E-3  40]); set(gca, 'FontSize', 14);
grid
subplot(2,3,4)
phase = squeeze(phase);
phase(1,:) = phase(1,:) - 720;
semilogx(squeeze(wout), squeeze(phase), 'LineWidth',2); xlabel("Frequency"); ylabel("Phase"); legend("Padé Approximant", "Exponential", 'Location','southwest');
grid
set(gca, 'FontSize', 14);

% Bode plot for tau3
s = tf('s');
f = [(1680 - 840*(s*tau3) + 180*(s*tau3)^2 - 20*(s*tau3)^3 + (s*tau3)^4)/(1680 + 840*(s*tau3) + 180*(s*tau3)^2 + 20*(s*tau3)^3 + (s*tau3)^4); exp(-s*tau3)];
[mag,phase,wout] = bode(f, {1E-1 40*pi}); % expects frequency in radians so 40*pi/pi
subplot(2,3,2)
loglog(squeeze(wout), squeeze(mag), 'LineWidth',2); ylabel("Magnitude"); title("$\tau_3$", "Interpreter","latex"); legend("Padé Approximant", "Exponential", 'Location','southwest');
set(gca, 'YLim',[1E-3  40]); set(gca, 'FontSize', 14);
grid
subplot(2,3,5)
phase = squeeze(phase);
phase(1,:) = phase(1,:) - 720;
semilogx(squeeze(wout), squeeze(phase), 'LineWidth',2); xlabel("Frequency"); ylabel("Phase"); legend("Padé Approximant", "Exponential", 'Location','southwest');
grid
set(gca, 'FontSize', 14);

% Bode plot for tau4
s = tf('s');
f = [(1680 - 840*(s*tau4) + 180*(s*tau4)^2 - 20*(s*tau4)^3 + (s*tau4)^4)/(1680 + 840*(s*tau4) + 180*(s*tau4)^2 + 20*(s*tau4)^3 + (s*tau4)^4); exp(-s*tau4)];
[mag,phase,wout] = bode(f, {1E-1 40*pi}); % expects frequency in radians so 40*pi/pi
subplot(2,3,3)
loglog(squeeze(wout), squeeze(mag), 'LineWidth',2); ylabel("Magnitude"); title("$\tau_4$", "Interpreter","latex"); legend("Padé Approximant", "Exponential", 'Location','southwest');
set(gca, 'YLim',[1E-3  40]); set(gca, 'FontSize', 14);
grid
subplot(2,3,6)
phase = squeeze(phase);
phase(1,:) = phase(1,:) - 720;
semilogx(squeeze(wout), squeeze(phase), 'LineWidth',2); xlabel("Frequency"); ylabel("Phase"); legend("Padé Approximant", "Exponential", 'Location','southwest');
grid
set(gca, 'FontSize', 14);

%% Bode plots - 5th order Pade
s = tf('s');
f = [(30240 - 15120*s*tau2 + 3360*(s*tau2)^2 - 420*(s*tau2)^3 + 30*(s*tau2)^4 - (s*tau2)^5)/(30240 + 15120*s*tau2 + 3360*(s*tau2)^2 + 420*(s*tau2)^3 + 30*(s*tau2)^4 + (s*tau2)^5); exp(-s*tau2)];
% f1 = tf([-1 12 -60 120], [120 60 12 1]);
[mag,phase,wout] = bode(f, {1E-1 20}); % expects frequency in radians so 40*pi/pi %1E3
figure
subplot(2,3,1)
loglog(squeeze(wout), squeeze(mag), 'LineWidth',2); ylabel("Magnitude"); title("$\tau_2$", "Interpreter","latex"); legend("Padé Approximant", "Exponential", 'Location','southwest');
set(gca, 'YLim',[1E-3  40]); set(gca, 'FontSize', 14);
grid
subplot(2,3,4)
phase = squeeze(phase);
phase(1,:) = phase(1,:) - 1080;
semilogx(squeeze(wout), squeeze(phase), 'LineWidth',2); xlabel("Frequency"); ylabel("Phase"); legend("Padé Approximant", "Exponential", 'Location','southwest');
grid
set(gca, 'FontSize', 14);

% Bode plot for tau3
s = tf('s');
f = [(30240 - 15120*s*tau3 + 3360*(s*tau3)^2 - 420*(s*tau3)^3 + 30*(s*tau3)^4 - (s*tau3)^5)/(30240 + 15120*s*tau3 + 3360*(s*tau3)^2 + 420*(s*tau3)^3 + 30*(s*tau3)^4 + (s*tau3)^5); exp(-s*tau3)];
[mag,phase,wout] = bode(f, {1E-1 20}); % expects frequency in radians so 40*pi/pi
subplot(2,3,2)
loglog(squeeze(wout), squeeze(mag), 'LineWidth',2); ylabel("Magnitude"); title("$\tau_3$", "Interpreter","latex"); legend("Padé Approximant", "Exponential", 'Location','southwest');
set(gca, 'YLim',[1E-3  40]); set(gca, 'FontSize', 14);
grid
subplot(2,3,5)
phase = squeeze(phase);
phase(1,:) = phase(1,:) - 1080;
semilogx(squeeze(wout), squeeze(phase), 'LineWidth',2); xlabel("Frequency"); ylabel("Phase"); legend("Padé Approximant", "Exponential", 'Location','southwest');
grid
set(gca, 'FontSize', 14);

% Bode plot for tau4
s = tf('s');
f = [(30240 - 15120*s*tau4 + 3360*(s*tau4)^2 - 420*(s*tau4)^3 + 30*(s*tau4)^4 - (s*tau4)^5)/(30240 + 15120*s*tau4 + 3360*(s*tau4)^2 + 420*(s*tau4)^3 + 30*(s*tau4)^4 + (s*tau4)^5); exp(-s*tau4)];
[mag,phase,wout] = bode(f, {1E-1 20}); % expects frequency in radians so 40*pi/pi
subplot(2,3,3)
loglog(squeeze(wout), squeeze(mag), 'LineWidth',2); ylabel("Magnitude"); title("$\tau_4$", "Interpreter","latex"); legend("Padé Approximant", "Exponential", 'Location','southwest');
set(gca, 'YLim',[1E-3  40]); set(gca, 'FontSize', 14);
grid
subplot(2,3,6)
phase = squeeze(phase);
phase(1,:) = phase(1,:) - 1080;
semilogx(squeeze(wout), squeeze(phase), 'LineWidth',2); xlabel("Frequency"); ylabel("Phase"); legend("Padé Approximant", "Exponential", 'Location','southwest');
grid
set(gca, 'FontSize', 14);

%% Bode plots - 5th order Pade - phase comparison
s = tf('s');
f = [(30240 - 15120*s*tau2 + 3360*(s*tau2)^2 - 420*(s*tau2)^3 + 30*(s*tau2)^4 - (s*tau2)^5)/(30240 + 15120*s*tau2 + 3360*(s*tau2)^2 + 420*(s*tau2)^3 + 30*(s*tau2)^4 + (s*tau2)^5); exp(-s*tau2)];
% f1 = tf([-1 12 -60 120], [120 60 12 1]);
[~,phase,wout] = bode(f, {1E-1 20}); % expects frequency in radians so 40*pi/pi %1E3
figure
subplot(3,1,1)
phase = squeeze(phase);
phase(1,:) = phase(1,:) - 1080;
semilogx(squeeze(wout), squeeze(phase), 'LineWidth',2); xlabel("Frequency"); ylabel("Phase"); legend("Padé Approximant", "Exponential", 'Location','southwest'); title("$\tau_2$", "Interpreter","latex");
grid
set(gca, 'FontSize', 14);


% Bode plot for tau3
s = tf('s');
f = [(30240 - 15120*s*tau3 + 3360*(s*tau3)^2 - 420*(s*tau3)^3 + 30*(s*tau3)^4 - (s*tau3)^5)/(30240 + 15120*s*tau3 + 3360*(s*tau3)^2 + 420*(s*tau3)^3 + 30*(s*tau3)^4 + (s*tau3)^5); exp(-s*tau3)];
[mag,phase,wout] = bode(f, {1E-1 20}); % expects frequency in radians so 40*pi/pi
subplot(3,1,2)
phase = squeeze(phase);
phase(1,:) = phase(1,:) - 1080;
semilogx(squeeze(wout), squeeze(phase), 'LineWidth',2); xlabel("Frequency"); ylabel("Phase"); legend("Padé Approximant", "Exponential", 'Location','southwest'); title("$\tau_3$", "Interpreter","latex");
grid
set(gca, 'FontSize', 14);

% Bode plot for tau4
s = tf('s');
f = [(30240 - 15120*s*tau4 + 3360*(s*tau4)^2 - 420*(s*tau4)^3 + 30*(s*tau4)^4 - (s*tau4)^5)/(30240 + 15120*s*tau4 + 3360*(s*tau4)^2 + 420*(s*tau4)^3 + 30*(s*tau4)^4 + (s*tau4)^5); exp(-s*tau4)];
[mag,phase,wout] = bode(f, {1E-1 20}); % expects frequency in radians so 40*pi/pi
subplot(3,1,3)
phase = squeeze(phase);
phase(1,:) = phase(1,:) - 1080;
semilogx(squeeze(wout), squeeze(phase), 'LineWidth',2); xlabel("Frequency"); ylabel("Phase"); legend("Padé Approximant", "Exponential", 'Location','southwest'); title("$\tau_4$", "Interpreter","latex");
grid
set(gca, 'FontSize', 14);

%% Low Pass Filter Disturbances
% https://ascelibrary.org/doi/epdf/10.1061/9780784482902.073
s = tf('s');
f = [1/(1+(s/40*pi))];
[mag,phase,wout] = bode(f, {55 5500}); % expects frequency in radians so 40*pi/pi
figure;
subplot(2,1,1)
loglog(squeeze(wout), squeeze(mag).^2, 'LineWidth',2); ylabel("Magnitude"); grid on; set(gca, 'FontSize', 14);
%phase = squeeze(phase);
%phase(1,:) = phase(1,:) - 720;
subplot(2,1,2)
semilogx(squeeze(wout), squeeze(phase), 'LineWidth',2); xlabel("Frequency"); ylabel("Phase"); grid; set(gca, 'FontSize', 14);
%legend("$J_{1Z(M)}^2$", "$J_{1Z(L)}^2$", "Interpreter", "latex");
%legend("Middle of train", "Leading edge of train", "Trailing edge of train", 'Location', 'southwest')
