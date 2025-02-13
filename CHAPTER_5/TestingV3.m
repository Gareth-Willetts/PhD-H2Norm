clear all; close all; clc
%% Following works well for obtaining symbolic TF
% Define symbols for components in the suspension system
time = [];
for i = 1:10
    tic
    syms ms mb mw ks kb kw cw z s cb cs x y t;
    syms x1  x2 %real
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
    [n, d] = numden(TFTF);
    time(i) = toc;
end
mean(time(10:end))


%% Evaluate symbolic TF so they're functions of x and y
symbolic = 0;
if symbolic == 1
    ms = 3500; mb = 250; mw = 350; ks = 141E3; kb = 1260E3; cb = 7100; kw = 8E9; cw = 670E3; %z = 200; cs = 8870;
    cd = coeffs(subs(d),s, 'All');
    cn = coeffs(subs(n),s,'All');
    cn = [zeros(1, length(cd) - length(cn) - 1) cn];
else
    ms = 3500; mb = 250; mw = 350; ks = 141E3; kb = 1260E3; kw = 8E9; cw = 670E3; %z = 200; cs = 8870;cb = 7100;
    x1 = 7616;
    x2 = 11844;
    cd = sym2poly(subs(d));
    cn = sym2poly(subs(n));
    cn = [zeros(1, length(cd) - length(cn) - 1) cn];
end


time = [];
%%
for j = 1:1000
    tic
    [n, pi, pi1] = H2norm_p_creation(cd);
    z0 = H2norm_z0_creation(cn,n,symbolic);
    [psi, epsilon, mu, gamma] = H2norm_intermediate_variables(pi, pi1, z0);
    [pi, pi1, z, mu, epsilon, psi, gamma] = H2norm_first_iteration(pi, pi1, z0, mu, epsilon, psi, gamma, n);
    [pi, pi1, z, mu, epsilon, psi, gamma] = H2norm_iterations(pi, pi1, z, mu, epsilon, psi, gamma, n);
    H2 = z/(2*cd(1)*pi1);
    time(j) = toc;
end

mean(time(10:end))

%%
time_norm = [];
for j = 1:1000
    tic
    norm(tf(cn,cd),2);
    time_norm(j) = toc;
end

mean(time_norm(10:end))

%%
semilogy(time(10:end));
hold on
semilogy(time_norm(10:end));
title("Semilog plot of computation times for 1000 H_2 norm calculations")
ylabel("log(Time)")
xlabel("Number of executions")
%%
function[n, pi, pi1] = H2norm_p_creation(cd)
    % Degree of transfer function = length of denominator - 1 (due to the
    % constant term)
    n = size(cd,2)-1;
    
    if mod(n,2) == 1
        % ODD
        pi1 = cd(2:2:end);
        pi = [cd(1:2:end) 0];
    else
        % EVEN
        pi = cd(1:2:end);
        pi1 = [cd(2:2:end) 0];
    end
end

function[z0] = H2norm_z0_creation(cn,n,symbolic)
    if symbolic
        syms t
        if mod(n,2) == 1
            % ODD
            % Define co and ce as the odd and even parts of the numerator polynomial,
            % respectively
            co = cn(2:2:end);
            ce = cn(1:2:end);
        
            % Square these polynomials, and add zeros to start and end of co2
            co2 = poly2sym(co,t)*poly2sym(co,t)*t;
            ce2 = poly2sym(ce,t)*poly2sym(ce,t);
            co2 = fliplr(coeffs(co2, t, 'All'));
            ce2 = fliplr(coeffs(ce2, t, 'All'));
            co2 = [0 co2 0];
        
            % Preallocate memory for z, change name of b to z0 (unnecessary)
            z0 = ce2 - co2;
        else
            % EVEN
            % Define co and ce as the odd and even parts of the numerator polynomial,
            % respectively
            % ce and co are the other way around in Tim's example?
            co = cn(1:2:end);
            ce = cn(2:2:end);
        
            % Square these polynomials, and add zeros to start and end of co2
            co2 = poly2sym(co,t)*poly2sym(co,t);
            ce2 = poly2sym(ce,t)*poly2sym(ce,t);
            co2 = coeffs(co2, t, 'All');
            ce2 = coeffs(ce2, t,'All');
        
            % CHECK NUMBER OF 0s
            %co2 = [0 co2 0];
        
            co2s = [co2, 0];
            %ce2 = [ce2, 0, 0];
            %co2s = co2;
            ns = size(co2s,2);
            for i = ns+1:n
                co2s = [0,co2s];
            end
            ns = size(ce2,2);
            for i = ns+1:n
                ce2 = [0,ce2];
            end
        
            % Preallocate memory for z0
            z0 = ce2 - co2s;
        end
    else
        if mod(n,2) == 1
            % ODD
            % Define co and ce as the odd and even parts of the numerator polynomial,
            % respectively
            co = cn(2:2:end);
            ce = cn(1:2:end);
        
            % Square these polynomials, and add zeros to start and end of co2
            co2 = conv(co, co);
            ce2 = conv(ce, ce);
            co2 = [0 co2 0];
        
            % Preallocate memory for b, which will be populated with the coefficients
            % of z_0
            b = zeros(1,length(cn));
        
            % Populates b with values of z_0
            for i = 1:length(cn)
                b(i) = ce2(i) - co2(i);
            end
        
            % Preallocate memory for z, change name of b to z0 (unnecessary)
            z0 = b;
        else
            % EVEN
            % Define co and ce as the odd and even parts of the numerator polynomial,
            % respectively
            ce = cn(2:2:end);
            co = cn(1:2:end);
        
            % Square these polynomials, and add zeros to start and end of co2
            co2 = conv(co, co);
            ce2 = conv(ce, ce);
            %co2 = [0 co2 0];
        
            co2s = [co2, 0];
            ns = size(co2s,2);
            for i = ns+1:n
                co2s = [0,co2s];
            end
            ns = size(ce2,2);
            for i = ns+1:n
                ce2 = [0,ce2];
            end
        
            % Preallocate memory for b, which will be populated with the coefficients
            % of z_0
            b = zeros(1,length(cn));
        
            % Populates b with values of z_0
            for i = 1:length(cn)
                b(i) = ce2(i) - co2s(i);
            end
        
            % Preallocate memory for z, change name of b to z0 (unnecessary)
            z0 = b;
        end
    end
end

function[psi, epsilon, mu, gamma] = H2norm_intermediate_variables(pi, pi1, z0)
    %% Variable definitions
    psi = pi(1);
    
    epsilon = z0(1);
    
    mu = pi1(1);
    
    gamma = pi(1);
end

function[pi, pi1, z, mu, epsilon, psi, gamma] = H2norm_first_iteration(pi, pi1, z0, mu, epsilon, psi, gamma, n)
    spi = [pi zeros(1, floor((n-1-1)/2))];
    z = psi*z0(2:end) - epsilon*spi(2:end);
    spi1 = [pi1 zeros(1, mod(n-1+1,2))];
    pi2 = mu*pi(2:end) - gamma*spi1(2:end);
    mu = pi2(1);
    epsilon = z(1);
    psi = pi1(1);
    gamma = pi1(1);

    pi = pi1;
    pi1 = pi2;
end

function[pi, pi1, z, mu, epsilon, psi, gamma] = H2norm_iterations(pi, pi1, z, mu, epsilon, psi, gamma, n)
    for i = 2:n-1
    % Recursive equations
        spi = [pi zeros(1, floor((n-i-1)/2))];
        z = psi*z(2:end) - epsilon*spi(2:end);
        spi1 = [pi1 zeros(1, mod(n-i+1,2))];
        pi2 = mu*pi(2:end) - gamma*spi1(2:end);

    % Variable definitions
        psi = pi1(1)/pi(1);        
        epsilon = z(1)/pi(1);
        mu = pi2(1)/pi(1);
        gamma = pi1(1)/pi(1);

        pi = pi1;
        pi1 = pi2;
    end
end