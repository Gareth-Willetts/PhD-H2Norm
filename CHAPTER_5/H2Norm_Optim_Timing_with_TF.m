clear all; close all; clc
%% Following works well for obtaining symbolic TF
% Define symbols for components in the suspension system
time_tot = zeros(1,100);
time_TF = zeros(1,100);
time_H2 = zeros(1,100);
time_optim = zeros(1,100);
for i = 1:100
    tic
    Train_Graph_TF
    time_TF(i) = toc;
    tic
    x=[x1;x2];
    
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
    [n, pi, pi1] = H2norm_p_creation(cd);
    %%
    z0 = H2norm_z0_creation(cn,n,symbolic);
    %%
    [psi, epsilon, mu, gamma] = H2norm_intermediate_variables(pi, pi1, z0);
    %%
    [pi, pi1, z, mu, epsilon, psi, gamma] = H2norm_first_iteration(pi, pi1, z0, mu, epsilon, psi, gamma, n);
    %%
    [pi, pi1, z, mu, epsilon, psi, gamma] = H2norm_iterations(pi, pi1, z, mu, epsilon, psi, gamma, n);
    %% Calculate H2 norm
    H2 = z/(2*cd(1)*pi1);
    time_H2(i) = toc;
    %x = 7616;
    %y = 11844;
    %H2nc2 = sqrt(double(subs(H2)))
    
    %% Optimisation?
    % tic
    % gradf = jacobian(H2, [x]); % column gradf
    % hessf = jacobian(gradf);
    % fh = matlabFunction(H2,gradf,hessf,'vars',{x});
    % options = optimoptions('fminunc', ...
    %     'SpecifyObjectiveGradient', true, ...
    %     'HessianFcn', 'objective', ...
    %     'Algorithm','trust-region', 'Display', 'none');%, ...
    %     %'Display','final');
    % [xfinal,fval,exitflag,output] = fminunc(fh,[7000; 11000],options);
    % time_optim(i) = toc;
    tic
    x0 = [7000;11000];
    fun = matlabFunction(H2, "Vars", {[x]});
    options = optimoptions('fminunc','SpecifyObjectiveGradient',false,...
        'Display', 'none');
    [x,fval] = fminunc(fun,x0,options);
    time_optim(i) = toc;
end
time_tot = time_TF + time_H2 + time_optim;
mean(time_tot(10:end))
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