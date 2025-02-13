clear all; close all; clc
syms s x
symbolic = 0;
time = [];
time_inbuilt = [];
mean_time_modified = [];
mean_time = [];
numIters = 10000;

for j = 5:45
    symb = (s+1)^j;
    cd = sym2poly(symb);
    cn = ones(1,length(cd)-1); %*x;
    for jj = 1:numIters
        tic
        [n, pi, pi1] = H2norm_p_creation(cd);
        z0 = H2norm_z0_creation(cn,n,symbolic);
        [psi, epsilon, mu, gamma] = H2norm_intermediate_variables(pi, pi1, z0);
        [pi, pi1, z, mu, epsilon, psi, gamma] = H2norm_first_iteration(pi, pi1, z0, mu, epsilon, psi, gamma, n);
        [pi, pi1, z, mu, epsilon, psi, gamma] = H2norm_iterations(pi, pi1, z, mu, epsilon, psi, gamma, n);
        %% Calculate H2 norm
        H2 = z/(2*cd(1)*pi1);
        time(j-4, jj) = toc;
    end
    mean_time(j-4) = mean(time(j-4,10:end));
end

for j = 5:45
   for jj = 1:numIters
        symb = (s+1)^(j);
        cd = sym2poly(symb);
        cn = ones(1,length(cd)-1);
        %%
        % Degree of transfer function = length of denominator - 1 (due to the
        % constant term)
        tic
        n = size(cd,2)-1;
        % This is correct
        if mod(n,2) == 1
            % ODD
            p2 = cd(2:2:end);
            p1 = [cd(1:2:end) 0];
        else
            % EVEN
            p1 = cd(1:2:end);
            p2 = [cd(2:2:end) 0];
        end
        
        %% Create z(s) from equation in the paper
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
        % Unchanged, so correct
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
       H2n2 = (p1(1) * gamma)/(4*pkhat(1)*cd(1));
       time(j-4,jj) = toc;
   end
   mean_time_modified(j-4) = mean(time(j-4,10:end));
end

%%
figure(1)
plot(5:45, mean_time_modified, 'g')
title("Polynomial Method (Modified) Timings Averaged Over 10000 Runs (first 10 removed)")
xlabel("n")
ylabel("Time (s)")
set(gca, "FontSize", 20)

figure(2)
plot(5:45, mean_time, 'b')
title("Polynomial Method Timings Averaged Over 10000 Runs (first 10 removed)")
xlabel("n")
ylabel("Time (s)")
set(gca, "FontSize", 20)

figure(3)
plot(5:45, log(mean_time_modified), 'g'); hold on;
plot(5:45, log(mean_time), 'b')
title("Different Method Timings Averaged Over 10000 Runs (first 10 removed)")
xlabel("n")
ylabel("log(Time (s))")
set(gca, "FontSize", 20)
legend("New Method (Modified)", "New Method")

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
            % Square these polynomials, and add zeros to start and end of co2
            co2 = poly2sym(co,t)*poly2sym(co,t);
            ce2 = poly2sym(ce,t)*poly2sym(ce,t);
            co2 = coeffs(co2, t, 'All');
            ce2 = coeffs(ce2, t,'All');
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