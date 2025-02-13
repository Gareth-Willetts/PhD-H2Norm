%% Define Polynomial
clear all; clc;
syms s
counter = 1;
j = 5:10:45;
mean_time = zeros(1,length(j));
H2n2_TC = zeros(1,length(j));
H2inbuilt = zeros(1,length(j));
for j = 5:10:45
   for jj = 1:1000
        symb = (s+1)^(j);
        cd = sym2poly(symb);
        cn = ones(1,length(cd)-1);      
        %%
        % Degree of transfer function = length of denominator - 1
        tic
        n = size(cd,2)-1;
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
            % Define co and ce as the odd and even parts of the numerator polynomial
            co = cn(2:2:end);
            ce = cn(1:2:end);
        
            % Square these polynomials, and add zeros to start and end of co2
            co2 = conv(co, co);
            ce2 = conv(ce, ce);
            co2 = [0 co2 0];
        
            % Preallocate memory which will be populated with the coefficients
            % of z_0
            z0 = zeros(1,length(cn));
        
            % Populates b with values of z_0
            for i = 1:length(cn)
                z0(i) = ce2(i) - co2(i);
            end
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
        
            % Preallocate memory which will be populated with the coefficients
            % of z_0
            z0 = zeros(1,length(cn));
        
            % Populates b with values of z_0
            for i = 1:length(cn)
                z0(i) = ce2(i) - co2s(i);
            end
        end
        %% Variable definitions
        nu = z0(end);
        gamma = nu;
        z0hat = z0/nu;
        nu = p2(end);
        gamma = gamma/nu;
        p2hat = p2/nu;
        
        %% Now z_1
        numZeros = length(z0hat) - length(p1);
        z1tilde = z0hat - [p1 zeros(1,numZeros)]*z0hat(1)/p1(1);
        nu = z1tilde(end);
        gamma = gamma*nu;
        z1hat = z1tilde(2:end)/nu;
        
        numZeros = length(p1) - length(p2hat);
        p3tilde = p1 - [p2hat zeros(1,numZeros)]*p1(1)/p2hat(1);
        nu = p3tilde(end);
        gamma = gamma/nu;
        p3hat = p3tilde(2:end)/nu;
        
        for k = 3:n
            if k == 3
                numZeros = length(z1hat) - length(p2hat);
                zk1tilde = z1hat - [p2hat zeros(1,numZeros)]*z1hat(1)/p2hat(1);
                nu = zk1tilde(end);
                gamma = gamma*nu;
                zk1hat = zk1tilde(2:end)/nu;
        
                numZeros = length(p2hat) - length(p3hat);
                pkplustilde = p2hat - [p3hat zeros(1,numZeros)]*p2hat(1)/p3hat(1);
                nu = pkplustilde(end);
                gamma = gamma/nu;
                pkplushat = pkplustilde(2:end)/nu;
        
        % Setting up variables for iterative calculations
                pkminushat = p3hat;
                pkhat = pkplushat;
            else
                numZeros = length(zk1hat) - length(pkminushat);
                zk1tilde = zk1hat - [pkminushat zeros(1,numZeros)]*zk1hat(1)/pkminushat(1);
                nu = zk1tilde(end);
                gamma = gamma*nu;
                zk1hat = zk1tilde(2:end)/nu;
        
                numZeros = length(pkminushat) - length(pkhat);
                pkplustilde = pkminushat - [pkhat zeros(1,numZeros)]*pkminushat(1)/pkhat(1);
                nu = pkplustilde(end);
                gamma = gamma/nu;
                pkplushat = pkplustilde(2:end)/nu;
                pkminushat = pkhat;
                pkhat = pkplushat;
            end
        end
       %% Calculate H2 norm
       H2n2_TC(counter) = (p1(1) * gamma)/(2*pkhat(1)*cd(1));
       t(counter,jj) = toc;
   end
   mean_time(counter) = mean(t(counter,10:end));
   H2inbuilt(counter) = norm(tf(cn,cd),2)^2;
   counter = counter + 1;
end

%%
figure(1)
plot(5:10:45, mean_time)
title("Polynomial Method Timings Averaged Over 1000 Runs (first 10 removed)")
xlabel("n")
ylabel("Time (s)")

% p = polyfit(5:80, mean_time, 2);
% y = polyval(p, 5:80);
% plot(mean_time); hold on; plot(y);
% Gives answer exactly double what the H2 norm is - missing a factor of 2
% somewhere