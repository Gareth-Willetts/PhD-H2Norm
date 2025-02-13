%% Define Polynomial
clear all; clc;
syms s
% digits(50)
time = [];
mean_time = [];
difference = [];

for j = 10:10:80
   for jj = 1:1000
        symb = (s+1)^(j);% * (s+1e-300);
        %symb = (s+100)*(s+23)*(s+2)*(s+1)*(s+0.023)*(s+574852000000);
        cd = sym2poly(symb);
        cn = ones(1,length(cd)-1);
        %cn(1) = 0;
        %cd = [306250000, 621118400000, 7050110278506400, 563527590388080000, 40239907967880000000, 90248344200000000000, 1421280000000000000000];       
        %cn = [0, 60436615680000, 729179561880000000, 90248344200000000000, 1421280000000000000000, 0];
        %cn = cn/cd(1);
        %cd = cd/cd(1);
        
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
       %normfunc = norm(tf(cn,cd),2)^2;
       H2n2 = (p1(1) * gamma)/(4*pkhat(1)*cd(1));
       %difference = [difference, H2n2-normfunc];
       %if H2n2 - normfunc > 0.01
       %    disp("ERROR");
       %end
       %time(j-4, jj) = toc;
       time(j-9*(j/10),jj) = toc;
   end
   %mean_time(j-4) = mean(time(j-4,10:end));
   mean_time(j-9*(j/10)) = mean(time(j-9*(j/10),10:end));
end

% p = polyfit(5:80, mean_time, 2);
% y = polyval(p, 5:80);
% plot(mean_time); hold on; plot(y);
% Gives answer exactly double what the H2 norm is - missing a factor of 2
% somewhere

%% Define Polynomial
syms s
digits(50)
% j = 50;
time = [];
mean_time_inbuilt = [];

for j = 10:10:80
   for jj = 1:100
       symb = (s+1)^(j-1) * (s+1e-300);
       cd = sym2poly(symb);
       cn = ones(1,length(cd)-1);
       tic
       % Calculate H2 norm
       norm(tf(cn,cd),2)^2;
       %time(j-4, jj) = toc;
       time(j-9*(j/10),jj) = toc;
   end
   %mean_time_inbuilt(j-4) = mean(time(j-4,10:end));
   mean_time_inbuilt(j-9*(j/10)) = mean(time(j-9*(j/10),10:end));
end


%%
figure(1)
plot(10:10:80, mean_time)
title("Polynomial Method Timings Averaged Over 1000 Runs (first 10 removed)")
xlabel("n")
ylabel("Time (s)")

figure(2)
plot(10:10:80, mean_time_inbuilt)
title("Inbuilt Method Timings Averaged Over 1000 Runs (first 10 removed)")
xlabel("n")
ylabel("Time (s)")

% p = polyfit(5:80, mean_time, 2);
% y = polyval(p, 5:80);
% plot(mean_time); hold on; plot(y);
% Gives answer exactly double what the H2 norm is - missing a factor of 2
% somewhere

% tic
% norm(tf(cn,cd),2)^2 %#ok<NOPTS>
% toc
% tic
% norm(tf(cn,cd),2)^2 %#ok<NOPTS>
% toc