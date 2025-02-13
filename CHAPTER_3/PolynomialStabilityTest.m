function isStable = PolynomialStabilityTest(cd, cn)
    % Degree of transfer function = length of denominator - 1 (due to the
    % constant term)
    n = size(cd,2)-1;
    isStable = 1;
    if mod(n,2) == 1
        % ODD
        pi1 = cd(2:2:end);
        pi = [cd(1:2:end) 0];
    else
        % EVEN
        pi = cd(1:2:end);
        pi1 = [cd(2:2:end) 0];
    end

    if pi(1) <= 0 || pi1(1) <= 0
        isStable = 0;
    else
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
    
            % Preallocate memory for b, which will be populated with the coefficients
            % of z_0
            z0 = zeros(1,length(cn));
    
            % Populates b with values of z_0
            for i = 1:length(cn)
                z0(i) = ce2(i) - co2s(i);
            end
        end
    %% Variable definitions
        psi = pi(1);
    
        epsilon = z0(1);
    
        mu = pi1(1);
    
        gamma = pi(1);
    %% Better memory allocation?
        for i = 1:n-1
        % Special case for i = 1
            if i == 1
                spi = [pi zeros(1, floor((n-i-1)/2))];
                z = psi*z0(2:end) - epsilon*spi(2:end);
                spi1 = [pi1 zeros(1, mod(n-i+1,2))];
                pi2 = mu*pi(2:end) - gamma*spi1(2:end);
                mu = pi2(1);
                if mu <= 0
                    isStable = 0;
                    break
                end
                epsilon = z(1);
                psi = pi1(1);
                gamma = pi1(1);
    
                pi = pi1;
                pi1 = pi2;
            else
            % Recursive equations
                spi = [pi zeros(1, floor((n-i-1)/2))];
                z = psi*z(2:end) - epsilon*spi(2:end);
                spi1 = [pi1 zeros(1, mod(n-i+1,2))];
                pi2 = mu*pi(2:end) - gamma*spi1(2:end);
                if pi2(1) <= 0
                    isStable = 0;
                    break
                else
                % Variable definitions
                    psi = pi1(1)/pi(1);        
                    epsilon = z(1)/pi(1);
                    mu = pi2(1)/pi(1);
                    gamma = pi1(1)/pi(1);
        
                    pi = pi1;
                    pi1 = pi2;
                end
            end
        end
    end
end