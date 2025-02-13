import numpy as np
from sympy.ntheory import binomial_coefficients
import control as ct

def H2norm_p_creation(cd):
    # Degree of transfer function
    n = np.int64(len(cd) - 1);
    
    if n%2 == 1:
        # ODD CASE
        # Define initial p polynomials
        pi = np.append(cd[::2], 0);
        pi1 = np.array(cd[1::2]);
    else:
        # EVEN CASE
        # Define initial p polynomials
        pi = np.array(cd[::2]);
        pi1 = np.append(cd[1::2], 0);
    
    return n, pi, pi1

def H2norm_z0_creation(cn, n, symbolic=0):
    # Create z(s) depending on degree of transfer function
    if n%2 == 1:
        # define co, ce polynomials
        #ce = np.poly1d(cn[::2]);
        #co = np.poly1d(cn[1::2]);
        ce = cn[::2];
        co = cn[1::2];
        # square them
        #co2 = (co**2).c
        #ce2 = (ce**2).c
        co2 = np.convolve(co, co);
        ce2 = np.convolve(ce, ce);
        co2 = np.append(co2, 0);
        co2 = np.insert(co2, 0, 0);
        
        # Create z_0
        z0 = np.array([], dtype = 'int64');
        for i in range(len(cn)):
            z0 = np.append(z0, np.int64(ce2[i] - co2[i]));
    else:
        # define co, ce polynomials
        #ce = np.poly1d(cn[1::2]);
        #co = np.poly1d(cn[::2]);
        ce = cn[1::2];
        co = cn[::2];
        # square them
        #co2 = (co**2).c
        #ce2 = (ce**2).c
        co2 = np.convolve(co, co);
        ce2 = np.convolve(ce, ce);
        #co2 = np.append(co2, 0);
        #ce2 = np.append(ce2, 0);
        co2s = np.append(co2, 0);
        ns = len(co2s);
        
        for i in range(ns+1,n):
            co2s = np.insert(co2, 0, 0);
        
        ns = len(ce2);
        for i in range(ns+1,n+1):
            ce2 = np.insert(ce2, 0, 0);
        
        z0 = np.array([], dtype = 'int64');
        for i in range(len(cn)):
            z0 = np.append(z0, np.int64(ce2[i] - co2s[i]));
        
    return z0


def H2norm_intermediate_variables(pi, pi1, z0):           
    # Declaring variables needed for recursive equations - could delete to free memory?
    psi = pi[0];
    epsilon = z0[0];
    mu = pi1[0];
    gamma = pi[0];
    return psi, epsilon, mu, gamma

def H2norm_iterations(pi, pi1, z0, mu, epsilon, psi, gamma, n):
    # RECURSIVE EQUATIONS
    for i in range(n-1):
        if i == 0:
            # this works
            spi = np.append(pi, np.int64(np.zeros(int(np.floor((n-i-2)/2)))));
            z = np.subtract([x*psi for x in z0[1::]], [y*epsilon for y in spi[1::]]);
            # this doesn't
            spi1 = np.append(pi1, np.int64(np.zeros(int(n-i)%2)));
            pi2 = np.subtract([x*mu for x in pi[1::]], [y*gamma for y in spi1[1::]]);
            # define for next iteration
            mu = pi2[0];
            epsilon = z[0];
            psi = pi1[0];
            gamma = pi1[0];
            prevLC = np.int64(1);
            
            pi = pi1;
            pi1 = pi2;
        else:
            spi = np.append(pi, np.int64(np.zeros(int(np.floor((n-i-2)/2)))));
            z = np.subtract([x*psi for x in z[1::]], [y*epsilon for y in spi[1::]]);
            z = np.int64(np.divide(z, prevLC));
            spi1 = np.append(pi1, np.int64(np.zeros(int(n-i)%2)));
            pi2 = np.subtract([x*mu for x in pi[1::]], [y*gamma for y in spi1[1::]]);
            pi2 = np.int64(np.divide(pi2, prevLC));
            
            prevLC = pi[0];
            psi = pi1[0];
            epsilon = z[0];
            mu = pi2[0];
            gamma = pi1[0];
            
            pi = pi1;
            pi1 = pi2;

    H2n = np.double(z)/(2*np.double(cd[0])*np.double(pi1))
    return H2n


if __name__ == '__main__':
    lowerN = np.int64(3)
    upperN = np.int64(8)
    numIters = 1 #1000000
    H2nNew = np.zeros(upperN-lowerN+1,dtype=float)
    H2NewTime = np.zeros(shape=(numIters,upperN-lowerN+1),dtype=float)
    H2nInbuiltTime = np.zeros(shape=(numIters,upperN-lowerN+1),dtype=float)
    H2nInbuilt = np.zeros(upperN-lowerN+1,dtype=float)
    H2nAllN = np.zeros(upperN-lowerN+1,dtype=float)
    H2nAllNTime = np.zeros(shape=(numIters,upperN-lowerN+1),dtype=float)
    
    counter = 0
    
    import time
    for n in np.arange(lowerN, upperN+1):
        for t in range(numIters):
            cd = binomial_coefficients(n)
            cd = np.array(list(cd.values()))
            cd = np.unique(cd)
            cdreverse = cd[::-1]
            if n%2 == 1:
                cd = np.append(cd, cdreverse)
                cn = np.ones(len(cd)-1)
            else:
                cd = np.append(cd, cdreverse[1::])
                cn = np.ones(len(cd)-1)
                
            cd = np.int64(cd)
            cn = np.int64(cn)
            
            t0 = time.time()
            n, pi, pi1 = H2norm_p_creation(cd);
            z0 = H2norm_z0_creation(cn, n);
            psi, epsilon, mu, gamma = H2norm_intermediate_variables(pi, pi1, z0);
            H2nNew[counter] = H2norm_iterations(pi, pi1, z0, mu, epsilon, psi, gamma, n);
            t1 = time.time()
            
            t0 = time.time()
            tf = ct.tf(cn,cd)
            ss = ct.tf2ss(tf)
            L = ct.lyap(np.transpose(ss.A), np.transpose(ss.C)*ss.C)
            H2nInbuilt[counter] = np.trace(np.transpose(ss.B)*L*ss.B)
            t1 = time.time()
            
        counter = counter+1
       
    import matplotlib.pyplot as plt
    
    H2Diff = np.subtract(H2nNew, H2nInbuilt)
    H2PropDiff = np.divide(H2Diff, H2nNew)
    n = np.arange(lowerN, upperN+1)
    plt.xlabel("n")
    plt.ylabel("Proportional Difference")
    plt.title("Proportional Difference in H2 Norm")
    plt.plot(n, H2PropDiff)
    plt.rc('font', size = 22)
    plt.show()
    