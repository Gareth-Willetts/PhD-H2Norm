import numpy as np
from sympy.ntheory import binomial_coefficients
import control as ct

def H2norm_p_creation(cd):
    # Degree of transfer function
    n = len(cd) - 1;
    
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
        z0 = []
        for i in range(len(cn)):
            z0 = np.append(z0, ce2[i] - co2[i]);
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
        
        z0 = [];
        for i in range(len(cn)):
            z0 = np.append(z0, ce2[i] - co2s[i]);
        
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
            spi = np.append(pi, np.zeros(int(np.floor((n-i-2)/2))));
            z = np.subtract([x*psi for x in z0[1::]], [y*epsilon for y in spi[1::]]);
            # this doesn't
            spi1 = np.append(pi1, np.zeros(int(n-i)%2));
            pi2 = np.subtract([x*mu for x in pi[1::]], [y*gamma for y in spi1[1::]]);
            # define for next iteration
            mu = pi2[0];
            epsilon = z[0];
            psi = pi1[0];
            gamma = pi1[0];
            
            pi = pi1;
            pi1 = pi2;
        else:
            spi = np.append(pi, np.zeros(int(np.floor((n-i-2)/2))));
            z = np.subtract([x*psi for x in z[1::]], [y*epsilon for y in spi[1::]]);
            spi1 = np.append(pi1, np.zeros(int(n-i)%2));
            pi2 = np.subtract([x*mu for x in pi[1::]], [y*gamma for y in spi1[1::]]);
            
            psi = np.divide(pi1[0],pi[0]);
            epsilon = np.divide(z[0],pi[0]);
            mu = np.divide(pi2[0],pi[0]);
            gamma = np.divide(pi1[0],pi[0]);
            
            pi = pi1;
            pi1 = pi2;

    H2n = z/(2*cd[0]*pi1)
    return H2n[0]

def H2norm_first_iterations_allN(z0,p1,p2):
    nu = np.mean(abs(z0));
    gamma = nu;
    z0hat = z0/nu;
    nu = np.mean(abs(p2));
    gamma = gamma/nu;
    p2hat = p2/nu;
    
    numZeros = len(z0hat) - len(p1);
    z1tilde = z0hat - np.append(p1, np.zeros(numZeros))*z0hat[0]/p1[0];
    nu = np.mean(abs(z1tilde[1::]));
    gamma = gamma*nu;
    z1hat = z1tilde[1::]/nu;
    
    numZeros = len(p1) - len(p2hat);
    p3tilde = p1 - np.append(p2hat, np.zeros(numZeros))*p1[0]/p2hat[0];
    nu = np.mean(abs(p3tilde[1::]));
    gamma = gamma/nu;
    p3hat = p3tilde[1::]/nu;
    
    return z1hat, p2hat, p3hat, nu, gamma

def H2norm_iterations_allN(z1hat, p2hat, p3hat, nu, gamma):
    for k in range(3, n+1):
        if k == 3:
            numZeros = len(z1hat) - len(p2hat);
            zk1tilde = z1hat - np.append(p2hat, np.zeros(numZeros))*z1hat[0]/p2hat[0];
            nu = np.mean(abs(zk1tilde[1::]));
            gamma = gamma*nu;
            zk1hat = zk1tilde[1::]/nu;
            
            numZeros = len(p2hat) - len(p3hat);
            pkplustilde = p2hat - np.append(p3hat, np.zeros(numZeros))*p2hat[0]/p3hat[0];
            nu = np.mean(abs(pkplustilde[1::]));
            gamma = gamma/nu;
            pkplushat = pkplustilde[1::]/nu;
            
            pkminushat = p3hat;
            pkhat = pkplushat;
        else:
            numZeros = len(zk1hat) - len(pkminushat);
            zk1tilde = zk1hat - np.append(pkminushat, np.zeros(numZeros))*zk1hat[0]/pkminushat[0];
            nu = np.mean(abs(zk1tilde[1::]));
            gamma = gamma*nu;
            zk1hat = zk1tilde[1::]/nu;
            
            numZeros = len(pkminushat) - len(pkhat);
            pkplustilde = pkminushat - np.append(pkhat, np.zeros(numZeros))*pkminushat[0]/pkhat[0];
            nu = np.mean(abs(pkplustilde[1::]));
            gamma = gamma/nu;
            
            pkplushat = pkplustilde[1::]/nu;      
            pkminushat = pkhat;
            pkhat = pkplushat;
            
    H2n = (p1[0] * gamma)/(4*pkhat[0]*cd[0])
    return H2n

if __name__ == '__main__':
    lowerN = 5
    upperN = 45
    numIters = 1 #1000000
    H2nNew = np.zeros(upperN-lowerN+1,dtype=float)
    H2nInbuilt = np.zeros(upperN-lowerN+1,dtype=float)
    H2nAllN = np.zeros(upperN-lowerN+1,dtype=float)
    
    for n in range(lowerN, upperN+1):
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
            
            n, pi, pi1 = H2norm_p_creation(cd);
            z0 = H2norm_z0_creation(cn, n);
            psi, epsilon, mu, gamma = H2norm_intermediate_variables(pi, pi1, z0);
            H2nNew[n-lowerN-1] = H2norm_iterations(pi, pi1, z0, mu, epsilon, psi, gamma, n);
            
            tf = ct.tf(cn,cd)
            ss = ct.tf2ss(tf)
            L = ct.lyap(np.transpose(ss.A), np.transpose(ss.C)*ss.C)
            H2nInbuilt[n-lowerN-1] = np.trace(np.transpose(ss.B)*L*ss.B)
        
            n, p1, p2 = H2norm_p_creation(cd);
            z0 = H2norm_z0_creation(cn, n);
            z1hat, p2hat, p3hat, nu, gamma = H2norm_first_iterations_allN(z0, p1, p2);
            
            
            H2nAllN[n-lowerN-1] = H2norm_iterations_allN(z1hat, p2hat, p3hat, nu, gamma);
       
    import matplotlib.pyplot as plt
    
    H2Diff = np.divide(np.subtract(H2nNew, H2nInbuilt), H2nNew)
    H2DiffAllN = np.divide(np.subtract(H2nAllN, H2nInbuilt), H2nAllN)
    # H2DiffNewMethods = np.subtract(H2nAllN, H2nNew)
    n = np.arange(lowerN, upperN+1)
    plt.xlabel("n")
    plt.ylabel("H2 norm")
    plt.plot(n, H2Diff)
    plt.plot(n, H2DiffAllN)
    # plt.plot(H2DiffNewMethods)
    plt.legend(["Symbolic Polynomial Method - Python's Inbuilt Method", "Numerical Polynomial Method - Python's Inbuilt method"])#, "Refined method - New method"])
    plt.rc('font', size = 22)
    plt.show()
    