import numpy as np

def findPotential(theta):

    k = [ 0.06450, 0.14975, 0.09458, 0.11251 ]
    n = [ 2, 3, 4, 5 ]
    delta = [ 0, np.pi, 0, 0 ]
    
    # Conversion factor from kcal/mol to J
    
    # conv = 6.9477e-21
    conv = 1
    
    # Initialization for V
    
    V = 0
    
    for i in range(0, 4):
    
        V = V + conv * k[i] * (1 + np.cos(n[i]*theta - delta[i]))
        
    return V