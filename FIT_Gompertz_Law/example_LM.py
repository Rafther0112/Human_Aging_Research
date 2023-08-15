# Code adapted from Gavin, H.P. (2020) The Levenberg-Marquardt method for 
# nonlinear least squares curve-fitting problems.
# https://people.duke.edu/~hpgavin/ce281/lm.pdf
    
import numpy as np
import matplotlib.pyplot as plt
import levenberg_marquardt as LM

def main(x,y,p_init):
    """
    
    Main function for performing Levenberg-Marquardt curve fitting.

    Parameters
    ----------
    x           : x-values of input data (m x 1), must be 2D array
    y           : y-values of input data (m x 1), must be 2D array
    p_init      : initial guess of parameters values (n x 1), must be 2D array
                  n = 4 in this example

    Returns
    -------
    p       : least-squares optimal estimate of the parameter values
    Chi_sq  : reduced Chi squared error criteria - should be close to 1
    sigma_p : asymptotic standard error of the parameters
    sigma_y : asymptotic standard error of the curve-fit
    corr    : correlation matrix of the parameters
    R_sq    : R-squared cofficient of multiple determination  
    cvg_hst : convergence history (col 1: function calls, col 2: reduced chi-sq,
              col 3 through n: parameter values). Row number corresponds to
              iteration number.

    """
    # close all plots
    plt.close('all')
    
    # minimize sum of weighted squared residuals with L-M least squares analysis
    p_fit,Chi_sq,sigma_p,sigma_y,corr,R_sq,cvg_hst = LM.lm(p_init,x,y)
    
    # plot results of L-M least squares analysis
    LM.make_lm_plots(x, y, cvg_hst)
    
    return p_fit,Chi_sq,sigma_p,sigma_y,corr,R_sq,cvg_hst
    
if __name__ == '__main__':
    
    x = 0
    y = 0
    p_init = [1,1]   

    p_fit,Chi_sq,sigma_p,sigma_y,corr,R_sq,cvg_hst = main(x,y,p_init)
