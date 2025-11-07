import numpy as np
import pandas as pd
from scipy.optimize import least_squares

def polyfit(time, mag, parameters):
    """
        Fit a combination of algebraic and trigonometric polynomials
        using scipy.optimize.least_squares instead of statsmodels.OLS.
    
        parameters = [deg_trend,
                      period1, deg1,
                      period2, deg2,
                      period3, deg3]
    """    
    t_mean = np.mean(time)
    temp_t = time - t_mean
    fit_arguments = []
    term_labels = []  # to remember what each column is
    
    # Algebraic polynomial part
    for i in range(parameters[0] + 1):
        fit_arguments.append(temp_t ** i)
        term_labels.append(f"poly__{i}")
        #print("Trend degree ", i, " added")
    
    # Up to 3 different periods
    for j in range(3):
        if parameters[j * 2 + 1] > 0.0 and parameters[j * 2 + 2] > 0:
            for i in range(1, parameters[j * 2 + 2] + 1):
                a = 2 * np.pi * i / parameters[j * 2 + 1] * temp_t
                fit_arguments.append(np.cos(a))
                term_labels.append(f"cos_{j}_{i}")
                fit_arguments.append(np.sin(a))
                term_labels.append(f"sin_{j}_{i}")

    #print(fit_arguments)
    #print("Total number of basis functions:", len(fit_arguments))
    X = np.column_stack(fit_arguments)
    
    cond = np.linalg.cond(X)
    print("Condition number:", cond)
    if cond > 1e10:
        print("**** WARNING!! The solution could be incorrect")
    
    # Define the residuals function for least_squares
    def residuals(params):
        return mag - X @ params
    
    # Initial guess for parameters
    p0 = np.zeros(X.shape[1])
    
    # Solve least squares
    result = least_squares(residuals, p0)
    coeffs = result.x
    
    fitted = X @ coeffs

    # --- PRINT ALL COEFFICIENTS ---
    print("\n=== FITTED COEFFICIENTS ===")
    for label, c in zip(term_labels, coeffs):
       print(f"Term {label}: {c}")
    print("============================\n")

    
    fit_result = pd.DataFrame({
        'Time': time, 
        'Mag': mag, 
        'Fit': fitted,
    })
    
    return fit_result
