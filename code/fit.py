import numpy as np
import pandas as pd
from scipy.optimize import minimize

# -------------------------
# Inner linear fit (with covariances)
# -------------------------
def polyfit_fixed_periods(time, mag, alg_poly, periods, degrees, return_X):
    # Centering and scaling time
    t_mean = np.mean(time)
    #temp_t = time - t_mean
    t_scale = np.max(np.abs(time - t_mean))
    #t_scale = 1.0
    temp_t = (time - t_mean) / t_scale # now |temp_t| <= 1
    
    fit_arguments = []

    # Algebraic polynomial
    for i in range(alg_poly + 1):
        fit_arguments.append(temp_t ** i)

    # Trigonometric part
    for j in range(len(periods)):
        period = periods[j] / t_scale
        degree = int(degrees[j])
        for i in range(1, degree + 1):
            a = 2 * np.pi * i / period * temp_t
            fit_arguments.append(np.cos(a))
            fit_arguments.append(np.sin(a))

    X = np.column_stack(fit_arguments)

    # Solve least squares
    coeffs, _, _, _ = np.linalg.lstsq(X, mag, rcond=None)
    fitted = X @ coeffs
    residuals = mag - fitted
    rss = np.sum(residuals**2)

    # Transform coefficients back
    # Algebraic polynomial
    for i in range(alg_poly + 1):
        coeffs[i] = coeffs[i] / (t_scale ** i)
        
    if return_X:
        for i in range(alg_poly + 1):
            X[:,i] = X[:,i] * (t_scale ** i)

    if return_X:
        return rss, fitted, coeffs, X
    else:
        return rss, fitted, coeffs

# -------------------------
# Utility: numerical Jacobian of residual vector wrt selected period variables
# -------------------------
def numerical_jacobian_residuals(time, mag, alg_poly, b_periods, degrees, period_indices, eps_rel=1e-6):
    # ChatGPT comment!
    """
    Compute Jacobian J (N x k) of residual vector r = mag - X@beta
    w.r.t. the vector of period variables (period_indices list).
    Use finite differences on periods. eps_rel relative step.
    """
    # base residuals and fitted
    _, fitted0, _ = polyfit_fixed_periods(time, mag, alg_poly, b_periods, degrees, False)
    #r0 = mag - fitted0
    N = len(mag)
    k = len(period_indices)
    J = np.zeros((N, k), dtype=float)

    for col, idx in enumerate(period_indices):
        p0 = b_periods[idx]
        # choose stepsize proportional to period
        h = max(abs(p0) * eps_rel, eps_rel)
        periods_pert = b_periods.copy()
        periods_pert[idx] = p0 + h
        # compute residuals at perturbed period
        _, fitted_p, _ = polyfit_fixed_periods(time, mag, alg_poly, periods_pert, degrees, False)
        #rp = mag - fitted_p
        # finite difference column
        #J[:, col] = (rp - r0) / h
        J[:, col] = (fitted0 - fitted_p) / h

    return J


def print_initial_periods(printLog, init_periods):
    # Initial periods
    printLog()
    if len(init_periods) < 1:
        return
    printLog("=== INITIAL PERIODS ============================================================")
    for i, pval in enumerate(init_periods):
        period_name = f"period_{i+1}"
        printLog(f"{period_name} = {pval:.12f}")
    printLog("================================================================================")
    printLog()
    
def print_optimized_periods(printLog, best_periods, se_periods, period_indices):
    printLog()
    if len(period_indices) < 1:
        return
    printLog("=== OPTIMIZED PERIODS ==========================================================")
    for i, ii in enumerate(period_indices):
        period_name = f"period_{ii+1}"
        if se_periods is not None:
            se = se_periods[i]
        else:
            se = np.nan
        printLog(f"{period_name} = {best_periods[ii]:.12f} ± {se:.12f} (1σ)")
    printLog("================================================================================")
    printLog()
    
def print_linear_coefficients_and_amplitudes(printLog, N, X, alg_poly, 
                                             initial_periods, degrees, 
                                             best_periods, best_period_indices, 
                                             coeffs, rss, sigma2, dof):
    # Errors of coefficients
    cov_coeffs = sigma2 * np.linalg.inv(X.T @ X)
    se_coeffs = np.sqrt(np.diag(cov_coeffs))

    p_lin = len(coeffs)
    term_labels = []
    for j in range(0, alg_poly + 1):
        term_labels.append(f'poly_{j}')
    for j in range(alg_poly + 1, p_lin, 2):
        n = (j - alg_poly - 1) // 2 + 1
        term_labels.append(f'cos_{n}')
        term_labels.append(f'sin_{n}')
    
    amplitudes = []
    se_amplitudes = []
    ampl_labels = []
    # Amplitudes of harmonics
    for j in range(alg_poly + 1, p_lin, 2):
        sub_cov_coeffs = cov_coeffs[j:j+2, j:j+2]
        amplitude = np.sqrt(coeffs[j]**2 + coeffs[j+1]**2)
        amplitudes.append(amplitude)
        #print(amplitude)
        g = np.array([coeffs[j]/amplitude, coeffs[j+1]/amplitude])
        #print(sub_cov_coeffs)
        #print(g)
        se_amplitude = np.sqrt(g @ sub_cov_coeffs @ g)
        se_amplitudes.append(se_amplitude)
        #print(se_amplitude)
        n = (j - alg_poly - 1) // 2 + 1
        ampl_labels.append(f'Ampl_{n}')
    
    printLog()
    printLog("=== LINEAR COEFFICIENTS ========================================================")
    for lbl, c, se in zip(term_labels, coeffs, se_coeffs):
        printLog(f"{lbl:12s}: {c:.8e} ± {se:.3e}")
    printLog("================================================================================")                        
    printLog()
    if len(initial_periods) > 0:
        period_infos = []
        for j in range(len(initial_periods)):
            period = initial_periods[j]
            if j in best_period_indices:
                period = best_periods[best_period_indices.index(j)]
            for d in range(degrees[j]):
                period_infos.append(str(d + 1) + " * " + str(period))
        #print(period_infos)
        
        printLog("=== AMPLITUDES =================================================================")
        for i, (lbl, a, se) in enumerate(zip(ampl_labels, amplitudes, se_amplitudes)):
            period_info = period_infos[i]
            printLog(f"{lbl:12s}: {a:.8e} ± {se:.3e}  Period: {period_info}")
        printLog("================================================================================")            
        printLog()
    n_panams = N - dof
    printLog(f"RSS [Sum(O-C)^2] = {rss:.6e}, sigma = {(sigma2**0.5):.6e}")
    printLog(f"Npoints = {N}, dof = {dof}, Nparams = {n_panams}")
    printLog(f"R.M.S. accuracy of the fit sigma[x_c] = {np.sqrt((n_panams / N / dof) * rss)}")
    printLog()
    printLog("================================================================================")
    printLog()

# -------------------------
# Outer optimizer with error estimates
# -------------------------
def optimize_periods_with_errors(time, mag,
                                 alg_poly,
                                 initial_periods,
                                 degrees,
                                 optimize_flags,
                                 method='Nelder-Mead',
                                 maxiter=2000, xtol=1e-8, ftol=1e-8,
                                 compute_bootstrap=False, n_bootstrap=200):
    """
    Outer nonlinear optimization of periods (returns error estimates).
    """

    log_message = ""
    
    def printLog(msg=""):
        print(msg)
        nonlocal log_message
        log_message += f"{msg}\n"

    # Initial periods
    print_initial_periods(printLog, initial_periods)

    # The period_indices list contains indices of poriods to be optimized
    period_indices = []    
    # Periods to optimize
    for i in range(len(optimize_flags)):
        if optimize_flags[i]:
            period_indices.append(i)

    if len(period_indices) == 0:
        # nothing to optimize -> just return linear fit results + cov
        rss, fitted, coeffs, X = polyfit_fixed_periods(time, mag, 
                                                       alg_poly, 
                                                       initial_periods, 
                                                       degrees, 
                                                       True)
        N = len(mag)
        p = len(coeffs)
        dof = max(N - p, 1)
        sigma2 = rss / dof
        print_linear_coefficients_and_amplitudes(printLog, N, X, alg_poly, initial_periods, degrees, [], [], coeffs, rss, sigma2, dof)
        return {
            'fitted': fitted, 
            'message': log_message
        }

    # initial guess vector for the optimizer
    x0 = np.array([initial_periods[i] for i in period_indices], dtype=float)

    def objective(period_values):
        # update param vector
        periods = initial_periods.copy()
        for idx, val in zip(period_indices, period_values):
            periods[idx] = float(val)
        rss, _, _ = polyfit_fixed_periods(time, mag, alg_poly, periods, degrees, False)
        return rss

    # run optimizer
    #print("Optimizer...")
    result = minimize(objective, x0, method=method,
                      options={'maxiter': maxiter, 'xatol': xtol, 'fatol': ftol, 'disp': True})

    #print(result)
    best_period_values = result.x
    # assemble full parameter vector with best periods plugged
    best_periods = initial_periods.copy()
    #print(best_periods)
    #print(best_period_values)
    for idx, val in zip(period_indices, best_period_values):
        best_periods[idx] = float(val)

    #print("Before final linear fit...")
    # final linear fit at best params (and get X, coeffs)
    rss, fitted, coeffs, X = polyfit_fixed_periods(time, mag, 
                                                   alg_poly,
                                                   best_periods,
                                                   degrees,
                                                   True)
    N = len(mag)
    p_lin = len(coeffs)
    p_tot = p_lin + len(period_indices)
    dof = max(N - p_tot, 1)   # conservative dof include nonlinear params
    sigma2 = rss / dof

    # Now compute Jacobian J of residuals wrt periods (finite differences)
    J = numerical_jacobian_residuals(time, mag, alg_poly, best_periods, degrees, period_indices)

    # This is ChatGPT comments!
    # Gauss-Newton approximation: Hessian ≈ 2 * J^T J, but covariance for periods:
    # Cov(periods) ≈ sigma2 * (J^T J)^{-1}
    JTJ = J.T @ J
    # handle possible singularity
    try:
        cov_periods = sigma2 * np.linalg.inv(JTJ)
        se_periods = np.sqrt(np.diag(cov_periods))
    except Exception:
        cov_periods = None
        se_periods = None
        printLog("Warning: J^T J is singular — cannot compute period covariance via Gauss-Newton.")

    # Print results
    print_optimized_periods(printLog, best_periods, se_periods, period_indices)

    print_linear_coefficients_and_amplitudes(printLog, N, X, alg_poly, initial_periods, degrees, best_periods, period_indices, coeffs, rss, sigma2, dof)

    output = {
        'fitted': fitted,
        'message': log_message
    }

    # Optional: bootstrap for period uncertainties (more robust)
    if compute_bootstrap:
        printLog("Running bootstrap...")
        print("(this may take a while)...")
        boot_periods = []
        rng = np.random.default_rng()
        for b in range(n_bootstrap):
            # sample residuals (residual bootstrap)
            # generate bootstrap mag = fitted + resampled residuals
            resampled_r = rng.choice((mag - fitted), size=N, replace=True)
            mag_b = fitted + resampled_r
            # re-optimize periods starting from best_period_values (fast)
            def obj_bpv(period_values):
                periods_b = initial_periods.copy()
                for idx, val in zip(period_indices, period_values):
                    periods_b[idx] = float(val)
                rss_b, _, _ = polyfit_fixed_periods(time, mag_b, alg_poly, periods_b, degrees, False)
                return rss_b
            #
            res_b = minimize(obj_bpv, best_period_values, method=method,
                             options={'maxiter': maxiter, 'xatol': xtol, 'fatol': ftol, 'disp': False})
            boot_periods.append(res_b.x)
        boot_periods = np.array(boot_periods)
        # compute bootstrap std dev for each optimized period
        boot_se = np.std(boot_periods, axis=0, ddof=1)
        printLog(f"Bootstrap periods standard errors: {boot_se}")
        #output['bootstrap_period_se'] = boot_se
        #output['bootstrap_periods'] = boot_periods
        output['message'] = log_message

    return output

# -------------------------
# Wrapper
# -------------------------
def polyfit(time, mag, alg_poly, periods, degrees, optimize_flags, **kwargs):
    normalized_periods = []
    normalized_degrees = []
    normalized_optimize_flags = []
    for period, degree, optimize in zip(periods, degrees, optimize_flags):
        if period > 0.0 and degree > 0:
            normalized_periods.append(period)
            normalized_degrees.append(degree)
            normalized_optimize_flags.append(optimize)
    
    out = optimize_periods_with_errors(time, mag, 
                                       alg_poly, 
                                       normalized_periods, 
                                       normalized_degrees, 
                                       normalized_optimize_flags,
                                       **kwargs)
    fit_result = pd.DataFrame({
        'Time': time,
        'Mag': mag,
        'Fit': out['fitted'],
    })
    return {'fit_result': fit_result, 'message': out['message']}
