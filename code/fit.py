import numpy as np
import pandas as pd
from scipy.optimize import minimize

# -------------------------
# Inner linear fit (with covariances)
# -------------------------
def polyfit_fixed_periods(time, mag, alg_poly, periods, degrees, return_X):
    t_mean = np.mean(time)
    temp_t = time - t_mean
    fit_arguments = []
    term_labels = []

    # Algebraic polynomial
    for i in range(alg_poly + 1):
        fit_arguments.append(temp_t ** i)
        term_labels.append(f"poly_{i}")

    # Up to 3 trigonometric series
    for j in range(len(periods)):
        period = periods[j]
        degree = int(degrees[j])
        for i in range(1, degree + 1):
            a = 2 * np.pi * i / period * temp_t
            fit_arguments.append(np.cos(a))
            term_labels.append(f"cos_{j+1}_{i}")
            fit_arguments.append(np.sin(a))
            term_labels.append(f"sin_{j+1}_{i}")

    X = np.column_stack(fit_arguments) if len(fit_arguments) > 0 else np.zeros((len(time), 0))

    # Solve least squares
    coeffs, residuals_sum, rank, s = np.linalg.lstsq(X, mag, rcond=None)
    fitted = X @ coeffs if X.size else np.zeros_like(mag)
    residuals = mag - fitted
    rss = np.sum(residuals**2)

    if return_X:
        return rss, fitted, coeffs, term_labels, X
    else:
        return rss, fitted, coeffs, term_labels

# -------------------------
# Utility: numerical Jacobian of residual vector wrt selected period variables
# -------------------------
def numerical_jacobian_residuals(time, mag, alg_poly, b_periods, degrees, period_indices, eps_rel=1e-6):
    """
    Compute Jacobian J (N x k) of residual vector r = mag - X@beta
    w.r.t. the vector of period variables (period_indices list).
    Use finite differences on periods. eps_rel relative step.
    """
    # base residuals and fitted
    rss0, fitted0, coeffs0, labels, X0 = polyfit_fixed_periods(time, mag, alg_poly, b_periods, degrees, True)
    r0 = mag - fitted0
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
        rss_p, fitted_p, _, _, Xp = polyfit_fixed_periods(time, mag, alg_poly, periods_pert, degrees, return_X=True)
        rp = mag - fitted_p
        # finite difference column
        J[:, col] = (rp - r0) / h

    return J, r0, X0, coeffs0, labels

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

    initial_parameters = list-like [deg_trend, p1, deg1, p2, deg2, p3, deg3]
    Only periods with period>0 and deg>0 are included for optimization.
    """

    period_indices = []    
    # Periods to optimize
    for i in range(len(optimize_flags)):
        if optimize_flags[i]:
            period_indices.append(i)

    if len(period_indices) == 0:
        # nothing to optimize -> just return linear fit results + cov
        rss, fitted, coeffs, term_labels, X = polyfit_fixed_periods(time, mag, 
                                                                    alg_poly, 
                                                                    initial_periods, 
                                                                    degrees, 
                                                                    True)
        N = len(mag)
        p = len(coeffs)
        dof = max(N - p, 1)
        sigma2 = rss / dof
        cov_coeffs = sigma2 * np.linalg.inv(X.T @ X) if p > 0 else np.zeros((0,0))
        se_coeffs = np.sqrt(np.diag(cov_coeffs)) if p > 0 else np.array([])
        return {
            'best_periods': initial_periods,
            'rss': rss, 
            'fitted': fitted, 
            'coeffs': coeffs, 
            'term_labels': term_labels,
            'cov_coeffs': cov_coeffs, 
            'se_coeffs': se_coeffs,
            'period_cov': None, 
            'se_periods': None
        }

    # initial guess vector for the optimizer
    x0 = np.array([initial_periods[i] for i in period_indices], dtype=float)

    def objective(period_values):
        # update param vector
        periods = initial_periods.copy()
        for idx, val in zip(period_indices, period_values):
            periods[idx] = float(val)
        rss, _, _, _ = polyfit_fixed_periods(time, mag, alg_poly, periods, degrees, False)
        return rss

    # run optimizer
    result = minimize(objective, x0, method=method,
                      options={'maxiter': maxiter, 'xatol': xtol, 'fatol': ftol, 'disp': True})

    best_period_values = result.x
    # assemble full parameter vector with best periods plugged
    best_periods = initial_periods.copy()
    for idx, val in zip(period_indices, best_period_values):
        best_periods[idx] = float(val)

    # final linear fit at best params (and get X, coeffs)
    rss, fitted, coeffs, term_labels, X = polyfit_fixed_periods(time, mag, 
                                                                alg_poly,
                                                                best_periods,
                                                                degrees,
                                                                True)
    N = len(mag)
    p_lin = len(coeffs)
    p_tot = p_lin + len(period_indices)
    dof = max(N - p_tot, 1)   # conservative dof include nonlinear params
    sigma2 = rss / dof

    # covariance of linear coefficients: sigma2 * (X^T X)^{-1}
    if p_lin > 0:
        XtX = X.T @ X
        cov_coeffs = sigma2 * np.linalg.inv(XtX)
        se_coeffs = np.sqrt(np.diag(cov_coeffs))
    else:
        cov_coeffs = np.zeros((0,0))
        se_coeffs = np.array([])

    # Now compute Jacobian J of residuals wrt periods (finite differences)
    J, r0, X0, coeffs0, labels = numerical_jacobian_residuals(time, mag, alg_poly, best_periods, degrees, period_indices)

    # Gauss-Newton approximation: Hessian ≈ 2 * J^T J, but covariance for periods:
    # Cov(periods) ≈ sigma2 * (J^T J)^{-1}
    JTJ = J.T @ J
    # handle possible singularity
    try:
        cov_periods = sigma2 * np.linalg.inv(JTJ)
        se_periods = np.sqrt(np.diag(cov_periods))
    except np.linalg.LinAlgError:
        cov_periods = None
        se_periods = None
        print("Warning: J^T J is singular — cannot compute period covariance via Gauss-Newton.")

    # Print results
    print("\n=== OPTIMIZED PERIODS ===")
    for local_i, global_idx in enumerate(period_indices):
        period_name = f"period_index_{global_idx}"  # maps to parameter slot 1,3,5
        pval = best_periods[global_idx]
        if se_periods is not None:
            print(f"{period_name} = {pval:.8f} ± {se_periods[local_i]:.8f} (1σ)")
        else:
            print(f"{period_name} = {pval:.8f} (no SE)")
    print("=========================\n")

    # Print linear coeffs with errors
    print("=== LINEAR COEFFICIENTS ===")
    for i, (lbl, c) in enumerate(zip(term_labels, coeffs)):
        se = se_coeffs[i] if i < len(se_coeffs) else np.nan
        print(f"{lbl:12s}: {c:.8e} ± {se:.3e}")
    print("===========================\n")
    print(f"RSS = {rss:.6e}, sigma2 (residual variance) = {sigma2:.6e}, dof = {dof}\n")

    output = {
        'best_periods': best_periods,
        'rss': rss,
        'fitted': fitted,
        'coeffs': coeffs,
        'term_labels': term_labels,
        'cov_coeffs': cov_coeffs,
        'se_coeffs': se_coeffs,
        'cov_periods': cov_periods,
        'se_periods': se_periods,
        'period_indices': period_indices,
        'J': J,
        'result': result
    }

    # Optional: bootstrap for period uncertainties (more robust)
    if compute_bootstrap:
        print("Running bootstrap (this may take a while)...")
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
                rss_b, _, _, _ = polyfit_fixed_periods(time, mag_b, alg_poly, periods_b, degrees, False)
                return rss_b
            res_b = minimize(obj_bpv, best_period_values, method=method,
                             options={'maxiter': maxiter, 'xatol': xtol, 'fatol': ftol, 'disp': False})
            boot_periods.append(res_b.x)
        boot_periods = np.array(boot_periods)
        # compute bootstrap std dev for each optimized period
        boot_se = np.std(boot_periods, axis=0, ddof=1)
        print("Bootstrap period 1..k standard errors:", boot_se)
        output['bootstrap_period_se'] = boot_se
        output['bootstrap_periods'] = boot_periods

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
    
    out = optimize_periods_with_errors(time, mag, alg_poly, normalized_periods, normalized_degrees, normalized_optimize_flags, **kwargs)
    fit_result = pd.DataFrame({
        'Time': time,
        'Mag': mag,
        'Fit': out['fitted'],
    })
    return fit_result
