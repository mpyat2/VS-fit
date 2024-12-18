import numpy as np
import statsmodels.api as sm
import pandas as pd

def polyfit(t, mag, parameters):
    #ndata = len(t)

    fit_time = []
    fit_mag = []

    min_t = min(t)
    max_t = max(t)
    #times = np.linspace(min(t), max(t), 1001)
    
    fit_arguments = []
    # algebraic polynomial
    for i in range(parameters[0] + 1):
        fit_arguments.append(t**i)
    # 1-st trigonometric polynomial
    if parameters[1] > 0.0 and parameters[2] > 0:
        for i in range(parameters[2] + 1):
            a = 2 * np.pi * i / parameters[1] * t
            fit_arguments.append(np.sin(a))
            fit_arguments.append(np.cos(a))
    # 2-nd trigonometric polynomial
    if parameters[3] > 0.0 and parameters[4] > 0:
        for i in range(parameters[4] + 1):
            a = 2 * np.pi * i / parameters[3] * t
            fit_arguments.append(np.sin(a))
            fit_arguments.append(np.cos(a))
    # 3-rd trigonometric polynomial
    if parameters[5] > 0.0 and parameters[6] > 0:
        for i in range(parameters[6] + 1):
            a = 2 * np.pi * i / parameters[5] * t
            fit_arguments.append(np.sin(a))
            fit_arguments.append(np.cos(a))

    #print(fit_arguments)
    #print(len(fit_arguments))
    mag_fit = sm.OLS(mag, np.column_stack(fit_arguments)).fit()
    #print(mag_fit)
    pred_ols = mag_fit.get_prediction()
    fit_result = pd.DataFrame({'Time': t, 
                               'Mag': mag, 
                               'Fit': mag_fit.fittedvalues,
                               'Fit_L': pred_ols.summary_frame()["obs_ci_lower"],
                               'Fit_U': pred_ols.summary_frame()["obs_ci_upper"]})
    return fit_result
