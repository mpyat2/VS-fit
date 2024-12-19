import numpy as np
import statsmodels.api as sm
#from statsmodels.stats.outliers_influence import summary_table
import pandas as pd

def polyfit(t, mag, parameters):
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
    #st, data, ss2 = summary_table(mag_fit, alpha=0.05)
    #fittedvalues = data[:, 2]
    #predict_mean_se  = data[:, 3]
    #predict_mean_ci_low, predict_mean_ci_upp = data[:, 4:6].T
    #predict_ci_low, predict_ci_upp = data[:, 6:8].T
    #fit_result = pd.DataFrame({'Time': t, 
    #                           'Mag': mag, 
    #                           'Fit': fittedvalues,
    #                           'Fit_L': predict_mean_ci_low,
    #                           'Fit_U': predict_mean_ci_upp
    #                           })

    pred_ols = mag_fit.get_prediction()    
    fit_result = pd.DataFrame({'Time': t, 
                               'Mag': mag, 
                               'Fit': mag_fit.fittedvalues,
#                               'Fit_L': pred_ols.summary_frame()["obs_ci_lower"],
#                               'Fit_U': pred_ols.summary_frame()["obs_ci_upper"]
                              })
    return fit_result
