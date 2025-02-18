import numpy as np
import statsmodels.api as sm
#from statsmodels.stats.outliers_influence import summary_table
import pandas as pd

def polyfit(time, mag, parameters):
    t_mean = np.mean(time)
    temp_t = time - t_mean
    fit_arguments = []
    # algebraic polynomial
    for i in range(parameters[0] + 1):
        fit_arguments.append(temp_t**i)
        print("Trend degree ", i, " added")
    # 1-st trigonometric polynomial
    if parameters[1] > 0.0 and parameters[2] > 0:
        for i in range(1, parameters[2] + 1):
            a = 2 * np.pi * i / parameters[1] * temp_t
            fit_arguments.append(np.sin(a))
            print("Poly degree1 sin ", i, " added")
            fit_arguments.append(np.cos(a))
            print("Poly degree1 cos ", i, " added")
    # 2-nd trigonometric polynomial
    if parameters[3] > 0.0 and parameters[4] > 0:
        for i in range(1, parameters[4] + 1):
            a = 2 * np.pi * i / parameters[3] * temp_t
            fit_arguments.append(np.sin(a))
            print("Poly degree2 sin ", i, " added")
            fit_arguments.append(np.cos(a))
            print("Poly degree2 cos ", i, " added")
    # 3-rd trigonometric polynomial
    if parameters[5] > 0.0 and parameters[6] > 0:
        for i in range(1, parameters[6] + 1):
            a = 2 * np.pi * i / parameters[5] * temp_t
            fit_arguments.append(np.sin(a))
            print("Poly degree3 sin ", i, " added")
            fit_arguments.append(np.cos(a))
            print("Poly degree3 cos ", i, " added")

    #print(fit_arguments)
    print(len(fit_arguments))
    mag_fit = sm.OLS(mag, np.column_stack(fit_arguments)).fit()

    pred_ols = mag_fit.get_prediction()    
    fit_result = pd.DataFrame({'Time': time, 
                               'Mag': mag, 
                               'Fit': mag_fit.fittedvalues,
#                               'Fit_L': pred_ols.summary_frame()["obs_ci_lower"],
#                               'Fit_U': pred_ols.summary_frame()["obs_ci_upper"]
                              })
    return fit_result
