import numpy as np
import statsmodels.api as sm
import pandas as pd

# Translated from R (Rprogs.r) (with a tiny fix)
# See https://www.aavso.org/software-directory, https://www.aavso.org/sites/default/files/software/Rcodes.zip
def dcdft(t, mag, lowfreq, hifreq, n_freq, mcv_mode):
    ndata = len(t)
    mag_var = np.var(mag)

    freq = []
    per = []
    power = []
    amp = []
    
    frequencies = np.linspace(lowfreq, hifreq, n_freq + 1)
    
    for nu in frequencies:
        freq.append(nu)
        per.append(1 / nu)
        # Calculate cos and sin (c1 and s1) of the time-sequence for the trial frequency
        a = 2 * np.pi * nu * t
        c1 = np.cos(a)
        s1 = np.sin(a)
        # Fit c1 and s1 to the signal (mag)
        mag_fit = sm.OLS(mag, sm.add_constant(np.column_stack((c1, s1)))).fit()
        # Get the amplitude of the trial frequency
        if len(mag_fit.params) < 3:
            raise Exception("OLS", "Error")
        a1 = mag_fit.params[1]
        a2 = mag_fit.params[2]
        amp_squared = a1 * a1 + a2 * a2
        amp.append(np.sqrt(amp_squared))
        # Get the power of the trial frequency
        pwr = np.var(mag_fit.fittedvalues)
        power.append(pwr)
    
    if mcv_mode:
        power = np.array(power) / mag_var
    else:
        power = np.array(power) * (ndata - 1) / mag_var / 2
        
    dcd = pd.DataFrame({'freq': freq, 'per': per, 'amp': amp, 'pow': power})
    #dcd = dcd.round(8)
    return dcd
