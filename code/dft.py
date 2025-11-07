import numpy as np
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
        # Design matrix
        X = np.column_stack((np.ones_like(t), c1, s1))
        # Least-squares solution
        params, residuals, rank, s = np.linalg.lstsq(X, mag, rcond=None)
        
        # Extract coefficients
        a0, a1, a2 = params
        
        amp.append(np.sqrt(a1**2 + a2**2))
        power.append(np.var(X @ params))
    
    if mcv_mode:
        power = np.array(power) / mag_var
    else:
        power = np.array(power) * (ndata - 1) / mag_var / 2
        
    dcd = pd.DataFrame({'freq': freq, 'per': per, 'amp': amp, 'pow': power})
    return dcd
