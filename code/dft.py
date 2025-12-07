import threading
import queue
import time as pytime
import numpy as np
import pandas as pd

stop_flag = {"stop": False, "running": False, "time": 0.0}   # mutable for safe thread sharing
result_queue = queue.Queue()  # worker → main communication

# Translated from R (Rprogs.r) (with a tiny fix)
# See https://www.aavso.org/software-directory, https://www.aavso.org/sites/default/files/software/Rcodes.zip
def dcdft(master, callback, time, mag, lowfreq, hifreq, n_intervals, mcv_mode=False):
    ndata = len(time)
    t_mean = np.mean(time)
    mag_var = np.var(mag)
    
    t = time - t_mean

    freq = []
    per = []
    power = []
    amp = []
    
    freq_step = (hifreq - lowfreq) / n_intervals
    for i in range(n_intervals + 1):
        if stop_flag["stop"]:
            return None
        nu = lowfreq + i * freq_step
        if nu <= 0:
            continue
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
        if (i + 1) % 1000 == 0:
            #print(f"{i + 1} of {n_intervals + 1} frequencies computed.")
            callback(master, None, f"{i + 1} of {n_intervals + 1} frequencies computed.", "progress")
        
    #print(f"Finished: {n_intervals + 1} frequencies computed.")
    
    if mcv_mode:
        power = np.array(power) / mag_var
    else:
        power = np.array(power) * (ndata - 1) / mag_var / 2
        
    dcd = pd.DataFrame({'freq': freq, 'per': per, 'amp': amp, 'pow': power})
    return dcd

def worker(master, callback, time, mag, lowfreq, hifreq, n_intervals):
    try:
        result = dcdft(
            master=master,
            callback=callback,
            time=time,
            mag=mag,
            lowfreq=lowfreq,
            hifreq=hifreq,
            n_intervals=n_intervals,
        )
        dcdft_result = {"data": result, "error": None}
    except Exception as e:
        dcdft_result = {"data": None, "error": str(e)}
    result_queue.put(dcdft_result)

def check_worker_result(master, callback):
    try:
        dcdft_result = result_queue.get_nowait()
    except queue.Empty:
        # Worker still running — keep polling
        master.after(100, lambda: check_worker_result(master, callback))
        return
    # Worker finished:
    stop_flag["running"] = False
    if dcdft_result is None:
        callback(master, None, "Unknown error", "error")
        return
    if dcdft_result["data"] is None:
        if dcdft_result["error"] is not None:
            callback(master, None, "Error: " + dcdft_result["error"], "error")
        else:
            callback(master, None, "DC DFT was stopped.", "stopped")
    else:
        msg = f"DC DFT calculation time {(pytime.time() - stop_flag['time']):.2f} s"
        callback(master, dcdft_result["data"], msg, "finished")

def stop_task():
    stop_flag["stop"] = True

def dcdft_async(master, callback, time, mag, lofreq, hifreq, n_intervals):
    stop_flag["stop"] = False  # reset stop flag
    threading.Thread(
        target=worker,
        args=(master, callback, time, mag, lofreq, hifreq, n_intervals),
        daemon=True
    ).start()
    stop_flag["time"] = pytime.time()    
    callback(master, None, None, "started")
    stop_flag["running"] = True
    master.after(100, lambda: check_worker_result(master, callback))
