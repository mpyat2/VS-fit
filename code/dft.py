import threading
import queue
import time as pytime
import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle

stop_flag = {"stop": False, "running": False, "time": 0.0}   # mutable for safe thread sharing
result_queue = queue.Queue()  # worker → main communication

def dcdft(master, callback, time, mag, lowfreq, hifreq, n_intervals, n_terms, chunk_size=1000):
    # time normalization to improve numerical stability of LombScargle
    t_span    = time.max() - time.min()

    frequencies = np.linspace(lowfreq, hifreq, n_intervals + 1)
    frequencies = frequencies[frequencies > 0]
    frequencies = frequencies * t_span  # rescale frequencies to normalized time units for LombScargle

    ls = LombScargle((time - time.min()) / t_span, mag, nterms=n_terms)

    freq_chunks = []
    power_chunks = []

    for i in range(0, len(frequencies), chunk_size):
        if stop_flag["stop"]:
            return None

        chunk_freqs = frequencies[i:i + chunk_size]
        power_chunks.append(ls.power(chunk_freqs))
        freq_chunks.append(chunk_freqs)

        done = min(i + chunk_size, len(frequencies))
        callback(master, None, f"{done} of {len(frequencies)} frequencies computed.", "progress")

    all_freq = np.concatenate(freq_chunks) / t_span  # rescale frequencies back to original units
    all_power = np.concatenate(power_chunks)

    dcd = pd.DataFrame({'freq': all_freq, 'per': 1.0 / all_freq, 'pow': all_power,})
    return dcd

def worker(master, callback, time, mag, lowfreq, hifreq, n_intervals, n_terms):
    try:
        result = dcdft(
            master=master,
            callback=callback,
            time=time,
            mag=mag,
            lowfreq=lowfreq,
            hifreq=hifreq,
            n_intervals=n_intervals,
            n_terms=n_terms
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
            callback(master, None, "DFT was stopped.", "stopped")
    else:
        msg = f"DFT calculation time {(pytime.time() - stop_flag['time']):.2f} s"
        callback(master, dcdft_result["data"], msg, "finished")

def stop_task():
    stop_flag["stop"] = True

def dcdft_async(master, callback, time, mag, lofreq, hifreq, n_intervals, n_terms):
    stop_flag["stop"] = False  # reset stop flag
    threading.Thread(
        target=worker,
        args=(master, callback, time, mag, lofreq, hifreq, n_intervals, n_terms),
        daemon=True
    ).start()
    stop_flag["time"] = pytime.time()    
    callback(master, None, None, "started")
    stop_flag["running"] = True
    master.after(100, lambda: check_worker_result(master, callback))
