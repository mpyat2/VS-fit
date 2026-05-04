from tkinter import Toplevel, Frame, Label, Entry, Button, StringVar, messagebox
import numpy as np
from utils import safe_eval

MAX_HARMONICS = 25

param_defined = False
param_lofreq = 0.0
param_hifreq = 1.0
param_n_intervals = 1000
param_n_terms = 1
params_initialized = False

def paramDialogDestroy(dialog):
    global param_defined
    param_defined = False
    #print("paramDialogDestroy: param_defined: ", param_defined)
    dialog.destroy()

def paramCheck(dialog, lofreq, hifreq, n_intervals, n_terms):
    global param_lofreq
    global param_hifreq
    global param_n_intervals
    global param_n_terms
    global param_defined
    param_defined = False
    param_lofreq = None
    param_hifreq = None
    param_n_intervals = None
    n_terms_s = n_terms.get()
    n_intervals_s = n_intervals.get()
    try:
        low, hi = checkLoHi(lofreq, hifreq)
        n = int(safe_eval(n_intervals_s))
        if n < 1:
            raise Exception('Error', 'Number of intervals must be > 0')
        harmonics = int(safe_eval(n_terms_s))
        if harmonics < 1 or harmonics > MAX_HARMONICS:
            raise Exception('Error', f'Number of harmonics must be between 1 and {MAX_HARMONICS}')
        param_lofreq = low
        param_hifreq = hi
        param_n_intervals = n
        param_n_terms = harmonics
        param_defined = True
        dialog.destroy()
        return
    except Exception as e:
        messagebox.showinfo(None, repr(e), parent=dialog)
        return

def checkLoHi(lofreq, hifreq):
    lof_s = lofreq.get()
    hif_s = hifreq.get()
    low = float(safe_eval(lof_s))
    if low < 0:
        raise Exception('Error', 'Low frequency must be >= 0')
    hi = float(safe_eval(hif_s))
    if hi <= 0:
        raise Exception('Error', 'High frequency must be > 0')
    if hi <= low:
        raise Exception('Error', 'High frequency must be > Low frequency')
    return low, hi

def recommendedInterval(dialog, lofreq, hifreq, n_intervals, n_terms, t):
    try:
        low, hi = checkLoHi(lofreq, hifreq)
        time_interval = max(t) - min(t)
        if time_interval > 0: # False also for 'nan'
            terms = int(safe_eval(n_terms.get()))
            if terms < 1 or terms > MAX_HARMONICS:
                raise Exception('Error', f'Number of harmonics must be between 1 and {MAX_HARMONICS}')
            recommended_freq_resolution = recommendedResolution(time_interval, terms)
            intervals = round((hi - low) / recommended_freq_resolution) + 1
        else:
            raise Exception('Error', 'Cannot compute recommended number of intervals: invalid time data')
    except Exception as e:
        messagebox.showinfo(None, repr(e), parent=dialog)
        return
    answer = messagebox.askyesno("Info", f"Recommended number of intervals: {intervals}\nUse this number?", parent=dialog)
    if answer:
        n_intervals.set(str(intervals))

def medianInterval(times):
    # Compute successive differences
    dt = np.diff(np.sort(times))
    # Keep only positive (or non-zero) intervals
    dt_nonzero = dt[dt != 0]
    # Compute the median
    return np.median(dt_nonzero)

def recommendedResolution(time_interval, terms):
    return 0.05 / time_interval / terms

def initializeDftParams(t):
    global params_initialized
    global param_lofreq
    global param_hifreq
    global param_n_intervals
    global param_n_terms
    if not params_initialized:
        params_initialized = True
        param_n_terms = 1
        interval = medianInterval(t)
        if interval > 0: # False also for 'nan'
            low = 0.0
            hi = 1.0 / interval / 2.0; # Approximate Nyquist
            if hi > 50.0:
                hi = 50.0
            time_interval = max(t) - min(t)
            if time_interval > 0: # False also for 'nan'
                recommended_freq_resolution = recommendedResolution(time_interval, param_n_terms)
                param_lofreq = low
                param_hifreq = hi
                param_n_intervals = round((hi - low) / recommended_freq_resolution) + 1


def dftParameters(master, t):
    initializeDftParams(t)

    global param_defined
    global param_lofreq
    global param_hifreq
    global param_n_intervals
    global param_n_terms
    global param_defined

    param_defined = False

    x = master.winfo_x()
    y = master.winfo_y()
    #print(x)
    #print(y)

    paramDialog = Toplevel(master)
    paramDialog.protocol("WM_DELETE_WINDOW", lambda: paramDialogDestroy(paramDialog))
    paramDialog.title("Parameters")
    paramDialog.geometry("+" + str(x+20) + "+" + str(y+20))

    frame = Frame(paramDialog)
    frame.pack(padx=10, pady=10)

    lofreq = StringVar(master=paramDialog)
    hifreq = StringVar(master=paramDialog)
    n_intervals = StringVar(master=paramDialog)
    n_terms = StringVar(master=paramDialog)

    label_lofreq = Label(frame, text="Low frequency")
    label_lofreq.grid(row=0, column=1)
    entry_lofreq = Entry(frame, textvariable = lofreq)
    entry_lofreq.insert(0, str(param_lofreq))
    entry_lofreq.grid(row=0, column=2)

    label_hifreq = Label(frame, text="High frequency")
    label_hifreq.grid(row=1, column=1)
    entry_hifreq = Entry(frame, textvariable = hifreq)
    entry_hifreq.insert(0, str(param_hifreq)) 
    entry_hifreq.grid(row=1, column=2)

    label_n_intervals = Label(frame, text="Intervals")
    label_n_intervals.grid(row=2, column=1)
    entry_n_intervals = Entry(frame, textvariable = n_intervals)
    entry_n_intervals.insert(0, str(param_n_intervals)) 
    entry_n_intervals.grid(row=2, column=2)

    buttonInfo = Button(frame, text="?", command=lambda: recommendedInterval(paramDialog, lofreq, hifreq, n_intervals, n_terms, t))
    buttonInfo.grid(row=2, column=3)
    
    label_n_terms = Label(frame, text="Harmonics")
    label_n_terms.grid(row=3, column=1)
    entry_n_terms = Entry(frame, textvariable = n_terms)
    entry_n_terms.insert(0, str(param_n_terms))
    entry_n_terms.grid(row=3, column=2)

    buttonOK = Button(frame, text="OK", command=lambda: paramCheck(paramDialog, lofreq, hifreq, n_intervals, n_terms))
    buttonOK.grid(row=4, column=1)
    buttonCancel = Button(frame, text="Cancel", command=lambda: paramDialogDestroy(paramDialog))
    buttonCancel.grid(row=4, column=2)
    paramDialog.transient(master)
    paramDialog.grab_set()
    paramDialog.focus_set()
    paramDialog.wait_window()
