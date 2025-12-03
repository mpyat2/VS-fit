from tkinter import Toplevel, Frame, Label, Entry, Button, StringVar, messagebox

param_defined = False
param_lofreq = 0.01
param_hifreq = 1
param_n_intervals = 1000
params_initialized = False

def paramDialogDestroy(dialog):
    global param_defined
    param_defined = False
    #print("paramDialogDestroy: param_defined: ", param_defined)
    dialog.destroy()

def paramCheck(dialog, lofreq, hifreq, n_intervals):
    global param_lofreq
    global param_hifreq
    global param_n_intervals
    global param_defined
    param_defined = False
    param_lofreq = None
    param_hifreq = None
    param_n_intervals = None
    lof_s = lofreq.get()
    hif_s = hifreq.get()
    n_intervals_s = n_intervals.get()
    #low = 0.0
    #hi = 0.0
    try:
        low = float(eval(lof_s, {}, {}))
        if low <= 0:
            raise Exception('Error', 'Low frequency must be > 0')
        hi = float(eval(hif_s, {}, {}))
        if hi <= 0:
            raise Exception('Error', 'High frequency must be > 0')
        if hi <= low:
            raise Exception('Error', 'High frequency must be > Low frequency')
        n = int(eval(n_intervals_s, {}, {}))
        if n < 1:
            raise Exception('Error', 'Number of intervals must be > 0')
        if n > 500000:
            raise Exception('Error', 'Number of intervals must be < 500000')
        #print(low)
        #print(hi)
        #print(n)
        param_lofreq = low
        param_hifreq = hi
        param_n_intervals = n
        param_defined = True
        dialog.destroy()
        return
    except Exception as e:
        messagebox.showinfo(None, repr(e), parent=dialog)
        return

def dftParameters(master):
    global param_defined
    global param_lofreq
    global param_hifreq
    global param_n_intervals
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
    entry_label_n = Entry(frame, textvariable = n_intervals)
    entry_label_n.insert(0, str(param_n_intervals)) 
    entry_label_n.grid(row=2, column=2)
    buttonOK = Button(frame, text="OK", command=lambda: paramCheck(paramDialog, lofreq, hifreq, n_intervals))
    buttonOK.grid(row=3, column=1)
    buttonCancel = Button(frame, text="Cancel", command=lambda: paramDialogDestroy(paramDialog))
    buttonCancel.grid(row=3, column=2)
    paramDialog.transient(master)
    paramDialog.grab_set()
    paramDialog.focus_set()
    paramDialog.wait_window()
