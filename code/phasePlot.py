import matplotlib.pyplot as plt
from tkinter import Toplevel, Frame, Label, Entry, Button, StringVar, messagebox
import numpy as np

param_period = None
param_epoch = None

def phaseDialogDestroy(dialog):
    #print("phaseDialogDestroy")
    dialog.destroy()

def phaseParamApply(phaseDialog, period, epoch, input_data):
    global param_period
    global param_epoch
    param_period = None
    param_epoch = None
    period_s = period.get()
    epoch_s = epoch.get()

    period_v = 0.0
    epoch_v = 0.0
    try:
        period_v = float(eval(period_s, {}, {}))
        if period_v <= 0:
            raise Exception('Error', 'Period must be > 0')
        epoch_v = float(eval(epoch_s, {}, {}))
        #if epoch_v <= 0:
        #    raise Exception('Error', 'Epoch > 0')
        #print(period_v)
        #print(epoch_v)
        
        param_period = period_v
        param_epoch = epoch_v
        fig = plt.figure(0)
        fig.clear()
        t = input_data['Time'].to_numpy()
        m = input_data['Mag'].to_numpy()
        std_phase = ((t - epoch_v) % period_v) / period_v
        std_phase2 = np.concatenate([std_phase, std_phase - 1.0])
        m2 = np.concatenate([m, m])
        
        plt.plot(std_phase2, m2, '.', color='royalblue')
        plt.ylim(max(input_data['Mag']), min(input_data['Mag']))
        plt.grid(True, linestyle='--', color='gray', alpha=0.3)
        plt.title('Phase Plot')
        plt.show(block=False)
        plt.pause(0.001)  # Forces redraw: needed in Spyder
    except Exception as e:
        messagebox.showinfo(None, repr(e), parent=phaseDialog)

def plotFolded(master, input_data):
    x = master.winfo_x()
    y = master.winfo_y()
    #print(x)
    #print(y)

    phaseDialog = Toplevel(master)
    phaseDialog.protocol("WM_DELETE_WINDOW", lambda: phaseDialogDestroy(phaseDialog))
    phaseDialog.title("Parameters")
    phaseDialog.geometry("240x100+" + str(x+60) + "+" + str(y+60))

    frame = Frame(phaseDialog)
    frame.pack(pady=10)

    period = StringVar()
    epoch = StringVar()

    label_period = Label(frame, text="Period")
    label_period.grid(row=0, column=1)
    entry_period = Entry(frame, textvariable = period)
    entry_period.insert(0, str(param_period) if param_period is not None else "") 
    entry_period.grid(row=0, column=2)

    label_epoch = Label(frame, text="Epoch")
    label_epoch.grid(row=3, column=1)
    entry_epoch = Entry(frame, textvariable = epoch)
    entry_epoch.insert(0, str(param_epoch) if param_epoch is not None else "") 
    entry_epoch.grid(row=3, column=2)
    buttonOK = Button(frame, text="Apply", command=lambda: phaseParamApply(phaseDialog, period, epoch, input_data))
    buttonOK.grid(row=5, column=1)
    buttonCancel = Button(frame, text="Close", command=lambda: phaseDialogDestroy(phaseDialog))
    buttonCancel.grid(row=5, column=2)
    phaseDialog.transient(master)
    phaseDialog.grab_set()
    phaseDialog.focus_set()
    phaseDialog.wait_window()
