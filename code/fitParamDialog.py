from tkinter import Toplevel, Frame, Label, Entry, Button, Checkbutton, StringVar, IntVar, messagebox

# Max number of periods to fit
MAX_PERIODS = 9

param_defined = False
param_algDegree = 0
param_periods = []
param_degrees = []
param_optFlags = []
param_bootstrapForErrors = 0

def paramDialogDestroy(dialog):
    global param_defined
    param_defined = False
    dialog.destroy()

def paramCheck(dialog,
               strVarAlgDeg,
               strVarPeriods,
               strVarDegrees,               
               optFlags,
               bootstrapErr):
    global param_defined    
    global param_algDegree
    global param_periods
    global param_degrees
    global param_optFlags
    global param_bootstrapForErrors
    
    param_defined = False
    param_algDegree = 0
    param_periods = []
    param_degrees = []
    param_optFlags = []
    param_bootstrapForErrors = 0
    
    n_of_periods = len(strVarPeriods)
    
    try:
        aD = int(eval(strVarAlgDeg.get().strip() or "0", {}, {}))
        if aD < 0:
            raise Exception('Error', 'Algebraic polynomial degree must be >= 0')
        param_algDegree = aD
            
        for i in range(n_of_periods):
            tD = int(eval(strVarDegrees[i].get().strip() or "0", {}, {}))
            if tD < 0 or tD > 9:
                raise Exception('Error', f'Trigonometric polynomial {i+1} degree must be between 0 and 9')
            tP = float(eval(strVarPeriods[i].get().strip() or "0", {}, {}))
            if tP < 0:
                raise Exception('Error', f'Trigonometric polynomial {i+1} period must be >= 0')
            tO = optFlags[i].get()
            param_periods.append(tP)
            param_degrees.append(tD)
            param_optFlags.append(tO)
        
        param_bootstrapForErrors = bootstrapErr.get()
        
        param_defined = True
        dialog.destroy()
        return
    except Exception as e:
        messagebox.showinfo(None, repr(e), parent=dialog)
        return

def createEntry(frame, labelText, textvariable, initValue, r1, c1, r2, c2):
    if labelText:
        Label(frame, text=labelText).grid(row=r1, column=c1)
    entry = Entry(frame, textvariable=textvariable)
    entry.insert(0, str(initValue)) 
    entry.grid(row=r2, column=c2)
    return entry

def createCheckBox(frame, labelText, intvariable, initValue, r1, c1, r2, c2):
    if labelText:
        Label(frame, text=labelText).grid(row=r1, column=c1)
    checkbox = Checkbutton(frame, variable=intvariable)
    checkbox.deselect() if initValue == 0 else checkbox.select()
    checkbox.grid(row=r2, column=c2)
    return checkbox

def fitParameters(master):
    global param_defined    
    global param_algDegree
    global param_periods
    global param_degrees
    global param_optFlags
    global param_bootstrapForErrors

    param_defined = False

    x = master.winfo_x()
    y = master.winfo_y()

    paramDialog = Toplevel(master)
    paramDialog.protocol("WM_DELETE_WINDOW", lambda: paramDialogDestroy(paramDialog))
    paramDialog.title("Parameters")
    paramDialog.geometry("+" + str(x+20) + "+" + str(y+20))

    frame = Frame(paramDialog)
    frame.pack(padx=10, pady=10)

    intVarBootstrapErr = IntVar(master=paramDialog)
    strVarAlgDeg = StringVar(master=paramDialog)
    strVarPeriods = []
    strVarDegrees = []
    intVarOptFlags = []
    for i in range(MAX_PERIODS):
        strVarPeriods.append(StringVar(master=paramDialog))
        strVarDegrees.append(StringVar(master=paramDialog))
        intVarOptFlags.append(IntVar(master=paramDialog))
        
    createEntry(frame, "Polynomial Degree", strVarAlgDeg, param_algDegree, 0, 1, 0, 2)
    
    for i in range(MAX_PERIODS):
        if i < len(param_periods):
            p = param_periods[i]
            d = param_degrees[i]
            o = param_optFlags[i]
        else:
            p = 0
            d = 0
            o = 0
        createEntry(frame, f"Trig. Polyn. Period {i+1}", strVarPeriods[i], p, i + 1, 1, i + 1, 2)
        createEntry(frame, "Degree", strVarDegrees[i], d, i + 1, 3, i + 1, 4)
        createCheckBox(frame, "Optimize", intVarOptFlags[i], o, i + 1, 5, i + 1, 6)
    
    label = Label(frame, text="Calculate period errors via bootstrap")
    label.grid(row=MAX_PERIODS + 1, column=4, columnspan=2)
    createCheckBox(frame, "",
                   intVarBootstrapErr, param_bootstrapForErrors, MAX_PERIODS + 1, 5, MAX_PERIODS + 1, 6)
    Label(frame, text="May take a while!").grid(row=MAX_PERIODS + 1, column=7)

    buttonOK = Button(frame, text="OK", 
                      command=lambda: paramCheck(paramDialog,
                                                 strVarAlgDeg,
                                                 strVarPeriods,
                                                 strVarDegrees,
                                                 intVarOptFlags,
                                                 intVarBootstrapErr))
    buttonOK.grid(row=MAX_PERIODS + 2, column=1)
    buttonCancel = Button(frame, text="Cancel", command=lambda: paramDialogDestroy(paramDialog))
    buttonCancel.grid(row=MAX_PERIODS + 2, column=4)
    paramDialog.transient(master)
    paramDialog.grab_set()
    paramDialog.focus_set()
    paramDialog.wait_window()
