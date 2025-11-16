from tkinter import Toplevel, Frame, Label, Entry, Button, Checkbutton, StringVar, IntVar, messagebox

param_defined = False
param_algeDegree = 0
param_trig1Period = 0.0
param_trig1Degree = 0
param_trig1Optimize = 0
param_trig2Period = 0.0
param_trig2Degree = 0
param_trig2Optimize = 0
param_trig3Period = 0.0
param_trig3Degree = 0
param_trig3Optimize = 0
param_bootstrapForErrors = 0

def paramDialogDestroy(dialog):
    global param_defined
    param_defined = False
    dialog.destroy()

def paramCheck(dialog, stringVars, optimize, bootstrapErr):
    global param_algeDegree
    global param_trig1Period
    global param_trig1Degree
    global param_trig1Optimize
    global param_trig2Period
    global param_trig2Degree
    global param_trig2Optimize
    global param_trig3Period
    global param_trig3Degree
    global param_trig3Optimize
    global param_bootstrapForErrors
    global param_defined
    param_defined = False
    param_algeDegree = None
    param_trig1Period = None
    param_trig1Degree = None
    param_trig1Optimize = 0
    param_trig2Period = None
    param_trig2Degree = None
    param_trig2Optimize = 0
    param_trig3Period = None
    param_trig3Degree = None
    param_trig3Optimize = 0
    param_bootstrapForErrors = 0
    
    algeDegree_s = stringVars[0].get()
    trig1Period_s = stringVars[1].get()
    trig1Degree_s = stringVars[2].get()
    trig2Period_s = stringVars[3].get()
    trig2Degree_s = stringVars[4].get()
    trig3Period_s = stringVars[5].get()
    trig3Degree_s = stringVars[6].get()
    trig1Optimize = optimize[0].get()
    trig2Optimize = optimize[1].get()
    trig3Optimize = optimize[2].get()
    bootstrapForErrors = bootstrapErr.get()
    try:
        aD = int(eval(algeDegree_s, {}, {}))
        if aD < 0:
            raise Exception('Error', 'Algebraic polynomial degree must be >= 0')
        t1P = float(eval(trig1Period_s, {}, {}))
        if t1P < 0:
            raise Exception('Error', 'Trigonometric polynomial 1 period must be >= 0')
        t1D = int(eval(trig1Degree_s, {}, {}))
        if t1D < 0 or t1D > 10:
            raise Exception('Error', 'Trigonometric polynomial 1 degree must be between 0 and 10')
        t2P = float(eval(trig2Period_s, {}, {}))
        if t2P < 0:
            raise Exception('Error', 'Trigonometric polynomial 2 period must be >= 0')
        t2D = int(eval(trig2Degree_s, {}, {}))
        if t2D < 0 or t2D > 10:
            raise Exception('Error', 'Trigonometric polynomial 2 degree must be between 0 and 10')
        t3P = float(eval(trig3Period_s, {}, {}))
        if t3P < 0:
            raise Exception('Error', 'Trigonometric polynomial 3 period must be >= 0')
        t3D = int(eval(trig3Degree_s, {}, {}))
        if t3D < 0 or t3D > 10:
            raise Exception('Error', 'Trigonometric polynomial 3 degree must be between 0 and 10')
        param_algeDegree = aD
        param_trig1Period = t1P
        param_trig1Degree = t1D
        param_trig2Period = t2P
        param_trig2Degree = t2D
        param_trig3Period = t3P
        param_trig3Degree = t3D
        param_trig1Optimize = trig1Optimize
        param_trig2Optimize = trig2Optimize
        param_trig3Optimize = trig3Optimize
        param_bootstrapForErrors = bootstrapForErrors
        param_defined = True
        dialog.destroy()
        return
    except Exception as e:
        messagebox.showinfo(None, repr(e), parent=dialog)
        return

def createEntry(frame, labelText, textvariable, initValue, r1, c1, r2, c2):
    Label(frame, text=labelText).grid(row=r1, column=c1)
    entry = Entry(frame, textvariable=textvariable)
    entry.insert(0, str(initValue)) 
    entry.grid(row=r2, column=c2)
    return entry

def createCheckBox(frame, intvariable, initValue, r1, c1):
    checkbox = Checkbutton(frame, variable=intvariable)
    checkbox.deselect() if initValue == 0 else checkbox.select()
    checkbox.grid(row=r1, column=c1)
    return checkbox

def fitParameters(master):
    global param_defined
    global param_algeDegree
    global param_trig1Period
    global param_trig1Degree
    global param_trig1Optimize
    global param_trig2Period
    global param_trig2Degree
    global param_trig2Optimize
    global param_trig3Period
    global param_trig3Degree
    global param_trig3Optimize
    global param_bootstrapForErrors

    param_defined = False

    x = master.winfo_x()
    y = master.winfo_y()
    #print(x)
    #print(y)

    paramDialog = Toplevel(master)
    paramDialog.protocol("WM_DELETE_WINDOW", lambda: paramDialogDestroy(paramDialog))
    paramDialog.title("Parameters")
    paramDialog.geometry("800x200+" + str(x+20) + "+" + str(y+20))

    frame = Frame(paramDialog)
    frame.pack(pady=10)

    algeDegree = StringVar(master=paramDialog)
    trig1Period = StringVar(master=paramDialog)
    trig1Degree = StringVar(master=paramDialog)
    trig1Optimize = IntVar(master=paramDialog)
    trig2Period = StringVar(master=paramDialog)
    trig2Degree = StringVar(master=paramDialog)
    trig2Optimize = IntVar(master=paramDialog)
    trig3Period = StringVar(master=paramDialog)
    trig3Degree = StringVar(master=paramDialog)
    trig3Optimize = IntVar(master=paramDialog)
    bootstrapErr = IntVar(master=paramDialog)

    createEntry(frame, "Polynomial Degree", algeDegree, param_algeDegree, 0, 1, 0, 2)
    createEntry(frame, "Trig. Polyn. Period, Degree, Optimize 1", trig1Period, param_trig1Period, 1, 1, 1, 2)
    createEntry(frame, "", trig1Degree, param_trig1Degree, 1, 3, 1, 4)
    createCheckBox(frame, trig1Optimize, param_trig1Optimize, 1, 5)
    createEntry(frame, "Trig. Polyn. Period, Degree, Optimize 2", trig2Period, param_trig2Period, 2, 1, 2, 2)
    createEntry(frame, "", trig2Degree, param_trig2Degree, 2, 3, 2, 4)
    createCheckBox(frame, trig2Optimize, param_trig2Optimize, 2, 5)
    createEntry(frame, "Trig. Polyn. Period, Degree, Optimize 3", trig3Period, param_trig3Period, 3, 1, 3, 2)
    createEntry(frame, "", trig3Degree, param_trig3Degree, 3, 3, 3, 4)
    createCheckBox(frame, trig3Optimize, param_trig3Optimize, 3, 5)
    Label(frame, text="Calculate period errors via bootstrap").grid(row=4, column=4)
    createCheckBox(frame, bootstrapErr, param_bootstrapForErrors, 4, 5)

    buttonOK = Button(frame, text="OK", 
                      command=lambda: paramCheck(paramDialog, 
                                                [algeDegree, 
                                                 trig1Period, trig1Degree,
                                                 trig2Period, trig2Degree,
                                                 trig3Period, trig3Degree],
                                                [trig1Optimize, 
                                                 trig2Optimize, 
                                                 trig3Optimize],
                                                 bootstrapErr))
    buttonOK.grid(row=5, column=1)
    buttonCancel = Button(frame, text="Cancel", command=lambda: paramDialogDestroy(paramDialog))
    buttonCancel.grid(row=5, column=4)
    paramDialog.transient(master)
    paramDialog.grab_set()
    paramDialog.focus_set()
    paramDialog.wait_window()
