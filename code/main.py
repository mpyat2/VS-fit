import os
import sys
import webbrowser
import threading
from tkinter import Tk, Frame, Label, Menu, Button, PhotoImage, filedialog, messagebox
from log_window import LogWindow
##
import matplotlib
matplotlib.use("TkAgg")
import plotWind
##
import dataio
import dft
import dftParamDialog
import fit
import fitParamDialog
import phasePlot

log_window = None

plotWind0 = None # Input data and approximation
plotWind1 = None # DFT result

input_data = None
dft_result = None
fit_result = None

def clear_log(master):
    try:
        log_window.clear()
    except Exception as e:
        messagebox.showinfo(None, "Error: " + str(e), parent=master)

def add_to_log(master, text):
    try:
        print(text)
        log_window.add_line(text)
    except Exception as e:
        messagebox.showinfo(None, "Error: " + str(e), parent=master)

def checkBackgroundTaskRunning(master):
    if dft.stop_flag["running"] or fit.stop_flag["running"]:
        messagebox.showinfo("Please wait", "Background task is still running. Please wait until it is finished.", parent=master)
        return True
    return False

def shutdown(master):
    # silent check
    if dft.stop_flag["running"] or fit.stop_flag["running"]:
        return
    try:
        master.quit() # stop mainloop if running: required in Spyder
    except:
        pass    
    master.destroy()

def doShutdown(master):
    if checkBackgroundTaskRunning(master): return
    shutdown(master)

def waitOverlay(master, name="waitOverlay", stop_call=None):
    overlay = Frame(master, name=name, bg="#cccccc")
    overlay.place(relx=0, rely=0, relwidth=1, relheight=1)
    lbl = None
    if stop_call is not None:
        btn = Button(overlay, text="Stop", command=lambda: stop_call())
        btn.place(relx=0.5, rely=0.0, anchor="n")
        lbl = Label(overlay, text="Processing...", bg="#cccccc")
        lbl.place(relx=0.5, rely=1.0, anchor="s")
    else:
        Label(overlay, text="Please wait...", bg="#c0c0c0").place(relx=0.5, rely=0.5, anchor="center")
    return overlay, lbl

def plotData(master):
    def plot_input(ax):
        ax.plot(input_data['Time'], input_data['Mag'], '.', color='royalblue')
        ax.set_ylim(max(input_data["Mag"]), min(input_data["Mag"]))
        ax.set_title('Input')
        ax.set_xlabel('Time')
        ax.set_ylabel('Magnitude')
        ax.grid(True, linestyle='--', color='gray', alpha=0.3)
        if fit_result is not None:
            ax.plot(fit_result['Time'], fit_result['Fit'], 'k.')
    global plotWind0
    if plotWind0 is None: 
        plotWind0 =  plotWind.PlotWindow(master, title="Input Data")
    plotWind0.show(plot_input)

def plotDftResult(master, plot_power, plot_frequency):
    def plot_dft(ax):
        global dft_result
        x_label = 'Period'
        y_label = 'Semi-amplitude'
        x_col = 'per'
        y_col = 'amp'
        if plot_power:
            y_col = 'pow'
            y_label = 'Power'
        if plot_frequency:
            x_col = 'freq'
            x_label = 'Frequency'
        ax.plot(dft_result[x_col], dft_result[y_col], color = 'k', linestyle='-')
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title('DC DFT')
        ax.grid(True, linestyle='--', color='gray', alpha=0.3)
    global plotWind1
    if plotWind1 is None: 
        plotWind1 =  plotWind.PlotWindow(master, title="Periodogram")
    plotWind1.show(plot_dft)

def doOpenFile(master):
    if checkBackgroundTaskRunning(master): return
    global input_data
    global dft_result
    global fit_result
    fileName = filedialog.askopenfilename(parent=master, filetypes=[('Data Files (*.dat *.txt *.csv *.tsv)', '*.dat *.txt *.csv *.tsv')])
    if fileName:
        try:
            global plotWind0
            global plotWind1
            if plotWind1 is not None: plotWind1.show(None)
            if plotWind0 is not None: plotWind0.show(None)
            dft_result = None; dftParamDialog.params_initialized = False
            input_data = None
            fit_result = None
            input_data = dataio.load_data(fileName)
        except Exception as e:
            input_data = None
            messagebox.showinfo(None, "Error: " + str(e), parent=master)
            return
        add_to_log(master, "")
        add_to_log(master, f"{fileName} loaded.")
        add_to_log(master, "")
        plotData(master)

def doSaveDftResult(master):
    if checkBackgroundTaskRunning(master): return
    global dft_result
    if dft_result is None:
        messagebox.showinfo("DC DFT", "No resulted data", parent=master)
        return;
    fileName = filedialog.asksaveasfilename(parent=master, filetypes=[('Tab-separated Files (*.tsv)', '*.tsv')])
    if fileName:
        try:
            dataio.save_result(fileName, dft_result)
        except Exception as e:
            messagebox.showinfo(None, "Error: " + str(e), parent=master)
            return

def doSaveFitResult(master):
    if checkBackgroundTaskRunning(master): return
    global fit_result
    if fit_result is None:
        messagebox.showinfo("Approximation", "No resulted data", parent=master)
        return;
    fileName = filedialog.asksaveasfilename(parent=master, filetypes=[('Tab-separated Files (*.tsv)', '*.tsv')])
    if fileName:
        try:
            dataio.save_result(fileName, fit_result)
        except Exception as e:
            messagebox.showinfo(None, "Error: " + str(e), parent=master)
            return

def doPlotData(master):
    if checkBackgroundTaskRunning(master): return
    global input_data
    if input_data is None:
        messagebox.showinfo("Plot", "No data file open", parent=master)
        return;
    plotData(master)

def doPlotFolded(master):
    if checkBackgroundTaskRunning(master): return
    global input_data
    global fit_result
    if input_data is None:
        messagebox.showinfo("Phase Plot", "No data file open", parent=master)
        return;
    global plotWind0
    if plotWind0 is None: 
        plotWind0 = plotWind.PlotWindow(master, title="Input Data")
    phasePlot.plotFolded(master, plotWind0, input_data, fit_result)

def doPlotDftResult(master, plot_power, plot_frequency):
    if checkBackgroundTaskRunning(master): return
    if dft_result is None:
        messagebox.showinfo("DC DFT", "No DC DFT result", parent=master)
        return;
    plotDftResult(master, plot_power, plot_frequency)

def dft_callback(master, result, msg, action):
    # thread-safe callback
    if threading.current_thread() is not threading.main_thread():
        master.after(0, dft_callback, master, result, msg, action) # arguments are frozen here!
        return
    if action == "started":
        dft_callback.overlay, dft_callback.progressLbl = waitOverlay(master, "waitOverlay_dft", dft.stop_task)
        master.update()
        return
    elif action == "progress":
        if hasattr(dft_callback, "progressLbl") and dft_callback.progressLbl is not None:
            dft_callback.progressLbl.config(text=msg)
        return
    else:
        if hasattr(dft_callback, "progressLbl") and dft_callback.progressLbl is not None:
            dft_callback.progressLbl = None
        if hasattr(dft_callback, "overlay") and dft_callback.overlay is not None:
            dft_callback.overlay.destroy()
            dft_callback.overlay = None
    global dft_result
    dft_result = result
    if msg is None:
        msg = "Unknown error"
    add_to_log(master, msg)
    if dft_result is not None:
        plotDftResult(master, True, True)
    else:
        messagebox.showinfo("DC DFT", msg, parent=master)
    
def doDCDFT(master):
    if checkBackgroundTaskRunning(master): return
    global input_data
    global dft_result

    if input_data is None:
        messagebox.showinfo("DC DFT", "No data file open", parent=master)
        return
        
    t = input_data['Time'].to_numpy()
    m = input_data['Mag'].to_numpy()
    
    dftParamDialog.dftParameters(master, t)
    if not dftParamDialog.param_defined:
        return
    global plotWind1
    if plotWind1 is not None:
        plotWind1.show(None)        
    dft_result = None
    add_to_log(master, "")
    add_to_log(master, "DC DFT started.")
    # Background task
    dft.dcdft_async(master, dft_callback, t, m, dftParamDialog.param_lofreq, dftParamDialog.param_hifreq, dftParamDialog.param_n_intervals)

def fit_callback(master, result, msg, action):
    # thread-safe callback
    if threading.current_thread() is not threading.main_thread():
        master.after(0, fit_callback, master, result, msg, action) # arguments are frozen here!
        return
    if action == "fit_started":
        fit_callback.static_overlay, _ = waitOverlay(master, "waitOverlay_fit1", None)
        master.update()
        return
    if action == "bootstrap_started":
        add_to_log(master, "Bootstrap started.")
        if hasattr(fit_callback, "static_overlay") and fit_callback.static_overlay is not None:
            fit_callback.static_overlay.destroy()
            fit_callback.static_overlay = None
        fit_callback.overlay, fit_callback.progressLbl = waitOverlay(master, "waitOverlay_fit_bootstrap", fit.stop_task)
        master.update()
        return
    elif action == "progress":
        if hasattr(fit_callback, "progressLbl") and fit_callback.progressLbl is not None:
            fit_callback.progressLbl.config(text=msg)
        return
    
    if hasattr(fit_callback, "static_overlay") and fit_callback.static_overlay is not None:
        fit_callback.static_overlay.destroy()
        fit_callback.static_overlay = None
    if hasattr(fit_callback, "progressLbl") and fit_callback.progressLbl is not None:
        fit_callback.progressLbl = None
    if hasattr(fit_callback, "overlay") and fit_callback.overlay is not None:
        fit_callback.overlay.destroy()
        fit_callback.overlay = None

    if action == "finished":
        global fit_result
        fit_result = result
        if msg is None:
            msg = "No message!"
        add_to_log(master, msg)
        if fit_result is not None:
            plotData(master)
        else:
            messagebox.showinfo("Fit", "No result!", parent=master)
        return
    if action == "bootstrap_finished":
        if msg is None:
            msg = "No message!"
        if result is not None:
            add_to_log(master, "")
            add_to_log(master, f"Bootstrap periods standard errors: {result}")
            add_to_log(master, "")
            add_to_log(master, msg)
        else:
            add_to_log(master, "")
            add_to_log(master, msg)
            messagebox.showinfo("Fit", "No bootstrap result!", parent=master)
        return
    if action == "error":
        if msg is None:
            msg = "Unknown error"
        add_to_log(master, msg)
        messagebox.showinfo("Fit", msg, parent=master)
        return
    if msg is None:
        msg = "No message!"
    add_to_log(master, msg)

def doPolyFit(master):
    if checkBackgroundTaskRunning(master): return
    global input_data
    global fit_result

    if input_data is None:
        messagebox.showinfo("Approximation", "No data file open", parent=master)
        return;

    fitParamDialog.fitParameters(master)
    if not fitParamDialog.param_defined:
        return
    fit_result = None
    plotData(master)
    t = input_data['Time'].to_numpy()
    m = input_data['Mag'].to_numpy()
    add_to_log(master, "")
    add_to_log(master, "PolyFit started.")
    # Partially background task
    fit.polyfit(master, fit_callback, 
                t, m,
                fitParamDialog.param_algDegree,
                fitParamDialog.param_periods,
                fitParamDialog.param_degrees,
                fitParamDialog.param_optFlags,
                fitParamDialog.param_bootstrapForErrors)

def doDetrend(master):
    # Replace input data with detrended one: like opening a new file
    if checkBackgroundTaskRunning(master): return    
    global input_data
    global dft_result
    global fit_result

    if fit_result is None:
        messagebox.showinfo("Detrend", "No fit result", parent=master)
        return
    
    if not messagebox.askyesno("Detrend", "Subtract the approximation from the input data?"):
        return
    
    global plotWind0
    global plotWind1
    if plotWind1 is not None: plotWind1.show(None)
    if plotWind0 is not None: plotWind0.show(None)
    input_data = dataio.create_frame(fit_result["Time"], fit_result["Mag"] - fit_result["Fit"])
    dft_result = None; #dftParamDialog.params_initialized = False
    fit_result = None
    add_to_log(master, "")
    add_to_log(master, "Input data replaced with detrended one.")
    add_to_log(master, "")
    plotData(master)
    

def doOpenHelp():
    try:
        webbrowser.open("https://github.com/mpyat2/VS-fit/blob/main/doc/GettingStarted.pdf")
    except Exception as e:
        messagebox.showerror("Error", f"Could not open help page:\n{e}")

# Ensure the 'root' at the top
def bring_to_front(root):
    root.lift()
    root.attributes('-topmost', True)
    root.after(100, lambda: root.attributes('-topmost', False))

##############################################################################

def main():
    print("App started")
    print()
    root = Tk()
    root.protocol("WM_DELETE_WINDOW", lambda: shutdown(root))
    root.title("V*-fit")
    menu = Menu(root)
    root.config(menu=menu)

    filemenu = Menu(menu, tearoff=False)
    menu.add_cascade(label='File', menu=filemenu)
    filemenu.add_command(label='Open Data File...', command=lambda: doOpenFile(root))
    saveresult = Menu(menu, tearoff=False)
    filemenu.add_cascade(label='Save Result...', menu=saveresult)
    saveresult.add_command(label='DCDFT Result...', command=lambda: doSaveDftResult(root))
    saveresult.add_command(label='Fit Result...', command=lambda: doSaveFitResult(root))
    filemenu.add_separator()
    filemenu.add_command(label='Exit', command=lambda: doShutdown(root))
    
    viewmenu = Menu(menu, tearoff=False)
    menu.add_cascade(label='View', menu=viewmenu)
    viewmenuInput = Menu(menu, tearoff=False)
    viewmenu.add_cascade(label='Plot Input', menu=viewmenuInput)
    viewmenuInput.add_command(label='Raw', command=lambda: doPlotData(root))
    viewmenuInput.add_command(label='Phase', command=lambda: doPlotFolded(root))
    viewmenuResult = Menu(menu, tearoff=False)
    viewmenu.add_cascade(label='Plot DFT Result', menu=viewmenuResult)
    viewmenuResult.add_command(label='Power(Frequency)', command=lambda: doPlotDftResult(root, True, True))
    viewmenuResult.add_command(label='Semi-amplitude(Frequency)', command=lambda: doPlotDftResult(root, False, True))
    viewmenuResult.add_command(label='Power(Period)', command=lambda: doPlotDftResult(root, True, False))
    viewmenuResult.add_command(label='Semi-amplitude(Period)', command=lambda: doPlotDftResult(root, False, False))
    viewmenu.add_separator()
    viewmenu.add_command(label='Clear log', command=lambda: clear_log(root))

    operationmenu = Menu(menu, tearoff=False)
    menu.add_cascade(label='Operations', menu=operationmenu)
    operationmenu.add_command(label='DC DFT (Ferraz-Mello)...', command=lambda: doDCDFT(root))
    operationmenu.add_command(label='Polynomial Fit...', command=lambda: doPolyFit(root))
    operationmenu.add_command(label='Detrend', command=lambda: doDetrend(root))

    helpmenu = Menu(menu, tearoff=False)
    menu.add_cascade(label='Help', menu=helpmenu)
    helpmenu.add_command(label='Getting started', command=lambda: doOpenHelp())
    helpmenu.add_separator()
    helpmenu.add_command(label='About...', 
                         command=lambda: messagebox.showinfo(
                             "About V*-fit",
                             "Light Curve fitting program by Maksym Yu. Pyatnytskyy\n\nhttps://github.com/mpyat2/VS-fit",
                             parent=root))

    try:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        img_path = os.path.join(script_dir, "icons")
        
        imgOpen = PhotoImage(master=root, file=os.path.join(img_path, 'Open.png'))
        btnOpen = Button(root, image=imgOpen, command=lambda: doOpenFile(root))
        btnOpen.image = imgOpen  # keep reference (? works without it)
        btnOpen.grid(row=0, column=0, padx=5, pady=5)
        
        imgDft = PhotoImage(master=root, file=os.path.join(img_path, 'DCDFT.png'))
        btnDft = Button(root, image=imgDft, command=lambda: doDCDFT(root))
        btnDft.image = imgDft  # keep reference (? works without it)
        btnDft.grid(row=0, column=1, padx=5, pady=5)
        
        imgApprox = PhotoImage(master=root, file=os.path.join(img_path, 'Approx.png'))
        btnApprox = Button(root, image=imgApprox, command=lambda: doPolyFit(root))
        btnApprox.image = imgApprox  # keep reference (? works without it)
        btnApprox.grid(row=0, column=2, padx=5, pady=5)
        
        imgDetrend = PhotoImage(master=root, file=os.path.join(img_path, 'Detrend.png'))
        btnDetrend = Button(root, image=imgDetrend, command=lambda: doDetrend(root))
        btnDetrend.image = imgDetrend  # keep reference (? works without it)
        btnDetrend.grid(row=0, column=3, padx=5, pady=5)
    except Exception as e:
        print(e)
    
    root.geometry("+40+40")
    bring_to_front(root)
    
    global log_window
    log_window = LogWindow(root, geometry="600x300+440+40")

    root.resizable(False, False)
    root.mainloop()

##############################################################################

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"Fatal Error: {e}")
    finally:
        print()
        #input("Press ENTER to continue: ")
        print("END")        
        sys.exit()
