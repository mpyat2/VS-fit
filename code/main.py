import os
import sys
from tkinter import Tk, Frame, Label, Menu, Button, PhotoImage, filedialog, messagebox
from log_window import LogWindow
##
import matplotlib
matplotlib.use("TkAgg")
import plotWind
##
import pandas as pd
import pathlib
import time
import dft
import dftParamDialog
import fit
import fitParamDialog
import phasePlot

## Data file spec (no header row, tab-separated)
NAMES = ['Time', 'Mag']
DTYPE = {'Time': 'float64', 'Mag': 'float64'}

log_window = None

plotWind0 = None # Input data
plotWind1 = None # DFT result
plotWind2 = None # Fit result

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
        log_window.add_line(text)
    except Exception as e:
        messagebox.showinfo(None, "Error: " + str(e), parent=master)

def shutdown(master):
    try:
        master.quit() # stop mainloop if running: required in Spyder
    except:
        pass    
    master.destroy()

def waitOverlay(master):
    overlay = Frame(master, bg="#cccccc")
    overlay.place(relx=0, rely=0, relwidth=1, relheight=1)
    Label(overlay, text="Please wait...", bg="#c0c0c0").place(relx=0.5, rely=0.5, anchor="center")
    return overlay

def plotData(master):
    def plot_input(ax):
        ax.plot(input_data['Time'], input_data['Mag'], '.', color='royalblue')
        ax.set_ylim(max(input_data["Mag"]), min(input_data["Mag"]))
        ax.set_title('Input')
        ax.set_xlabel('Time')
        ax.set_ylabel('Magnitude')
        ax.grid(True, linestyle='--', color='gray', alpha=0.3)
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

def plotFitResult(master):
    def plot_fit(ax):
        global fit_result    
        ax.plot(fit_result['Time'], fit_result['Mag'], '.', color='royalblue')
        ax.plot(fit_result['Time'], fit_result['Fit'], 'k.')
        y_min = min(min(fit_result['Mag']), min(fit_result['Fit']))
        y_max = max(max(fit_result['Mag']), max(fit_result['Fit']))
        ax.set_ylim(y_max, y_min)
        ax.set_title('Approximation')
        ax.set_xlabel('Time')
        ax.set_ylabel('Magnitude')
        ax.grid(True, linestyle='--', color='gray', alpha=0.3)
    global plotWind2        
    if plotWind2 is None: 
        plotWind2 =  plotWind.PlotWindow(master, title="Approximation")
    plotWind2.show(plot_fit)

def openFile(master):
    global input_data
    global dft_result
    global fit_result
    fileName = filedialog.askopenfilename(parent=master, filetypes=[('Data Files (*.dat *.txt *.csv *.tsv)', '*.dat *.txt *.csv *.tsv')])
    if fileName:
        try:
            global plotWind0
            global plotWind1
            global plotWind2
            if plotWind2 is not None: plotWind2.show(None)
            if plotWind1 is not None: plotWind1.show(None)
            if plotWind0 is not None: plotWind0.show(None)
            dft_result = None; dftParamDialog.params_initialized = False
            input_data = None
            fit_result = None
            input_data = pd.read_csv(fileName, 
                                     comment='#', skip_blank_lines=True,
                                     names=NAMES, dtype=DTYPE, header=None, usecols=[0, 1],
                                     sep=r'\s+', engine='python')
        except Exception as e:
            input_data = None
            messagebox.showinfo(None, "Error: " + str(e), parent=master)
            return
        add_to_log(master, "")
        add_to_log(master, f"{fileName} loaded.")
        add_to_log(master, "")
        plotData(master)

def save_result(fileName, data):
    filePath = pathlib.Path(fileName)
    if filePath.suffix == "":
        fileName += ".tsv"
    with open(fileName, 'w', newline='') as f:
        # write header line with '#'
        f.write('#' + '\t'.join(data.columns) + '\n')
        data.to_csv(f, index=False, header=False, sep='\t')

def saveDftResult(master):
    global dft_result
    if dft_result is None:
        messagebox.showinfo("DC DFT", "No resulted data", parent=master)
        return;
    fileName = filedialog.asksaveasfilename(parent=master, filetypes=[('Tab-separated Files (*.tsv)', '*.tsv')])
    #print("Selected file: ", fileName);
    if fileName:
        try:
            save_result(fileName, dft_result)
        except Exception as e:
            messagebox.showinfo(None, "Error: " + str(e), parent=master)
            return

def saveFitResult(master):
    global fit_result
    if fit_result is None:
        messagebox.showinfo("Approximation", "No resulted data", parent=master)
        return;
    fileName = filedialog.asksaveasfilename(parent=master, filetypes=[('Tab-separated Files (*.tsv)', '*.tsv')])
    if fileName :
        try:
            save_result(fileName, fit_result)
        except Exception as e:
            messagebox.showinfo(None, "Error: " + str(e), parent=master)
            return

def doPlotData(master):
    global input_data
    if input_data is None:
        messagebox.showinfo("Plot", "No data file open", parent=master)
        return;
    plotData(master)

def doPlotFolded(master):
    global input_data
    if input_data is None:
        messagebox.showinfo("Phase Plot", "No data file open", parent=master)
        return;
    global plotWind0
    if plotWind0 is None: 
        plotWind0 = plotWind.PlotWindow(master, title="Input Data")
    phasePlot.plotFolded(master, plotWind0, input_data)

def doPlotDftResult(master, plot_power, plot_frequency):
    if dft_result is None:
        messagebox.showinfo("DC DFT", "No DC DFT result", parent=master)
        return;
    plotDftResult(master, plot_power, plot_frequency)

def doDCDFT(master):
    global input_data
    global dft_result

    if input_data is None:
        messagebox.showinfo("DC DFT", "No data file open", parent=master)
        return;
        
    t = input_data['Time'].to_numpy()
    m = input_data['Mag'].to_numpy()
    
    if not dftParamDialog.params_initialized:
        dftParamDialog.params_initialized = True
        interval = dft.median_interval(t)
        if interval > 0: # False also for 'nan'
            hifreq = 1.0 / interval / 2.0; # Approximate Nyquist
            if hifreq > 50.0: hifreq = 50.0
            time_interval = max(t) - min(t)
            recommended_freq_resolution = 0.05 / time_interval
            lofreq = recommended_freq_resolution
            if lofreq < 1.0: lofreq = round(lofreq, 6)
            if lofreq < 0.000001: lofreq = 0.000001
            if hifreq <= lofreq:
                hifreq = lofreq + 100 * recommended_freq_resolution
            dftParamDialog.param_lofreq = lofreq
            dftParamDialog.param_hifreq = hifreq
            dftParamDialog.param_n_intervals = round((hifreq - lofreq) / recommended_freq_resolution) + 1
    
    dftParamDialog.dftParameters(master)
    if not dftParamDialog.param_defined:
        return
    global plotWind1
    if plotWind1 is not None:
        plotWind1.show(None)        
    dft_result = None
    add_to_log(master, "DC DFT started.")
    try:
        master.focus_force()
        master.config(cursor="wait")
        overlay = waitOverlay(master)
        master.update()
        try:
            t0 = time.time()
            dft_result = dft.dcdft(t, m, 
                                   lowfreq=dftParamDialog.param_lofreq, 
                                   hifreq=dftParamDialog.param_hifreq, 
                                   n_intervals=dftParamDialog.param_n_intervals, 
                                   mcv_mode=False)
            msg = f"DC DFT calculation time {(time.time() - t0):.2f} s"
            print(msg)
            add_to_log(master, msg)
        finally:
            overlay.destroy()
            master.config(cursor="")
    except Exception as e:
        messagebox.showinfo(None, "Error: " + str(e), parent=master)
        return
    plotDftResult(master, True, True)

def doPolyFit(master):
    global input_data
    global fit_result

    if input_data is None:
        messagebox.showinfo("Approximation", "No data file open", parent=master)
        return;
    fitParamDialog.fitParameters(master)
    if not fitParamDialog.param_defined:
        return
    global plotWind2
    if plotWind2 is not None:
        plotWind2.show(None)        
    fit_result = None
    t = input_data['Time'].to_numpy()
    m = input_data['Mag'].to_numpy()
    add_to_log(master, "PolyFit started.")
    try:
        master.focus_force()
        master.config(cursor="wait")
        overlay = waitOverlay(master)
        master.update()
        try:
            t0 = time.time()
            out = fit.polyfit(t, m,
                              fitParamDialog.param_algDegree,
                              fitParamDialog.param_periods,
                              fitParamDialog.param_degrees,
                              fitParamDialog.param_optFlags,
                              compute_bootstrap=fitParamDialog.param_bootstrapForErrors)
            fit_result = out['fit_result']
            msg = f"PolyFit calculation time {(time.time() - t0):.2f} s"
            print(msg)
            add_to_log(master, msg)
            add_to_log(master, out['message'])
        finally:
            overlay.destroy()
            master.config(cursor="")
    except Exception as e:
        messagebox.showinfo(None, "Error: " + str(e), parent=master)
        return
    plotFitResult(master)


def doDetrend(master):
    # Replace input data with detrended one: like opening a new file
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
    global plotWind2
    if plotWind2 is not None: plotWind2.show(None)
    if plotWind1 is not None: plotWind1.show(None)
    if plotWind0 is not None: plotWind0.show(None)
    input_data = pd.DataFrame({
        "Time": fit_result["Time"],
        "Mag": fit_result["Mag"] - fit_result["Fit"]
    })
    dft_result = None; dftParamDialog.params_initialized = False
    fit_result = None
    add_to_log(master, "")
    add_to_log(master, "Input data replaced with detrended one.")
    add_to_log(master, "")
    plotData(master)
    

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
    filemenu.add_command(label='Open Data File...', command=lambda: openFile(root))
    saveresult = Menu(menu, tearoff=False)
    filemenu.add_cascade(label='Save Result...', menu=saveresult)
    saveresult.add_command(label='DCDFT Result...', command=lambda: saveDftResult(root))
    saveresult.add_command(label='Fit Result...', command=lambda: saveFitResult(root))
    filemenu.add_separator()
    filemenu.add_command(label='Exit', command=lambda: shutdown(root))
    
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
    helpmenu.add_command(label='About...', 
                         command=lambda: messagebox.showinfo(
                             "About V*-fit",
                             "Light Curve fitting program by Maksym Yu. Pyatnytskyy\n\nhttps://github.com/mpyat2/VS-fit",
                             parent=root))

    try:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        img_path = os.path.join(script_dir, "icons")
        
        imgOpen = PhotoImage(master=root, file=os.path.join(img_path, 'Open.png'))
        btnOpen = Button(root, image=imgOpen, command=lambda: openFile(root))
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
