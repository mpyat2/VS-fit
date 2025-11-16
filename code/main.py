import sys
from tkinter import Tk, Menu, filedialog, messagebox
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

def plotDftResult(master, plot_power):
    def plot_dft(ax):
        global dft_result
        if plot_power:
            ax.plot(dft_result['freq'], dft_result['pow'], color = 'k', linestyle='-')
            ax.set_xlabel('Frequency')
            ax.set_ylabel('Power')
        else:
            ax.plot(dft_result['freq'], dft_result['amp'], color = 'k', linestyle='-')
            ax.set_xlabel('Frequency')
            ax.set_ylabel('Semi-amplitude')
        ax.set_title('DCDFT')
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
    fileName = filedialog.askopenfilename(parent=master, filetypes=[('Data Files (*.dat *.txt *.csv *.tsv)', '*.dat *.txt *.csv *.tsv')])
    if fileName:
        try:
            global plotWind0
            global plotWind1
            global plotWind2
            if plotWind2 is not None: plotWind2.show(None)
            if plotWind1 is not None: plotWind1.show(None)
            if plotWind0 is not None: plotWind0.show(None)
            dft_result = None
            input_data = None
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
        messagebox.showinfo(None, "No resulted data", parent=master)
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
        messagebox.showinfo(None, "No resulted data", parent=master)
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
        messagebox.showinfo(None, "No data file open", parent=master)
        return;
    plotData(master)

def doPlotFolded(master):
    global input_data
    if input_data is None:
        messagebox.showinfo(None, "No data file open", parent=master)
        return;
    global plotWind0
    if plotWind0 is None: 
        plotWind0 = plotWind.PlotWindow(master, title="Input Data")
    phasePlot.plotFolded(master, plotWind0, input_data)

def doPlotDftResult(master, plot_power):
    if dft_result is None:
        messagebox.showinfo(None, "No resulted data", parent=master)
        return;
    plotDftResult(master, plot_power)

def doDCDFT(master):
    global input_data
    global dft_result

    if input_data is None:
        messagebox.showinfo(None, "No data file open", parent=master)
        return;
    dftParamDialog.dftParameters(master)
    if not dftParamDialog.param_defined:
        return
    global plotWind1
    if plotWind1 is not None:
        plotWind1.show(None)        
    dft_result = None
    t = input_data['Time'].to_numpy()
    m = input_data['Mag'].to_numpy()
    add_to_log(master, "DCDFT started.")
    try:
        master.config(cursor="watch")
        master.update()
        try:
            t0 = time.time()
            dft_result = dft.dcdft(t, m, 
                                   lowfreq=dftParamDialog.param_lofreq, 
                                   hifreq=dftParamDialog.param_hifreq, 
                                   n_freq=dftParamDialog.param_n_intervals, 
                                   mcv_mode=False)
            msg = f"DCDFT calculation time {(time.time() - t0):.2f} s"
            print(msg)
            add_to_log(master, msg)
        finally:
            master.config(cursor="")
    except Exception as e:
        messagebox.showinfo(None, "Error: " + str(e), parent=master)
        return
    plotDftResult(master, True)

def doPolyFit(master):
    global input_data
    global fit_result

    if input_data is None:
        messagebox.showinfo(None, "No data file open", parent=master)
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
        master.config(cursor="watch")
        master.update()
        try:
            t0 = time.time()
            out = fit.polyfit(t, m,
                              fitParamDialog.param_algeDegree,
                              [fitParamDialog.param_trig1Period,
                                fitParamDialog.param_trig2Period,
                                fitParamDialog.param_trig3Period],
                              [fitParamDialog.param_trig1Degree,                                     
                                fitParamDialog.param_trig2Degree,
                                fitParamDialog.param_trig3Degree],
                              [fitParamDialog.param_trig1Optimize,
                                fitParamDialog.param_trig2Optimize,
                                fitParamDialog.param_trig3Optimize],
                              compute_bootstrap=fitParamDialog.param_bootstrapForErrors)
            fit_result = out['fit_result']
            msg = f"PolyFit calculation time {(time.time() - t0):.2f} s"
            print(msg)
            add_to_log(master, msg)
            add_to_log(master, out['message'])
            
        finally:
            master.config(cursor="")
    except Exception as e:
        messagebox.showinfo(None, "Error: " + str(e), parent=master)
        return
    plotFitResult(master)


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
    viewmenuResult.add_command(label='Power', command=lambda: doPlotDftResult(root, True))
    viewmenuResult.add_command(label='Semi-amplitude', command=lambda: doPlotDftResult(root, False))
    viewmenu.add_separator()
    viewmenu.add_command(label='Clear log', command=lambda: clear_log(root))

    operationmenu = Menu(menu, tearoff=False)
    menu.add_cascade(label='Operations', menu=operationmenu)
    operationmenu.add_command(label='DCDFT (Ferraz-Mello)', command=lambda: doDCDFT(root))
    operationmenu.add_command(label='Polynomial Fit', command=lambda: doPolyFit(root))

    root.geometry("320x100+40+40")

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
