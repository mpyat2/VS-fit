import sys
from tkinter import Tk, Menu, mainloop, filedialog, messagebox
from log_window import LogWindow
import matplotlib.pyplot as plt
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

input_data = None
dft_result = None
fit_result = None
log_window = None

def shutdown(master):
    plt.close('all')
    master.destroy()

def plotData():
    fig = plt.figure(0)
    fig.clear()
    global input_data
    plt.plot(input_data['Time'], input_data['Mag'], 'r.')
    plt.ylim(max(input_data['Mag']), min(input_data['Mag']))
    plt.title('Input')
    plt.show()
    plt.pause(0.001)  # Forces redraw: needed in Spyder

def plotDftResult(plot_power):
    fig = plt.figure(1)
    fig.clear()
    global dft_result
    if plot_power:
        plt.plot(dft_result['freq'], dft_result['pow'], color = 'k', linestyle='-')
        #plt.xlabel('Frequency', fontsize=15)
        #plt.ylabel('Power', fontsize=15)
    else:
        plt.plot(dft_result['freq'], dft_result['amp'], color = 'k', linestyle='-')
        #plt.xlabel('Frequency', fontsize=15)
        #plt.ylabel('Semi-amplitude', fontsize=15)
    #plt.tick_params(axis='both', which='major', labelsize=15)
    plt.title('DCDFT')
    plt.show()
    plt.pause(0.001)  # Forces redraw: needed in Spyder

def plotFitResult():
    fig = plt.figure(2)
    fig.clear()
    global fit_result    
    plt.plot(fit_result['Time'], fit_result['Mag'], 'r.')
    plt.plot(fit_result['Time'], fit_result['Fit'], 'k.')
    #plt.plot(fit_result['Time'], fit_result['Fit_L'], 'b--')
    #plt.plot(fit_result['Time'], fit_result['Fit_U'], 'b--')
    y_min = min(min(fit_result['Mag']), min(fit_result['Fit']))
    y_max = max(max(fit_result['Mag']), max(fit_result['Fit']))
    plt.ylim(y_max, y_min)
    #plt.xlabel('Time', fontsize=15)
    #plt.ylabel('Magnitude', fontsize=15)
    #plt.tick_params(axis='both', which='major', labelsize=15)
    plt.title('Approximation')
    plt.show()
    plt.pause(0.001)  # Forces redraw: needed in Spyder

def openFile(master):
    global input_data
    global dft_result
    fileName = filedialog.askopenfilename(filetypes=[('Data Files (*.dat *.txt *.csv *.tsv)', '*.dat *.txt *.csv *.tsv')])
    #print("Selected file: ", fileName);
    if fileName:
        try:
            if plt.fignum_exists(2):
                plt.figure(2).clear()
            if plt.fignum_exists(1):
                plt.figure(1).clear()
            if plt.fignum_exists(0):
                plt.figure(0).clear()
            dft_result = None
            input_data = None
            input_data = pd.read_csv(fileName, names=NAMES, dtype=DTYPE, sep='\t', header=None, usecols=[0, 1])
        except Exception as e:
            input_data = None
            messagebox.showinfo(None, "Error: " + str(e))
            return
        log_window.add_line(f"{fileName} loaded.")
        plotData()
    else:
        input_data = None

def saveResult(master):
    global dft_result
    if dft_result is None:
        messagebox.showinfo(None, "No resulted data")
        return;
    fileName = filedialog.asksaveasfilename(filetypes=[('Tab-separated Files (*.tsv)', '*.tsv')])
    #print("Selected file: ", fileName);
    if fileName:
        try:
            filePath = pathlib.Path(fileName)
            if filePath.suffix == "":
                fileName += ".tsv"
            dft_result.to_csv(fileName, index=False, sep='\t')
        except Exception as e:
            messagebox.showinfo(None, "Error: " + str(e))
            return

def saveFitResult(master):
    global fit_result
    if fit_result is None:
        messagebox.showinfo(None, "No resulted data")
        return;
    fileName = filedialog.asksaveasfilename(filetypes=[('Tab-separated Files (*.tsv)', '*.tsv')])
    if fileName :
        try:
            filePath = pathlib.Path(fileName)
            if filePath.suffix == "":
                fileName += ".tsv"
            fit_result.to_csv(fileName, index=False, sep='\t')
        except Exception as e:
            messagebox.showinfo(None, "Error: " + str(e))
            return

def doPlotData(master):
    global input_data
    if input_data is None:
        messagebox.showinfo(None, "No data file open")
        return;
    plotData()

def doPlotFolded(master):
    global input_data
    if input_data is None:
        messagebox.showinfo(None, "No data file open")
        return;
    phasePlot.plotFolded(master, input_data)

def doPlotDftResult(master, plot_power):
    if dft_result is None:
        messagebox.showinfo(None, "No resulted data")
        return;
    plotDftResult(plot_power)

def doDCDFT(master):
    global input_data
    global dft_result

    if input_data is None:
        messagebox.showinfo(None, "No data file open")
        return;
    dftParamDialog.dftParameters(master)
    if not dftParamDialog.param_defined:
        return
    fig = plt.figure(1)
    fig.clear()
    dft_result = None
    t = input_data['Time'].to_numpy()
    m = input_data['Mag'].to_numpy()
    log_window.add_line("DCDFT started.")
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
            log_window.add_line(msg)
        finally:
            master.config(cursor="")
    except Exception as e:
        messagebox.showinfo(None, "Error: " + str(e))
        return
    plotDftResult(True)

def doPolyFit(master):
    global input_data
    global fit_result

    if input_data is None:
        messagebox.showinfo(None, "No data file open")
        return;
    fitParamDialog.fitParameters(master)
    if not fitParamDialog.param_defined:
        return
    fig = plt.figure(2)
    fig.clear()
    fit_result = None
    t = input_data['Time'].to_numpy()
    m = input_data['Mag'].to_numpy()
    log_window.add_line("PolyFit started.")
    try:
        master.config(cursor="watch")
        master.update()
        try:
            t0 = time.time()
            fit_result = fit.polyfit(t, m,
                                     fitParamDialog.param_algeDegree,
                                    [fitParamDialog.param_trig1Period,
                                     fitParamDialog.param_trig2Period,
                                     fitParamDialog.param_trig3Period],
                                    [fitParamDialog.param_trig1Degree,                                     
                                     fitParamDialog.param_trig2Degree,
                                     fitParamDialog.param_trig3Degree],
                                    [fitParamDialog.param_trig1Optimize,
                                     fitParamDialog.param_trig2Optimize,
                                     fitParamDialog.param_trig3Optimize])
            msg = f"DCDFT calculation time {(time.time() - t0):.2f} s"
            print(msg)
            log_window.add_line(msg)
        finally:
            master.config(cursor="")
    except Exception as e:
        messagebox.showinfo(None, "Error: " + str(e))
        return
    plotFitResult()


# Ensure the 'root' at the top
def bring_to_front(root):
    root.lift()
    root.attributes('-topmost', True)
    root.after(100, lambda: root.attributes('-topmost', False))

def clear_log(root):
    log_window.clear()

##############################################################################

def main():
    print("App started")
    print()
    root = Tk()
    root.protocol("WM_DELETE_WINDOW", lambda: shutdown(root))
    root.title("V*-mini")
    menu = Menu(root)
    root.config(menu=menu)

    filemenu = Menu(menu, tearoff=False)
    menu.add_cascade(label='File', menu=filemenu)
    filemenu.add_command(label='Open Data File...', command=lambda: openFile(root))
    saveresult = Menu(menu, tearoff=False)
    filemenu.add_cascade(label='Save Result...', menu=saveresult)
    saveresult.add_command(label='DCDFT Result...', command=lambda: saveResult(root))
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

    root.geometry("320x100+100+100")

    bring_to_front(root)
    
    global log_window
    log_window = LogWindow(root)
    
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
