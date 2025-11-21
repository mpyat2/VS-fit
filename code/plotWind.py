# Created by ChatGPT (with modifications)
import tkinter as tk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

class PlotWindow:
    def __init__(self, root, title="Plot Window", size="800x600"):
        self.root = root
        self.title = title
        self.size = size

        self.win = None
        #self.fig = None
        self.canvas = None
        self.toolbar = None

    def show(self, plot_func):
        # Create window if needed
        if self.win is None or not self.win.winfo_exists():
            self.win = tk.Toplevel(self.root)
            self.win.title(self.title)
            self.win.geometry(self.size)

        # Destroy toolbar BEFORE canvas
        if self.toolbar is not None:
            self.toolbar.destroy()
            self.toolbar = None

        # Destroy old widgets
        if self.canvas is not None:
            #print("self.canvas.get_tk_widget().destroy()")
            self.canvas.get_tk_widget().destroy()
            self.canvas = None
        
        # Always create a new figure (do NOT reuse old one to avoid toolbar errors)
        fig = Figure()
            
        # Create axes
        ax = fig.add_subplot(111)

        # Plot logic
        if plot_func is not None:
            plot_func(ax)

        # Embed into Tk
        self.canvas = FigureCanvasTkAgg(fig, master=self.win)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        # Add toolbar
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.win)
        self.toolbar.update()
        self.toolbar.pack(fill="x")

        self.win.lift()
        self.win.focus_force()

        # return axes if you want to customize further
        return ax
