# Created by ChatGPT, Claude (with modifications)
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

#    def _on_right_click(self, event):
#        # event.x/y are canvas pixel coords; use inaxes to get data coords
#        if event.inaxes is None:
#            return
#        x = event.xdata
#        if x is None:
#            return
#        self.root.clipboard_clear()
#        self.root.clipboard_append(str(x))
#        self.root.update()  # keeps clipboard alive after window closes

    def _on_right_click(self, event):
        if event.inaxes is None:
            return
        x = event.xdata
        if x is None:
            return
        self.root.clipboard_clear()
        self.root.clipboard_append(str(x))
        self.root.update()
        self._show_copied_hint(f"X = {x}  copied to clipboard")

    def _show_copied_hint(self, text, duration=2000):
        label = tk.Label(
            self.win, text=text,
            bg="lightyellow", fg="black",
            relief="solid", borderwidth=1,
            font=("TkDefaultFont", 9)
        )
        label.place(relx=0.5, rely=0.02, anchor="n")  # top-center of the plot window

        # def blink(count=6):
        #     if count <= 0:
        #         label.destroy()
        #         return
        #     current = label.cget("fg")
        #     label.config(fg="lightyellow" if current == "black" else "black")
        #     self.win.after(200, lambda: blink(count - 1))

        def hide_label():
            label.destroy()
        
        self.win.after(duration, hide_label)

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

        # Bind right-click via matplotlib (works correctly with zoom/pan)
        self.canvas.mpl_connect("button_press_event", lambda e: self._on_right_click(e) if e.button == 3 else None)

        # Add toolbar
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.win)
        self.toolbar.update()
        self.toolbar.pack(fill="x")

        self.win.lift()
        self.win.focus_force()

        # return axes if you want to customize further
        return ax
