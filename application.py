import math
import os
import glob
import platform
import random
import time
import shutil
import datetime
import warnings

warnings.filterwarnings("ignore")

from multiprocessing import Process, Pipe, Queue
from threading import Thread

import tkinter as tk
from tkinter import ttk
from tkinter.colorchooser import askcolor

from settings import *
from settings import Settings

from simulation_wrapper import  run_simulation

#--------- Do analizy  -----------------------------
import matplotlib.ticker as ticker
from   matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from   matplotlib import figure
from   matplotlib import pyplot as plt
from   matplotlib import gridspec
import matplotlib.image as pltimage
plt.style.use("ggplot")

import seaborn     as sbn
import pandas      as pd
import numpy       as np
import scipy.stats as sci

from  statsmodels.graphics.tsaplots import plot_acf




#------------------------------------------------


#Słownik do zmiany opisu osi Y histogramu na podstawie jego rodzaju
ylabel = {
            "count": "Liczba symulacji",
            "frequency": "Częstotliwość",
            "proportion": "Proporcje",
            "density":  "Prawdopodobieństwo",
            "percent":  "Procent ogółu"
}

#Opcje dla rozwiajalnego menu Typy histogramu
hist_options = ("count", "density")

#Mapwanie nazw histogramu do ładnych opisów na osi Y histogramy
hist_type_to_name = {
    "count": "Ilość obserwacji",
    "density": "Prawdopodobieństwo obserwacji"
}

#Mapowanie opisów osi Y histogramu do jego nazwy w scipy
hist_name_to_type = {
    "Ilość obserwacji": "count",
    "Prawdopodobieństwo obserwacji": "density"
}

#'johnsonsb','johnsonsu
#nazwy rozkładów gęstości dls scipy.stats
pdf_name = ['alpha',
        'anglit','arcsine','beta','betaprime','bradford','burr','burr12',
        'cauchy','chi','chi2','cosine','dgamma','dweibull','erlang','expon','exponnorm',
        'exponweib','exponpow','f','fatiguelife','fisk','foldcauchy','foldnorm',
        'genlogistic','genpareto','gennorm','genexpon','genextreme','gausshyper',
        'gamma','gengamma','genhalflogistic','gilbrat','gompertz','gumbel_r','gumbel_l',
        'halfcauchy','halflogistic','halfnorm','halfgennorm','hypsecant','invgamma','invgauss',
        'invweibull','kstwobign','laplace','levy','levy_l','logistic',
        'loggamma','loglaplace','lognorm','lomax','maxwell','mielke','nakagami','ncx2','ncf','nct',
        'norm','pareto','pearson3','powerlaw','powerlognorm','powernorm','rdist','reciprocal','rayleigh',
        'rice','recipinvgauss','semicircular','t','triang','truncexpon','truncnorm','tukeylambda','uniform',
        'vonmises','vonmises_line','wald','weibull_min','weibull_max']




class appGUI:
    '''
    Main class, starting point for whole simulation programm , Creates GUI,
    and mediate betwwen, simulation process, data collector and pygame. It is also
    responsible for presenting result of simulation

    '''

    def __init__(self) -> None:
        self.settings = Settings.load_from_file('settings.json')
        self.process: Process = None
        self.thread: Thread = None

    
        self.data_collected: list = None 

        self.simulation_ends_sucesfuly: bool = False
        self.simulation_running: bool = False
        self.best_fitting_pdf_iter = []
        self.best_fitting_pdf_ipd = []

        self._widget_states = []

        self.make_ui()
        self.fill_entries_from_setting()                                                
        self.on_dry_run_checked_actions()  
        self.on_bin_ckecked()      

        self.pd_series = []
        self.cykle_series = []    
        self.dataframe_all= pd.DataFrame()                                




    @staticmethod
    def hex_to_rgb(value) -> tuple[int, int, int]:
        '''
        Converts hex string #RRGGBB to tuple
        (R, G, B)
        '''
        value = value.lstrip('#')
        lv = len(value)
        return tuple(int(value[i:i+lv//3], 16) for i in range(0, lv, lv//3))



    @staticmethod
    def rgb_to_hex(rgb) -> str:
        '''
        Converts RGB tuple to string hex representation
        #RRGGBB
        '''
        return '#%02x%02x%02x' % tuple(rgb)



#------------------------------------------------------------------------------------------
#                                       ustawienia
# -----------------------------------------------------------------------------------------

    def update_settings(self) -> None:
        """
        Updates settings object with values and states from 
        entries and widgets
        """

        s = self.settings

        s.GRID_SIZE             = (self.grid_width_in_tails.get(), self.grid_height_in_tails.get())
        s.GRID_SIZE_PIXELS      = (self.pygame_window_width.get(),self.pygame_window_height.get())
        s.P1_PROBABILITY        = float(self.P1_probability_value.get())
        s.P2_PROBABILITY        = float(self.P2_probability_value.get())
        s.HEAL_CYKLES           = self.heal_cykles.get()
        s.FPS                   = self.simulation_framerate.get()
        s.ITERATIONS            = self.iteration_limit.get()
        s.DRYRUN                = self.dry_run.get()
        s.SIMULATION_CYKLES     = self.simcykles.get()



    def fill_entries_from_setting(self) -> None:
        '''
        Takes values from settings obj and atribute them to widgets
        '''
        self.grid_width_in_tails.set(self.settings.GRID_SIZE[0])
        self.grid_height_in_tails.set(self.settings.GRID_SIZE[1])
        self.pygame_window_width.set(self.settings.GRID_SIZE_PIXELS[0])
        self.pygame_window_height.set(self.settings.GRID_SIZE_PIXELS[1])
        self.P1_probability_value.set(self.settings.P1_PROBABILITY)
        self.P2_probability_value.set(self.settings.P2_PROBABILITY)
        self.heal_cykles.set(self.settings.HEAL_CYKLES)
        self.simulation_framerate.set(self.settings.FPS)
        self.iteration_limit.set(self.settings.ITERATIONS)
        self.dry_run.set(1 if self.settings.DRYRUN else 0)
        self.simcykles.set(self.settings.SIMULATION_CYKLES)



    def make_ui(self):

        self.window = tk.Tk()
        self.window.title("Badanie dynamiki ewolucji defektów w dwuwymiarowych modelowych układach za pomocą metody Monte Carlo")
        self.window.geometry("1270x800+150+20")
        self.window.resizable(0,0)
        self.window.protocol("WM_DELETE_WINDOW", self.on_close)

        self.grid_width_in_tails  = tk.IntVar()       #szerokosc siatki
        self.grid_height_in_tails = tk.IntVar()       #wysokosc siatki
        self.pygame_window_width  = tk.IntVar()       # szerokość okna w px
        self.pygame_window_height = tk.IntVar()       # wysokość okna symulacji w px
        self.P1_probability_value = tk.StringVar()    # prawdopodobieństwo P1
        self.P2_probability_value = tk.StringVar()    # prawdopodobieństwo P2
        self.heal_cykles          = tk.IntVar()       # cykle leczenia
        self.simulation_framerate = tk.IntVar()       # framerate dla pygame
        self.dry_run              = tk.BooleanVar()   # czy pokazywać okno symulacji
        self.iteration_limit      = tk.IntVar()       # ile razy powtórzyć symulacje
        self.simcykles            = tk.IntVar()       # Po ilu cyklach kończy się symulacja, jeżeli 0 to symulacja trwa aż do momentu znalezienia ścieżki perkolacyjnej
        self.combo                = tk.StringVar()
        self.kde                  = tk.IntVar(value=0)
        self.progressbar          = tk.IntVar(value=0)
        self.progress_label       = tk.StringVar(value="0 % [0/0]")
        self.pdf                  = tk.IntVar(value = 0)
        self.eta_label            = tk.StringVar(value="ETA:")
        self.bins                 = tk.IntVar(value=1)
        self.bin_check            = tk.BooleanVar(value=False)
        self.bin_slider_value     = tk.StringVar(value="Liczba przedziałów: 1")
        self.toolbar_check        = tk.BooleanVar(value=False)
        self.meanline             = tk.BooleanVar(value=False)
        self.boxplot_show         = tk.BooleanVar(value=True)
        self.boxplot_show_fliers  = tk.BooleanVar(value=False)


        

        self.frame1 = ttk.LabelFrame(self.window, text="Symulacja", relief=tk.GROOVE)
        self.frame1.grid(column=0, row = 0, columnspan=8, sticky=tk.N + tk.W + tk.S + tk.E)

        self.frame2 = ttk.LabelFrame(self.window, text="Opcje", relief=tk.GROOVE)
        self.frame2.grid(column=0, row=1, columnspan=8,sticky=tk.N + tk.W + tk.S + tk.E)

        self.frame4 = ttk.LabelFrame(self.window, text="Kolory", relief=tk.GROOVE)
        self.frame4.grid(column=0, row=2, columnspan=8,sticky=tk.N + tk.W + tk.S + tk.E)

        self.frame5 = ttk.Labelframe(self.window, relief=tk.GROOVE)
        self.frame5.grid(column=0, row=4, rowspan = 3, columnspan=8, sticky='nswe')

        self.frame3 = ttk.LabelFrame(self.window, text="Kontrola", relief=tk.GROOVE)
        self.frame3.grid(column=0, row=3,columnspan=8,sticky=tk.W + tk.E)

        self.notebook = ttk.Notebook(master=self.window)
        self.notebook.grid(column=9, row=0, rowspan=6, columnspan=15, sticky='nwse')

        self.matplotframe = ttk.Frame(self.notebook)
        self.matplotframe.grid(column=0, row=0, rowspan=6, columnspan=15, sticky='nswe')

        self.infopage= tk.Text(self.notebook, height=3, width=10)
        self.infopage.grid(column=0, row=0)
    
        
        self.notebook.add(self.matplotframe, text="Wykresy")
        self.notebook.add(self.infopage, text="Statystyka")


        self.type_of_plots = tk.Frame(self.window)
        self.type_of_plots.grid(column=9, row = 6, sticky=tk.N + tk.W + tk.S + tk.E)

        progress = ttk.Progressbar(master=self.frame1, orient="horizontal", mode='determinate', variable=self.progressbar, length=350)
        progress.grid(column=0, row=0, columnspan=4, padx=20)
        progress_label = tk.Label(master=self.frame1, textvariable=self.progress_label, font=("Arial", 10))
        progress_label.grid(column=0, row = 1, rowspan=3, sticky='w', padx=20)
        eta_label = tk.Label(master=self.frame1, textvariable=self.eta_label, font=("Arial", 10))
        eta_label.grid(column=3, row = 1, rowspan=3, sticky='w')

        self.pygame_width_label = ttk.Label(self.frame2, text="Szerokość okna symulacji (px)")
        self.pygame_width_label.grid(column=0, row=0,sticky=tk.E, pady=5, padx=10)
        self.pygame_width_entry = ttk.Entry(self.frame2, textvariable=self.pygame_window_width, width=6)
        self.pygame_width_entry.grid(column=1, row=0)

        self.pygame_height_label = ttk.Label(self.frame2, text="Wysokość okna symulacji (px)")
        self.pygame_height_label.grid(column=0, row=1, sticky=tk.E, pady=5, padx=10)
        self.pygame_height_entry = ttk.Entry(self.frame2, textvariable=self.pygame_window_height, width=7)
        self.pygame_height_entry.grid(column=1, row=1)

        self.fps_label = ttk.Label(self.frame2, text="FPS")
        self.fps_label.grid(column=0, row=2, sticky=tk.E, pady=5, padx=10)
        self.fps_entry = ttk.Entry(self.frame2, textvariable=self.simulation_framerate, width=7)
        self.fps_entry.grid(column=1, row=2)

        self.nx_label = ttk.Label(self.frame2, text="Ilość komórek w osi X")
        self.nx_label.grid(column=0, row=3, sticky=tk.E, pady=5, padx=10)
        self.nx_entry = ttk.Entry(self.frame2, textvariable=self.grid_width_in_tails, width=7)
        self.nx_entry.grid(column=1, row=3)

        self.ny_label = ttk.Label(self.frame2, text="Ilość komórek w osi Y")
        self.ny_label.grid(column=0, row=4, sticky=tk.E, pady=5, padx=10)
        self.ny_entry = ttk.Entry(self.frame2, textvariable=self.grid_height_in_tails, width=7)
        self.ny_entry.grid(column=1, row=4)

        p1_label = ttk.Label(self.frame2, text="Prawdopodobieństwo p1")
        p1_label.grid(column=0, row=6, sticky=tk.E, pady=5, padx=10)
        p1_entry = ttk.Entry(self.frame2, textvariable=self.P1_probability_value, width=7)
        p1_entry.grid(column=1, row=6)

        p2_label = ttk.Label(self.frame2, text="Prawdopodobieństwo p2")
        p2_label.grid(column=0, row=7, sticky=tk.E, pady=5, padx=10)
        p2_entry = ttk.Entry(self.frame2, textvariable=self.P2_probability_value, width=7)
        p2_entry.grid(column=1, row=7)

        cykles_label = ttk.Label(self.frame2, text="Liczba cykli leczenia (N)")
        cykles_label.grid(column=0, row=8, sticky=tk.E, pady=5, padx=10)
        cykles_entry = ttk.Entry(self.frame2, textvariable=self.heal_cykles, width=7)
        cykles_entry.grid(column=1, row=8)

        self.iteration_label = ttk.Label(self.frame2, text="Liczba powtórzeń (I)")
        self.iteration_label.grid(column=0, row = 9, sticky=tk.E, pady=5, padx=10)
        self.iteration_entry = ttk.Entry(self.frame2, textvariable=self.iteration_limit, width=7)
        self.iteration_entry.grid(column=1, row=9)

        simcykles_label = ttk.Label(self.frame2, text="Limit kroków symulacji [T] (0 = inf.)")
        simcykles_label.grid(column = 0 , row = 10, sticky=tk.E,pady=5, padx=10)
        simcykles_entry = ttk.Entry(self.frame2, textvariable=self.simcykles, width=7)
        simcykles_entry.grid(column=1, row=10, sticky=tk.W, pady=5, padx=10)

        self.start_button = ttk.Button(self.frame3, text="START", command=self.on_process_start)
        self.start_button.grid(column=0, row=0)

        self.stop_button =ttk.Button(self.frame3, text="STOP", command=self.on_process_end, name="stop_button")
        self.stop_button.grid(column=1, row=0)

        self.close_button = ttk.Button(self.frame3, text="ZAKOŃCZ", command=self.on_close)
        self.close_button.grid(column=2, row=0)

        self.label_dry_run = ttk.Label(self.frame3, text="Bez wizualizacji")
        self.label_dry_run.grid(column=6, row=0, sticky=tk.W, padx=5)
        self.checkdryrun = ttk.Checkbutton(self.frame3, variable=self.dry_run, command=self.on_dry_run_checked_actions)
        self.checkdryrun.grid(column=7, row=0, sticky=tk.W, padx=5)


        self.frame4.columnconfigure(0, weight=2)
        self.frame4.columnconfigure(2, weight=2)
        self.frame4.columnconfigure(4, weight=2)
        self.frame4.columnconfigure(6, weight=2)
        lp1color = ttk.Label(self.frame4, text="P1")
        lp1color.grid(column=0, row=0, sticky='e')
        Cp1color = tk.Canvas(self.frame4, width=20, height=20, bg=self.rgb_to_hex(self.settings.P1_COLOR), relief=tk.SUNKEN)
        Cp1color.bind('<Button-1>', lambda event: self.changecolor(event, P1_COLOR))
        Cp1color.grid(column=1, row=0, padx=5)

        label_p2_color = ttk.Label(self.frame4, text="P2")
        label_p2_color.grid(column=2, row=0,sticky=tk.E)
        canvas_p2_color = tk.Canvas(self.frame4, width=20, height=20, bg=self.rgb_to_hex(self.settings.P2_COLOR), relief=tk.SUNKEN)
        canvas_p2_color.bind('<Button-1>', lambda event: self.changecolor(event, P2_COLOR))
        canvas_p2_color.grid(column=3, row=0, padx=5)

        lDcolor = ttk.Label(self.frame4, text="D")
        lDcolor.grid(column=4, row=0,sticky='e')
        CDcolor = tk.Canvas(self.frame4, width=20, height=20, bg=self.rgb_to_hex(self.settings.D_COLOR), relief=tk.SUNKEN)
        CDcolor.bind('<Button-1>', lambda event: self.changecolor(event, D_COLOR))
        CDcolor.grid(column=5, row=0, padx=5)

        lPDcolor = ttk.Label(self.frame4, text="PD")
        lPDcolor.grid(column=6, row=0,sticky=tk.E)
        CPDcolor = tk.Canvas(self.frame4, width=20, height=20, bg=self.rgb_to_hex(self.settings.PD_COLOR), relief=tk.SUNKEN)
        CPDcolor.bind('<Button-1>', lambda event: self.changecolor(event, PD_COLOR))
        CPDcolor.grid(column=7, row=0, padx=5)

        lGcolor = ttk.Label(self.frame4, text="P")
        lGcolor.grid(column=0, row=1,sticky=tk.E)
        CGcolor = tk.Canvas(self.frame4, width=20, height=20, bg=self.rgb_to_hex(self.settings.DP_COLOR), relief=tk.SUNKEN)
        CGcolor.bind('<Button-1>', lambda event: self.changecolor(event, DP_COLOR))
        CGcolor.grid(column=1, row=1, padx=5, sticky='w')

        lPcolor = ttk.Label(self.frame4, text="PS")
        lPcolor.grid(column=2, row=1,sticky=tk.E)
        CPcolor = tk.Canvas(self.frame4, width=20, height=20, bg=self.rgb_to_hex(self.settings.FP_COLOR), relief=tk.SUNKEN)
        CPcolor.bind('<Button-1>', lambda event: self.changecolor(event, FP_COLOR))
        CPcolor.grid(column=3, row=1, padx=5)

        lBcolor = ttk.Label(self.frame4, text="Kolor tła")
        lBcolor.grid(column=4, row=1,sticky=tk.E)
        CBcolor = tk.Canvas(self.frame4, width=20, height=20, bg=self.rgb_to_hex(self.settings.BACKGROUND_COLOR), relief=tk.SUNKEN)
        CBcolor.bind('<Button-1>', lambda event: self.changecolor(event, BACKGROUND_COLOR))
        CBcolor.grid(column=5, row=1, padx=5)

        self.fig = figure.Figure(figsize=(7,5), dpi=120)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.matplotframe)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.button_histplot_iterations = tk.Button(master=self.type_of_plots, text = "Hist. czasu", command=self.plot_hist_time)
        self.button_histplot_iterations.grid(column=0, row=0, pady=10, sticky='ew')
        self.button_barplot_iterations = tk.Button(master = self.type_of_plots, text = "Regresja liniowa", command=self.linear_regresion_plot)
        self.button_barplot_iterations.grid(column=1, row= 0, sticky="ew")
        self.button_histplot_PDpercentage = ttk.Button(master = self.type_of_plots, text="Hist. komórek PD", command=self.plot_hist_PD)
        self.button_histplot_PDpercentage.grid(column=2, row=0, sticky='ew')
        self.button_show_screencapture = ttk.Button(master=self.type_of_plots, text="Zapisz wyniki", command=self.save_data_from_current_run)
        self.button_show_screencapture.grid(column=5, row=1, sticky='ew')
        self.button_count_plot = ttk.Button(master=self.type_of_plots, text="Kierunki", command=self.direction_count)
        self.button_count_plot.grid(row=0, column=3)
        self.button_show_screencapture = ttk.Button(master=self.type_of_plots, text="Wynik wizualizacji", command=self.show_captured_image)
        self.button_show_screencapture.grid(column=5, row=0, sticky='ew')
      
        check_kde = ttk.Checkbutton(master=self.type_of_plots, variable=self.kde, onvalue=1, offvalue=0, text="Estymacja metodą KDE")
        check_kde.grid(column=7,row=0, sticky= 'w')

        check_kde = ttk.Checkbutton(master=self.type_of_plots, variable=self.pdf, onvalue=1, offvalue=0, text="Pokaż estymowany rozkład gęstości")
        check_kde.grid(column=7,row=1, sticky='w')

        check_toolbar = ttk.Checkbutton(master=self.type_of_plots, onvalue=1, offvalue=0, command=self.toggle_toolbar, variable=self.toolbar_check, text="Pokaż toolbar")
        check_toolbar.grid(column=7, row=2, sticky = 'w')
        self.plot_toolbar = NavigationToolbar2Tk(self.canvas, self.matplotframe, pack_toolbar=False)

        komunikacja = tk.Label(master=self.frame5, text="Praca inżynierska:\n autor: Marcin Ciachowski \n kontakt: m.ciachowski@gmail.com",
                              font=("Arial", 15)  )
        komunikacja.pack(expand=True, fill=tk.BOTH, side=tk.TOP, anchor=tk.CENTER)


        label_listbox= ttk.Label(master = self.type_of_plots, text="Typ histogramu")
        label_listbox.grid(column=0, row = 1, padx=5)   
        combobox = ttk.Combobox(master=self.type_of_plots, textvariable=self.combo, state="readonly")
        combobox.grid(column=1, row = 1, padx = 5, sticky='ew', columnspan=2)
        combobox["values"] = [hist_type_to_name[opt] for opt in hist_options]
        combobox.current(0)

        button_test = ttk.Button(self.type_of_plots, text="Testuj PRNG", command=self.prng_testing)
        button_test.grid(column=3, row=1, sticky="we")

        self.check_bins= ttk.Checkbutton(self.type_of_plots, onvalue=True, offvalue=False, variable=self.bin_check, text="Ręczny podział", command=self.on_bin_ckecked)
        self.check_bins.grid(column=0, row=2)

        bins_show_value = ttk.Label(self.type_of_plots, textvariable=self.bin_slider_value)
        bins_show_value.grid(column=3, row = 2, columnspan=3, sticky='w')
        self.bins_scale = ttk.Scale(self.type_of_plots, orient=tk.HORIZONTAL, length=300, variable=self.bins, from_=1, to=50,
            command= lambda x: self.bin_slider_value.set("Liczba przedziałów: " + str(self.bins.get())))
        self.bins_scale.grid(column=1, row=2, sticky = 'w', columnspan=2)

        boxplot_button_cykle = ttk.Button(master=self.type_of_plots, text="Boxplot", command=self.plot_boxplot_cykle)
        boxplot_button_cykle.grid(column=0, row=3, sticky='we')
        
        boxplot_button_pd = ttk.Button(master=self.type_of_plots, text="Boxplot_PD", command=self.plot_boxplot_pd)
        boxplot_button_pd.grid(column=1, row=3, sticky="we")

        boxplot_button_grid = ttk.Button(master=self.type_of_plots, text="Siatka cykle", command=self.plot_boxplot_grid_cykle)
        boxplot_button_grid.grid(column=2, row=3, sticky="we")

        boxplot_button_grid = ttk.Button(master=self.type_of_plots, text="Siatka PD", command=self.plot_boxplot_grid_pd)
        boxplot_button_grid.grid(column=3, row=3, sticky="we")

        check_medianline = ttk.Checkbutton(master=self.type_of_plots, onvalue=1, offvalue=0, variable=self.meanline, text="Połącz mediany")
        check_medianline.grid(column=4,row=3, columnspan=2)

        check_boxplot = ttk.Checkbutton(master=self.type_of_plots, onvalue=1, offvalue=0, variable=self.boxplot_show, text="Pokaż boxploty")
        check_boxplot.grid(column=6,row=3, columnspan=2, sticky='w')

        check_boxplot_fliers = ttk.Checkbutton(master=self.type_of_plots, onvalue=1, offvalue=0, variable=self.boxplot_show_fliers, text="Pokaż outlinerów")
        check_boxplot_fliers.grid(column=0,row=4, columnspan=2, sticky='w')
        
        



#----------------------------------------------------------------------------------------------------------------
#                                           ANALIZA DANYCH
#----------------------------------------------------------------------------------------------------------------

    def two_subplots(self, x: bool) -> list[plt.Axes]:

        for axes in self.fig.axes:
            self.fig.delaxes(axes)

        axes: plt.Axes
        if x:
            spec = gridspec.GridSpec(ncols=2, nrows=1, width_ratios=[6,1], hspace=0.1)
            self.fig.add_subplot(spec[0])
            self.fig.add_subplot(spec[1])

        else:
            spec = gridspec.GridSpec(ncols=1, nrows=1)
            self.fig.add_subplot(spec[0])

        axes = self.fig.axes
        return axes



    def analayze(self):
        """
        Extract data from Collector object to pandas DataFrame
        """
        column_names = ["X", "Y","P1_PROB", "P2_PROB", "N", "PD", "Cykle", "Time", "S/s", "NS", "WE"]
        converted_data = []

        for i in range(len(self.data_collected)):
            info, run, times = self.data_collected[i]

            size = info[4] * info[5]
            data = (
                info[4],               #X
                info[5],               #Y
                info[1],               #P1_prob
                info[2],               #P2_prob
                info[3],               #N
                run[3] * 100 / size,   #PD%
                run[4],                #cykle
                times[0],              #Time
                #times[1],
                times[2],              #Sim/s
                run[5],                #NS
                run[6]                 #WE
            )
            converted_data.append(data)

        self.dataframe = pd.DataFrame(converted_data, columns=column_names)

        self.dataframe_all = pd.concat([self.dataframe_all, self.dataframe])



    def plot_boxplot_cykle(self):
        ax = self.two_subplots(False)[0]
        ax.clear()
        df = self.dataframe_all
        uniq_n = df["N"].unique().tolist()
        if self.boxplot_show.get():
            sbn.boxplot(x ="N", y = "Cykle", data=df, ax=ax, showmeans=False, showfliers=False if not self.boxplot_show_fliers.get() else True)
        if self.meanline.get():
            sbn.pointplot(data=df.groupby("N", as_index=False).median(), x="N", y="Cykle", ax= ax, color="black")
        
        ax.get_yaxis().set_minor_locator(ticker.AutoMinorLocator())
        ax.get_yaxis().set_major_locator(ticker.AutoLocator())
        ax.set_xticklabels(sorted(uniq_n))
        ax.grid(visible=True, which="both", axis="both")
        self.canvas.draw()


    def plot_boxplot_pd(self):
        ax = self.two_subplots(False)[0]
        ax.clear()
        df = self.dataframe_all
        uniq_n = df["N"].unique().tolist()
        if self.boxplot_show.get():
            sbn.boxplot(x ="N", y = "PD", data=df, ax=ax, showmeans=False, showfliers=False if not self.boxplot_show_fliers.get() else True)
        if self.meanline.get():
            sbn.pointplot(data=df.groupby("N", as_index=False).median(), x="N", y="PD", ax= ax, color="black")

        ax.get_yaxis().set_minor_locator(ticker.AutoMinorLocator())
        ax.get_yaxis().set_major_locator(ticker.AutoLocator())
        ax.set_xticklabels(sorted(uniq_n))
        ax.grid(visible=True, which="both", axis="both")
        self.canvas.draw()


    def plot_boxplot_grid_cykle(self):
        ax = self.two_subplots(False)[0]
        ax.clear()
        df = self.dataframe_all
        uniq_x = df["X"].unique().tolist()
        if self.boxplot_show.get():
            sbn.boxplot(x ="X", y = "Cykle", data=df, ax=ax, showmeans=False, showfliers=False if not self.boxplot_show_fliers.get() else True)
        if self.meanline.get():
            sbn.pointplot(data=df.groupby("X", as_index=False).median(), x="X", y="Cykle", ax= ax, color="black")
        
        ax.get_yaxis().set_minor_locator(ticker.AutoMinorLocator())
        ax.get_yaxis().set_major_locator(ticker.AutoLocator())
        ax.set_xticklabels(sorted(uniq_x))
        ax.grid(visible=True, which="both", axis="both")
        self.canvas.draw()


    def plot_boxplot_grid_pd(self):
        ax = self.two_subplots(False)[0]
        ax.clear()
        df = self.dataframe_all
        uniq_x = df["X"].unique().tolist()
        if self.boxplot_show.get():
            sbn.boxplot(x ="X", y = "PD", data=df, ax=ax, showmeans=True, showfliers=False if not self.boxplot_show_fliers.get() else True)
        if self.meanline.get():
            sbn.pointplot(data=df.groupby("X", as_index=False).median(), x="X", y="PD", ax= ax, color="black")

        ax.get_yaxis().set_minor_locator(ticker.AutoMinorLocator())
        ax.get_yaxis().set_major_locator(ticker.AutoLocator())

        ax.set_xticklabels(sorted(uniq_x))
        ax.grid(visible=True, which="both", axis="both")
        self.canvas.draw()




    def direction_count(self):
        if hasattr(self, 'dataframe'):
            df = self.dataframe
            rt = df[['NS', 'WE']].replace({True: 1, False: 0})
            

            ax = self.two_subplots(False)[0]
            ax.clear()
            v = rt[['NS', 'WE']].value_counts()
            v.index = ["NS", "WE", "NS i WE"]

            ax_df = v.plot.barh(ax=ax, alpha=0.5, edgecolor="black")
            ax_df.tick_params(axis="x", which="minor", bottom="off")
            ax_df.set_title("Porównanie liczebności kierunków perkolacji")
            ax_df.bar_label(ax_df.containers[0])
            ax_df.set_xlabel("Liczba perkolacji")
            ax_df.set_ylabel("kierunek perkolacji")
            
            self.canvas.draw()
        else:
            print("Dataframe niezainicjowany. Przeprowadź serie symulacji żeby zebrać dane")

    

    def prng_uniformity_testing(self):
        from collections import namedtuple

        def chi2(classes: list[float|int], prob_size: int, intervals: int) -> list[float]:
            expected = prob_size / intervals
            return sum([ ((x - expected)**2 / expected) for x in classes])

        Chi2Test = namedtuple("chi_test", ["size", "chi2_cp"])

        #(wielkosc próby, punkt krytyczny ch^2 dla tej wielkosci proby)
        test1 = Chi2Test(100, 16.9190)
        test2 = Chi2Test(1000, 44.9853)
        test3 = Chi2Test(10000, 123.225)
        tests = [test1, test2, test3]
        mean = 0

        ITERATIONS = 100

        string = ("\nWeryfikacja rozkładu jednostajnego PRNG języka Python")
        string += '\n' + ("za pomocą testu Chi^2")
        string += '\n' + ("--------------------------------------------------------------------")
        test_time = time.time()
        for item in tests:
            result = []
            start = time.time()
            for i in range(ITERATIONS):
                #przygotuj serie liczb z generatora pseudoloswego
                sample = pd.Series([random.random() for _ in  range(item.size)])
                #oblicz ilosc przedzialow histogramu
                bins = math.floor(math.sqrt(item.size)) 
                #stwórz histogram licznosci
                hist = pd.cut(sample, bins=bins).value_counts().sort_index().tolist()
                #wykonaj test chi^2 i zapisz wynik
                c = chi2(hist, item.size, bins)
                result.append(True if c < item.chi2_cp else False)

            suma   = sum(result)
            passed = suma / ITERATIONS *100
            mean  += passed
            s1 = f'Wielkość próby: {item.size:5d}'
            s2 = f' Ilość iteracji:{ITERATIONS:5d}'
            s3 = f' Test przeszło: {passed:.1f} % prób ({suma}/{ITERATIONS})'
            s4 = f' Czas testu: {time.time() - start: .2f} s'
            string += '\n' + s1
            string += '\n' + s2
            string += '\n' + s3
            string += '\n' + s4

        s5 = f'Średnia z 3 próbek: {mean / 3:.2f} % testów zakończyło się pomyślnie'
        s6 = f'Czas testu: {time.time() - test_time: .3f} s'
        string += '\n' + "--------------------------------------------------------------------"
        string += '\n' + s5
        string += '\n' + s6

        return string


    def statisticks(self):
        if hasattr(self, 'dataframe'):
            pd = self.dataframe
            s = f'{pd[["Cykle", "PD"]].describe()}\n'

            title = f'-------------------------------------\n'
            p1 = self.P1_probability_value.get()
            p2 = self.P2_probability_value.get()
            n = self.heal_cykles.get()
            t = self.simcykles.get()
            y = self.grid_height_in_tails.get()
            x= self.grid_width_in_tails.get()
            subtitle = f'Params: X({x})x Y({y}); P1({p1}); P2({p2}); N({n}); T({t})\n'
            self.infopage.insert(tk.END, title + subtitle + s )


    def prng_testing(self):

        uniformity_test = self.prng_uniformity_testing()
        self.infopage.insert(tk.END, uniformity_test)
        
        n = 1000
        ax = self.two_subplots(False)[0]
        ax.clear()

        s = np.array([random.randint(0, 100) for _ in range(n)])
        plot_acf(s, ax=ax, lags=100)

        ax.set_ylim(-0.1, 0.1)
        ax.set_title("Autokorelacja dla " + str(n) + " liczb z genratora PRNG")
        ax.set_ylabel("Współczynnik korelacji")
        ax.set_xlabel("Przesunięcie (lag)")
        self.canvas.draw()
        


    def fit_test(self):
        """
        Finds best fitting probabilistic density function (PDF)

        Based on the code found on:
        https://medium.com/@amirarsalan.rajabi/distribution-fitting-with-python-scipy-bb70a42c0aed
        https://stackoverflow.com/questions/6620471/fitting-empirical-distribution-to-theoretical-ones-with-scipy-python
        """
        df = self.dataframe
        p_value = 1
        statistic = 0

        seria_iter = df["Cykle"]
        seria_pd =   df["PD"]
        result_tab_iter = []
        result_tab_pd = []

        self.progressbar.set(0)
        self.progress_label.set("Dopasowywanie PDF [0/" + str(len(pdf_name)))
  
        for i, dist in enumerate(pdf_name):
            p = math.floor(i * 100 / len(pdf_name))
            self.progress_label.set('Dopasowywanie funkcji gęstości ('+ dist + ')')
            self.eta_label.set(str(i) + "/" + str(len(pdf_name)))
            self.progressbar.set(p)
            pdf = getattr(sci, dist)
        
            parameters_iter = pdf.fit(seria_iter)
            parameters_pd   = pdf.fit(seria_pd)

            result_iter = sci.kstest(seria_iter, dist, args=parameters_iter)
            result_tab_iter.append((dist, result_iter[statistic], result_iter[p_value], pdf, parameters_iter))

            result_pd = sci.kstest(seria_pd, dist, args=parameters_pd)
            result_tab_pd.append((dist, result_pd[statistic], result_pd[p_value], pdf, parameters_pd))

        result_tab_iter.sort(key = lambda x: float(x[2]), reverse = True)
        result_tab_pd.sort(key   = lambda x: float(x[2]), reverse = True)

        for i in range(2):
            item_pd = result_tab_pd[i]
            item_iter = result_tab_iter[i]

            iter_name = item_iter[0]
            pd_name = item_pd[0]

            pd_arg, iter_arg     = item_pd[4][:-2], item_iter[4][:-2]
            pd_loc, iter_loc     = item_pd[4][-2],  item_iter[4][-2]
            pd_scale, iter_scale = item_pd[4][-1],  item_iter[4][-1]

            dist_pd = item_pd[3]
            dist_iter = item_iter[3]

            start_pd = dist_pd.ppf(0.0001, *pd_arg, loc=pd_loc, scale=pd_scale)
            end_pd = dist_pd.ppf(0.999, *pd_arg, loc=pd_loc, scale=pd_scale)

            x = np.linspace(start_pd,end_pd,10000)
            y = dist_pd.pdf(x, loc=pd_loc, scale=pd_scale, *pd_arg)
            pdf_function = pd.Series(y, x)
            self.best_fitting_pdf_pd.append((pd_name, pdf_function))

            start_iter = dist_iter.ppf(0.0001, *iter_arg, loc=iter_loc, scale=iter_scale)
            end_iter   = dist_iter.ppf(0.999, *iter_arg, loc=iter_loc, scale=iter_scale)

            x2 = np.linspace(start_iter,end_iter,10000)
            y2 = dist_iter.pdf(x2, loc=iter_loc, scale=iter_scale, *iter_arg)
            pdf_function2 = pd.Series(y2, x2)
            self.best_fitting_pdf_iter.append((iter_name,pdf_function2))



    def plot_subtitle(self):
        """
        Creates plot subtitle containing simulation parameters
        """

        op1 = "Parametry symulacji: "
        op2 = f'P1: {self.P1_probability_value.get()}, '
        op3 = f'P2: {self.P2_probability_value.get()}'
        op4 = f'  Siatka ({self.grid_width_in_tails.get()} x {self.grid_height_in_tails.get()})'
        op5 = f' N = {self.heal_cykles.get()}, T = {self.simcykles.get() if self.simcykles.get() > 0 else "inf"}'

        return  op1 + op2 + op3 + op4 + op5

        

    def plot_hist_time(self, show = True, kde_show = False, pdf_show = False, stat = None, discrete = False):

        df = self.dataframe["Cykle"]
        # ax = self.fig.gca()
        # ax.clear()
        axes = self.two_subplots(True)
        boxplot_ax = axes[1]
        hist_ax = axes[0]

        if self.bin_check.get():
            print(self.bins.get())
            bins = self.bins.get()
            p = sbn.histplot(
                df,
                ax=hist_ax, 
                kde = True if self.kde.get() or kde_show else False,
                stat = hist_name_to_type[self.combo.get()] if not stat else stat,
                alpha=0.5,
                edgecolor="black",
                legend=True,
                discrete=discrete,
                bins= bins
            )
        else:
            p = sbn.histplot(
                df,
                ax=hist_ax, 
                kde = True if self.kde.get() or kde_show else False,
                stat = hist_name_to_type[self.combo.get()] if not stat else stat,
                alpha=0.5,
                edgecolor="black",
                legend=True,
                discrete=discrete
            )


        boxplot = sbn.boxplot(df, ax=boxplot_ax)

        if not stat:
            p.set(xlabel='$t_p$', 
                 ylabel=ylabel[hist_name_to_type[self.combo.get()]])
        elif stat =="density":
             p.set(xlabel='$t_p$', 
                 ylabel="Prawdopodobieństwo")

        elif stat =="count":
             p.set(xlabel='$t_p$', 
                 ylabel="Liczba symulacji")

        hist_ax.get_xaxis().set_minor_locator(ticker.AutoMinorLocator())
        hist_ax.get_yaxis().set_minor_locator(ticker.AutoMinorLocator())
        hist_ax.get_xaxis().set_major_locator(ticker.AutoLocator())
        hist_ax.get_yaxis().set_major_locator(ticker.AutoLocator())
        hist_ax.grid(visible=True, which="both", axis="both")

        boxplot_ax.get_xaxis().set_minor_locator(ticker.AutoMinorLocator())
        boxplot_ax.get_yaxis().set_minor_locator(ticker.AutoMinorLocator())
        boxplot_ax.get_xaxis().set_major_locator(ticker.AutoLocator())
        boxplot_ax.get_yaxis().set_major_locator(ticker.AutoLocator())
        boxplot_ax.grid(visible=True, which="both", axis="both")
        self.fig.suptitle(self.plot_subtitle(), fontsize=10)
        #hist_ax.set_title(self.plot_subtitle(), fontsize=8)

        if self.pdf.get() or pdf_show:
            for i, it in enumerate(self.best_fitting_pdf_iter):
                name, pdf = it
                hist_ax.plot(pdf, label=name + "(best)" if i == 0 else name)
                hist_ax.legend()
        if show:
            self.canvas.draw()



    def linear_regresion_plot(self, show = True):
        df = self.dataframe
        ax = self.two_subplots(False)[0]
        ax.clear()

        slope, intercept , p_value, r_value, std_err = sci.stats.linregress(df['Cykle'], df['PD'])
        
        sbn.regplot(data=df, x="Cykle", y="PD", ax=ax,
                    line_kws={'label':"y={0:.1f}x+{1:.1f}".format(slope,intercept)})
        #self.fig.suptitle("Wykres regresji liniowej", fontsize=15)
        ax.set_title(self.plot_subtitle(), fontsize=10)

        ax.get_xaxis().set_minor_locator(ticker.AutoMinorLocator())
        ax.get_yaxis().set_minor_locator(ticker.AutoMinorLocator())
        ax.get_xaxis().set_major_locator(ticker.AutoLocator())
        ax.get_yaxis().set_major_locator(ticker.AutoLocator())
        ax.grid(visible=True, which="both", axis="both")
        ax.set_xlabel('$t_p$', fontsize=14)
        ax.set_ylabel("Udział komórek PD [%]", fontsize=12)
        ax.legend()

        if show:
            self.canvas.draw()
            s =f'Parametry regresji Liniowej:\n Nachylenie (slope): {slope}\n Wyraz wolny (intercept): {intercept}\n'
            s2 = f'p_value: {p_value} \n r_value: {r_value}, \n bład standardowy {std_err}'
            self.infopage.insert(tk.END, s + s2)



    def plot_hist_PD(self,show = True,
                                kde_show = False,
                                pdf_show = False,
                                stat=None,
                                discrete = False):
        df = self.dataframe["PD"]

        hist_ax, box_ax = self.two_subplots(True)

        if self.bin_check.get():
            bins = self.bins.get()
            p = sbn.histplot(
                df,
                ax=hist_ax, 
                kde = True if self.kde.get() or kde_show else False,
                stat = hist_name_to_type[self.combo.get()] if not stat else stat,
                alpha=0.5,
                edgecolor="black",
                color='purple',
                legend=True,
                discrete=discrete,
                bins = bins
            )
        else:
            p = sbn.histplot(
                df,
                ax=hist_ax, 
                kde = True if self.kde.get() or kde_show else False,
                stat = hist_name_to_type[self.combo.get()] if not stat else stat,
                alpha=0.5,
                edgecolor="black",
                color='purple',
                legend=True,
                discrete=discrete
            )

        if not stat:
            p.set(xlabel='PD [%]', 
                 ylabel=ylabel[hist_name_to_type[self.combo.get()]])
        elif stat =="density":
             p.set(xlabel='PD [%]', 
                 ylabel="Prawdopodobieństwo")

        elif stat =="count":
             p.set(xlabel='PD [%]', 
                 ylabel="Liczba symulacji")


        # ax.set_xlabel('$t_p$', fontsize=14)
        # ax.set_ylabel("Udział komórek PD [%]", fontsize=12)

        # ax.set_xlabel("Udział komórek PD [%]")
        # ax.set_ylabel(ylabel[hist_name_to_type[self.combo.get()]])

        hist_ax.get_xaxis().set_minor_locator(ticker.AutoMinorLocator())
        hist_ax.get_yaxis().set_minor_locator(ticker.AutoMinorLocator())
        hist_ax.get_xaxis().set_major_locator(ticker.AutoLocator())
        hist_ax.get_yaxis().set_major_locator(ticker.AutoLocator())
        hist_ax.grid(visible=True, which="both", axis="both")
        self.fig.suptitle(self.plot_subtitle(), fontsize=10)
        #ax.set_title(self.plot_subtitle(), fontsize=8)

        if self.pdf.get() or pdf_show:
            for i, it in enumerate(self.best_fitting_pdf_pd):
                name, pdf = it
                hist_ax.plot(pdf, label=name + "(best)" if i == 0 else name)
                hist_ax.legend()


        box_ax.get_xaxis().set_minor_locator(ticker.AutoMinorLocator())
        box_ax.get_yaxis().set_minor_locator(ticker.AutoMinorLocator())
        box_ax.get_xaxis().set_major_locator(ticker.AutoLocator())
        box_ax.get_yaxis().set_major_locator(ticker.AutoLocator())
        box_ax.grid(visible=True, which="both", axis="both")
        sbn.boxplot(df, ax=box_ax)

        if show:
            self.canvas.draw()




    def show_captured_image(self, show = True):
        '''
        Wyświetla ostatnio zapisany wynik symulacji z wizualizacją.
        
        Parameters:
        ----------
        show:bool (True)
            Określa czy wyświetlić obraz na płótnie. Jeżeli nie obraz jest
            tylko wczytywany ale niewyświetlany. 
            Domyślnie = True
        '''

        list_of_files = glob.glob('*.png')
        latest_file   = max(list_of_files, key=os.path.getctime)
        img           = pltimage.imread(latest_file)

        a = self.fig.gca()

        #Wyłącz siatke, wyłącz pokazywanie osi
        a.grid(visible=False)
        a.clear()
        a.set_axis_off()
        a.imshow(img)
        #wyświetl  tytuły wykresu
        self.fig.suptitle(" ", fontsize=5)
        a.set_title(self.plot_subtitle(), fontsize= 10)

        if show:
            self.canvas.draw()



    def save_data_from_current_run(self):
        '''
        Zapisuję wyekstrahowane dane do pliku CSV, obrazy wszystkich wykresów
        do pliku ZIP
        '''

        #format pliku z wynikami badania:
        #badanie_YYYYMMDD-HMS.zip
        #trzeba go zapisać w folderze badania w katalogu uruchowminiowym programu
        #niestety w linuksie używa się slasha (/) a w windowsie backslasha (\)
        if platform.system() == "Linux":
            date = time.strftime("%Y%m%d-%H%M%S")
            name = f'badanie_{date}'
            path = (os.path.realpath(os.path.dirname(__file__)))
            experiment_dir = f'{path}/Badania'
            new_exam = f'{experiment_dir}/{name}'
            badanie= new_exam + f'/{name}'

        elif platform.system() == "Windows":
            pass


        os.mkdir(new_exam)
        self.dataframe.to_csv(badanie + '.csv', index=True, sep=';')
        if not self.dry_run.get():
            self.show_captured_image(False)
            self.fig.savefig(badanie + '_IMAGE.png', dpi=200, format="png", bbox_inches = 'tight')
            self.plot_hist_PD(False)
            self.fig.savefig(badanie + '_HISTPD.png', dpi=200, format="png", bbox_inches = 'tight')
            self.plot_hist_time(False)
            self.fig.savefig(badanie+ '_HISTCYKLE.png', dpi=200, format="png", bbox_inches = 'tight')
            
        else:
            self.plot_hist_PD(False, stat="count")
            self.fig.savefig(badanie + '_HISTPD.png', dpi=200, format="png", bbox_inches = 'tight')
            self.linear_regresion_plot(False)
            self.fig.savefig(badanie + '_LINREG.png', dpi=200, format="png", bbox_inches = 'tight')
            self.plot_hist_time(False, stat="count")
            self.fig.savefig(badanie+ '_HISTCYKLE.png', dpi=200, format="png", bbox_inches = 'tight')
            self.plot_hist_PD(False, True, stat="count")
            self.fig.savefig(badanie+ '_HISTPD+KDE.png', dpi=200, format="png", bbox_inches = 'tight')
            self.plot_hist_time(False, True, stat="count")
            self.fig.savefig(badanie + '_HISTCYKLE+KDE.png', dpi=200, format="png", bbox_inches = 'tight')
            self.plot_hist_PD(False, False, True, "density")
            self.fig.savefig(badanie+ '_HISTPD+PDF.png', dpi=200, format="png", bbox_inches = 'tight')
            self.plot_hist_time(False, False, True, "density")
            self.fig.savefig(badanie + '_HISTCYKLE+PDF.png', dpi=200, format="png", bbox_inches = 'tight')
            self.plot_hist_PD(False, True, True, "density")
            self.fig.savefig(badanie+ '_HISTPD+KDE+PDF.png', dpi=200, format="png", bbox_inches = 'tight')
            self.plot_hist_time(False, True,True, "density")
            self.fig.savefig(badanie + '_HISTCYKLE+KDE+PDF.png', dpi=200, format="png", bbox_inches = 'tight')


        op2 = f'P1({self.P1_probability_value.get()})'
        op3 = f'P2({self.P2_probability_value.get()})'
        op4 = f'G({self.grid_width_in_tails.get()}x{self.grid_height_in_tails.get()})'
        op5 = f'N({self.heal_cykles.get()})T({self.simcykles.get() if self.simcykles.get() > 0 else "inf"})'

        arch_name = experiment_dir + '/' + name + "--" + op2 + op3 + op4 + op5

        #Pusty plik. W nazwie pliku zakodowane są parmetry symulacji
        with open(new_exam + '/' + op2 + op3 + op4 + op5, 'w') as f:
            pass

        #Stwórz archiwum ze skompresowaną zawartością katalogu z badaniami
        shutil.make_archive(arch_name, format='zip', root_dir=new_exam)



            
        
        
#------------------------       KONIEC AALIZY DANYCH     ----------------------------------------------------

#--------------------------------------------------------------------------------------------
#                                    MULIPROCESSING
#--------------------------------------------------------------------------------------------
    def thread_worker(self, pipe: Pipe, process:Process, queue:Queue):
        """
        Wątek sprawdza cyklicznie czy symulacja się zakończyła ( czy istnieje process
        który ją prowadzi), aktualizuję progressbar, sprawdza czy 
        proces wysłął dane. Wątek jest konieczny aby nie blokować GUI"""

        #zmienna która odbiera stan symulacji
        d: tuple
        run_thread = True

        labstart = f'0 % [0 / {self.iteration_limit.get()}]'
        self.progress_label.set(labstart)
        self.best_fitting_pdf_iter = []
        self.best_fitting_pdf_pd = []

        try:
            self.simulation_running = True
            self.interface_lock(True)
            while process and process.is_alive() or run_thread:
                try:
                   d = queue.get(block = False)
                except:
                    pass
                else:
                    self.progressbar.set(d[2])
                    progress = f'{d[2]} % [{d[0]}/{d[1]} ({d[3]:.2f} [Symulacji/s])]'
                    eta = math.floor(d[4])
                    self.progress_label.set(progress)
                    self.eta_label.set('ETA: ' + str(datetime.timedelta(seconds=eta)))

                if pipe.poll():
                    # print("PROCESS SENT DATA")
                    self.data_collected: list = pipe.recv()
                    if len(self.data_collected) > 0:
                        self.simulation_ends_sucesfuly = True
                    else:
                        self.simulation_ends_sucesfuly = False
                    
                    # print("PROCESS ENDED")

                    self.analayze()
                    # if self.dry_run.get():
                    #     self.fit_test()
                    #     self.unlock_plot_buttons(True)
                    # else:
                    #     self.unlock_plot_buttons(False)

                    self.statisticks()
                    # zakończ wątek

                    run_thread = False
            else:
                # print("PROCESS TERMINATED BY USER")
                self.simulation_ends_sucesfuly = False
        except EOFError:
             print("EOF")
            
        finally:
            pipe.close()
            self.simulation_running = False
            self.interface_lock(False)
            self.progressbar.set(100)
            self.progress_label.set("100 %")

                

    def on_process_start(self):
        """
        proces uruchamia symulacje, a wątek sprawdza cyklicznie czy symulacja sie skończyła, 
        sprawdzanie jest w œatku
        bo inaczej zablokowałbym całe gui czekając na Pipe.recv()
        """
        
        self.update_settings()
        self.simulation_ends_sucesfuly = False
            
        parent_conn, child_conn = Pipe()
        queue = Queue()
        self.process = Process(target=run_simulation, args=(self.settings, child_conn, queue))
        self.thread  = Thread (target=self.thread_worker, args=(parent_conn, self.process, queue))

        self.process.start()
        child_conn.close()
        self.thread.start()



    def on_process_end(self):
        '''
        Kończy process jeżeli wciąż istnieje
        '''

        if self.process and self.process.is_alive():
            self.process.terminate()


#------------------------------------ MULTIPROCESSNG END --------------------------------------------

#----------------------------------------------------------------------------------------------------
#                                MANIPULACJA INTERFEJSEM
#----------------------------------------------------------------------------------------------------

    def on_bin_ckecked(self, *args):
        if  self.bin_check.get():
            self.bins_scale.configure(state='normal')
        else:
            self.bins_scale.configure(state="disable")


    def toggle_toolbar(self):
        if self.toolbar_check.get():
            self.window.geometry("1270x800")
            self.plot_toolbar.pack()
        else:
            self.window.geometry("1270x800")
            self.plot_toolbar.pack_forget()



    def unlock_plot_buttons(self, unlock):
        state = "disable"
        if unlock:
            state = "normal"
  
        self.button_barplot_iterations.configure(state = state)
        self.button_histplot_iterations.configure(state = state)
        self.button_histplot_PDpercentage.configure(state = state)
    
    
    
    def changecolor(self, event, name):
        color = askcolor(initialcolor=event.widget['bg'])
        print(color, True if color else False)

        if color[0] and color[1]:
            event.widget['bg'] = color[1]
            if name == P1_COLOR:
                self.settings.P1_COLOR = color[0]

            elif name == P2_COLOR:
                self.settings.P2_COLOR = color[0]

            elif name == D_COLOR:
                self.settings.D_COLOR = color[0]

            elif name == PD_COLOR:
                self.settings.PD_COLOR = color[0]

            elif name == FP_COLOR:
                self.settings.FP_COLOR = color[0]

            elif name == DP_COLOR:
                self.settings.DP_COLOR = color[0]

            elif name == BACKGROUND_COLOR:
                self.settings.BACKGROUND_COLOR = color[0]



    def interface_lock(self, lock = True):
        state = "disable"
        self._widget_states = []

        if not lock:
            state = "normal"

        for child in self.frame1.winfo_children():
            if not isinstance(child, ttk.Progressbar) and not isinstance(child, tk.Label):
                child.configure(state = state)

        for child in self.frame2.winfo_children():
            child.configure(state = state)

        for child in self.frame3.winfo_children():
            child.configure(state = state)

        for child in self.frame4.winfo_children():
            child.configure(state = state)

        for child in self.type_of_plots.winfo_children():
            child.configure(state = state)

        
        if self.simulation_running:
            self.stop_button.configure(state = "normal")
        else:
            self.stop_button.configure(state = "disable")
            if self.dry_run.get():
                self.on_dry_run_checked_actions()


        self.on_bin_ckecked()




    def on_dry_run_checked_actions(self):
    
        state = 'normal'
        if self.dry_run.get():
            state = 'disable'
        else:
            self.iteration_limit.set(1)

        self.pygame_height_entry.configure(state=state)
        self.pygame_height_label.configure(state=state)
        self.pygame_width_entry.configure(state=state)
        self.pygame_width_label.configure(state=state)
        self.pygame_height_label.configure(state=state)
        self.fps_entry.configure(state=state)
        self.fps_label.configure(state=state)
        self.button_show_screencapture.configure(state=state)
        self.iteration_label.configure(state= "disable" if not self.dry_run.get() else "normal")
        self.iteration_entry.configure(state= "disable" if not self.dry_run.get() else "normal")
        self.button_barplot_iterations.configure(state= "disable" if not self.dry_run.get() else "normal")
        self.button_histplot_PDpercentage.configure(state= "disable" if not self.dry_run.get() else "normal")
        self.button_histplot_iterations.configure(state= "disable" if not self.dry_run.get() else "normal")
        


    def run(self):
        self.window.mainloop()



    def on_close(self):
        """
        Actions that have to be taken before 
        we can close progrram
        """
        self.on_process_end()
        self.update_settings()
        self.settings.save_to_file('settings.json')
        self.window.destroy()




if __name__ == "__main__":
    app = appGUI()
    app.run()