"""
A PyQt5 gui for exploring all the expected fluxes I've made 
"""

# Gui stuff
from PyQt5.QtWidgets import QMainWindow, QWidget, QApplication, QFileDialog
from PyQt5 import QtCore, QtGui, QtWidgets

# embedding matplotlib stuff
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
import matplotlib.cm as cm


import sys # exit
from math import pi # everyone needs pi
import numpy as np
import pickle # used to load in the fluxes
import os
from glob import glob

from cascade.utils import SterileParams, gen_filename

def app_zero(obj):
    return np.concatenate(( np.array([0]), obj))

class base_gui(object):
    """
    implements the UI for the window

    This was generated in Qt Designer, so it's kinda messy and gross 
    """
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(800, 600)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.centralwidget)
        self.verticalLayout.setObjectName("verticalLayout")
        
        # I added this, let's put a matplotlib figure in here!
        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, MainWindow)
        self.verticalLayout.addWidget(self.toolbar)
        self.verticalLayout.addWidget(self.canvas)       

        self.formLayout = QtWidgets.QFormLayout()
        self.formLayout.setFormAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTop|QtCore.Qt.AlignTrailing)
        self.formLayout.setObjectName("formLayout")
        self.msq_slider = QtWidgets.QSlider(self.centralwidget)
        self.msq_slider.setMaximum(40)
        self.msq_slider.setOrientation(QtCore.Qt.Horizontal)
        self.msq_slider.setObjectName("msq_slider")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.msq_slider)
        self.msq_lbl = QtWidgets.QLabel(self.centralwidget)
        self.msq_lbl.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.msq_lbl.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.msq_lbl.setObjectName("msq_lbl")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.msq_lbl)

        self.muon_lbl = QtWidgets.QLabel(self.centralwidget)
        self.muon_lbl.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.muon_lbl.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.muon_lbl.setObjectName("muon_lbl")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.muon_lbl)
        self.muon_slider = QtWidgets.QSlider(self.centralwidget)
        self.muon_slider.setOrientation(QtCore.Qt.Horizontal)
        self.muon_slider.setMaximum(90)
        self.muon_slider.setObjectName("muon_combo")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.muon_slider)


        self.tau_lbl = QtWidgets.QLabel(self.centralwidget)
        self.tau_lbl.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.tau_lbl.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.tau_lbl.setObjectName("tau_lbl")
        self.formLayout.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.tau_lbl)
        self.tau_slider = QtWidgets.QSlider(self.centralwidget)
        self.tau_slider.setOrientation(QtCore.Qt.Horizontal)
        self.tau_slider.setObjectName("tau_slider")
        self.tau_slider.setMaximum(90)
        self.formLayout.setWidget(2, QtWidgets.QFormLayout.FieldRole, self.tau_slider)

        self.width_lbl = QtWidgets.QLabel(self.centralwidget)
        self.width_lbl.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.width_lbl.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.width_lbl.setObjectName("width_lbl")
        self.formLayout.setWidget(3, QtWidgets.QFormLayout.LabelRole, self.width_lbl)
        self.width_slider = QtWidgets.QSlider(self.centralwidget)
        self.width_slider.setOrientation(QtCore.Qt.Horizontal)
        self.width_slider.setObjectName("width_slider")
        self.width_slider.setMinimum(1)
        self.width_slider.setMaximum(100)
        self.formLayout.setWidget(3, QtWidgets.QFormLayout.FieldRole, self.width_slider)
        
        self.verticalLayout.addLayout(self.formLayout)
        MainWindow.setCentralWidget(self.centralwidget)
      
       	self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 22))
        self.menubar.setObjectName("menubar")
        self.menuFile = QtWidgets.QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")
        self.menuFlavors = QtWidgets.QMenu(self.menubar)
        self.menuFlavors.setObjectName("menuFlavors")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.actionQuit = QtWidgets.QAction(MainWindow)
        self.actionQuit.setObjectName("actionQuit")
        self.menuFile.addAction(self.actionQuit)
         
        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

        self.msq_slider.setTickInterval(1)
        self.tau_slider.setTickInterval(1)
        self.width_slider.setTickInterval(1)
        self.width_slider.setValue(40)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "Flux Explorer"))
        self.msq_lbl.setText(_translate("MainWindow", "Delta M:  0.00"))
        self.muon_lbl.setText(_translate("MainWindow", "Theta-24: 0.00"))
        self.tau_lbl.setText(_translate("MainWindow", "Theta-34: 0.00"))
        self.width_lbl.setText(_translate("MainWindow", "Width: 0.10"))
        self.menuFile.setTitle(_translate("MainWindow", "File"))
        self.actionQuit.setText(_translate("MainWindow", "Quit"))



class main_window(QMainWindow):
    """
    Main function that creates the window 
    """
    def __init__(self, parent=None):
        QWidget.__init__(self, parent)
        # set up the UI
        self.ui = base_gui()
        self.ui.setupUi(self)

        self.folder = "/home/bsmithers/software/data/hg_sib/expected_fluxes_reco"
        self.filenames = "expected_flux_smearedwell.dat"

        # whenever the sliders value change, we update the plots 
        # may want to change this to (mouse released) if this takes a while 

        self.ui.tau_slider.valueChanged.connect(self.update_plot)
        self.ui.muon_slider.valueChanged.connect(self.update_plot)
        self.ui.msq_slider.valueChanged.connect(self.update_plot)
        self.ui.width_slider.valueChanged.connect(self.update_width)

        self.core_b = -0.98
        self.mantle_b = -0.83


        self.tau_a = 0.0
        self.muon_a = 0.0
        self.msq = 0.0

        # copied from the make_dag python script
        # this just ensures that the names are all right 
        n_grid = 90
        self.theta13s = app_zero(np.arcsin(np.sqrt(np.logspace(-3,0, n_grid)))/2)
        self.theta23s = app_zero(np.arcsin(np.sqrt(np.logspace(-3,0, n_grid)))/2)
        self.msqs = app_zero(np.logspace(-2,2,40))

        self._n_fluxes = len(self.theta13s)*len(self.msqs)*len(self.theta23s)

        

        # load the null flux!
        null_dict = self.load(SterileParams())

        self.en = null_dict["e_edges"]
        self.an = np.linspace(-1,1,11)

        self.width = 0.50

        self.update_plot()

    def reload_null(self):
        pass

    def update_labels(self):
        self.muon_a = self.theta13s[self.ui.muon_slider.value()]
        self.tau_a = self.theta23s[self.ui.tau_slider.value()]
        self.msq  = self.msqs[self.ui.msq_slider.value()]

        self.ui.muon_lbl.setText("Theta-24: {:.4f}, {:.2f}".format(self.muon_a, self.muon_a*180/pi))
        self.ui.tau_lbl.setText("Theta-34: {:.4f}, {:.2f}".format(self.tau_a, self.tau_a*180/pi))
        self.ui.msq_lbl.setText("Delta M: {:.4f}".format(self.msq))

    def load(self, param):
        filename = gen_filename(self.folder, self.filenames, param)
        f = open(filename, 'rb')
        data = pickle.load(f)
        f.close()
        return data

    def update_plot(self):
        self.update_labels()

        null = SterileParams()
        this_s = SterileParams(theta13=self.muon_a, theta23=self.tau_a, msq2=self.msq)

        null_dat = self.load(null)
        ster_dat = self.load(this_s)

        self.ui.figure.clear()
        ax = self.ui.figure.add_subplot(111)

        pmesh = ax.pcolormesh(self.an, self.en, ster_dat["event_rate"]/null_dat["event_rate"], vmin=1.0-self.width, vmax=1.0+self.width, cmap="coolwarm")
        ax.set_yscale('log')
        ax.set_ylim([1e2,1e6])
        ax.set_xlim([-1,0.2])
        self.cbar = self.ui.figure.colorbar(pmesh, ax=ax)
        self.cbar.set_label("Difference")

        plt.vlines(self.core_b,ymin=10**2, ymax=10**6, colors="white", ls="-")
        plt.vlines(self.mantle_b,ymin=10**2, ymax=10**6, colors="white", ls="--")

        self.ui.canvas.draw()


    def update_width(self):
        pass

app = QApplication(sys.argv)
app_instance = main_window()
if __name__=="__main__":
    app_instance.show()
    sys.exit(app.exec_())
