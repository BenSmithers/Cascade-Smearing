# -*- coding: utf-8 -*-
"""
Ben Smithers

This is the cool script! It opens up a GUI allowing the user to explore theta_14 and theta_34 space 
"""
# Gui stuff 
from PyQt5.QtWidgets import QMainWindow, QWidget, QApplication
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

from cascade.utils import config # base configuration file, we need this to get the location of hte data

# these functions are used for the interpolating 
from cascade.utils import get_loc, bilinear_interp
from cascade.utils import SterileParams, gen_filename
from cascade.utils import bhist 
from cascade.utils import Data

from cascade.cross_section_test import get_total_flux as get_xs
from cascade.nus_utils import get_flavor, get_neut, get_curr 

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
        self.electron_slider = QtWidgets.QSlider(self.centralwidget)
        self.electron_slider.setMaximum(100)
        self.electron_slider.setOrientation(QtCore.Qt.Horizontal)
        self.electron_slider.setObjectName("electron_slider")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.electron_slider)
        self.electron_lbl = QtWidgets.QLabel(self.centralwidget)
        self.electron_lbl.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.electron_lbl.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.electron_lbl.setObjectName("electron_lbl")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.electron_lbl)

        self.muon_lbl = QtWidgets.QLabel(self.centralwidget)
        self.muon_lbl.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.muon_lbl.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.muon_lbl.setObjectName("muon_lbl")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.muon_lbl)
        self.muon_combo = QtWidgets.QComboBox(self.centralwidget)
        self.muon_combo.setObjectName("muon_combo")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.muon_combo)


        self.tau_lbl = QtWidgets.QLabel(self.centralwidget)
        self.tau_lbl.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.tau_lbl.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.tau_lbl.setObjectName("tau_lbl")
        self.formLayout.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.tau_lbl)
        self.tau_slider = QtWidgets.QSlider(self.centralwidget)
        self.tau_slider.setOrientation(QtCore.Qt.Horizontal)
        self.tau_slider.setObjectName("tau_slider")
        self.tau_slider.setMaximum(100)
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

        self.recoBox = QtWidgets.QCheckBox(self.centralwidget)
        self.recoBox.setObjectName("recoBox")
        self.recoBox.setChecked(False)
        self.verticalLayout.addLayout(self.formLayout)
        self.verticalLayout.addWidget(self.recoBox)
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

        self.actionElectron = QtWidgets.QAction(MainWindow)
        self.actionElectron.setCheckable(True)
        self.actionElectron.setChecked(True)
        self.actionElectron.setObjectName("actionElectron")
        self.actionNuEBar = QtWidgets.QAction(MainWindow)
        self.actionNuEBar.setCheckable(True)
        self.actionNuEBar.setChecked(True)
        self.actionNuEBar.setObjectName("actionNuEBar")
        self.actionNuMu = QtWidgets.QAction(MainWindow)
        self.actionNuMu.setCheckable(True)
        self.actionNuMu.setChecked(True)
        self.actionNuMu.setObjectName("actionNuMu")
        self.actionNuMuBar = QtWidgets.QAction(MainWindow)
        self.actionNuMuBar.setCheckable(True)
        self.actionNuMuBar.setChecked(True)
        self.actionNuMuBar.setObjectName("actionNuMuBar")
        self.actionNuTau = QtWidgets.QAction(MainWindow)
        self.actionNuTau.setCheckable(True)
        self.actionNuTau.setChecked(True)
        self.actionNuTau.setObjectName("actionNuTau")
        self.actionNuTauBar = QtWidgets.QAction(MainWindow)
        self.actionNuTauBar.setCheckable(True)
        self.actionNuTauBar.setChecked(True)
        self.actionNuTauBar.setObjectName("actionNuTauBar")
        self.menuFile.addAction(self.actionQuit)
        self.menuFlavors.addAction(self.actionElectron)
        self.menuFlavors.addAction(self.actionNuEBar)
        self.menuFlavors.addAction(self.actionNuMu)
        self.menuFlavors.addAction(self.actionNuMuBar)
        self.menuFlavors.addAction(self.actionNuTau)
        self.menuFlavors.addAction(self.actionNuTauBar)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuFlavors.menuAction())
 
        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

        self.electron_slider.setTickInterval(1)
        self.tau_slider.setTickInterval(1)
        self.width_slider.setTickInterval(1)
        self.width_slider.setValue(10)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "Flux Explorer"))
        self.electron_lbl.setText(_translate("MainWindow", "Theta 1-4:  0.00"))
        self.muon_lbl.setText(_translate("MainWindow", "Theta 2-4: 0.00"))
        self.tau_lbl.setText(_translate("MainWindow", "Theta 3-4: 0.00"))
        self.width_lbl.setText(_translate("MainWindow", "Width: 0.10"))
        self.recoBox.setText(_translate("MainWnidow", "Reconstructed Fluxes"))
        self.menuFile.setTitle(_translate("MainWindow", "File"))
        self.actionQuit.setText(_translate("MainWindow", "Quit"))

        self.menuFlavors.setTitle(_translate("MainWindow", "Neutrinos"))
        self.actionQuit.setText(_translate("MainWindow", "Quit"))
        self.actionElectron.setText(_translate("MainWindow", "NuE"))
        self.actionNuEBar.setText(_translate("MainWindow", "NuEBar"))
        self.actionNuMu.setText(_translate("MainWindow", "NuMu"))
        self.actionNuMuBar.setText(_translate("MainWindow", "NuMuBar"))
        self.actionNuTau.setText(_translate("MainWindow", "NuTau"))
        self.actionNuTauBar.setText(_translate("MainWindow", "NuTauBar"))


class main_window(QMainWindow):
    """
    Main function that creates the window 
    """
    def __init__(self, parent=None):
        QWidget.__init__(self, parent)
        # set up the UI
        self.ui = base_gui()
        self.ui.setupUi(self)

        # whenever the sliders value change, we update the plots 
        # may want to change this to (mouse released) if this takes a while 
        self.ui.tau_slider.valueChanged.connect(self.update_plot)
        self.ui.electron_slider.valueChanged.connect(self.update_plot)
        self.ui.width_slider.valueChanged.connect(self.update_width)


        #whenever the checkactions change, call this other function
        self.ui.actionElectron.triggered.connect(self.checked_changed)
        self.ui.actionNuEBar.triggered.connect(self.checked_changed)
        self.ui.actionNuMu.triggered.connect(self.checked_changed)
        self.ui.actionNuMuBar.triggered.connect(self.checked_changed)
        self.ui.actionNuTau.triggered.connect(self.checked_changed)
        self.ui.actionNuTauBar.triggered.connect(self.checked_changed)    
        self.ui.recoBox.clicked.connect(self.checked_changed)

        self.core_b = -0.98
        self.mantle_b = -0.83

        self.tau_angle = 0.0
        self.electron_angle = 0.0

        # copied from the make_dag python script
        # this just ensures that the names are all right 
        n_grid = 20
        self.theta03s = np.linspace(0, pi, n_grid) #el

        self.theta13s = [0.0,0.05, 0.1, 0.160875, 0.2]
        #self.theta13s = [0.160875]
        self.thetamu =self.theta13s[3]
        self.theta23s = np.linspace(0, 40*pi/180., n_grid) #tau
        self.msq = 4.47
        
        for value in self.theta13s:
            self.ui.muon_combo.addItem(str(value))
        
        self.ui.muon_combo.setCurrentIndex( self.theta13s.index(self.thetamu) )
        self.ui.muon_combo.currentIndexChanged.connect(self.update_plot)

        # load the null flux!
        self.reload_null()

        self.width = 0.10

        self.update_plot()

    def update_width(self):
        self.width = float(self.ui.width_slider.value())/100.
        self.ui.width_lbl.setText("Width: {:.2f}".format(self.width))
        self.cbar.mappable.set_clim(vmin=1-self.width, vmax=1+self.width)
        self.ui.canvas.draw()

    def checked_changed(self):
        """
        This is called whenever one of the checkboxes are checked/unchecked 
        """
        #self.ui.tau_slider.setEnabled(self.ui.recoBox.isChecked())
        #self.ui.electron_slider.setEnabled(self.ui.recoBox.isChecked())

        self.reload_null()
        self.update_plot()

    def check_key(self, key):
        """
        Passed a key, checks whether the relevant action box is checked 

        Returns True
        """
        # only return the NC key flux when looking at the raw flux
        # this prevents us from double-counting electron and tau fluxes
        if ("CC" in key) and (not self.ui.recoBox.isChecked()):
            value = True
        else:
            value = True

        if "E_nu_" in key:
            return self.ui.actionElectron.isChecked() and value
        elif "E_nuBar_" in key:
            return self.ui.actionNuEBar.isChecked() and value
        elif "Mu_nu_" in key:
            return self.ui.actionNuMu.isChecked() and value
        elif "Mu_nuBar_" in key:
            return self.ui.actionNuMuBar.isChecked() and value
        elif "Tau_nu_" in key:
            return self.ui.actionNuTau.isChecked() and value
        elif "Tau_nuBar_" in key:
            return self.ui.actionNuTauBar.isChecked() and value
        else:
            raise ValueError("Not sure what to do with key {}".format(key))

    def apply_xs(self, flux,key):
        neut = get_neut(key)
        curr = get_curr(key)
        flav = get_flavor(key)

        for i_e in range(len(self.e_reco)):
            xs = get_xs(self.e_reco[i_e], flav, neut, curr)
            flux[i_e] *= xs

        return(flux)

    def reload_null(self):
        """
        This function reloads the null flux 
        """

        sp = SterileParams(0., 0., 0., 0.)

        if self.ui.recoBox.isChecked():
            which = config["recon_flux"]+".dat"
        else:
            which = config["nu_flux_downsize"]+".dat"

        f = open(gen_filename(config["datapath"], which, sp), 'rb')
        all_data  = pickle.load(f)
        f.close()
       
        if self.ui.recoBox.isChecked():
            self.e_reco = np.array(bhist([all_data["e_reco"]]).centers)
            self.a_reco = np.array(bhist([all_data["a_reco"]]).centers)
        else:
            self.e_reco = np.array(all_data["e_true"])
            self.a_reco = np.array(all_data["a_true"])

        all_data = all_data["flux"]
                
        self.flux_null = np.zeros(shape = (len(self.e_reco), len(self.a_reco)))
        for key in all_data.keys():
            if self.check_key(key):
                self.flux_null += np.array(all_data[key]) if self.ui.recoBox.isChecked() else self.apply_xs(np.array(all_data[key]), key)
        
    def update_plot(self):
        """
        This should be called whenever we choose new angles and need to update the plots
        """
        self.update_angles()

        self.ui.figure.clear()
        ax = self.ui.figure.add_subplot(111)
        
        flux = self.get_interp_flux()
        
        limiter = False

        pmesh = ax.pcolormesh(self.a_reco, self.e_reco/(1e9), flux/self.flux_null, cmap=cm.coolwarm, vmin=None if limiter else 1.0-self.width, vmax=None if limiter else 1.+self.width)
        ax.set_yscale('log')
        ax.set_ylim([10**2, 10**6])
        plot_labels = "Reco" if self.ui.recoBox.isChecked() else "True"
        ax.set_ylabel("{} Energy [GeV]".format(plot_labels), size=14)
        ax.set_xlabel(r"{} $\cos\theta$".format(plot_labels),size=14)
        self.cbar = self.ui.figure.colorbar(pmesh, ax=ax)
        self.cbar.set_label("Sterile Flux / Null Flux")

        plt.vlines(self.core_b,ymin=10**2, ymax=10**6, colors="white", ls="-")
        plt.vlines(self.mantle_b,ymin=10**2, ymax=10**6, colors="white", ls="--")

        self.ui.canvas.draw()

    def _load_flux_file(self, filename):
        """
        simplified little thing. Just loads in the fluxes we want 
        """
        f = open(filename, 'rb')
        all_data = pickle.load(f)["flux"]
        f.close()

        if False:
            dirname, filename = os.path.split(filename)
            temp = Data(filename)
            all_data = temp.fluxes

        flux = np.zeros(shape=(len(self.e_reco), len(self.a_reco)))
        for key in all_data.keys():
            if self.check_key(key):
                flux += np.array(all_data[key]) if self.ui.recoBox.isChecked() else self.apply_xs(np.array(all_data[key]),key)

        return(flux)
 

    def get_interp_flux(self):
        """
        This creates the interpolated flux by accessing the currently set angles, finding the four neighboring fluxes, and then performing a bilinear interpolation 
        """
    
        # this gets the indices of the two mixing angle values neighboring the intermediate one we hav enow
        i_x1, i_x2 = get_loc(self.electron_angle, self.theta03s)
        i_y1, i_y2 = get_loc(self.tau_angle, self.theta23s)

        # now let's build the parameter objects using those neighboring points we have 
        param_11 = SterileParams(self.theta03s[i_x1], self.thetamu, self.theta23s[i_y1],self.msq)
        param_12 = SterileParams(self.theta03s[i_x1], self.thetamu, self.theta23s[i_y2],self.msq)
        param_21 = SterileParams(self.theta03s[i_x2], self.thetamu, self.theta23s[i_y1],self.msq)
        param_22 = SterileParams(self.theta03s[i_x2], self.thetamu, self.theta23s[i_y2],self.msq)

        which = (config["recon_flux"] if self.ui.recoBox.isChecked() else config["nu_flux_downsize"]) + ".dat"

        # using those indices, we generate the names of the flux files and load
        flux_11 = self._load_flux_file(gen_filename(config["datapath"], which, param_11))
        flux_12 = self._load_flux_file(gen_filename(config["datapath"], which, param_12))
        flux_21 = self._load_flux_file(gen_filename(config["datapath"], which, param_21))
        flux_22 = self._load_flux_file(gen_filename(config["datapath"], which, param_22))

        # these are useful intermediates used for my bilinear interpolation function 
        p0 = (self.electron_angle, self.tau_angle)
        p1 = (self.theta03s[i_x1], self.theta23s[i_y1])
        p2 = (self.theta03s[i_x2], self.theta23s[i_y2])

        return bilinear_interp(p0, p1, p2, flux_11, flux_12, flux_21, flux_22)

    def update_angles(self):
        """
        Grab the angles from the sliders, update the slider text
        """
        self.thetamu = float(self.ui.muon_combo.currentText())
        self.tau_angle = float(self.ui.tau_slider.value())*(25./100)*pi/180
        self.electron_angle = float(self.ui.electron_slider.value())*pi/100
        
        self.ui.muon_lbl.setText("Theta-24: {:.4f}, {:.0f} deg".format(self.thetamu, self.thetamu*180./pi))
        self.ui.electron_lbl.setText("Theta-14: {:.4f}, {:.0f} deg".format(self.electron_angle, self.electron_angle*180/pi))
        self.ui.tau_lbl.setText("Theta-34: {:.4f}, {:.0f} deg".format(self.tau_angle, self.tau_angle*180/pi))


app = QApplication(sys.argv)
app_instance = main_window()
if __name__=="__main__":
    app_instance.show()
    sys.exit(app.exec_())
