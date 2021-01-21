# -*- coding: utf-8 -*-
"""
Ben Smithers
"""

from PyQt5.QtWidgets import QMainWindow, QWidget, QApplication
from PyQt5 import QtCore, QtGui, QtWidgets

import sys
from math import pi

# embedding matplotlib stuff

from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt

import numpy as np

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
        self.tau_lbl = QtWidgets.QLabel(self.centralwidget)
        self.tau_lbl.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.tau_lbl.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.tau_lbl.setObjectName("tau_lbl")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.tau_lbl)
        self.tau_slider = QtWidgets.QSlider(self.centralwidget)
        self.tau_slider.setOrientation(QtCore.Qt.Horizontal)
        self.tau_slider.setObjectName("tau_slider")
        self.tau_slider.setMaximum(100)
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.tau_slider)
        self.verticalLayout.addLayout(self.formLayout)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 22))
        self.menubar.setObjectName("menubar")
        self.menuFile = QtWidgets.QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.actionQuit = QtWidgets.QAction(MainWindow)
        self.actionQuit.setObjectName("actionQuit")
        self.menuFile.addAction(self.actionQuit)
        self.menubar.addAction(self.menuFile.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

        self.electron_slider.setTickInterval(1)
        self.tau_slider.setTickInterval(1)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.electron_lbl.setText(_translate("MainWindow", "Theta e-s:  0.00"))
        self.tau_lbl.setText(_translate("MainWindow", "Theta tau-s: 0.00"))
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

        # whenever the sliders value change, we update the plots 
        # may want to change this to (mouse released) if this takes a while 
        self.ui.tau_slider.valueChanged.connect(self.update_plot)
        self.ui.electron_slider.valueChanged.connect(self.update_plot)

        self.tau_angle = 0
        self.electron_angle = 0

        self.update_angles()
        
    def update_plot(self):
        """
        This should be called whenever we choose new angles and need to update the plots
        """
        self.update_angles()

        self.ui.figure.clear()
        ax = self.ui.figure.add_subplot(111)
        
        xs = np.linspace(0,3,50)
        ys = xs*self.tau_angle + self.electron_angle

        ax.set_xlim([0,3])
        ax.set_ylim([0,4])

        ax.plot(xs,ys)
        self.ui.canvas.draw()

    def update_angles(self):
        """
        Grab the angles from the sliders, update the slider text
        """
        self.tau_angle = float(self.ui.tau_slider.value())*pi/100.
        self.electron_angle = float(self.ui.electron_slider.value())*pi/100.

        self.ui.electron_lbl.setText("Theta e-s: {:.4f}".format(self.electron_angle))
        self.ui.tau_lbl.setText("Theta tau-s: {:.4f}".format(self.tau_angle))


app = QApplication(sys.argv)
app_instance = main_window()
if __name__=="__main__":
    app_instance.show()
    sys.exit(app.exec_())
