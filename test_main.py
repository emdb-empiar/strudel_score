from recordtype import recordtype
from PyQt5 import uic
import shutil
import os
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QMessageBox, QVBoxLayout, QHBoxLayout, QGridLayout, QFileDialog, QLabel, QPushButton,\
    QMainWindow, QCheckBox, QMenuBar
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.ticker import FormatStrFormatter
import pandas as pd
import numpy as np

# from threed_strudel.utils import bioUtils as bu
# from threed_strudel.chop.chopMap import ChopMap
# from threed_strudel.parse.mapParser import MapParser
# from threed_strudel import nomenclature


CHIMERA_MODE = True
try:
    from chimerax.core.tools import ToolInstance
    from chimerax.core.commands import Command
    from chimerax.map.volume import Volume
    from chimerax.atomic.structure import AtomicStructure


except ImportError:
    CHIMERA_MODE = False

class ToolBase:
    def __init__(self, session, tool_name):
        self.tool_name = tool_name


if CHIMERA_MODE:
    Base = ToolInstance
else:
    Base = ToolBase


class StrudelScore(Base):
    SESSION_ENDURING = True
    SESSION_SAVE = False         # No local_session saving for now
    display_name = "Strudel Tool"

    def __init__(self, session, tool_name):
        Base.__init__(self, session, tool_name)
        self.local_session = session
        self.menu_bar = None
        self.ui_area = QtWidgets.QWidget()
        self.init_ui()

    def init_ui(self):
        print("Creating main")
        self.mw = QtWidgets.QWidget()
        self.basedir = os.path.dirname(os.path.abspath(__file__))
        # self.mw = QtWidgets.QWidget()
        uifile = os.path.join(self.basedir, 'strudel_score', 'src','ui', 'test_w.ui')
        print("Load ui")
        self.ui = uic.loadUi(uifile, self.mw)
        self.ui.setWindowTitle('Strudel score')
        # self.ui.show()

        self.menu_bar = self.create_menubar()
        vbox = QVBoxLayout()
        self.ui_area.setLayout(vbox)
        vbox.addWidget(self.menu_bar)
        # vbox.addWidget(self.ui)
        self.menu_bar.setNativeMenuBar(False)
        self.ui_area.show()

    def create_menubar(self):
        menu_bar = QMenuBar()
        file_menu = menu_bar.addMenu('File')
        edit_menu = menu_bar.addMenu('Edit')
        # exit_action = QAction('Exit', window)
        # exit_action.triggered.connect(exit)
        # file_menu.addAction(exit_action)
        # redo_action = QAction('Redo', window)
        # redo_action.triggered.connect(redoClicked)
        # edit_menu.addAction(redo_action)
        menu_bar.setNativeMenuBar(False)
        return menu_bar

if __name__ == '__main__':
    import sys
    app = QtWidgets.QApplication(sys.argv)
    session = ''
    tool_name = 'Strudel Score'
    window = StrudelScore(session, tool_name)
    app.exec_()