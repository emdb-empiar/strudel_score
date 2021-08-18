"""
tool.py

Strudel score visualisation tool

Copyright [2013] EMBL - European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the
"License"); you may not use this file except in
compliance with the License. You may obtain a copy of
the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on
an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
"""

__author__ = 'Andrei Istrate'
__email__ = 'andrei@ebi.ac.uk'
__date__ = '2020-05-28'


from recordtype import recordtype
from PyQt5 import uic
import shutil
import os
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QMessageBox, QVBoxLayout, QHBoxLayout, QGridLayout, QFileDialog, QLabel, QPushButton,\
    QWidget, QCheckBox, QFrame, QMenu
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.ticker import FormatStrFormatter
import pandas as pd
import numpy as np

from threed_strudel.utils import bio_utils as bu
from threed_strudel.chop.chop_map import ChopMap
from threed_strudel.parse.map_parser import MapParser
from threed_strudel import nomenclature
try:
    from .functions import find_y_label_coordinates
except ImportError:
    from functions import find_y_label_coordinates

CHIMERA_MODE = True
try:
    from chimerax.core.tools import ToolInstance
    from chimerax.core.commands import Command
    from chimerax.map.volume import Volume
    from chimerax.atomic.structure import AtomicStructure

except ImportError:
    CHIMERA_MODE = False

SD_LEVEL = 4


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
    display_name = "Strudel Score"

    def __init__(self, session, tool_name):
        Base.__init__(self, session, tool_name)
        self.local_session = session
        if CHIMERA_MODE:
            self.rc = Command(self.local_session)
            from chimerax.ui.gui import MainToolWindow
            self.tw = MainToolWindow(self)
            self.ui_tw_area = self.tw.ui_area
        else:
            self.ui_tw_area = QtWidgets.QWidget()
            # self.ui_tw_area = QMainWindow()
        self.parent = self.ui_tw_area

        self.basedir = os.path.dirname(os.path.abspath(__file__))
        self.config_path = os.path.join(self.basedir, 'config.txt')
        self.libs_path = None
        self.ui = None
        self.close_to_ligand = {}
        self.axvl = None
        self.axvtxt = None
        self.ax_residues = []
        self.tmp = None
        self.P = Parameters()
        self.k = DictKeys()
        self.line_h = 18
        self.files_lbls = []
        self.clear_indices = []

        # Layouts
        main_layout = QVBoxLayout(self.parent)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(1)
        self.parent.setLayout(main_layout)
        if CHIMERA_MODE:
            self.tw.manage(placement=None)

        self.left_plots_vbox = None
        self.mean_cc_grid = None
        # Figures
        self.right_plot_fig = None
        self.right_ax = None
        self.right_annot = None
        self.right_plot_canvas = None
        # Widgets
        self.chain_selector = None
        self.left_axis_selector = None
        self.lib_selector = None
        self.current_chain_data = None

        mbar = self.create_menubar1(self.parent)
        mbar.setFixedHeight(25)
        main_layout.addWidget(mbar)
        dw = self.load_designer_widget()
        main_layout.addWidget(dw)
        main_layout.setStretchFactor(mbar, 0)
        # main_layout.addStretch(1)
        # self.init_ui()

        self.active_lib_path = None
        self.map_path = None
        self.map_chim = None
        self.map_obj = None
        self.model_obj = None
        self.model_chim = None
        self.model_path = None
        self.segments_path = None
        self.sel_res_atomspec = None
        self.all_rotamers_id = []
        self.chop = ChopMap()
        self.tmp = '/tmp/strudel'
        if not os.path.exists(self.tmp):
            os.makedirs(self.tmp)

        self.csv_paths = {}
        self.executing = False
        self.df = None
        self.all_chains = ['No chains found']
        self.all_libs = ['No data loaded']
        self.y_axis_ranges = ['Auto', '0 to 1', '0.3 to 1', '0.5 to 1', '0.6 to 1', '0.7 to 1', '0.8 to 1']
        self.right_y_axis_ranges = ['Auto', '0.5 to 1', '0.6 to 1', '0.7 to 1', '0.8 to 1', '0.9 to 1']
        self.res_per_row = ['100', '90', '80', '70', '60', '50']
        self.metrics = ['Correlation', 'Z-score']
        self.per_chain_mean_cc = None
        self.active_chain = None
        self.active_residue = None
        self.active_residue_data = None
        self.active_res_map_obj = None
        self.active_res_model_obj = None
        self.open_motifs = []
        # type - amino acid residue type, 3 letter lower case,
        # mtp (model type) - can be same, top, other
        self.Motif = recordtype('Motif', ['type', 'map_id', 'model_id', 'show', 'mtp'])

        self.inject_data()
        # self.ui.show()
        if self.df is not None:
            # self.draw_right_plot()
            self.draw_left_plot_vline()
        # self.ui_tw_area.show()
        if CHIMERA_MODE:
            session.strudel = self

    def load_designer_widget(self):
        mw = QWidget()
        uifile = os.path.join(self.basedir, 'ui', 'test_w.ui')
        # print("Load ui")
        self.ui = ui = uic.loadUi(uifile, mw)
        # self.ui.setWindowTitle('Strudel score')
        self.chain_selector = ui.chain_comboBox
        self.resized = QtCore.pyqtSignal()

        self.ui.metric_label.setHidden(True)
        self.ui.metric_comboBox.setHidden(True)

        self.left_axis_selector = ui.y_axis_range_comboBox
        self.right_axis_selector = ui.r_y_axis_range_comboBox
        self.lib_selector = ui.lib_comboBox
        self.mean_cc_grid = ui.mean_score_gridLayout
        self.mean_cc_grid.setContentsMargins(0, 0, 0, 0)
        self.create_left_plots_container()
        self.right_plot_fig = Figure()
        self.right_ax = self.right_plot_fig.add_subplot(111)
        self.right_plot_canvas = FigureCanvas(self.right_plot_fig)
        self.right_plot_fig.subplots_adjust(left=0.25, top=0.95, bottom=0.08)
        self.create_right_plot_container()
        w = ui.width()
        ui.plots_splitter.setSizes([int(w * .8), int(w * .2)])
        self.set_labels_style()
        mw.resizeEvent = self.resized_actions

        ui.same_map_checkBox.stateChanged.connect(self.update_motifs_view)
        ui.same_model_checkBox.stateChanged.connect(self.update_motifs_view)
        ui.top_map_checkBox.stateChanged.connect(self.update_motifs_view)
        ui.top_model_checkBox.stateChanged.connect(self.update_motifs_view)
        return mw

    def resized_actions(self, event):
        if self.df is not None:
            self.draw_right_plot()

    # def init_ui(self):
    #     print("Creating main")
    #     self.mw = QtWidgets.QMainWindow()
    #     # self.mw = QtWidgets.QWidget()
    #     uifile = os.path.join(self.basedir, 'ui', 'qt_ui2_b.ui')
    #     print("Load ui")
    #     self.ui = uic.loadUi(uifile, self.mw)
    #     self.ui.setWindowTitle('Strudel score')
    #     self.setup_menubar()
    #     self.chain_selector = self.ui.chain_comboBox
    #     self.resized = QtCore.pyqtSignal()
    #
    #     # self.ui.metric_label.setHidden(True)
    #     # self.ui.metric_comboBox.setHidden(True)
    #
    #     self.left_axis_selector = self.ui.y_axis_range_comboBox
    #     self.right_axis_selector = self.ui.r_y_axis_range_comboBox
    #     self.lib_selector = self.ui.lib_comboBox
    #     self.mean_cc_grid = self.ui.mean_score_gridLayout
    #     self.mean_cc_grid.setContentsMargins(0, 0, 0, 0)
    #     self.create_left_plots_container()
    #     self.right_plot_fig = Figure()
    #     self.right_ax = self.right_plot_fig.add_subplot(111)
    #     self.right_plot_canvas = FigureCanvas(self.right_plot_fig)
    #     self.right_plot_fig.subplots_adjust(left=0.25, top=0.95, bottom=0.08)
    #     self.create_right_plot_container()
    #     w = self.ui.width()
    #     self.ui.plots_splitter.setSizes([int(w * .8), int(w * .2)])
    #     self.set_labels_style()
    #     self.mw.resizeEvent = self.resized_actions
    #
    #     self.ui.same_map_checkBox.stateChanged.connect(self.update_motifs_view)
    #     self.ui.same_model_checkBox.stateChanged.connect(self.update_motifs_view)
    #     self.ui.top_map_checkBox.stateChanged.connect(self.update_motifs_view)
    #     self.ui.top_model_checkBox.stateChanged.connect(self.update_motifs_view)

    def integrate(self):
        parent = self.ui_tw_area
        layout = QtWidgets.QHBoxLayout()
        # print("Insert main in main")
        layout.addWidget(self.mw)
        layout.setStretchFactor(self.mw, 0)
        layout.setContentsMargins(0, 0, 0, 0)
        parent.setLayout(layout)
        self.tw.manage(placement=None)
        # print("Integrated")

    def inject_data(self):
        # print('Inject data called')
        self.fill_selectors()
        self.update_chain_cc_lbls()
        self.fill_files_data()
        if self.df is not None:
            # self.find_ligand_close_contacts()
            self.update_left_plots()
            self.draw_right_plot()
            self.draw_left_plot_vline()

    def fill_files_data(self):
        for lbl in self.files_lbls:
            lbl.setParent(None)
        self.files_lbls = []
        if self.map_path is not None:
            map_file = os.path.basename(self.map_path)
            lbl1 = QLabel(map_file)
            self.files_lbls.append(lbl1)
            self.ui.files_gridLayout.addWidget(lbl1, 0, 1)
            self.ui.files_gridLayout.setColumnStretch(1, 1)

        if self.model_path is not None:
            model_file = os.path.basename(self.model_path)
            lbl2 = QLabel(model_file)
            self.files_lbls.append(lbl2)
            self.ui.files_gridLayout.addWidget(lbl2, 1, 1)

    def setup_menubar(self):
        menu_bar = self.ui.menuBar()
        self.ui.actionOpen_Project.setShortcut("Ctrl+O")
        self.ui.actionOpen_Project.triggered.connect(self.open_project)
        self.ui.actionSet_motif_library.triggered.connect(self.set_library)
        self.ui.actionChop_residue.triggered.connect(self.chop_selected_residue)
        self.ui.actionChop_residue_environement_sensitive.triggered.connect(self.chop_selected_residue_env_sensitive)
        self.ui.actionLoad_rotamers.triggered.connect(self.load_all_rotamers)
        self.ui.actionClear_rotamers.triggered.connect(self.clear_rotamers)
        menu_bar.setNativeMenuBar(False)

    def create_menubar1(self, parent):
        mbar = QFrame(parent)
        mbar.setStyleSheet('QFrame { background-color: None; }')
        layout = QHBoxLayout(mbar)
        # layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(10)

        file_menu_entries = (
            ('Open Project', self.open_project),
            ('separator', None),
            ('Set Motif Library', self.set_library)
        )

        fb = QPushButton('File', mbar)
        button_style = 'QPushButton { border: none; } QPushButton::menu-indicator { image: none; }'
        fb.setStyleSheet(button_style)
        layout.addWidget(fb)
        fmenu = QMenu(fb)
        fb.setMenu(fmenu)
        for text, func in file_menu_entries:
            if text == 'separator':
                fmenu.addSeparator()
            else:
                fmenu.addAction(text, func)
        # layout.addStretch(1)

        action_menu_entries = (
            ('Chop residue', self.chop_selected_residue),
            ('Chop residue environement sensitive', self.chop_selected_residue_env_sensitive),
            ('separator', None),
            ('Load rotamers', self.load_all_rotamers),
            ('Clear rotamers', self.clear_rotamers)
        )

        ab = QPushButton('Actions', mbar)
        ab.setStyleSheet(button_style)
        layout.addWidget(ab)
        amenu = QMenu(ab)
        ab.setMenu(amenu)
        for text, func in action_menu_entries:
            if text == 'separator':
                amenu.addSeparator()
            else:
                amenu.addAction(text, func)
        layout.addStretch(1)
        return mbar

    def set_library(self, lib_path=None):
        if lib_path is None:
            lib_path = str(QFileDialog.getExistingDirectory())
        else:
            lib_path = os.path.abspath(os.path.expanduser(os.path.expandvars(lib_path)))

        files = os.listdir(lib_path)
        files = [f for f in files if f.startswith('motifs_')]
        if len(files) == 0:
            self.show_mesage("Error", f"Could not find strudel libraries in the directory {lib_path}\n"
                             f"Please select a directory which contain folders named motifs_***")
            return

        record = f'lib_path={lib_path}\n'
        lines = [record]
        try:
            with open(self.config_path, 'r') as f:
                lines = f.readlines()
                for i, line in enumerate(lines):
                    if 'lib_path' in line:
                        lines[i] = record
        except FileNotFoundError:
            pass
        with open(self.config_path, 'w') as f:
            for line in lines:
                f.write(line)
        self.libs_path = lib_path
        self.show_mesage('Info', f'Strudel library set to: {lib_path}')

    def read_config(self):
        path = None
        try:
            with open(self.config_path, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    if 'lib_path' in line:
                        path = line.split('=')[-1].strip()
        except FileNotFoundError:
            pass
        if path is None:
            self.show_nolib_text()
        else:
            self.libs_path = path

    def show_nolib_text(self):
        text = 'Please set motifs library path!\n\n' \
                  'File/Set Motif Library'
        title = "No motif library set"
        self.show_mesage(title, text)

    def show_mesage(self, title, text):
        msg = QMessageBox()
        msg.setWindowTitle(title)
        msg.setText(text)
        msg.setIcon(QMessageBox.Warning)
        x = msg.exec_()

    def download_library(self):
        pass

    def open_project(self, set_path=None):
        proj_path = None
        if set_path is None:
            set_path = str(QFileDialog.getExistingDirectory())
        set_path = os.path.abspath(os.path.expanduser(os.path.expandvars(set_path)))
        try:
            folders = os.listdir(set_path)
        except FileNotFoundError:
            return
        for file in folders:
            if file.startswith('vs_'):
                proj_path = set_path
        if proj_path is None:
            for file in folders:
                if file.startswith('out'):
                    proj_path = os.path.join(set_path, file)
                    break

        if proj_path:
            self.csv_paths = {}
            self.set_paths(proj_path)
            self.set_data()
            self.load_map_model()
            self.inject_data()
            self.clear_tmp()
            # self.ui.show()
            self.read_config()
            self.set_active_lib()

        else:
            self.show_mesage('Fail', 'The specified directory does not contain Strudel validation results files!')

    def set_active_lib(self):
        avail_libs = os.listdir(self.libs_path)
        avail_libs = [l for l in avail_libs if not l.startswith('.')]
        lib_path = None
        for l in avail_libs:
            if l.endswith(self.lib_selector.currentText()):
                lib_path = os.path.join(self.libs_path, l, 'motifs')
        self.active_lib_path = lib_path

    def load_map_model(self):
        if CHIMERA_MODE:
            self.run_x('close all')
            # self.map_chim = self.run_x(f'open {self.map_path}')
            # self.run_x(f'volume #{self.map_chim.id_string} sdLevel {SD_LEVEL + 2} style mesh step 1')
            # self.model_chim = self.run_x(f'open {self.model_path}')
            self.model_chim = self.open_model(self.model_path)
            self.map_chim = self.open_model(self.map_path)
            self.run_x('clipper assoc #2 to #1')
            self.run_x('size stickRadius 0.03')
            from chimerax.std_commands import cofr, camera
            cofr.cofr(self.local_session, 'centerOfView', show_pivot=False)
            # self.run_x('style stick')
            # self.run_x('size stickRadius 0.03')
            # self.run_x('~ribbon')
            # self.run_x('show atoms')
        self.model_obj = bu.load_structure(self.model_path)
        self.map_obj = MapParser(self.map_path)

    def find_ligand_close_contacts(self, d=2):
        ligands = []
        for chain in self.model_obj[0]:
            tmp = []
            for res in chain:
                r_name = res.get_resname().upper()
                if r_name not in nomenclature.AA_RESIDUES_LIST:
                    ligands.append(res)

        for chain in self.model_obj[0]:
            tmp = []
            for res in chain:
                r_name = res.get_resname().upper()
                if r_name in nomenclature.AA_RESIDUES_LIST:
                    for atom in res:
                        for res2 in ligands:
                            for atom2 in res2:
                                dist = atom - atom2
                                if dist < d:
                                    tmp.append(res.get_id()[1])
            tmp = list(set(tmp))
            self.close_to_ligand[chain.get_id()] = tmp

    def clear_tmp(self):
        try:
            shutil.rmtree(self.tmp)
        except FileNotFoundError:
            pass
        os.makedirs(self.tmp)

    def run_x(self, command):
        if CHIMERA_MODE:
            # print(f'Running: {command}')
            try:
                # returns a list of lists
                return self.rc.run(command)
            except TypeError:
                return None
        else:
            self.no_chimera_message()

    def open_model(self, path):
        try:
            model_obj = self.local_session.open_command.open_data(path)[0][0]
            if type(model_obj) == AtomicStructure or type(model_obj) == Volume:
                self.local_session.models.add([model_obj])
                return model_obj
        except AttributeError:
            model_obj = self.run_x(f'open {path}')[0][0]
            if type(model_obj) == AtomicStructure or type(model_obj) == Volume:
                return model_obj
            else:
                self.show_mesage('Error', 'Could not get model object reference after opening')

    def no_chimera_message(self):
        text = 'Functionality not available in standalone mode.' \
               '\nPlease run as ChimeraX Tool'
        self.show_mesage('Functionality unavailable', text)

    def fill_left_plot_axix_range_selector(self):

        combo_box = self.ui.y_axis_range_comboBox
        combo_box.addItems(self.y_axis_ranges)
        if self.df is not None:
            combo_box.currentIndexChanged.connect(self.switch_axis_range)

    def switch_axis_range(self):
        if self.df is not None:
            self.update_left_plots()

    def fill_right_plot_axix_range_selector(self):
        combo_box = self.ui.r_y_axis_range_comboBox
        combo_box.addItems(self.right_y_axis_ranges)
        if self.df is not None:
            combo_box.currentIndexChanged.connect(self.switch_right_axis_range)

    def switch_right_axis_range(self):
        if self.df is not None:
            self.draw_right_plot()

    def update_chain_cc_lbls(self):
        h = self.line_h
        if self.df is not None:
            for i in reversed(range(self.mean_cc_grid.count())):
                try:
                    self.mean_cc_grid.itemAt(i).widget().setParent(None)
                except AttributeError:
                    self.mean_cc_grid.removeItem(self.mean_cc_grid.itemAt(i))

            self.mean_cc_grid.setSpacing(0)
            i = 0
            j = 0

            for chain, value in self.per_chain_mean_cc.items():
                text = f'<font style="{self.P.grey_text} font-size: {self.P.cc_text_size}px">Chain {chain}: </font>' \
                       f'<font style="font-size: {self.P.cc_text_size}px">{value[0]:.3f} ({value[1]:.1f}% red)</font>'
                item = QLabel(text)
                self.mean_cc_grid.addWidget(item, j, i)
                i += 1
                if i == self.P.per_chain_cc_columns:
                    j += 1
                    i = 0

            if j > 2:
                h = self.line_h * 3
            elif j > 0:
                h = (j + 1) * self.line_h
        self.ui.mean_score_scrollArea.setFixedHeight(h)

    def fill_selectors(self):
        # Lib selector
        try:
            self.ui.lib_comboBox.currentIndexChanged.disconnect(self.switch_lib)
        except TypeError:
            pass
        self.ui.lib_comboBox.clear()
        self.ui.lib_comboBox.addItems(self.all_libs)
        self.ui.lib_comboBox.currentIndexChanged.connect(self.switch_lib)
        # Chain selector
        try:
            self.ui.chain_comboBox.currentIndexChanged.disconnect(self.switch_chain)
        except TypeError:
            pass
        self.ui.chain_comboBox.clear()
        self.ui.chain_comboBox.addItems(self.all_chains)

        self.ui.chain_comboBox.currentIndexChanged.connect(self.switch_chain)
        # Y range selector
        combo_box = self.ui.y_axis_range_comboBox
        self.set_combobox_style1(combo_box)
        if combo_box.currentText() not in self.y_axis_ranges:
            combo_box.addItems(self.y_axis_ranges)
        if self.df is not None:
            combo_box.currentIndexChanged.connect(self.switch_axis_range)

        combo_box = self.ui.r_y_axis_range_comboBox
        self.set_combobox_style1(combo_box)
        if combo_box.currentText() not in self.y_axis_ranges:
            combo_box.addItems(self.right_y_axis_ranges)
        if self.df is not None:
            combo_box.currentIndexChanged.connect(self.switch_right_axis_range)

        combo_box = self.ui.res_per_row_comboBox
        self.set_combobox_style1(combo_box)
        if combo_box.currentText() not in self.res_per_row:
            combo_box.addItems(self.res_per_row)
        combo_box.setCurrentText(self.P.default_res_per_row)
        if self.df is not None:
            combo_box.currentIndexChanged.connect(self.update_left_plots)

        combo_box = self.ui.metric_comboBox
        if combo_box.currentText() not in self.metrics:
            combo_box.addItems(self.metrics)
        combo_box.setCurrentText(self.P.default_metrics)
        if self.df is not None:
            combo_box.currentIndexChanged.connect(self.update_left_plots)

        # self.color_selectors([self.ui.res_per_row_comboBox])

    def set_combobox_style1(self, combo_box):
        """
        Custom QComboBox for plots area to make sure the text is visible
        when dark color scheme is used
        :param combo_box:
        :return:
        """
        # combo_box.setObjectName("myParentWidget")
        # combo_box.setStyleSheet('QWidget#myParentWidget { background-color: red}')
        combo_box.setStyleSheet('QComboBox { color: black}')

    def do_nothing(self):
        pass

    def switch_lib(self):
        self.set_data()
        self.set_active_lib()
        self.update_chain_cc_lbls()
        self.switch_chain()
        self.clear_residue_data()
        self.active_residue = self.current_chain_data[self.k.RES_NR].iloc[0]
        self.draw_left_plot_vline()
        self.ui.show()

    def switch_chain(self):
        self.active_chain = self.chain_selector.currentText()
        self.current_chain_data = self.get_chain_data(self.df, self.active_chain)
        self.update_left_plots()
        self.clear_residue_data()
        self.active_residue = self.current_chain_data[self.k.RES_NR].iloc[0]
        self.draw_right_plot()
        self.draw_left_plot_vline()

    def clear_residue_data(self):
        self.open_motifs = []
        self.active_res_map_obj = None
        self.active_res_model_obj = None
        if CHIMERA_MODE and len(self.clear_indices) > 0:
            # for index in self.clear_indices:
            #     self.run_x(f'close #{index}')
            self.run_x(f'close #{",".join(self.clear_indices)}')
            self.clear_indices = []

    def create_left_plots_container(self):
        ui = self.ui
        grid = QGridLayout()
        y_label = QLabel()
        y_label.setPicture(self.create_left_axis_y_label())
        # y_label.show()
        x_label = QLabel()
        txt_lst = self.P.x_axis_text
        x_label_text = f'<font style="{self.P.grey_text} {self.P.left_plot_axis_font}">{txt_lst[0]} </font>' \
                       f'<font style="color: black; {self.P.left_plot_axis_font}">{txt_lst[1]}</font>'
        x_label.setText(x_label_text)
        x_label.setAlignment(QtCore.Qt.AlignHCenter)
        grid.setSpacing(5)
        grid.addWidget(y_label, 0, 0)
        scroll_area = ui.lef_plot_scrollArea
        grid.addWidget(scroll_area, 0, 1)
        grid.addWidget(x_label, 1, 1)
        grid.setContentsMargins(12, 6, 0, 12)
        lf = ui.left_plot_frame
        ui.left_plot_frame_verticalLayout.addLayout(grid)
        self.left_plots_vbox = self.ui.left_plot_scrollarea_verticalLayout
        p = lf.palette()
        lf.setAutoFillBackground(True)
        p.setColor(lf.backgroundRole(), QtCore.Qt.white)


        lf.setPalette(p)

        # lf.setObjectName("myParentWidget")
        # lf.setStyleSheet('QWidget#myParentWidget { background-color: white}')

    def create_right_plot_container(self):
        p = self.P.r_tick_params
        self.right_ax.tick_params(labelsize=p['labelsize'], length=p['length'], pad=p['pad'], width=p['width'])
        [self.right_ax.spines[i].set_linewidth(p['width']) for i in self.right_ax.spines.keys()]

        left_arrow = self.ui.left_pushButton
        right_arrow = self.ui.right_pushButton
        left_red = self.ui.left_red_pushButton
        right_red = self.ui.right_red_pushButton

        self.ui.right_plot_show_pushButton.clicked.connect(self.show_selected_residue)
        self.ui.right_plot_show_pushButton.setHidden(True)
        left_arrow.clicked.connect(self.navigate_residue_left)
        right_arrow.clicked.connect(self.navigate_residue_right)
        left_red.clicked.connect(self.navigate_red_residue_left)
        right_red.clicked.connect(self.navigate_red_residue_right)
        self.right_arrow = right_arrow
        vlayout = self.ui.right_plot_frame_verticalLayout
        vlayout.insertWidget(1, self.right_plot_canvas)
        vlayout.setContentsMargins(0, 0, 0, 0)

        frame = self.ui.right_plot_frame

        p = frame.palette()
        frame.setAutoFillBackground(True)
        p.setColor(frame.backgroundRole(), QtCore.Qt.white)
        frame.setPalette(p)
        #
        # frame.setObjectName("myParentWidget")
        # frame.setStyleSheet('QWidget#myParentWidget { background-color: None}')


        # frame.show()

    def navigate_residue(self, direction):
        data = self.current_chain_data
        if data is None:
            return

        previous_active = self.active_residue

        red_res = data[data[self.k.TOP_SAME_TYPE_DIFF] > self.P.red_difference][self.k.RES_NR].tolist()
        avail_res = data[self.k.RES_NR].tolist()
        index = avail_res.index(self.active_residue)

        if direction == 'left':
            if index > 0:
                self.active_residue = avail_res[index-1]
        elif direction == 'right':
            if index < len(avail_res)-1:
                self.active_residue = avail_res[index + 1]
        elif direction == 'red_left':
            red_res = reversed(red_res)
            if index > 0:
                for res in red_res:
                    if res < self.active_residue:
                        self.active_residue = res
                        break
        elif direction == 'red_right':
            if index < len(avail_res):
                for res in red_res:
                    if res > self.active_residue:
                        self.active_residue = res
                        break
        else:
            return
        if previous_active != self.active_residue:
            self.clear_residue_data()
            self.draw_right_plot()
            if CHIMERA_MODE:
                self.view_residue()
                self.draw_right_plot()
            self.draw_left_plot_vline()

    def show_selected_residue(self):
        self.clear_residue_data()
        if CHIMERA_MODE:
            self.view_residue()

    def navigate_residue_left(self):
        self.navigate_residue('left')

    def navigate_residue_right(self):
        self.navigate_residue('right')

    def navigate_red_residue_left(self):
        self.navigate_residue('red_left')

    def navigate_red_residue_right(self):
        self.navigate_residue('red_right')

    def update_left_plots(self):
        y_range = self.left_axis_selector.currentText()
        self.current_chain_data = self.get_chain_data(self.df, self.active_chain)
        res = int(self.ui.res_per_row_comboBox.currentText())
        metric = self.ui.metric_comboBox.currentText()
        all_left = self.draw_left_plots(self.current_chain_data,
                                        max_res=res,
                                        red_diff=self.P.red_difference,
                                        blue_threshold=self.P.blue_threshold,
                                        y_range=y_range,
                                        metric=metric)
        # Clear left plots
        for i in reversed(range(self.left_plots_vbox.count())):
            self.left_plots_vbox.itemAt(i).widget().setParent(None)

        for plot in all_left:
            plot.setMinimumHeight(self.P.left_plot_min_h)
            self.left_plots_vbox.addWidget(plot)

    def draw_left_plots(self, data, max_res=60, red_diff=0.05, blue_threshold=0.7, y_range='Auto', metric='Correlation'):
        all_plots = []
        nr_list = data[self.k.RES_NR].tolist()
        corr_cc_list = data[self.k.SAME_TYPE_CC].tolist()
        corr_cc_list = np.array(corr_cc_list)
        corr_cc_list[np.isnan(corr_cc_list)] = 1
        txt_size = self.P.l_lbl_size

        x_min = nr_list[0]
        x_max = nr_list[-1]
        y_max = 1.05
        if y_range == '0 to 1':
            y_min = 0
        elif y_range == '0.3 to 1':
            y_min = 0.3
        elif y_range == '0.5 to 1':
            y_min = 0.5
        elif y_range == '0.6 to 1':
            y_min = 0.6
        elif y_range == '0.7 to 1':
            y_min = 0.7
        elif y_range == '0.8 to 1':
            y_min = 0.8
        else:
            y_min = min(corr_cc_list) - 0.2

        if y_min is None:
            y_min = 0

        nr_elem = x_max - x_min
        nr = nr_elem // max_res
        if nr_elem % max_res != 0:
            nr += 1
        self.ax_residues = []
        for i in range(nr):
            ax_lim = (round(y_min, 2), round(y_max, 2))
            x_start = x_min + i * max_res
            x_end = x_start + max_res
            fig = Figure()
            ax = fig.add_subplot(111)
            ax.set_xlim(x_start - 0.5, x_end)
            if metric == 'Correlation':
                ax.set_ylim(ax_lim)
            else:
                ax.set_ylim([-3, 3])
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')
            ax.yaxis.set_label_position("right")
            p = self.P.l_tick_params
            ax.tick_params(labelsize=p['labelsize'], length=p['length'], pad=p['pad'], width=p['width'])
            [ax.spines[i].set_linewidth(p['width']) for i in ax.spines.keys()]
            ax.set_xticks([i for i in range(x_start, x_end + 1, p['xgap'])])

            green_row_data = data[(x_start <= data[self.k.RES_NR]) & (data[self.k.RES_NR] < x_end)]
            green_x = green_row_data[self.k.RES_NR].tolist()
            self.ax_residues.append([ax, fig, green_x])

            if metric == 'Correlation':
                green_y = green_row_data[self.k.SAME_TYPE_CC].tolist()
                red_row_data = green_row_data[green_row_data[self.k.TOP_SAME_TYPE_DIFF] > red_diff]
                red_x = red_row_data[self.k.RES_NR].tolist()
                red_y = red_row_data[self.k.TOP_CC].tolist()

                blue_row_data = green_row_data[green_row_data[self.k.SAME_TYPE_CC] < blue_threshold]
                blue_x = blue_row_data[self.k.RES_NR].tolist()
                blue_y = blue_row_data[self.k.SAME_TYPE_CC].tolist()
                top_key = self.k.TOP_CC
                same_key = self.k.SAME_TYPE_CC

            elif metric == 'Z-score':
                green_y = green_row_data[self.k.SAME_TYPE_Z].tolist()
                red_row_data = green_row_data[green_row_data[self.k.TOP_SAME_TYPE_DIFF_Z] > red_diff]
                red_x = red_row_data[self.k.RES_NR].tolist()
                red_y = red_row_data[self.k.TOP_Z].tolist()

                blue_row_data = green_row_data[green_row_data[self.k.SAME_TYPE_CC] < blue_threshold]
                blue_x = blue_row_data[self.k.RES_NR].tolist()
                blue_y = blue_row_data[self.k.SAME_TYPE_Z].tolist()
                top_key = self.k.TOP_Z
                same_key = self.k.SAME_TYPE_Z
            else:
                top_key = self.k.TOP_CC
                same_key = self.k.SAME_TYPE_CC
                green_y = []
                red_y = []
                red_x = []
                blue_x = []
                blue_y = []
                red_row_data = []
                blue_row_data = []

            # z_green_y = []

            # if metric == 'Z-score':
            #     for index, row in green_row_data.iterrows():
            #         all_corr = []
            #         for name, value in row.iteritems():
            #             if name.endswith(self.k.TYPE_TOP_CC):
            #                 all_corr.append(value)
            #         all_corr = np.array(all_corr)
            #         aver = all_corr.sum() / len(all_corr)
            #         z = (row[self.k.SAME_TYPE_CC] - aver) / all_corr.std()
            #         z_green_y.append(z)

            # ax.plot(green_x, green_y, 'o', c=self.P.l_same_type_color,
            #         markersize=self.P.l_marker_size,
            #         picker=self.P.l_picker_size)

            ax.plot(green_x, green_y, 'o', c=self.P.same_type_c,
                    markersize=self.P.l_marker_size,
                    picker=True, pickradius=self.P.l_picker_size)
            # ax.plot(red_x, red_y, 'o', c=self.P.l_top_type_color,
            #         markersize=self.P.l_marker_size,
            #         picker=self.P.l_picker_size)
            ax.plot(red_x, red_y, 'o', c=self.P.outlier_c,
                    markersize=self.P.l_marker_size,
                    picker=True, pickradius=self.P.l_picker_size)
            # ax.plot(blue_x, blue_y, 'o', c=self.P.l_same_type_lowcorr_color,
            #         markersize=self.P.l_marker_size,
            #         picker=self.P.l_picker_size)
            ax.plot(blue_x, blue_y, 'o', c=self.P.very_low_c,
                    markersize=self.P.l_marker_size,
                    picker=True, pickradius=self.P.l_picker_size)
            try:
                if self.close_to_ligand[self.active_chain]:
                    rl_x = [r for r in self.close_to_ligand[self.active_chain] if r in green_x]
                    ax.plot(rl_x, [ax_lim[0]+0.05 for i in range(len(rl_x))], '*', c='b')
            except KeyError:
                pass
            # if metric == 'Correlation':

            for index, row in red_row_data.iterrows():
                ax.text(row[self.k.RES_NR], row[top_key], row[self.k.M_TOP_TYPE].upper(),
                        color=self.P.outlier_c,
                        size=txt_size, horizontalalignment='left', verticalalignment='bottom', rotation=50)

                ax.text(row[self.k.RES_NR], row[same_key], row[self.k.RES_TYPE].upper(),
                        color=self.P.l_lbl_same_type_color,
                        size=txt_size, horizontalalignment='right', verticalalignment='top', rotation=50)
            for index, row in blue_row_data.iterrows():
                ax.text(row[self.k.RES_NR], row[same_key], row[self.k.RES_TYPE].upper(),
                        color=self.P.l_lbl_same_type_lowcorr_color,
                        size=txt_size, horizontalalignment='right', verticalalignment='top', rotation=50)

            fig.subplots_adjust(left=0.04, right=0.98, bottom=0.2, top=0.9)

            canvas = FigureCanvas(fig)
            # fig.savefig('/Users/andrei/tmp/valid/6272_100/out/vs_motifs_2.0-2.5/fig/fig' + str(i) + '.png', dpi=200)
            fig.canvas.mpl_connect('pick_event', self.update_right_plot)
            all_plots.append(canvas)
        return all_plots

    def draw_left_plot_vline(self):
        if not self.ax_residues:
            return
        try:
            self.axvl.remove()
            self.axvtxt.remove()
        except AttributeError:
            pass
        for pair in self.ax_residues:
            if self.active_residue in pair[2]:
                ax = pair[0]
                self.axvl = ax.axvline(self.active_residue, color=self.P.vline_color,
                                          linewidth=self.P.vline_width, ls='--')
                self.axvtxt = ax.text(self.active_residue, 1.06, str(self.active_residue),
                                             color=self.P.vline_color,
                                             size=self.P.l_lbl_size, horizontalalignment='center')
            pair[1].canvas.draw_idle()

    def update_right_plot(self, event):
        thisline = event.artist
        xdata = thisline.get_xdata()
        ind = event.ind
        self.clear_residue_data()
        self.active_residue = xdata[ind][0]
        self.draw_right_plot()
        self.draw_left_plot_vline()
        self.view_residue()
        self.draw_right_plot()

    def draw_right_plot(self, max_res_display=20, y_max=1):
        metric = self.ui.metric_comboBox.currentText()
        width = self.right_plot_canvas.width()
        height = self.right_plot_canvas.height()
        if height == 0:
            return
        self.right_ax.clear()

        r_nr = self.active_residue
        r = self.current_chain_data[self.current_chain_data[self.k.RES_NR] == r_nr]
        self.active_residue_data = r
        same_type = r[self.k.RES_TYPE].tolist()[0].upper()

        if metric == 'Correlation':
            key = self.k.TYPE_TOP_CC
        else:
            key = self.k.TYPE_TOP_Z
        y_data = []
        for name, value in r.iteritems():
            if name.endswith(key):
                y_data.append((value.tolist()[0], name.split('_')[0].upper()))


        y_data.sort(key=lambda v: v[0], reverse=True)
        y_data = y_data[:max_res_display]
        y_lst = [v[0] for v in y_data]
        y_txt = [v[1] for v in y_data]

        y_range = self.right_axis_selector.currentText()
        if y_range == '0.5 to 1':
            y_min = 0.5
        elif y_range == '0.6 to 1':
            y_min = 0.6
        elif y_range == '0.7 to 1':
            y_min = 0.7
        elif y_range == '0.8 to 1':
            y_min = 0.8
        elif y_range == '0.9 to 1':
            y_min = 0.9
        else:
            y_min = min(y_lst) - 0.05

        text_size = self.P.r_lbl_size
        # print(metric)
        if metric == 'Z-score':
            y_min = -3
            y_max = 3
        #     tmp = np.array(y_lst)
        #     mean = tmp.sum() / len(tmp)
        #     std = tmp.std()
        #     y_lst = [(v - mean) / std for v in tmp]
        #     y_min = -2
        #     y_max = 2
        # y_lbls = [i for i in y_lst if y_min < i < y_range[1]]
        y_lbls, text_size = find_y_label_coordinates(y_lst, [y_min, y_max], height, text_size, spacing=1)
        text_w = text_size / width * 4

        # Consider using dash option for labels
        for i in range(len(y_lbls)):
            if same_type != y_txt[i]:
                self.right_ax.plot([r_nr - 1.3 * text_w, r_nr], [y_lbls[i], y_lst[i]], '--', linewidth=0.3,
                                   c=self.P.outlier_c)
                self.right_ax.plot(r_nr, y_lst[i], 'o', markersize=2.5, c=self.P.outlier_c)
                self.right_ax.text(r_nr - 2.5 * text_w, y_lbls[i], y_txt[i], color=self.P.outlier_c, size=text_size,
                                   picker=True, verticalalignment='center', horizontalalignment='left')
            else:
                self.right_ax.plot([r_nr - 1.3 * text_w, r_nr], [y_lbls[i], y_lst[i]], '--', linewidth=0.3,
                                   c=self.P.same_type_c)
                self.right_ax.plot(r_nr, y_lst[i], 'o', markersize=2.5, c=self.P.same_type_c)
                self.right_ax.text(r_nr - 2.5 * text_w, y_lbls[i], y_txt[i], color=self.P.same_type_c, size=text_size,
                                   picker=True, verticalalignment='center', horizontalalignment='left')
            for motif in self.open_motifs:
                if y_txt[i].lower() == motif.type and motif.show:
                    self.right_ax.text(r_nr - 3 * text_w, y_lbls[i], '*', size=text_size, color=self.P.vline_color,
                                       verticalalignment='center')
        self.right_annot = self.right_ax.annotate("", xy=(0, 0), xytext=(20, 20), textcoords="offset points",
                            bbox=dict(boxstyle="round", fc="w"),
                            arrowprops=dict(arrowstyle="->"))
        self.right_annot.set_visible(False)

        self.right_ax.set_ylim(y_min, y_max)
        self.right_ax.set_xlim(r_nr - 0.7, r_nr + 0.3)
        self.right_ax.set_xticks([r_nr])
        self.right_ax.set_xlabel(self.P.x_axis_text[1], size=self.P.r_ax_title_lbl_size)
        self.right_ax.set_ylabel(self.P.y_axis_text[1], size=self.P.r_ax_title_lbl_size)
        self.right_ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        self.right_plot_fig.subplots_adjust(left=0.25, top=0.98, bottom=0.08)
        self.right_plot_fig.canvas.mpl_connect('pick_event', self.switch_motif)
        # self.right_plot_fig.canvas.mpl_connect('motion_notify_event', self.right_plot_label_hover)
        # self.right_plot_fig.canvas.mpl_connect('axes_enter_event', self.enter_axes)
        # self.right_plot_fig.canvas.mpl_connect('axes_leave_event', self.leave_axes)
        self.right_plot_fig.canvas.draw_idle()

    # def right_plot_label_hover(self, event):
    #     # thisline = event.artist
    #     # xdata = thisline.get_xdata()
    #     thisline = event.artist
    #     res_type = thisline.get_text().lower()
    #     print(f'Hover over {res_type}')
    #
    #     def update_annot(ind):
    #
    #         # pos = sc.get_offsets()[ind["ind"][0]]
    #         pos = [1,1]
    #         self.right_annot.xy = pos
    #         text = 'some text'
    #         self.right_annot.set_text(text)
    #         # self.right_annot.get_bbox_patch().set_facecolor(cmap(norm(c[ind["ind"][0]])))
    #         # self.right_annot.get_bbox_patch().set_alpha(0.4)
    #
    #     vis = self.right_annot.get_visible()
    #     if event.inaxes == self.right_ax:
    #         cont, ind = self.right_ax.contains(event)
    #         if cont:
    #             update_annot(ind)
    #             self.right_annot.set_visible(True)
    #             self.right_plot_fig.canvas.draw_idle()
    #         else:
    #             if vis:
    #                 self.right_annot.set_visible(False)
    #                 self.right_plot_fig.canvas.draw_idle()
    #
    # def enter_axes(self, event):
    #     print('enter_axes', event.inaxes)
    #     event.inaxes.patch.set_facecolor('yellow')
    #     event.canvas.draw()
    #
    # def leave_axes(self, event):
    #     print('leave_axes', event.inaxes)
    #     event.inaxes.patch.set_facecolor('white')
    #     event.canvas.draw()


    def switch_motif(self, event):
        thisline = event.artist
        res_type = thisline.get_text().lower()
        if not self.open_motifs and CHIMERA_MODE:
            self.view_residue()

        if res_type not in [m.type for m in self.open_motifs]:
            r = self.active_residue_data
            motif_name = r[res_type + self.k.TYPE_TOP_NAME].tolist()[0]

            try:
                matrix = r[res_type + self.k.TYPE_TOP_MATRIX].tolist()[0]
            except:
                matrix = None
            map_path = os.path.join(self.active_lib_path, motif_name + '.mrc')
            model_path = os.path.join(self.active_lib_path, motif_name + '.cif')

            if os.path.exists(map_path) and os.path.exists(model_path):
                map_obj, model_obj = self.display_motif(map_path, model_path, matrix)
                self.open_motifs.append(self.Motif(type=res_type, map_id=map_obj.id_string, model_id=model_obj.id_string,
                                                   show=False, mtp='other'))

            else:
                text = f'Missing files in {self.active_lib_path}\n' \
                       f'{motif_name + ".mrc"},  {motif_name + ".cif"}'
                self.show_mesage('Motif library missing', text)

        other = True
        for motif in self.open_motifs:
            if motif.type == res_type and motif.mtp in ['same', 'top']:
                motif.show = not motif.show
                other = False
                if motif.mtp == 'same':
                    self.ui.same_map_checkBox.setChecked(motif.show)
                    # self.ui.same_model_checkBox.setChecked(motif.show)
                elif motif.mtp == 'top':
                    self.ui.top_map_checkBox.setChecked(motif.show)
        del_ind = []
        if other:
            for i, motif in enumerate(self.open_motifs):
                if motif.mtp == 'other':
                    if motif.type == res_type:
                        motif.show = not motif.show
                    else:
                        self.run_x(f'close #{motif.map_id},{motif.model_id}')
                        del_ind.append(i)
            for i in del_ind:
                del self.open_motifs[i]
        self.update_motifs_view()

    def open_chopped_residue(self, map_path, model_path, color='#B2B2B2'):
        res_map = self.open_model(map_path)
        self.clear_indices.append(res_map.id_string)
        self.run_x(f'volume #{res_map.id_string} sdLevel {SD_LEVEL}')
        res_model = self.open_model(model_path)
        self.clear_indices.append(res_model.id_string)
        self.run_x(f'color #{res_map.id_string} {color}')
        self.run_x(f'hide #{res_model.id_string} models')
        self.run_x(f'hide #{res_map.id_string} models')
        self.run_x(f'transparency #{res_map.id_string} 70 s')
        return res_map, res_model

    def display_motif(self, map_path, model_path, matrix, map_color='#B2B2B2', model_color='#B2B2B2'):

        r_nr = int(self.active_residue)
        c = self.active_chain

        map_obj = self.open_model(map_path)
        self.clear_indices.append(map_obj.id_string)
        self.run_x(f'volume #{map_obj.id_string} sdLevel {SD_LEVEL} color {map_color}')
        model_obj = self.open_model(model_path)
        self.clear_indices.append(model_obj.id_string)
        self.run_x(f'color #{model_obj.id_string} {model_color}')
        # self.run_x(f'color #{model_obj.id_string} {map_color}')
        if matrix:
            self.run_x(f'view matrix mod #{model_obj.id_string},{matrix}')
            self.run_x(f'view matrix mod #{map_obj.id_string},{matrix}')
        else:
            if self.active_res_model_obj is None or self.active_res_map_obj is None:
                res = self.model_obj[0][c][r_nr]
                ma_path, mo_path = self.chop_residue(res)
                self.active_res_map_obj, self.active_res_model_obj = self.open_chopped_residue(ma_path, mo_path)
            self.run_x(f'align #{model_obj.id_string}@c,ca,n to #{self.active_res_model_obj.id_string}/{c}:{r_nr}@c,ca,n')
            self.run_x(f'view position #{map_obj.id_string} sameAsModels #{model_obj.id_string}')
            self.run_x(f'fitmap #{map_obj.id_string} inMap #{self.active_res_map_obj.id_string} metric correlation')
            self.run_x(f'view position #{model_obj.id_string} sameAsModels #{map_obj.id_string}')
        self.run_x(f'transparency #{map_obj.id_string} 70 s')
        return map_obj, model_obj

    def create_left_axis_y_label(self):
        pic = QtGui.QPicture()
        height = self.P.left_plots_H
        pic.setBoundingRect(QtCore.QRect(0, int(-height / 2), 20, int(height - 120)))
        y_label = QtGui.QPainter(pic)
        y_label.rotate(-90)
        document = QtGui.QTextDocument()
        rect = QtCore.QRectF(0, -height / 2.0, height, height * 2)
        document.setTextWidth(rect.width())
        txt_lst = self.P.y_axis_text
        document.setHtml(f'<font style="{self.P.grey_text} {self.P.left_plot_axis_font}">{txt_lst[0]} </font>'
                         f'<font style="color: black; {self.P.left_plot_axis_font}">{txt_lst[1]}</font>')
        document.drawContents(y_label, rect)
        y_label.end()
        return pic

    def set_paths(self, project_dir):
        folders = os.listdir(project_dir)

        tmp = os.path.join(project_dir, 'segments')
        if os.path.exists(tmp):
            self.segments_path = tmp
        for folder in folders:
            if folder.startswith('vs_'):
                path = os.path.join(project_dir, folder, 'scores_top.csv')
                if os.path.exists(path):
                    lib = folder.split('_')[-1]
                    self.csv_paths[lib] = path
        self.all_libs = list(sorted(self.csv_paths.keys()))
        inp = os.path.join(project_dir, 'input')
        if os.path.exists(inp):
            in_files = os.listdir(inp)
            for ext in ['.cif', '.mmcif', '.pdb', '.ent']:
                for f in in_files:
                    if f.endswith(ext):
                        self.model_path = os.path.join(inp, f)
                        break
            for ext in ['.map', '.mrc']:
                for f in in_files:
                    if f.endswith(ext):
                        self.map_path = os.path.join(inp, f)
                        break

    def set_data(self):
        if self.csv_paths:
            try:
                csv_path = self.csv_paths[self.lib_selector.currentText()]
            except (KeyError, AttributeError) as e:
                csv_path = self.csv_paths[sorted(self.csv_paths.keys())[0]]
            self.df = self.read_data(csv_path)
            self.calculate_z_scores()
            self.all_chains = self.get_chain_list(self.df)
            self.active_chain = self.all_chains[0]
            self.current_chain_data = self.get_chain_data(self.df, self.active_chain)
            self.active_residue = self.current_chain_data[self.k.RES_NR].iloc[0]
            self.right_plot_min_max = self.get_min_max_correlation()
            self.per_chain_mean_cc = self.get_per_chain_mean_cc_and_red_perc(self.P.red_difference)
            self.update_chain_cc_lbls()

    def read_data(self, csv_path):
        df = pd.read_csv(csv_path)
        df[self.k.TOP_SAME_TYPE_DIFF] = (df[self.k.TOP_CC] - df[self.k.SAME_TYPE_CC]) / df[self.k.TOP_CC]
        df[self.k.CHAIN] = df[self.k.CHAIN].apply(str)
        return df

    def calculate_z_scores(self):

        df = self.df
        nr_rows = df.shape[0]
        zeros = [0.0 for _ in range(nr_rows)]
        df[self.k.SAME_TYPE_Z] = zeros
        df[self.k.TYPE_TOP_Z] = zeros
        df[self.k.TOP_Z] = zeros
        df[self.k.STD] = zeros

        for index, row in df.iterrows():
            # print('Index', index)
            all_corr = []
            for name, value in row.iteritems():
                if name.endswith(self.k.TYPE_TOP_CC):
                    all_corr.append(value)
            # print(all_corr)
            all_corr = np.array(all_corr)

            std = all_corr.std()
            aver = all_corr.sum() / len(all_corr)
            df.at[index, self.k.MEAN] = aver
            df.at[index, self.k.STD] = std
            # print('std', std)
            # print('aver', aver)
            for name, value in row.iteritems():
                if name.endswith('_rscc'):
                    df.at[index, name[:-4]+'z'] = (value - aver) / std
                    # print(name, value)
                    # print(index, name[:-4]+'z', (value - aver) / std)
            # break
        df[self.k.TOP_SAME_TYPE_DIFF_Z] = (df[self.k.TOP_Z] - df[self.k.SAME_TYPE_Z]) / df[self.k.TOP_Z] * df[self.k.STD]
        # df.to_csv('/Volumes/data/test.csv')
    def get_chain_list(self, df):
        chains = df[self.k.CHAIN].unique()
        # chains = sorted(chains)
        # chains = [str(c) for c in chains]
        return sorted(chains)

    def get_chain_data(self, df, chain):
        out_df = df[df[self.k.CHAIN] == chain]
        return out_df

    def get_min_max_correlation(self):
        all_max = []
        all_min = []
        for name, value in self.df.iteritems():
            if name.endswith(self.k.TYPE_TOP_CC):
                all_max.append(value.max())
                all_min.append(value.min())
        return min(all_min), max(all_max)

    def get_per_chain_mean_cc_and_red_perc(self, max_diff):
        cc_dict = {}
        for chain in self.all_chains:
            chain_data = self.get_chain_data(self.df, chain)
            mean_same_type_cc = chain_data[self.k.SAME_TYPE_CC].mean()
            red_row_data = chain_data[chain_data[self.k.TOP_SAME_TYPE_DIFF] > max_diff]
            cc_dict[chain] = [round(mean_same_type_cc, 5), len(red_row_data) / len(chain_data) * 100]
        return cc_dict

    def chop_selected_residue(self):
        if not CHIMERA_MODE:
            self.no_chimera_message()
            return
        if self.active_residue_data is None:
            return
        r_nr = int(self.active_residue)
        c = self.active_chain
        res = self.model_obj[0][c][r_nr]
        res_map_path, res_model_path = self.chop_residue(res)
        res_map = self.open_model(res_map_path)
        self.clear_indices.append(res_map.id_string)
        self.run_x(f'volume #{res_map.id_string} sdLevel {SD_LEVEL}')
        res_model = self.open_model(res_model_path)
        self.clear_indices.append(res_model.id_string)
        self.run_x(f'size #{res_model.id_string} stickRadius 0.055')
        self.run_x(f'color #{res_map.id_string} #B2B2B2')
        self.run_x(f'transparency #{res_map.id_string} 70 s')
        self.run_x(f'hide #{res_model.id_string} models')

    def chop_selected_residue_env_sensitive(self):
        if not CHIMERA_MODE:
            self.no_chimera_message()
            return
        if self.active_residue_data is None:
            # print('No data yet')
            return
        r_nr = int(self.active_residue)
        c = self.active_chain
        res = self.model_obj[0][c][r_nr]
        res_map_path, res_model_path = self.chop_residue_environ_sensitive(res)
        res_map = self.open_model(res_map_path)
        self.clear_indices.append(res_map.id_string)
        self.run_x(f'volume #{res_map.id_string} sdLevel {SD_LEVEL}')
        self.run_x(f'color #{res_map.id_string} #B9BF6C')
        self.run_x(f'transparency #{res_map.id_string} 70 s')

    def load_all_rotamers(self):
        if not CHIMERA_MODE:
            self.no_chimera_message()
            return
        if self.active_residue_data is None:
            return
        self.clear_rotamers()
        r_nr = int(self.active_residue)
        c = self.active_chain
        res = self.model_obj[0][c][r_nr]
        res_type = res.get_resname().lower()
        files = os.listdir(self.active_lib_path)
        model_files = [file for file in files if file.startswith(res_type) and file.endswith('cif')]

        model_rmsds = []
        for model_file in model_files:
            path = os.path.join(self.active_lib_path, model_file)
            model = bu.load_structure(path)
            aligned_model = bu.superimpose_n_ca_c(res, model)
            rmsd = bu.calc_residue_static_pairwise_rmsd(res, aligned_model)
            model_rmsds.append((model_file, rmsd))
        model_rmsds = sorted(model_rmsds, key=lambda x: x[1])

        for model, _ in model_rmsds:
            model_path = os.path.join(self.active_lib_path, model)
            res_model = self.open_model(model_path)
            self.run_x(f'align #{res_model.id_string}@c,ca,n to #{self.model_chim.id_string}/{c}:{r_nr}@c,ca,n')
            self.run_x(f'size #{res_model.id_string} stickRadius 0.055')
            map_path = os.path.join(self.active_lib_path, model.split('.')[0] + '.mrc')
            res_map = self.open_model(map_path)
            self.run_x(f'volume #{res_map.id_string} sdLevel {SD_LEVEL}')
            self.run_x(f'view position #{res_map.id_string} sameAsModels #{res_model.id_string}')
            self.all_rotamers_id.append(res_map.id_string)
            self.all_rotamers_id.append(res_model.id_string)
            self.run_x(f'transparency #{res_map.id_string} 70 s')
            self.run_x(f'color #{res_model.id_string} rgb({",".join([str(i) for i in list(res_map.surfaces[0].color[:3])])})')

    def clear_rotamers(self):
        if not CHIMERA_MODE:
            self.no_chimera_message()
            return
        if self.all_rotamers_id:
            self.run_x(f'close #{",".join(self.all_rotamers_id)}')

    def chop_residue(self, res):
        res_type = res.get_resname().lower()
        chain_id = res.parent.id
        res_nr = res.id[1]
        name_prefix = '{}-{}-{}'.format(res_type, res_nr, chain_id)
        out_map_path = os.path.join(self.tmp, name_prefix + '.mrc')
        out_res_path = os.path.join(self.tmp, name_prefix + '.cif')
        if not os.path.exists(out_map_path) and not os.path.exists(out_res_path):
            cube_map_obj, shifts = self.chop.chop_cube(res, self.map_obj, 4, zero_origin=False)
            cube_map_obj.grid_resample_emda(0.25)
            cube_map_obj.write_map(os.path.join(self.tmp, name_prefix + 'cube.mrc'))
            struct = bu.residues2structure(res)
            bu.shift_coord(shifts, struct)
            bu.save_model(struct, out_res_path)
            for residue in struct.get_residues():
                bu.del_main_chain(residue)
            fin_map = self.chop.chop_soft_radius(struct, cube_map_obj, hard_radius=2, soft_radius=1)
            fin_map.write_map(out_map_path)
        return out_map_path, out_res_path

    def chop_residue_environ_sensitive(self, res):
        res_type = res.get_resname().lower()
        chain_id = res.parent.id
        res_nr = res.id[1]
        name_prefix = '{}-{}-{}_env'.format(res_type, res_nr, chain_id)
        out_map_path = os.path.join(self.tmp, name_prefix + '.mrc')
        out_res_path = os.path.join(self.tmp, name_prefix + '.cif')
        if not os.path.exists(out_map_path) and not os.path.exists(out_res_path):
            cube_map_obj, shifts = self.chop.chop_cube(res, self.map_obj, 4, zero_origin=False)
            cube_map_obj.grid_resample_emda(0.25)
            cube_map_obj.write_map(os.path.join(self.tmp, name_prefix + 'cube.mrc'))
            struct = bu.residues2structure(res)
            bu.shift_coord(shifts, struct)
            bu.save_model(struct, out_res_path)
            for residue in struct.get_residues():
                bu.del_main_chain(residue)
            fin_map = self.chop.chop_soft_radius_watershed(struct, cube_map_obj, self.model_obj, shifts, radius=2, soft_radius=1)
            fin_map.write_map(out_map_path)
        return out_map_path, out_res_path


    @staticmethod
    def first(trigger, trigger_data):
        print('trigger =', trigger)
        print('  trigger_data =', trigger_data)

    def view_residue(self):
        if not CHIMERA_MODE:
            self.no_chimera_message()
            return

        if self.active_residue_data is None:
            return

        r_nr = int(self.active_residue)
        c = self.active_chain
        # self.run_x(f'clipper isolate #{self.model_chim.id_string}/{c}:{r_nr}@ca focus true')
        # self.run_x('clipper spotlight')
        rr = None
        for ch in self.model_chim.chains:
            if ch.chain_id == c.upper():
                for r in ch.residues:
                    if r is None:
                        continue
                    if r.number == r_nr:
                        rr = r

        # self._new_camera_position(rr)
        from chimerax.isolde.navigate import get_stepper
        get_stepper(self.model_chim).step_to(rr)
        # self.run_x('wait 100')
        # self.local_session.triggers.activate_trigger('model position changed', self.model_chim)
        # import time
        # from chimerax.core import triggerset
        # ts = triggerset.TriggerSet()
        # ts.add_trigger('conrad')
        # h1 = ts.add_handler('conrad', self.first)
        # ts.activate_trigger('conrad', 1)
        self.model_chim.session.models[0].update()
        # delete the previous residue label and create a new one
        if self.sel_res_atomspec is not None:
            self.run_x(f'~label #{self.sel_res_atomspec}')
        self.sel_res_atomspec = f'{self.model_chim.id_string}/{c}:{r_nr}'
        self.run_x(f'label #{self.sel_res_atomspec} height 0.4')
        self.run_x('wait')
        # time.sleep(3)

        if self.active_lib_path is None:
            text = f"Sorry no motif library found for {self.lib_selector.currentText()} resolution range\n" \
                   f"in {self.libs_path}"
            self.show_mesage('Motif library missing', text)

        r = self.active_residue_data
        same_type_name = r[self.k.SAME_TYPE_NAME].tolist()[0]
        top_motif_name = r[self.k.M_TOP_NAME].tolist()[0]

        same_type = same_type_name.split('_')[0].lower()
        top_type = top_motif_name.split('_')[0].lower()

        try:
            same_type_matrix = r[self.k.SAME_TYPE_MATRIX].tolist()[0]
            top_motif_matrix = r[self.k.M_TOP_MATRIX].tolist()[0]
        except:
            same_type_matrix, top_motif_matrix = None, None

        if same_type_name == "ala_rotamer_1":
            same_type_name = 'ala_rotamer_single-conf'
        if top_motif_name == "ala_rotamer_1":
            top_motif_name = 'ala_rotamer_single-conf'


        # self.run_x('close #3-100')
        self.clear_residue_data()
        # self.run_x('show #1,2 models')
        # self.run_x(f'color #1 #00FFFF')


        if self.active_lib_path is not None:
            try:
                lib_same_map_path = os.path.join(self.active_lib_path, same_type_name + '.mrc')
                lib_same_model_path = os.path.join(self.active_lib_path, same_type_name + '.cif')
            except TypeError:
                pass
            else:
                if os.path.exists(lib_same_map_path) and os.path.exists(lib_same_model_path):
                    same_map, same_model = self.display_motif(lib_same_map_path, lib_same_model_path, same_type_matrix,
                                                              map_color=self.P.same_map_color,
                                                              model_color=self.P.same_model_color)
                    same_map.display = self.ui.same_map_checkBox.isChecked()
                    same_model.display = self.ui.same_model_checkBox.isChecked()
                    self.open_motifs.append(self.Motif(type=same_type, map_id=same_map.id_string, model_id=same_model.id_string,
                                                       show=self.ui.same_map_checkBox.isChecked(), mtp='same'))
                else:
                    text = f'Missing files in {self.active_lib_path}\n' \
                           f'{same_type_name + ".mrc"},  {same_type_name + ".cif"}'
                    self.show_mesage('Motif library missing', text)

            try:
                lib_top_map_path = os.path.join(self.active_lib_path, top_motif_name + '.mrc')
                lib_top_model_path = os.path.join(self.active_lib_path, top_motif_name + '.cif')
            except TypeError:
                pass
            else:
                if top_motif_name != same_type_name:
                    if os.path.exists(lib_top_map_path) and os.path.exists(lib_top_model_path):
                        top_map, top_model = self.display_motif(lib_top_map_path, lib_top_model_path,
                                                                top_motif_matrix,
                                                                map_color=self.P.top_map_color,
                                                                model_color=self.P.top_model_color)
                        top_map.display = self.ui.top_map_checkBox.isChecked()
                        top_model.display = self.ui.top_model_checkBox.isChecked()
                        self.open_motifs.append(self.Motif(type=top_type, map_id=top_map.id_string, model_id=top_model.id_string,
                                                           show=self.ui.top_map_checkBox.isChecked(), mtp='top'))

                    else:
                        text = f'Missing files in {self.active_lib_path}\n' \
                               f'{top_motif_name + ".mrc"},  {top_motif_name + ".cif"}'
                        self.show_mesage('Motif library missing', text)

        # self.run_x(f'view #{self.model_chim.id_string}/{c}:{r_nr}@ca pad 0.7')

        # self.run_x('size stickRadius 0.03')
        # self.run_x(f'zone #{self.model_chim.id_string}/{c}:{r_nr}@ca surfaceDistance 8 label false')
        # self.run_x('clip near -6')
        # self.run_x('clip far 6')

        # self.update_motifs_view()

    def update_motifs_view(self):
        if self.df is None:
            return
        same_ma = self.ui.same_map_checkBox.isChecked()
        same_mo = self.ui.same_model_checkBox.isChecked()
        top_ma = self.ui.top_map_checkBox.isChecked()
        top_mo = self.ui.top_model_checkBox.isChecked()

        for motif in self.open_motifs:
            if motif.mtp == 'same':
                motif.show = same_ma
                if same_ma:
                    self.run_x(f'show #{motif.map_id} models')
                else:
                    self.run_x(f'hide #{motif.map_id} models')
                if same_mo:
                    self.run_x(f'show #{motif.model_id} models')
                else:
                    self.run_x(f'hide #{motif.model_id} models')
            if motif.mtp == 'top':
                motif.show = top_ma
                if top_ma:
                    self.run_x(f'show #{motif.map_id} models')
                else:
                    self.run_x(f'hide #{motif.map_id} models')
                if top_mo:
                    self.run_x(f'show #{motif.model_id} models')
                else:
                    self.run_x(f'hide #{motif.model_id} models')
            if motif.mtp == 'other':
                self.run_x(f'hide #{motif.model_id} models')
                if not motif.show:
                    self.run_x(f'hide #{motif.map_id} models')
                else:
                    self.run_x(f'show #{motif.map_id} models')
        self.draw_right_plot()

    def set_labels_style(self):
        ui = self.ui
        for label in [ui.y_axis_range_label_1,
                      ui.y_axis_range_label_2,
                      ui.motif_label_1,
                      ui.motif_label_2,
                      ui.map_label,
                      ui.model_label,
                      ui.lib_label,
                      ui.chain_label,
                      ui.mean_score_label,
                      ui.in_row_label,
                      ui.metric_label]:
            label.setStyleSheet(f"QLabel {self.P.label_titles_style}")

class DictKeys:
    RES_TYPE = 'residue_type'
    RES_NR = 'residue_nr'
    CHAIN = 'chain'
    RES_MAP = 'residue_map'
    RES_MODEL = 'residue_model'
    M_TOP_TYPE = 'top_motif_type'
    M_TOP_NAME = 'top_motif_name'
    M_TOP_MATRIX = 'top_motif_matrix'
    TOP_CC = 'top_rscc'
    TOP_Z = 'top_z'
    SAME_TYPE_NAME = 'same_type_motif_name'
    SAME_TYPE_CC = 'same_type_motif_rscc'
    SAME_TYPE_Z = 'same_type_motif_z'
    SAME_TYPE_MATRIX = 'same_type_motif_matrix'
    # FULL NAME EXAMPLE "ala_cc" suffix for top correlation of specific residue type library motifs
    TYPE_TOP_CC = '_top_motif_rscc'
    TYPE_TOP_Z = '_top_motif_z'

    # FULL NAME EXAMPLE "ala_motif_name"
    TYPE_TOP_NAME = '_top_motif_name'
    TYPE_TOP_MATRIX = '_top_motif_matrix'
    SCORES = 'scores'
    ALL_TOP_MOTIFS = 'all_top_motifs'
    TOP_SAME_TYPE_DIFF = 'top_same_dif'
    TOP_SAME_TYPE_DIFF_Z = 'top_same_dif_z'
    COMPLETE_DATA = 'complete_data'
    MEAN = 'mean'
    STD = 'std'


class Parameters:
    def __init__(self):
        # Tolerance between maximum CC and claimed residue type CC
        # calculated by: (TOP_CC - SAME_TYPE_CC) / TOP_CC
        self.red_difference = 0.05
        # left plots default residues per row
        self.default_res_per_row = '60'
        # If claimed residue type RSCC is lower than it is coloured blue
        self.blue_threshold = 0.7
        self.default_metrics = 'Correlation'

        self.line_height = 28
        self.elem_spacing = 20
        self.selector_frame_size = 60
        self.main_window_width = 1100
        self.main_window_height = 1000
        self.left_plot_min_h = 100

        self.left_plots_H = 300# - self.elem_spacing *10
        # Plots colors
        self.outlier_c = '#CC3311'
        self.same_type_c = '#009998'
        self.very_low_c = '#0077BB'
        # self.very_low_c = 'b'

        # Left plots parameters
        # self.l_same_type_color = 'g'
        # self.l_top_type_color = 'r'
        # self.l_same_type_lowcorr_color = 'b'
        self.l_marker_size = 2.5
        self.l_picker_size = 2.5
        self.l_lbl_size = 6
        self.l_lbl_same_type_color = 'k'
        self.l_lbl_top_type_color = 'r'
        self.l_lbl_same_type_lowcorr_color = 'b'
        # self.l_ax_lbl_size = 4
        self.l_ax_txt_size = 10
        self.l_tick_params = {'labelsize': 5, 'length': 1.5, 'pad': 1, 'width': 0.4, 'xgap': 5}
        self.vline_color = '#6A6D70'
        self.vline_width = 0.7
        # Right plot parameters
        self.r_marker_size = 2.5
        self.r_picker_size = 2.5
        self.r_lbl_size = 8
        self.r_ax_lbl_size = 5
        self.r_ax_txt_size = 5
        self.r_ax_title_lbl_size = 8
        self.per_chain_cc_columns = 5
        self.r_tick_params = {'labelsize': 5, 'length': 1.5, 'pad': 1, 'width': 0.4}
        self.y_axis_text = ['Y axis:', 'Strudel score']
        self.x_axis_text = ['X axis:', 'Residue number']
        self.cc_text_size = 11
        self.grey_text = 'color: rgb(106, 109, 112);'
        self.left_plot_axis_font = 'font-family: Helvetica; font-size: 14px; font-weight: bold;'
        self.label_titles_style = '{color: rgb(34, 93, 153); font-size: 12px;}'
        self.horisontal_line_color = 'color: rgb(139,139,139);'

        # models display
        self.same_map_color = '#26FF53'
        self.same_model_color = '#26FF53'
        self.top_map_color = '#FF3837'
        self.top_model_color = '#FF3837'
        self.other_map_color = '#B2B2B2'
        self.other_model_color = '#B2B2B2'


if __name__ == '__main__':
    import sys
    app = QtWidgets.QApplication(sys.argv)
    session = ''
    tool_name = 'Strudel Score'
    window = StrudelScore(session, tool_name)
    app.exec_()
