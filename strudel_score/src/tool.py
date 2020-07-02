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



CHIMERA_MODE = True
try:
    from chimerax.core.tools import ToolInstance
    from chimerax.core.commands import Command
except ImportError:
    CHIMERA_MODE = False
from PyQt5 import uic
from math import fabs
import shutil
import os
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QMessageBox, QVBoxLayout, QHBoxLayout, QGridLayout, QFileDialog, QLabel, QPushButton, QMainWindow
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.ticker import FormatStrFormatter
import pandas as pd
import numpy as np

from strudel.utils import bioUtils as bu
from strudel.chop.chopMap import ChopMap
from strudel.parse.mapParser import MapParser
from strudel import nomenclature

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
    display_name = "Strudel Tool"

    def __init__(self, session, tool_name):
        Base.__init__(self, session, tool_name)
        self.local_session = session
        if CHIMERA_MODE:
            from chimerax.ui.gui import MainToolWindow
            self.tool_window = MainToolWindow(self)
            self.ui_area = self.tool_window.ui_area
        else:
            self.ui_area = QMainWindow()

        self.basedir = os.path.dirname(os.path.abspath(__file__))
        self.config_path = os.path.join(self.basedir, 'config.txt')
        self.libs_path = None

        self.close_to_ligand = {}
        self.axvl = None
        self.axvtxt = None
        self.ax_residues = []
        self.tmp = None
        self.P = Parameters()
        self.k = DictKeys()
        self.line_h = 18
        self.files_lbls = []

        # Layouts
        self.left_plots_vbox = None
        self.mean_cc_grid = None
        # Figures
        self.right_plot_fig = None
        self.right_ax = None
        self.right_plot_canvas = None
        # Widgets
        self.chain_selector = None
        self.left_axis_selector = None
        self.lib_selector = None
        self.current_chain_data = None
        self.ui = None
        self.init_ui()
        if CHIMERA_MODE:
            self.integrate()

        self.active_lib_path = None
        self.map_path = None
        self.map_chim = None
        self.map_obj = None
        self.model_obj = None
        self.model_chim = None
        self.model_path = None
        self.segments_path = None
        self.all_rotamers_id = []
        self.chop = ChopMap()
        self.tmp = '/tmp/strudel'
        if not os.path.exists(self.tmp):
            os.makedirs(self.tmp)

        self.csv_paths = {}

        self.df = None
        self.all_chains = ['No chains found']
        self.all_libs = ['No data loaded']
        self.y_axis_ranges = ['Expanded', '0 to 1', '0.3 to 1', '0.5 to 1', '0.6 to 1', '0.7 to 1', '0.8 to 1']
        self.per_chain_mean_cc = None
        self.active_chain = None
        self.active_residue = None
        self.active_residue_data = None
        self.right_plot_size = None

        self.inject_data()
        self.ui.show()
        if self.df is not None:
            self.draw_right_plot()
            self.draw_left_plot_vline()
        self.ui.show()

    def init_ui(self):
        self.mw = QtWidgets.QMainWindow()
        uifile = os.path.join(self.basedir, 'ui', 'qt_ui.ui')
        self.ui = uic.loadUi(uifile, self.mw)
        self.ui.setWindowTitle('Strudel score')
        self.setup_menubar()
        self.chain_selector = self.ui.chain_comboBox
        self.left_axis_selector = self.ui.y_axis_range_comboBox
        self.lib_selector = self.ui.lib_comboBox
        self.mean_cc_grid = self.ui.mean_score_gridLayout
        self.mean_cc_grid.setContentsMargins(0, 0, 0, 0)
        self.create_left_plots_container()
        self.right_plot_fig = Figure()
        self.right_ax = self.right_plot_fig.add_subplot(111)
        self.right_plot_canvas = FigureCanvas(self.right_plot_fig)
        self.right_plot_fig.subplots_adjust(left=0.25, top=0.95, bottom=0.08)
        self.create_right_plot_container()
        w = self.ui.width()
        self.ui.plots_splitter.setSizes([int(w * .8), int(w * .2)])

    def integrate(self):
        parent = self.ui_area
        layout = QtWidgets.QHBoxLayout()
        layout.addWidget(self.mw)
        layout.setStretchFactor(self.mw, 0)
        layout.setContentsMargins(0, 0, 0, 0)
        parent.setLayout(layout)
        self.tool_window.manage(placement=None)

    def inject_data(self):
        print('Inject data called')
        self.fill_selectors()
        self.update_chain_cc_lbls()
        self.fill_files_data()
        if self.df is not None:
            # self.find_ligand_close_contacts()
            self.update_left_plots()
            self.draw_right_plot()
            self.draw_left_plot_vline()

    def fill_files_data(self):
        self.ui.map_label.setText(f'<font style="{self.P.titles_style}">Map:</font>')
        self.ui.model_label.setText(f'<font style="{self.P.titles_style}">Model:</font>')
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

    def set_library(self):
        lib_path = str(QFileDialog.getExistingDirectory())
        self.libs_path = lib_path
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

    def read_config(self):
        path = None
        try:
            with open(self.config_path, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    if 'lib_path' in line:
                        path = line.split('=')[-1].strip()
                        print(path)
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

    def open_project(self):
        proj_path = None
        set_path = str(QFileDialog.getExistingDirectory())
        folders = os.listdir(set_path)
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
            self.right_plot_size = (self.right_plot_canvas.width(), self.right_plot_canvas.height())
            self.set_data()
            self.load_map_model()
            self.inject_data()
            self.clear_tmp()
            self.ui.show()
            self.read_config()
            self.set_active_lib()

        else:
            self.show_mesage('Fail', 'The specified directory does not contain Strudel validation results files!')

    def set_active_lib(self):
        avail_libs = os.listdir(self.libs_path)
        lib_path = None
        for l in avail_libs:
            if l.endswith(self.lib_selector.currentText()):
                lib_path = os.path.join(self.libs_path, l, 'motifs')
        self.active_lib_path = lib_path

    def load_map_model(self):
        if CHIMERA_MODE:
            self.run_x('close all')
            self.map_chim = self.run_x(f'open {self.map_path}')
            self.run_x(f'volume #{self.map_chim.id_string} sdLevel {SD_LEVEL + 2} style mesh step 1')
            self.model_chim = self.run_x(f'open {self.model_path}')

            self.run_x('style stick')
            self.run_x('size stickRadius 0.03')
            self.run_x('~ribbon')
            self.run_x('show atoms')
        self.model_obj = bu.load_structure(self.model_path)
        self.map_obj = MapParser(self.map_path)

    def find_ligand_close_contacts(self, d=2):
        ligands = []
        for chain in self.model_obj[0]:
            tmp = []
            for res in chain:
                r_name = res.get_resname().upper()
                if r_name not in nomenclature.AMINOACID_RESIDUES:
                    ligands.append(res)

        for chain in self.model_obj[0]:
            tmp = []
            for res in chain:
                r_name = res.get_resname().upper()
                if r_name in nomenclature.AMINOACID_RESIDUES:
                    for atom in res:
                        for res2 in ligands:
                            for atom2 in res2:
                                dist = atom - atom2
                                if dist < d:
                                    tmp.append(res.get_id()[1])
            tmp = list(set(tmp))
            self.close_to_ligand[chain.get_id()] = tmp

    def clear_tmp(self):
        shutil.rmtree(self.tmp)
        os.makedirs(self.tmp)

    def run_x(self, command):
        if CHIMERA_MODE:
            rc = Command(self.local_session)
            print(f'Running: {command}')
            try:
                return rc.run(command)[0][0]
            except TypeError:
                return None
        else:
            self.no_chimera_message()

    def no_chimera_message(self):
        text = 'Functionality not available in standalone mode.' \
               '\nPlease run as ChimeraX Tool'
        self.show_mesage('Functionality unavailable0', text)

    def fill_left_plot_axix_range_selector(self):
        label = self.ui.y_axis_range_label
        label.setText(f'<font style="{self.P.titles_style}">Y axis range</font>')
        combo_box = self.ui.y_axis_range_comboBox
        combo_box.addItems(self.y_axis_ranges)
        if self.df is not None:
            combo_box.currentIndexChanged.connect(self.switch_axis_range)

    def switch_axis_range(self):
        if self.df is not None:
            self.update_left_plots(self.ui.y_axis_range_comboBox.currentText())

    def update_chain_cc_lbls(self):
        self.ui.mean_score_label.setText(f'<font style="{self.P.titles_style}">Mean per chain strudel score</font>')
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
                       f'<font style="font-size: {self.P.cc_text_size}px">{value[0]:.5f} ({value[1]:.1f}% red)</font>'
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
        print("fill selectors called")
        # Lib selector
        try:
            self.ui.lib_comboBox.currentIndexChanged.disconnect(self.switch_lib)
        except TypeError:
            pass
        self.ui.lib_comboBox.clear()
        self.ui.lib_comboBox.addItems(self.all_libs)
        self.ui.lib_comboBox.currentIndexChanged.connect(self.switch_lib)
        self.ui.lib_label.setText(f'<font style="{self.P.titles_style}">Select motif library</font>')
        # Chain selector
        try:
            self.ui.chain_comboBox.currentIndexChanged.disconnect(self.switch_chain)
        except TypeError:
            pass
        self.ui.chain_comboBox.clear()
        print('all chains (fill selectors)', self.all_chains)
        self.ui.chain_comboBox.addItems(self.all_chains)

        self.ui.chain_comboBox.currentIndexChanged.connect(self.switch_chain)
        self.ui.chain_label.setText(f'<font style="{self.P.titles_style}">Select polypeptide chain</font>')
        # Y range selector
        combo_box = self.ui.y_axis_range_comboBox
        if combo_box.currentText() not in self.y_axis_ranges:
            label = self.ui.y_axis_range_label
            label.setText(f'<font style="{self.P.titles_style}">Y axis range</font>')
            combo_box.addItems(self.y_axis_ranges)
        if self.df is not None:
            combo_box.currentIndexChanged.connect(self.switch_axis_range)

    def do_nothing(self):
        pass

    def switch_lib(self):
        self.set_data()
        self.set_active_lib()
        self.update_chain_cc_lbls()
        self.switch_chain()
        self.active_residue = self.current_chain_data[self.k.RES_NR].iloc[0]
        self.draw_left_plot_vline()
        self.ui.show()

    def switch_chain(self):
        self.active_chain = self.chain_selector.currentText()
        self.current_chain_data = self.get_chain_data(self.df, self.active_chain)
        self.update_left_plots(self.left_axis_selector.currentText())
        self.active_residue = self.current_chain_data[self.k.RES_NR].iloc[0]
        self.draw_right_plot()
        self.draw_left_plot_vline()

    def create_left_plots_container(self):
        ui = self.ui
        grid = QGridLayout()
        y_label = ui.y_left_plot_label
        y_label.setPicture(self.create_left_axis_y_label())
        y_label.show()

        scroll_area = ui.lef_plot_scrollArea

        x_label = ui.x_left_plot_label
        txt_lst = self.P.x_axis_text
        x_label_text = f'<font style="{self.P.grey_text} {self.P.left_plot_axis_font}">{txt_lst[0]} </font>' \
                       f'<font style="{self.P.left_plot_axis_font}">{txt_lst[1]}</font>'
        x_label.setText(x_label_text)
        x_label.setAlignment(QtCore.Qt.AlignHCenter)
        grid.setSpacing(5)
        grid.addWidget(y_label, 0, 0)
        grid.addWidget(scroll_area, 0, 1)
        grid.addWidget(x_label, 1, 1)
        grid.setContentsMargins(12, 6, 0, 12)
        lf = ui.left_plot_frame
        ui.left_plot_frame.setLayout(grid)
        self.left_plots_vbox = self.ui.left_plot_scrollarea_verticalLayout
        p = lf.palette()
        lf.setAutoFillBackground(True)
        p.setColor(lf.backgroundRole(), QtCore.Qt.white)
        lf.setPalette(p)

    def create_right_plot_container(self):
        p = self.P.r_tick_params
        self.right_ax.tick_params(labelsize=p['labelsize'], length=p['length'], pad=p['pad'], width=p['width'])
        [self.right_ax.spines[i].set_linewidth(p['width']) for i in self.right_ax.spines.keys()]
        frame = self.ui.right_plot_frame
        vlayout = QVBoxLayout()
        hlayout = QHBoxLayout()
        left_arrow = QPushButton('<')
        right_arrow = QPushButton('>')
        hlayout.addStretch()
        hlayout.addWidget(left_arrow)
        hlayout.addWidget(self.ui.right_plot_pushButton)
        hlayout.addWidget(right_arrow)
        hlayout.addStretch()
        vlayout.addWidget(self.right_plot_canvas)
        vlayout.addLayout(hlayout)
        self.ui.right_plot_pushButton.clicked.connect(self.display_in_chimera)
        left_arrow.clicked.connect(self.navigate_residue_left)
        right_arrow.clicked.connect(self.navigate_residue_right)
        frame.setLayout(vlayout)
        vlayout.setSpacing(0)
        vlayout.setContentsMargins(0, 0, 0, 0)

        p = frame.palette()
        frame.setAutoFillBackground(True)
        p.setColor(frame.backgroundRole(), QtCore.Qt.white)
        frame.setPalette(p)
        frame.show()

    def navigate_residue(self, direction):
        if self.current_chain_data is None:
            return
        avail_res = self.current_chain_data[self.k.RES_NR].tolist()
        index = avail_res.index(self.active_residue)
        if direction == 'left':
            if index > 0:
                self.active_residue = avail_res[index-1]
            self.draw_right_plot()
            self.draw_left_plot_vline()
        if direction == 'right':
            if index < len(avail_res):
                self.active_residue = avail_res[index + 1]
            self.draw_right_plot()
            self.draw_left_plot_vline()
        else:
            pass

    def navigate_residue_left(self):
        self.navigate_residue('left')

    def navigate_residue_right(self):
        self.navigate_residue('right')

    def update_left_plots(self, y_range='Expanded'):
        self.current_chain_data = self.get_chain_data(self.df, self.active_chain)
        all_left = self.all_left_plots(self.current_chain_data, y_range=y_range)
        # Clear left plots
        for i in reversed(range(self.left_plots_vbox.count())):
            self.left_plots_vbox.itemAt(i).widget().setParent(None)

        for plot in all_left:
            plot.setMinimumHeight(self.P.left_plot_min_h)
            self.left_plots_vbox.addWidget(plot)

    def all_left_plots(self, data, max_res=80, max_diff=0.05, min_corr=0.7, y_range='Expanded'):
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
            # y_max = min(max(top_cc_list) + 0.1, y_max)
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
            ax.set_ylim(ax_lim)
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')
            ax.yaxis.set_label_position("right")
            p = self.P.l_tick_params
            ax.tick_params(labelsize=p['labelsize'], length=p['length'], pad=p['pad'], width=p['width'])
            [ax.spines[i].set_linewidth(p['width']) for i in ax.spines.keys()]
            ax.set_xticks([i for i in range(x_start, x_end + 1, p['xgap'])])

            green_row_data = data[(x_start <= data[self.k.RES_NR]) & (data[self.k.RES_NR] < x_end)]
            green_x = green_row_data[self.k.RES_NR].tolist()
            green_y = green_row_data[self.k.SAME_TYPE_CC].tolist()
            self.ax_residues.append([ax, fig, green_x])

            red_row_data = green_row_data[green_row_data[self.k.TOP_SAME_TYPE_DIFF] > max_diff]
            red_x = red_row_data[self.k.RES_NR].tolist()
            red_y = red_row_data[self.k.TOP_CC].tolist()

            blue_row_data = green_row_data[green_row_data[self.k.SAME_TYPE_CC] < min_corr]
            blue_x = blue_row_data[self.k.RES_NR].tolist()
            blue_y = blue_row_data[self.k.SAME_TYPE_CC].tolist()

            ax.plot(green_x, green_y, 'o', c=self.P.l_same_type_color,
                    markersize=self.P.l_marker_size,
                    picker=self.P.l_picker_size)
            ax.plot(red_x, red_y, 'o', c=self.P.l_top_type_color,
                    markersize=self.P.l_marker_size,
                    picker=self.P.l_picker_size)
            ax.plot(blue_x, blue_y, 'o', c=self.P.l_same_type_lowcorr_color,
                    markersize=self.P.l_marker_size,
                    picker=self.P.l_picker_size)
            try:
                if self.close_to_ligand[self.active_chain]:
                    rl_x = [r for r in self.close_to_ligand[self.active_chain] if r in green_x]
                    ax.plot(rl_x, [ax_lim[0]+0.05 for i in range(len(rl_x))], '*', c='b')
            except KeyError:
                pass

            for index, row in red_row_data.iterrows():
                ax.text(row[self.k.RES_NR], row[self.k.TOP_CC], row[self.k.M_TOP_TYPE].upper(),
                        color=self.P.l_top_type_color,
                        size=txt_size, horizontalalignment='left', verticalalignment='bottom', rotation=50)

                ax.text(row[self.k.RES_NR], row[self.k.SAME_TYPE_CC], row[self.k.RES_TYPE].upper(),
                        color=self.P.l_lbl_same_type_color,
                        size=txt_size, horizontalalignment='right', verticalalignment='top', rotation=50)
            for index, row in blue_row_data.iterrows():
                ax.text(row[self.k.RES_NR], row[self.k.SAME_TYPE_CC], row[self.k.RES_TYPE].upper(),
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
        ydata = thisline.get_ydata()
        ind = event.ind
        points = tuple(zip(xdata[ind], ydata[ind]))
        self.active_residue = xdata[ind][0]
        self.draw_right_plot()
        self.draw_left_plot_vline()

    def draw_right_plot(self, max_res_display=20):
        if self.right_plot_size is None:
            width = self.right_plot_canvas.width()
            height = self.right_plot_canvas.height()
        else:
            width = self.right_plot_size[0]
            height = self.right_plot_size[1]
            self.right_plot_size = None

        self.right_ax.clear()
        r_nr = self.active_residue
        text_size = self.P.r_lbl_size
        xlim = 1
        r = self.current_chain_data[self.current_chain_data[self.k.RES_NR] == r_nr]
        self.active_residue_data = r
        same_type = r[self.k.RES_TYPE].tolist()[0].upper()
        y_data = []
        for name, value in r.iteritems():
            if name.endswith(self.k.TYPE_TOP_CC):
                y_data.append((value.tolist()[0], name.split('_')[0].upper()))

        y_data.sort(key=lambda v: v[0], reverse=True)
        y_data = y_data[:max_res_display]
        y_lst = [v[0] for v in y_data]
        y_txt = [v[1] for v in y_data]
        y_delta = max(y_lst) - min(y_lst)

        text_h = text_size * y_delta / height * 3
        text_w = text_size * xlim / width * 4
        space = text_h * 0.8
        y_lbls = self.find_y_label_coordinates(y_lst, space)
        # Consider using dash option for labels
        for i in range(len(y_lst)):
            if same_type == y_txt[i]:
                self.right_ax.plot([r_nr - 1.3 * text_w, r_nr], [y_lbls[i], y_lst[i]], '--', linewidth=0.3, c='g')
                self.right_ax.plot(r_nr, y_lst[i], 'o', markersize=2.5, c='g')
                self.right_ax.text(r_nr - 2.5 * text_w, y_lbls[i] - text_h / 4, y_txt[i], color='g', size=text_size)
            else:
                self.right_ax.plot([r_nr - 1.3 * text_w, r_nr], [y_lbls[i], y_lst[i]], '--', linewidth=0.3, c='r')
                self.right_ax.plot(r_nr, y_lst[i], 'o', markersize=2.5, c='r')
                self.right_ax.text(r_nr - 2.5 * text_w, y_lbls[i] - text_h / 4, y_txt[i], color='r', size=text_size)
        self.right_ax.set_xlim(r_nr - 0.7, r_nr + 0.3)
        self.right_ax.set_xticks([r_nr])
        self.right_ax.set_xlabel(self.P.x_axis_text[1], size=self.P.r_ax_title_lbl_size)
        self.right_ax.set_ylabel(self.P.y_axis_text[1], size=self.P.r_ax_title_lbl_size)
        self.right_ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        self.right_plot_fig.subplots_adjust(left=0.25, top=0.98, bottom=0.08)
        self.right_plot_fig.canvas.draw_idle()

    @staticmethod
    def find_y_label_coordinates(lst, min_d):
        """
        Takes a list of numbers and makes the distance between any adjacent
        number greater than min_d. If the output values form clusters
        :param lst:
        :param min_d:
        :return:
        """
        # Create clusters of points
        clusters = []
        cluster = [lst[0]]
        for i, val in enumerate(lst[:-1]):
            if fabs(val - lst[i + 1]) < min_d:
                cluster.append(lst[i + 1])
            else:
                clusters.append(cluster)
                cluster = [lst[i + 1]]
        clusters.append(cluster)

        # Merge clusters which overlap
        if len(clusters) > 2:
            tmp = clusters
            while True:
                merged = []
                for i, c1 in enumerate(tmp[:-1]):
                    c2 = tmp[i + 1]
                    d_c1 = min(c1) + (max(c1) - min(c1)) / 2
                    d_c2 = min(c2) + (max(c2) - min(c2)) / 2
                    d = fabs(d_c1 - d_c2)

                    if d < (len(c1) / 2 + len(c2) / 2) * min_d:
                        merged.append(c1 + c2)
                        merged = merged + tmp[i + 2:]
                        break
                    else:
                        merged.append(c1)
                else:
                    merged.append(tmp[-1])
                if merged == tmp:
                    break
                else:
                    tmp = merged
        else:
            merged = clusters

        # Calculate label positions
        out = []
        for cluster in merged:
            ln = len(cluster)
            if ln > 1:
                mdl_v = min(cluster) + (max(cluster) - min(cluster)) / 2
                for i, v in enumerate(cluster):
                    if v > mdl_v:
                        val = mdl_v + (ln / 2 - i - 0.5) * min_d
                    else:
                        val = mdl_v - (i - ln / 2 + 0.5) * min_d
                    out.append(val)
            else:
                out.append(cluster[0])
        return out

    def create_left_axis_y_label(self):
        pic = QtGui.QPicture()
        height = self.P.left_plots_H
        pic.setBoundingRect(QtCore.QRect(0, -height / 2, 20, height - 120))
        y_label = QtGui.QPainter(pic)
        y_label.rotate(-90)
        document = QtGui.QTextDocument()
        rect = QtCore.QRectF(0, -height / 2.0, height, height * 2)
        document.setTextWidth(rect.width())
        txt_lst = self.P.y_axis_text
        document.setHtml(f'<font style="{self.P.grey_text} {self.P.left_plot_axis_font}">{txt_lst[0]} </font>'
                         f'<font style="{self.P.left_plot_axis_font}">{txt_lst[1]}</font>')
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
            self.all_chains = self.get_chain_list(self.df)
            print('All chains (set data)', self.all_chains)
            self.active_chain = self.all_chains[0]
            print('active chain (set data)', self.active_chain)
            self.current_chain_data = self.get_chain_data(self.df, self.active_chain)
            self.active_residue = self.current_chain_data[self.k.RES_NR].iloc[0]
            self.right_plot_min_max = self.get_min_max_correlation()
            self.per_chain_mean_cc = self.get_per_chain_mean_cc_and_red_perc()
            self.update_chain_cc_lbls()

    def read_data(self, csv_path):
        df = pd.read_csv(csv_path)
        df[self.k.TOP_SAME_TYPE_DIFF] = (df[self.k.TOP_CC] - df[self.k.SAME_TYPE_CC]) / df[self.k.TOP_CC]
        df[self.k.CHAIN] = df[self.k.CHAIN].apply(str)
        return df

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

    def get_per_chain_mean_cc_and_red_perc(self, max_diff=0.05):
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
            print('No data yet')
            return
        r_nr = int(self.active_residue)
        c = self.active_chain
        res = self.model_obj[0][c][r_nr]
        res_map_path, res_model_path = self.chop_residue(res)
        res_map = self.run_x(f'open {res_map_path}')
        self.run_x(f'volume #{res_map.id_string} sdLevel {SD_LEVEL}')
        res_model = self.run_x(f'open {res_model_path}')
        self.run_x(f'size #{res_model.id_string} stickRadius 0.055')
        self.run_x(f'color #{res_map.id_string} #B2B2B2')
        self.run_x(f'align #{res_model.id_string}@c,ca,n to #{self.model_chim.id_string}/{c}:{r_nr}@c,ca,n')
        self.run_x(f'view position #{res_map.id_string} sameAsModels #{res_model.id_string}')
        self.run_x(f'transparency #{res_map.id_string} 70 s')
        self.run_x(f'hide #{res_model.id_string} models')

    def chop_selected_residue_env_sensitive(self):
        if not CHIMERA_MODE:
            self.no_chimera_message()
            return
        if self.active_residue_data is None:
            print('No data yet')
            return
        r_nr = int(self.active_residue)
        c = self.active_chain
        res = self.model_obj[0][c][r_nr]
        res_map_path, res_model_path = self.chop_residue_environ_sensitive(res)
        res_map = self.run_x(f'open {res_map_path}')
        self.run_x(f'volume #{res_map.id_string} sdLevel {SD_LEVEL}')
        # res_model = self.run_x(f'open {res_model_path}')
        self.run_x(f'color #{res_map.id_string} #B9BF6C')
        self.run_x(f'transparency #{res_map.id_string} 70 s')

    def load_all_rotamers(self):
        if not CHIMERA_MODE:
            self.no_chimera_message()
            return
        if self.active_residue_data is None:
            print('No data yet')
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
            res_model = self.run_x(f'open {model_path}')
            self.run_x(f'align #{res_model.id_string}@c,ca,n to #{self.model_chim.id_string}/{c}:{r_nr}@c,ca,n')
            self.run_x(f'size #{res_model.id_string} stickRadius 0.055')
            map_path = os.path.join(self.active_lib_path, model.split('.')[0] + '.mrc')
            res_map = self.run_x(f'open {map_path}')
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
            cube_map_obj, shifts = self.chop.chop_cube_list(res, self.map_obj, 4)
            cube_map_obj.grid_resample_emda(0.25)
            cube_map_obj.write_map(os.path.join(self.tmp, name_prefix + 'cube.mrc'))
            struct = bu.residues2structure(res)
            bu.shift_coord(shifts, struct)
            bu.save_model(struct, out_res_path)
            for residue in struct.get_residues():
                bu.del_main_chain(residue)
            fin_map = self.chop.chop_soft_radius(struct, cube_map_obj, radius=2, soft_radius=1)
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
            cube_map_obj, shifts = self.chop.chop_cube_list(res, self.map_obj, 4, zero_origin=False)
            cube_map_obj.grid_resample_emda(0.25)
            cube_map_obj.write_map(os.path.join(self.tmp, name_prefix + 'cube.mrc'))
            struct = bu.residues2structure(res)
            bu.shift_coord(shifts, struct)
            bu.save_model(struct, out_res_path)
            for residue in struct.get_residues():
                bu.del_main_chain(residue)
            fin_map = self.chop.chop_soft_radius4(struct, cube_map_obj, self.model_obj, shifts, radius=2, soft_radius=1)
            fin_map.write_map(out_map_path)
        return out_map_path, out_res_path

    def display_in_chimera(self):
        if not CHIMERA_MODE:
            self.no_chimera_message()
            return

        if self.active_residue_data is None:
            return

        if self.active_lib_path is None:
            text = f"Sorry no motif library found for {self.lib_selector.currentText()} resolution range\n" \
                   f"in {self.libs_path}"
            self.show_mesage('Motif library missing', text)

        r = self.active_residue_data
        same_type_name = r[self.k.SAME_TYPE_NAME].tolist()[0]
        top_motif_name = r[self.k.M_TOP_NAME].tolist()[0]
        try:
            same_type_matrix = r[self.k.SAME_TYPE_MATRIX].tolist()[0]
            top_motif_matrix = r[self.k.M_TOP_MATRIX].tolist()[0]
        except:
            same_type_matrix, top_motif_matrix = None, None

        if same_type_name == "ala_rotamer_1":
            same_type_name = 'ala_rotamer_single-conf'
        if top_motif_name == "ala_rotamer_1":
            top_motif_name = 'ala_rotamer_single-conf'

        # superimpose
        r_nr = int(self.active_residue)
        c = self.active_chain
        self.run_x('close #3-10')
        self.run_x('show #1,2 models')
        self.run_x(f'color #1 #00FFFF')
        self.run_x(f'label #{self.model_chim.id_string}/{c}:{r_nr} height 0.4')

        if self.active_lib_path is not None:
            lib_same_map_path = os.path.join(self.active_lib_path, same_type_name + '.mrc')
            lib_same_model_path = os.path.join(self.active_lib_path, same_type_name + '.cif')

            if not same_type_matrix or not top_motif_matrix:
                res = self.model_obj[0][c][r_nr]
                res_map_path, res_model_path = self.chop_residue(res)
                res_map = self.run_x(f'open {res_map_path}')
                self.run_x(f'volume #{res_map.id_string} sdLevel {SD_LEVEL}')
                res_model = self.run_x(f'open {res_model_path}')
                self.run_x(f'color #{res_map.id_string} #B2B2B2')
                self.run_x(f'align #{res_model.id_string}@c,ca,n to #{self.model_chim.id_string}/{c}:{r_nr}@c,ca,n')
                self.run_x(f'view position #{res_map.id_string} sameAsModels #{res_model.id_string}')


            if os.path.exists(lib_same_map_path) and os.path.exists(lib_same_model_path):
                same_map = self.run_x(f'open {lib_same_map_path}')
                self.run_x(f'volume #{same_map.id_string} sdLevel {SD_LEVEL} color #26FF53')
                same_model = self.run_x(f'open {lib_same_model_path}')
                self.run_x(f'color #{same_model.id_string} #26FF53')
                if same_type_matrix:
                    self.run_x(f'view matrix mod #{same_model.id_string},{same_type_matrix}')
                    self.run_x(f'view matrix mod #{same_map.id_string},{same_type_matrix}')
                else:
                    self.run_x(f'align #{same_model.id_string}@c,ca,n to #{res_model.id_string}/{c}:{r_nr}@c,ca,n')
                    self.run_x(f'view position #{same_map.id_string} sameAsModels #{same_model.id_string}')
                    self.run_x(f'fitmap #{same_map.id_string} inMap #{res_map.id_string} metric correlation')
                    self.run_x(f'view position #{same_model.id_string} sameAsModels #{same_map.id_string}')

                self.run_x(f'transparency #{same_map.id_string} 70 s')
                self.run_x(f'hide #{same_model.id_string} models')
                self.run_x(f'hide #{same_map.id_string} models')
            else:
                text = f'Missing files in {self.active_lib_path}\n' \
                       f'{same_type_name + ".mrc"},  {same_type_name + ".cif"}'
                self.show_mesage('Motif library missing', text)

            lib_top_map_path = os.path.join(self.active_lib_path, top_motif_name + '.mrc')
            lib_top_model_path = os.path.join(self.active_lib_path, top_motif_name + '.cif')
            if top_motif_name != same_type_name:
                if os.path.exists(lib_top_map_path) and os.path.exists(lib_top_model_path):
                    top_map = self.run_x(f'open {lib_top_map_path}')
                    self.run_x(f'volume #{top_map.id_string} sdLevel {SD_LEVEL} color #FF3837')
                    top_model = self.run_x(f'open {lib_top_model_path}')
                    self.run_x(f'color #{top_model.id_string} #FF3837')
                    if top_motif_matrix:
                        self.run_x(f'view matrix mod #{top_model.id_string},{top_motif_matrix}')
                        self.run_x(f'view matrix mod #{top_map.id_string},{top_motif_matrix}')
                    else:
                        self.run_x(f'align #{top_model.id_string}@c,ca,n to #{res_model.id_string}/{c}:{r_nr}@c,ca,n')
                        self.run_x(f'view position #{top_map.id_string} sameAsModels #{top_model.id_string}')

                        self.run_x(f'fitmap #{top_map.id_string} inMap #{res_map.id_string} metric correlation')
                        self.run_x(f'view position #{top_model.id_string} sameAsModels #{top_map.id_string}')
                    self.run_x(f'transparency #{top_map.id_string} 70 s')
                    self.run_x(f'hide #{top_model.id_string} models')
                    self.run_x(f'hide #{top_map.id_string} models')
                else:
                    text = f'Missing files in {self.active_lib_path}\n' \
                           f'{top_motif_name + ".mrc"},  {top_motif_name + ".cif"}'
                    self.show_mesage('Motif library missing', text)

        # self.run_x(f'hide #{res_map.id_string} models')
        self.run_x(f'view #{self.model_chim.id_string}/{c}:{r_nr}@ca pad 0.7')
        # self.run_x('cofr #{}@ca'.format(res_model.id_string))
        self.run_x('size stickRadius 0.03')
        try:
            self.run_x(f'transparency #{res_map.id_string} 70 s')
            self.run_x(f'hide #{res_model.id_string} models')
        except:
            pass
        self.run_x(f'zone #{self.model_chim.id_string}/{c}:{r_nr}@ca surfaceDistance 8 label false')
        self.run_x('clip near -6')
        self.run_x('clip far 6')


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
    SAME_TYPE_NAME = 'same_type_motif_name'
    SAME_TYPE_CC = 'same_type_motif_rscc'
    SAME_TYPE_MATRIX = 'same_type_motif_matrix'
    # FULL NAME EXAMPLE "ala_cc" suffix for top correlation of specific residue type library motifs
    TYPE_TOP_CC = '_top_motif_rscc'
    # FULL NAME EXAMPLE "ala_motif_name"
    TYPE_TOP_NAME = '_top_motif_name'
    TYPE_TOP_MATRIX = '_top_motif_matrix'
    SCORES = 'scores'
    ALL_TOP_MOTIFS = 'all_top_motifs'
    TOP_SAME_TYPE_DIFF = 'top_same_dif'
    COMPLETE_DATA = 'complete_data'


class Parameters:
    def __init__(self):
        self.line_height = 28
        self.elem_spacing = 20
        self.selector_frame_size = 60
        self.main_window_width = 1100
        self.main_window_height = 800
        self.left_plot_min_h = 100

        self.left_plots_H = 300# - self.elem_spacing *10
        # Left plots parameters
        self.l_same_type_color = 'g'
        self.l_top_type_color = 'r'
        self.l_same_type_lowcorr_color = 'b'
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
        self.r_lbl_size = 7
        self.r_ax_lbl_size = 4
        self.r_ax_txt_size = 5
        self.r_ax_title_lbl_size = 8
        self.per_chain_cc_columns = 5
        self.r_tick_params = {'labelsize': 5, 'length': 1.5, 'pad': 1, 'width': 0.4}
        self.y_axis_text = ['Y axis:', 'Strudel score']
        self.x_axis_text = ['X axis:', 'Residue number']
        self.cc_text_size = 11
        self.grey_text = 'color: rgb(106, 109, 112);'
        self.left_plot_axis_font = 'font-family: Helvetica; font-size: 14px; font-weight: bold;'
        # self.titles_style = 'color: rgb(0, 102, 204);'
        self.titles_style = 'color: rgb(34, 93, 153);'
        self.horisontal_line_color = 'color: rgb(139,139,139);'


if __name__ == '__main__':
    import sys
    app = QtWidgets.QApplication(sys.argv)
    session = ''
    tool_name = 'Strudel Score'
    window = StrudelScore(session, tool_name)
    app.exec_()
