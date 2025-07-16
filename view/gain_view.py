from PyQt5.QtWidgets import (
    QWidget, QLabel, QLineEdit, QPushButton, QFileDialog,
    QHBoxLayout, QVBoxLayout, QSizePolicy
)
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import numpy as np
from astropy.io import fits

class FileInput(QWidget):
    def __init__(self, label_text):
        super().__init__()
        self.label = QLabel(label_text)
        self.label.setFixedWidth(110)

        self.line_edit = QLineEdit()
        self.line_edit.setReadOnly(True)
        self.line_edit.setMinimumWidth(200)

        self.browse_button = QPushButton("Browse")
        self.browse_button.clicked.connect(self.browse)

        layout = QHBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(5)
        layout.addWidget(self.label)
        layout.addWidget(self.line_edit)
        layout.addWidget(self.browse_button)
        self.setLayout(layout)

    def browse(self):
        path, _ = QFileDialog.getOpenFileName(self)
        if path:
            self.line_edit.setText(path)

    def get_path(self):
        return self.line_edit.text()

class PlotCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super().__init__(fig)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

class GainCalculator(QWidget):
    def __init__(self):
        super().__init__()

        self.flat1_input_field = FileInput("Flat File1:")
        self.flat2_input_field = FileInput("Flat File2:")
        self.bias1_input_field = FileInput("Bias File1:")
        self.bias2_input_field = FileInput("Bias File2:")

        self.calculate_button = QPushButton("Calculate")
        self.calculate_button.setFixedWidth(100)

        top_bar = QHBoxLayout()
        top_bar.setContentsMargins(5, 5, 5, 5)
        top_bar.setSpacing(10)
        top_bar.addWidget(self.flat1_input_field)
        top_bar.addWidget(self.flat2_input_field)
        top_bar.addWidget(self.bias1_input_field)
        top_bar.addWidget(self.bias2_input_field)
        top_bar.addWidget(self.calculate_button)

        self.canvas = PlotCanvas(self, width=6, height=5, dpi=100)

        self.main_layout = QVBoxLayout()
        self.main_layout.addLayout(top_bar)
        self.main_layout.addWidget(self.canvas)

        self.setLayout(self.main_layout)
        self.setWindowTitle("Gain Calculator")
        self.resize(1000, 700)
        self.show()

    def update_plot(self,gain):
        self.canvas.figure.clf()
        toolbar = NavigationToolbar(self.canvas, self)
        self.main_layout.addWidget(toolbar)
        axes = self.canvas.figure.subplots()
        axes.imshow(gain,origin="lower",cmap="gray")
        self.canvas.figure.tight_layout()
        self.canvas.figure.subplots_adjust(hspace=0, wspace=0)
        self.canvas.draw()