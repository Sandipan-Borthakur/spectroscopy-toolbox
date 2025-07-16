from PyQt5.QtWidgets import (
    QWidget, QLabel, QLineEdit, QPushButton, QFileDialog,
    QHBoxLayout, QVBoxLayout, QSizePolicy
)
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import numpy as np

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

class ResolvingPowerCalculator(QWidget):
    def __init__(self):
        super().__init__()

        self.lamp_input_field = FileInput("Lamp File:")
        self.linelist_input_field = FileInput("Linelist File:")

        self.ndim_label = QLabel("Dimensions:")
        self.ndim_input_field = QLineEdit()
        self.ndim_input_field.setFixedWidth(50)

        self.calculate_button = QPushButton("Calculate")
        self.calculate_button.setFixedWidth(100)

        top_bar = QHBoxLayout()
        top_bar.setContentsMargins(5, 5, 5, 5)
        top_bar.setSpacing(10)
        top_bar.addWidget(self.lamp_input_field)
        top_bar.addWidget(self.linelist_input_field)
        top_bar.addWidget(self.ndim_label)
        top_bar.addWidget(self.ndim_input_field)
        top_bar.addWidget(self.calculate_button)

        self.canvas = PlotCanvas(self, width=6, height=5, dpi=100)

        main_layout = QVBoxLayout()
        main_layout.addLayout(top_bar)
        main_layout.addWidget(self.canvas)

        self.setLayout(main_layout)
        self.setWindowTitle("Resolving Power Calculator")
        self.resize(1000, 700)
        self.show()

    def update_plot(self, fullvel, fullnormflux, fullnormfluxfit, fullresolving_power):
        nrows, ncols, _ = np.shape(fullvel)
        self.canvas.figure.clf()
        axes = self.canvas.figure.subplots(nrows, ncols,sharex=True,sharey=True)

        if nrows == 1:
            axes = np.array([axes])
        if ncols == 1:
            axes = np.array([[ax] for ax in axes])

        for i in range(nrows):
            for j in range(ncols):
                ax = axes[i][j]
                ax.plot(fullvel[i, j], fullnormflux[i, j], 'k.', alpha=0.5)
                ax.plot(fullvel[i, j], fullnormfluxfit[i, j], 'r', label=f"R={fullresolving_power[i, j][0]:.0f}")
                ax.legend(fontsize=6)
                ax.tick_params(labelsize=6)

        self.canvas.figure.tight_layout()
        self.canvas.figure.subplots_adjust(hspace=0, wspace=0)
        self.canvas.draw()
