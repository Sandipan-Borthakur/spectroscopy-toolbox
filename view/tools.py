# views/image_view.py

import numpy as np
from PyQt5.QtWidgets import QWidget, QVBoxLayout
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.widgets import Slider

class ImageView(QWidget):
    def __init__(self, data):
        super().__init__()

        self.data = data
        self.vmin_init = np.mean(data) - np.std(data)
        self.vmax_init = np.mean(data) + np.std(data)

        # Create a matplotlib figure
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.ax = self.figure.add_subplot(111)

        # Add the image
        self.im = self.ax.imshow(data, cmap='gray', vmin=self.vmin_init, vmax=self.vmax_init,
                                 aspect="auto", origin="lower")
        self.figure.colorbar(self.im, ax=self.ax)

        # Adjust layout for sliders
        self.figure.subplots_adjust(bottom=0.25)

        # Create slider axes within the figure (not PyQt)
        self.ax_vmin = self.figure.add_axes([0.2, 0.1, 0.65, 0.03])
        self.ax_vmax = self.figure.add_axes([0.2, 0.05, 0.65, 0.03])

        self.slider_vmin = Slider(self.ax_vmin, 'Min',
                                  np.mean(data) - 5 * np.std(data),
                                  np.mean(data), valinit=self.vmin_init)
        self.slider_vmax = Slider(self.ax_vmax, 'Max',
                                  np.mean(data),
                                  np.mean(data) + 5 * np.std(data), valinit=self.vmax_init)

        # Connect sliders
        self.slider_vmin.on_changed(self.update_contrast)
        self.slider_vmax.on_changed(self.update_contrast)

        # Setup layout
        layout = QVBoxLayout()
        layout.addWidget(self.canvas)
        self.setLayout(layout)

    def update_contrast(self, val):
        self.im.set_clim(vmin=self.slider_vmin.val, vmax=self.slider_vmax.val)
        self.canvas.draw_idle()
