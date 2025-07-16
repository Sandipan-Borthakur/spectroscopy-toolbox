from PyQt5 import uic
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget,QVBoxLayout
import sys
from controller.readoutnoise_controller import ReadoutNoiseController
from controller.gain_controller import GainController
from controller.resolving_power_controller import ResolvingPowerController
from controller.wavelength_solution_controller import WavelengthSolutionController

def test():
    print("Yes")

class MyApp(QMainWindow):
    def __init__(self):
        super().__init__()
        uic.loadUi("Spectroscopy-Toolbox.ui", self)

        self.layout = self.centralwidget.layout()
        if self.layout is None:
            self.layout = QVBoxLayout(self.centralwidget)
            self.layout.addWidget(self.textBrowser)

        options = [self.actionResolving_Power_Calculator,
        self.actionReadout_Noise_Calculator,
        self.actionGain_Calculator,
        self.actionSwath_Width_Calculator,
        self.actionWavelength_Solution_Generator]

        optionclasses = [
            ResolvingPowerController,
            ReadoutNoiseController,
            GainController,
            ReadoutNoiseController,
            WavelengthSolutionController
        ]
        for option,optionclass in zip(options,optionclasses):
            option.triggered.connect(lambda _,oc=optionclass: self.open_options(oc()))

        self.actionSave.triggered.connect(test)
        self.actionExit.triggered.connect(self.close)
        self.actionAbout.triggered.connect(self.show_initial_view)
        self.actionTools_help.triggered.connect(test)

    def open_options(self,OptionClass):
        self.option_controller = OptionClass
        while self.layout.count():
            child = self.layout.takeAt(0)
            if child.widget():
                child.widget().setParent(None)
        self.layout.addWidget(self.option_controller.view)

    def show_initial_view(self):
        self.clear_layout()
        self.layout.addWidget(self.textBrowser)

    def clear_layout(self):
        while self.layout.count():
            child = self.layout.takeAt(0)
            if child.widget():
                child.widget().setParent(None)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MyApp()
    window.show()
    sys.exit(app.exec_())