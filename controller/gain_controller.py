from model.gain_calculator import gain_calculator
from view.gain_view import GainCalculator

class GainController:
    def __init__(self):
        self.view = GainCalculator()
        self.view.calculate_button.clicked.connect(self.calculate)

    def calculate(self):
        try:
            gain = gain_calculator([self.view.flat1_input_field.get_path(),self.view.flat2_input_field.get_path()],
                                                                                          [self.view.bias1_input_field.get_path(),self.view.bias2_input_field.get_path()])

            self.view.update_plot(gain)
        except ValueError:
            self.view.result_label.setText("Invalid input")