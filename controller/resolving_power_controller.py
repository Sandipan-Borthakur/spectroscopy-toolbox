from model.resolving_power_calculator import resolving_power_calculator
from view.resolving_power_view import ResolvingPowerCalculator

class ResolvingPowerController:
    def __init__(self):
        self.view = ResolvingPowerCalculator()
        self.view.calculate_button.clicked.connect(self.calculate)

    def calculate(self):
        try:
            fullvel,fullnormflux,fullnormfluxfit,fullresolving_power = resolving_power_calculator(self.view.lamp_input_field.get_path(),self.view.linelist_input_field.get_path(),int(self.view.ndim_input_field.text()))

            self.view.update_plot(fullvel,fullnormflux,fullnormfluxfit,fullresolving_power)
        except ValueError:
            self.view.result_label.setText("Invalid input")