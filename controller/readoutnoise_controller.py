from model.readoutnoise_calculator import readoutnoise_calculator
from view.readoutnoise_view import ReadoutNoiseCalculator

class ReadoutNoiseController:
    def __init__(self):
        self.view = ReadoutNoiseCalculator()
        self.view.calculate_button.clicked.connect(self.calculate)

    def calculate(self):
        input_val = [self.view.bias1_input_field.get_path(),self.view.bias2_input_field.get_path()]
        try:
            result = readoutnoise_calculator(input_val)
            self.view.result_label.setText(f"Result: {result}")
        except ValueError:
            self.view.result_label.setText("Invalid input")