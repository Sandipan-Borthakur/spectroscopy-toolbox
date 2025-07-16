from PyQt5.QtWidgets import QWidget, QLabel, QLineEdit, QPushButton, QVBoxLayout, QHBoxLayout, QFileDialog

class FileInput(QWidget):
    def __init__(self, label_text):
        super().__init__()
        self.label = QLabel(label_text)
        self.label.setFixedWidth(150)

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

class ReadoutNoiseCalculator(QWidget):
    def __init__(self):
        super().__init__()

        self.bias1_input_field = FileInput("Bias File1:")
        self.bias1_input_field.setFixedWidth(500)
        self.bias2_input_field = FileInput("Bias File2 (Optional):")
        self.bias2_input_field.setFixedWidth(500)
        self.calculate_button = QPushButton("Calculate")
        self.calculate_button.setFixedWidth(100)

        self.result_label = QLabel("Readout noise, R = ")

        top_layout = QHBoxLayout()
        top_layout.addWidget(self.bias1_input_field)
        top_layout.addWidget(self.bias2_input_field)

        main_layout = QVBoxLayout()
        main_layout.addLayout(top_layout)
        main_layout.addWidget(self.calculate_button)
        main_layout.addWidget(self.result_label)

        self.setLayout(main_layout)