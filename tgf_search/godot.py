"""Child class of Detector customized to handle data from GODOT."""
from detector import Detector
from scintillator import Scintillator


class Godot(Detector):
    def __init__(self, unit, first_sec, modes, print_feedback=False):
        super().__init__(unit, first_sec, modes, print_feedback)

        self.long_event_scint_list = ['NaI', 'LP']
        self.calibration_params = {'bin_range': 15008.0, 'bin_size': 16, 'band_starts': [38, 94],
                                   'band_ends': [75, 125], 'template_bin_plot_edge': 200}
        self.default_data_loc = '/media/AllDetectorData/Detectors/GODOT/Data'
        self.location = self.get_location(self.default_data_loc[:-5])
        self.import_path = f'{self.default_data_loc}/{self.date_str}'
        self.regex = self.regex = lambda eRC: f'eRC{eRC}_lm*_{self.date_str}_*'
        self.scintillators = {'NaI': Scintillator('NaI', '1490'), 'LP': Scintillator('LP', '1491')}
        self.scint_list = list(self.scintillators.keys())

        self.check_processed()
        self.check_gui()
