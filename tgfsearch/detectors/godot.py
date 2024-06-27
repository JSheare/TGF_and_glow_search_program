"""Child class of Detector customized to handle data from GODOT."""
from tgfsearch.detectors.detector import Detector
from tgfsearch.detectors.scintillator import Scintillator


class Godot(Detector):
    def __init__(self, unit, date_str, print_feedback=False):
        super().__init__(unit, date_str, print_feedback)

        self.long_event_scint_list = ['NaI', 'LP']
        self.calibration_params = {'bin_range': 15008.0, 'bin_size': 16, 'template_bin_plot_edge': 200}
        self.default_data_loc = '/media/tgfdata/Detectors/GODOT/Data'
        self.location = self.get_location(self.default_data_loc[:-5])
        self._import_loc = f'{self.default_data_loc}/{self.date_str}'
        self._scintillators = {'NaI': Scintillator('NaI', '1490'), 'LP': Scintillator('LP', '1491')}
        self.scint_list = list(self._scintillators.keys())

    def file_form(self, eRC):
        return f'eRC{eRC}_*_{self.date_str}_*'
