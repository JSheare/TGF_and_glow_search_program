"""Child class of Detector customized to handle data from GODOT."""
from tgfsearch.detectors.detector import Detector
from tgfsearch.detectors.scintillator import Scintillator


class Godot(Detector):
    def __init__(self, unit, date_str, print_feedback=False):
        super().__init__(unit, date_str, print_feedback)

        self.spectra_params = {'bin_range': 15008.0, 'bin_size': 16}
        self._import_loc = f'{self.default_data_loc}/GODOT/Data/{self.date_str}'
        self._scintillators = {'NaI': Scintillator('NaI', '1490'), 'LP': Scintillator('LP', '1491')}
        self.scint_list = list(self._scintillators.keys())

    def file_form(self, eRC):
        return f'eRC{eRC}_*_{self.date_str}_*'
