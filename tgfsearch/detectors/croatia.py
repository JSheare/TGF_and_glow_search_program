"""Child class of Detector customized to handle data from the Croatia instrument."""
from tgfsearch.detectors.detector import Detector
from tgfsearch.detectors.scintillator import Scintillator


class Croatia(Detector):
    def __init__(self, unit, date_str, print_feedback=False):
        super().__init__(unit, date_str, print_feedback)

        self.spectra_params = {'bin_range': 15008.0, 'bin_size': 16}
        self._import_loc = f'{self.default_data_loc}/SANTIS/Data/{self.date_str}'
        self._scintillators = {'MP': Scintillator('MP', '4193'), 'LP': Scintillator('LP', '2549')}
        self.scint_list = list(self._scintillators.keys())
