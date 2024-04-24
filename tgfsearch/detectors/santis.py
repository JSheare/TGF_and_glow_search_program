"""Child class of Detector customized to handle data from the Santis instrument."""
from tgfsearch.detectors.detector import Detector
from tgfsearch.detectors.scintillator import Scintillator


class Santis(Detector):
    def __init__(self, unit, date_str, modes, print_feedback=False):
        super().__init__(unit, date_str, modes, print_feedback)

        self.long_event_scint_list = ['LP']
        self.calibration_params = {'bin_range': 15008.0, 'bin_size': 16, 'template_bin_plot_edge': 400}
        self.default_data_loc = '/media/tgfdata/Detectors/SANTIS/Data'
        self.location = self.get_location(self.default_data_loc[:-5])
        self.import_loc = f'{self.default_data_loc}/{self.date_str}'
        self._scintillators = {'LP': Scintillator('LP', '2549')}
        self.scint_list = list(self._scintillators.keys())

        self.check_processed()
        self.check_custom()

    def file_form(self, eRC):
        return f'eRC{eRC}*_*_{self.date_str}_*'
