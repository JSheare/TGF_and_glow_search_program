"""Child class of Detector customized to handle data from the Santis instrument."""
from tgf_search.detectors.detector import Detector
from tgf_search.detectors.scintillator import Scintillator


class Santis(Detector):
    def __init__(self, unit, first_sec, modes, print_feedback=False):
        super().__init__(unit, first_sec, modes, print_feedback)

        self.long_event_scint_list = ['LP']
        self.calibration_params = {'bin_range': 15008.0, 'bin_size': 16, 'band_starts': [38, 94],
                                   'band_ends': [75, 125], 'template_bin_plot_edge': 400}
        self.default_data_loc = '/media/AllDetectorData/Detectors/SANTIS/Data'
        self.location = self.get_location(self.default_data_loc[:-5])
        self.import_path = f'{self.default_data_loc}/{self.date_str}'
        self.regex = self.regex = lambda eRC: f'eRC{eRC}*_lm_{self.date_str}_*'
        self.scintillators = {'LP': Scintillator('LP', '2549')}
        self.scint_list = list(self.scintillators.keys())

        self.check_processed()
        self.check_gui()
