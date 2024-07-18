"""Child class of Detector customized to handle data from THOR."""
from tgfsearch.detectors.detector import Detector
from tgfsearch.detectors.scintillator import Scintillator


class Thor(Detector):
    def __init__(self, unit, date_str, print_feedback=False):
        super().__init__(unit, date_str, print_feedback)

        self.calibration_params = {'bin_range': 65535.0, 'bin_size': 1, 'template_bin_plot_edge': 8000}
        self.default_data_loc = '/media/tgfdata/Detectors/THOR'
        self.location = self.get_location(self.default_data_loc)
        self._import_loc = f'{self.default_data_loc}/{self.unit}/Data/{self.date_str}'
        serial_nums = self._get_eRC()
        self._scintillators = {'NaI': Scintillator('NaI', serial_nums[0]), 'SP': Scintillator('SP', serial_nums[1]),
                               'MP': Scintillator('MP', serial_nums[2]), 'LP': Scintillator('LP', serial_nums[3])}
        self.scint_list = list(self._scintillators.keys())

    def file_form(self, eRC):
        return f'eRC{eRC}*_*_{self.date_str}_*'

    def _get_eRC(self):
        """Returns a list of all THOR eRC serial numbers for the instantiated THOR unit."""
        # All lists are in this form: NaI, small plastic, medium plastic, large plastic
        unit = int(self.unit[4:])
        if unit == 1:
            return ['4179', '4194', '4189', '4195']
        elif unit == 2:
            return ['4182', '4172', '4167', '4187']
        elif unit == 3:
            return ['4169', '4175', '4174', '4185']
        elif unit == 4:
            return ['4177', '4191', '4192', '4181']
        elif unit == 5:
            # THOR5 MP was replaced on Nov. 9th 2022
            return ['4188', '4190', '4169' if int(self.date_str) >= 221109 else '4178', '4173']
        elif unit == 6:
            return ['4186', '4176', '4183', '4180']

    def is_named(self, name):
        return True if 'THOR' in name.upper() else False
