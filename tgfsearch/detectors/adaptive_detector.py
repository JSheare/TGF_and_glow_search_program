"""Child class of Detector made to handle data for instruments with arbitrary scintillator configurations."""

from tgfsearch.detectors.detector import Detector
from tgfsearch.detectors.scintillator import Scintillator


class AdaptiveDetector(Detector):
    def __init__(self, date_str, print_feedback=False):
        super().__init__('ADAPTIVE', date_str, print_feedback)
        self.default_scintillator = ''

        self.spectra_params = {'bin_range': 65535.0, 'bin_size': 1}
        self._scintillators = {}
        self._has_identity = False

    def __bool__(self):
        if self._has_identity:
            return super().__bool__()

        return False

    def _reset_identity(self):
        """Resets everything back to its default state (no identity)."""
        self.clear()
        self._results_loc = self._results_loc.replace(f'Results/{self.unit}', 'Results/ADAPTIVE')
        self.unit = 'ADAPTIVE'
        self.deployment = self._get_deployment()
        self._scintillators = {}
        self.scint_list = []
        self._has_identity = False

    def _infer_identity(self):
        """Uses files in import loc to determine instrument type and scintillator configuration."""
        all_files = self._get_serial_num_filelist('*')
        if len(all_files) == 0:
            raise FileNotFoundError('no data files to infer instrument configuration from.')

        # Attempting to infer instrument name based on the names of parent directories
        parent_dirs = all_files[0].replace('\\', '/').split('/')[:-1]
        if len(parent_dirs) >= 3:
            self.unit = parent_dirs[-3].upper()
            # Regenerating deployment info and results directory with the new name
            self.deployment = self._get_deployment()
            self._results_loc = self._results_loc.replace('Results/ADAPTIVE', f'Results/{self.unit}')

        # Attempting to infer scintillator configuration based on the data files present
        # Walking each file to determine its corresponding scintillator
        for file in all_files:
            index = 0
            for i in range(len(file) - 1, 0, -1):
                if file[i] == '/' or file[i] == '\\':
                    break

                index = i

            # Making a new Scintillator if one doesn't exist already
            match file[index + 7: index + 10]:
                case 'lpl':
                    scintillator = 'LP'
                case 'ipl':
                    scintillator = 'IP'
                case 'spl':
                    scintillator = 'SP'
                case 'nai':
                    scintillator = 'NaI'
                case 'mpl':
                    scintillator = 'MP'
                case _:
                    raise ValueError('unsupported scintillator type')

            if scintillator not in self._scintillators:
                eRC = file[index + 3: index + 7]
                self._scintillators[scintillator] = Scintillator(scintillator, eRC)

        self.scint_list = list(self._scintillators.keys())

        # Assigning the default scintillator based on the following priority (greatest to least)
        for scint in ['LP', 'IP', 'MP', 'NaI', 'SP']:
            if scint in self._scintillators:
                self.default_scintillator = scint
                break

        self._has_identity = True

    def set_import_loc(self, loc):
        self._reset_identity()
        super().set_import_loc(loc)
        self._infer_identity()

    def is_named(self, name):
        return name.upper() in self.unit

    def get_clone(self):
        clone = type(self)(self.date_str, print_feedback=self.print_feedback)
        clone.unit = self.unit
        clone.deployment = clone._get_deployment()
        clone.default_scintillator = self.default_scintillator
        clone._import_loc = self._import_loc
        clone._results_loc = self._results_loc
        for scintillator in self._scintillators:
            clone._scintillators[scintillator] = Scintillator(
                self._scintillators[scintillator].name,
                self._scintillators[scintillator].eRC)

        clone.scint_list = list(clone._scintillators.keys())
        clone.processed = self.processed
        clone._has_identity = self._has_identity
        return clone
