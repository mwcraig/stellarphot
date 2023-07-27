# Class to representa file in AAVSO extended format
from dataclasses import dataclass, field
from pathlib import Path
from importlib import resources
import yaml

from astropy.table import Table, Column

__all__ = [
    "AAVSOExtendedFileFormat"
]


@dataclass
class AAVSOExtendedFileFormat:
    obscode : str
    _type: str = 'EXTENDED'
    _software: str = 'stellarphot'
    delim: str = ","
    date_format: str = 'JD'
    date: Column = field(default_factory=Column)
    obstype: str = 'CCD'
    magnitude: Column = field(default_factory=Column)
    magerr: Column = field(default_factory=Column)
    filter: Column = field(default_factory=Column)
    kmag: Column = field(default_factory=Column)
    cmag: Column = field(default_factory=Column)
    ensemble: bool = False
    group : int  = ""
    trans : str = "NO"
    chart : str = ""
    notes : str = ""
    starid : str = ""
    cname : str = ""
    kname : str = ""
    airmass : Column = field(default_factory=Column)
    mtype : str = "STD"

    @property
    def type(self):
        return self._type

    @property
    def software(self):
        return self._software

    def set_data_column(self, data, column_map, star_id=None):
        if star_id is not None:
            use_data = data[data['id'] == star_id]
        else:
            use_data = data
        for source_column, destination_column in column_map.items():
            try:
                setattr(self, destination_column, use_data[source_column])
            except AttributeError:
                raise AttributeError(f"Column {destination_column} not allowed in AAVSO extended format")
            except KeyError:
                raise KeyError(f"Column {source_column} not found in data")

    def set_column_from_single_value(self, value, column):
        setattr(self, column, Column([value] * len(self.variable_mag), name=column))

    def set_data(self, destination, data,
                 time='time', mag='mag', error='err', airmass='airmass'):
        """
        Set data for each thingy
        """
        # Standardize column order

        standard_data = {
            time: data[time],
            airmass: data[airmass],
            mag: data[mag],
            error: data[error],
        }
        standard_data = Table(standard_data)

        if destination == "variable":
            self.variable_data = standard_data
        elif destination == "check":
            self.check_data = standard_data
            del self.check_data[airmass]
        elif destination == "comparison":
            if self.ensemble:
                raise ValueError("Cannot set comparison data for ensemble magnitudeas")
            self.comparison_data = standard_data
            del self.comparison_data[airmass]

    def write(self, file):
        p = Path(file)

        # Make a table
        table_description = resources.read_text('stellarphot.io', 'aavso_submission_schema.yml')
        table_structure = yaml.safe_load(table_description)
        for key in table_structure['comments'].keys():
            print(f"{key}: {hasattr(self, key.lower())}")
        table_dict = {}
        length = 0
        for key in table_structure['data'].keys():
            print(f"{key}: {hasattr(self, key.lower())}")
            if isinstance((item := getattr(self, key.lower())), Column):
                if len(item) == 0:
                    item = "na"
                table_dict[key] = item
                length = len(item) if len(item) > length else length
            else:
                table_dict[key] = item
        print(length)
        for k, v in table_dict.items():
            if len(v) != length:
                print(f"{k=} {v=} {len(v)=}")

                table_dict[k] = Column([v] * length, name=k)

        table = Table(table_dict)
        table.meta['comments'] = []
        for key in table_structure['comments'].keys():
            if key == "DATE":
                value = self.date_format
            else:
                value = getattr(self, key.lower())
            table.meta['comments'].append(f"{key} = {value}")

        print(table)
        return table