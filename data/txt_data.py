
import re

import numpy as np
from astropy import units as u

class TxtTable:
    def __init__(self,filepath: str):
        self.file = filepath
        self.column={}
        self._read_data()
        self._read_header()
        self._column_info()
        
    def _read_header(self, startwith='#'):
        header = []
        with open(self.file,'r') as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith(startwith):
                    header.append(line)
                else:
                    break
        self.header = header
    def _column_info(self):
        if 'header' in dir(self):
            pass
        else:
            self.read_header()
        for i in self.header:
            find = re.search(r'column \d+: (.*)', i,re.IGNORECASE)
            if find:
                col_n = int(re.search(r'column \d+',find.group(),re.IGNORECASE).group().split()[-1])
                col_name = re.search(r': (.*) ',find.group()).group()[2:-1]
                col_units = re.search(r'\((.*?)\)',find.group()).group()[1:-1]
                self.column[col_name] = self.data[:,col_n-1]*u.Unit(col_units)
    def _read_data(self):
        self.data = np.loadtxt(self.file)
    