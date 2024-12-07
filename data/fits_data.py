

from astropy.io import fits


class FitsImage:
    def __init__(self,filepath):
        self.file = filepath
        data, header = fits.getdata(self.file, header=True)
        self.data =data
        self.header = header