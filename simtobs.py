
import configparser
import psutil


from .grabsim import sim_data
from .runskirt import skirt


class Sim2Obs:
    def __init__(self,config: str):
        self._config = configparser.RawConfigParser(comment_prefixes = ('#',';','---'),inline_comment_prefixes=('#'),
                                                   interpolation=configparser.ExtendedInterpolation())
        self._config.optionxform = lambda option: option        # it prevents option key from being lowercase,
        self._config.read(config)
        
        print(f'Obtaining parameters from file \'{config}\'\n')
        
        self.systeminfo()
        self.sim = sim_data.Simdata(self._config)
        self.skirt = skirt.Skirt(self._config)
        #self.__write_ski()
        #self.__write_config()
    
    def pre_init(self,**kwargs):
        self.sim.load_init()
        self.skirt.ski_structure()
        return
    
    def pre_skirt(self,**kwargs):
        self.sim.to_file(**kwargs)
        self.__write_ski()
        self.__write_config()
        return
        
    def __write_config(self):
        output_path = self._config.get('IO','OUTPUT_PATH')
        print(f'Saving config file: {output_path}/sim2obs.ini')
        with open(f'{output_path}/sim2obs.ini','w') as configfile:
            self._config.write(configfile)
        return
    def __write_ski(self):
        self.skirt.write_ski()
        self.skirt.set_parameter()
        return
        
    def systeminfo(self):
        print(f"System Info: ")
        print(f"---Available CPU (total: {psutil.cpu_count()}): {(100-psutil.cpu_percent(interval=0.2)):.2f}%")
        output_path = self._config.get('IO','OUTPUT_PATH')
        checkusage = psutil.disk_usage(path=output_path)
        print(f'---Available storage ({(100-checkusage.percent):.2f}%): {sim_data.readableFileSize(checkusage.free)}/{sim_data.readableFileSize(checkusage.total)}')
        checkmemory = psutil.virtual_memory()
        print(f'---Available memory ({(100-checkmemory.percent):.2f}%): {sim_data.readableFileSize(checkmemory.available)}/{sim_data.readableFileSize(checkmemory.total)}\n')
        return