
import os
import time

import numpy as np 
import pynbody

from . import read_sim_io

from .util import make_header, make_fmt, readableFileSize

class Simdata:
    
    def __init__(self, config : str):
        
        self.config = config

        file_operator = self.config.get('SimIO','FILE_OPERATOR')
        print(f'Using simulation operator: {file_operator} \n')
        try:
            operator_func = eval(f'read_sim_io.{file_operator}')
            pass
        except:
            raise KeyError(f'No operator {file_operator} in read_sim_io')


        self.file_operator = operator_func(self.config)
        
        self.file_list = {}
        for i in self.config.sections():
            if (i[:6] == 'Source') or (i[:6] == 'Medium'):
                self.file_list[i] = self.config.get(i,'columns').split(',')
        #self.to_file()

    
    def to_file(self,**kwargs):
        
        for i in kwargs:
            if i in self.config['IO']:
                self.config.set('IO',i,str(kwargs[i]))
            if i in self.config['SimIO']:
                self.config.set('SimIO',i,str(kwargs[i]))
                
        output_path = self.config.get('IO','OUTPUT_PATH')
        print(f'Output path: {output_path}')
        folder = os.path.exists(output_path)
        if not folder:
            os.makedirs(output_path)
            print(f'Create folder: {output_path}')
        print(f'Files containing the data extracted from simulation:')
        totalsize = []
        totaltime = []
        for i in self.file_list:
            time1 = time.time()
            print(f'- Saving {i}.txt, columns: {self.file_list[i]}')
            file_header = make_header(i,self.file_list[i])
            file_fmt = make_fmt(i,self.file_list[i])
            if i[:6] == 'Source':
                info = np.column_stack([self.file_operator.get_star(j, particle=i) for j in self.file_list[i]])
            elif i[:6] == 'Medium':
                info = np.column_stack([self.file_operator.get_gas(j, particle=i) for j in self.file_list[i]])
            else:
                print(f'No match {self.file_list[i]}: saving none')
                info = []            
            np.savetxt(f'{output_path}/{i}.txt', info, header=file_header,fmt = file_fmt)
            totalsize.append(os.path.getsize(f'{output_path}/{i}.txt'))
            time2 = time.time()
            totaltime.append(time2-time1)
            print(f'------ size: {readableFileSize(totalsize[-1])}, time: {(time2-time1):.1f} s')
        print(f'-- Total size: {readableFileSize(np.sum(totalsize))}, total time: {np.sum(totaltime):.1f} s \n')
    
    def load_init(self):
        self.file_operator.load_init()
        return
    
    

