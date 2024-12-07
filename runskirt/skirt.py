
from .skirt_element import Element_ski





class Skirt:
    def __init__(self,config):
        self._config = config
        self._skirt = Element_ski()
        
        self.set_element()
        #self.ski_structure()
        
        simulation_mode = self._config.get('SimIO','ski_simulationMode'.lower(),fallback = 'ExtinctionOnly')

    
    
    
    def write_ski(self,filename = None, xml = None, filepath = None):
        if filepath:
            output_path = filepath
        else:
            output_path = self._config.get('IO','OUTPUT_PATH')
        if filename:
            pass
        else:
            filename = 'sim2obs.ski'
        print(f'Creating {output_path}/{filename} file:')
        self._skirt.write_ski(filename=f'{output_path}/{filename}',xmlelemenmt=xml)
    
    
    def ski_structure(self):
        iterdict = {i: iter(self._skirt._default[i]) for i in self._skirt._default}
        print('\nThe Structure of the .ski file:')
        print(f'MonteCarloSimulation')
        print(f"-random")
        print(f"-units: {next(iterdict['Units'])}")
        print(f"-cosmology: {next(iterdict['Cosmology'])}")
        print(f"-sourceSystem: {self._skirt._default['sources']}")
        for i in iterdict['sources']:
            print(f"---{i}")
            for j in next(iterdict['sources_options']).split(','):
                print(f"-----{j}: {next(iterdict[j])}")
        
        
        print(f"-mediumSystem: {self._skirt._default['mediums']}")
        for i in iterdict['mediums']:
            print(f"---{i}")
            for j in next(iterdict['mediums_options']).split(','):
                print(f"-----{j}: {next(iterdict[j])}")
        print(f"-instrumentSystem: {self._skirt._default['Instrument']}")
        print(f"-probeSystem")
        return

    def set_element(self,config = None):
        
        if config:
            thisconfig=config
        else:
            thisconfig = self._config
        
        self._skirt._default['Units'] = [thisconfig.get('units','type',fallback = 'ExtragalacticUnits')]
        self._skirt._default['Cosmology'] = [thisconfig.get('cosmology','type',fallback = 'FlatUniverseCosmology')]
        
        self._set_sourceSystem(config = thisconfig)
        self._set_mediumSystem(config = thisconfig)
        self._set_InstrumentSystem(config = thisconfig)
        
        
    def _set_sourceSystem(self,config=None):
        if config:
            thisconfig=config
        else:
            thisconfig = self._config
        self._skirt._default['sources'] = []
        self._skirt._default['sources_options'] = []
        self._skirt._default['smoothingKernel'] = []
        self._skirt._default['wavelengthBiasDistribution'] = []
        self._skirt._default['sedFamily'] = []
        for i in thisconfig.sections():
            if (i[:6] == 'Source'):
                self._skirt._default['sources'].append(thisconfig.get(i,'type'))
                self._skirt._default['sources_options'].append(thisconfig.get(i,'options'))
                for j in self._skirt._default['sources_options'][-1].split(','):
                    if j not in self._skirt._default:
                        self._skirt._default[j] = []
                    self._skirt._default[j].append(thisconfig.get(i,j))
                

    def _set_mediumSystem(self,config=None):
        if config:
            thisconfig=config
        else:
            thisconfig = self._config
        self._skirt._default['mediumSystem_options'] = [thisconfig.get('mediumSystem','options')]
        if 'radiationFieldOptions' in self._skirt._default['mediumSystem_options'][-1]:
             self._skirt._default['radiationFieldWLG'] = [thisconfig.get('mediumSystem','radiationFieldWLG')]
        if 'grid' in self._skirt._default['mediumSystem_options'][-1]:
            self._skirt._default['grid'] = [thisconfig.get('mediumSystem','grid')]
            self._skirt._default['policy'] = [thisconfig.get('mediumSystem','policy')]
            
        self._skirt._default['mediums'] = []
        self._skirt._default['mediums_options'] = []
        if 'smoothingKernel' not in self._skirt._default:
            self._skirt._default['smoothingKernel'] = []
        self._skirt._default['materialMix'] = []
        for i in thisconfig.sections():
            if (i[:6] == 'Medium'):
                self._skirt._default['mediums'].append(thisconfig.get(i,'type'))
                self._skirt._default['mediums_options'].append(thisconfig.get(i,'options'))
                for j in self._skirt._default['mediums_options'][-1].split(','):
                    if j not in self._skirt._default:
                        self._skirt._default[j] = []
                    self._skirt._default[j].append(thisconfig.get(i,j))
                
            
            
            
    def _set_InstrumentSystem(self,config=None):
        if config:
            thisconfig=config
        else:
            thisconfig = self._config
        self._skirt._default['defaultWavelengthGrid'] = [thisconfig.get('InstrumentSystem','defaultWavelengthGrid')]
        
        self._skirt._default['Instrument'] = []
        self._skirt._default['wavelengthGrid'] = []
        for i in thisconfig['InstrumentSystem']:
            if ('instruments' in i) and ('_' not in i):
                self._skirt._default['Instrument'].append(thisconfig.get('InstrumentSystem',i))
                if i+'_wavelengthGrid' in thisconfig['InstrumentSystem']:
                    self._skirt._default['wavelengthGrid'].append(thisconfig.get('InstrumentSystem',i+'_wavelengthGrid'))
                
    def set_parameter(self, file=None):
        import xml.etree.ElementTree as ET
        if file:
            h= ET.parse(file)
        else:
            output_path = self._config.get('IO','OUTPUT_PATH')
            filename = 'sim2obs.ski'
            file = f'{output_path}/{filename}'
            h= ET.parse(file)
        hroot = h.getroot()
        
        recordSource = {}
        recordMedium = {}
        recordInstrument = {}
        for i in self._config:
            if i == 'MonteCarloSimulation':
                self._pa_MonteCarloSimulation(hroot,i)
            if i == 'random':
                self._pa_Random(hroot,i)
            if i == 'units':
                self._pa_Units(hroot,i)
            if i == 'cosmology':
                self._pa_Cosmology(hroot,i)
                
            if i == 'sourceSystem':
                self._pa_SourceSystem(hroot,i)
            if i[:6] == 'Source':
                print(f'\nSet {i}:')
                SourceType = self._config[i]['type']
                if SourceType in recordSource:
                    recordSource[SourceType] = recordSource[SourceType] + 1
                else:
                    recordSource[SourceType] = 0
                    
                if len(hroot.findall(f'.//{SourceType}')) > recordSource[SourceType]:
                    for j in self._config[i]:
                        if j[:5] == 'type_':
                            pa = j.split('_')[1]
                            if len(hroot.findall(f'.//{SourceType}[@{pa}]')) > recordSource[SourceType]:
                                hroot.findall(f'.//{SourceType}[@{pa}]')[recordSource[SourceType]].attrib[pa] = self._config[i][j]
                                print(f'---{pa} = {self._config[i][j]}')
                            continue
                        if ('_' in j) and (j.split('_')[0] in self._config[i]['options']):
                            ta = j.split('_')[0]
                            pa = j.split('_')[1]
                            if hroot.findall(f'.//{SourceType}')[recordSource[SourceType]].findall(f'.//{ta}//*[@{pa}]'):
                                hroot.findall(f'.//{SourceType}')[recordSource[SourceType]].findall(f'.//{ta}//*[@{pa}]')[0].attrib[pa] = self._config[i][j]
                                print(f'---{ta}: {pa} = {self._config[i][j]}')
                            continue
                else:
                    print(f'! Some wrong, do not find enough {SourceType} in .ski file')
                    
            if i == 'mediumSystem':
                self._pa_MediumSystem(hroot,i)
                    
            if i[:6] == 'Medium':
                print(f'\nSet {i}:')
                MediumType = self._config[i]['type']
                if MediumType in recordMedium:
                    recordMedium[MediumType] = recordMedium[MediumType] + 1
                else:
                    recordMedium[MediumType] = 0
                    
                if len(hroot.findall(f'.//{MediumType}')) > recordMedium[MediumType]:
                    for j in self._config[i]:
                        if j[:5] == 'type_':
                            pa = j.split('_')[1]
                            if len(hroot.findall(f'.//{MediumType}[@{pa}]')) > recordMedium[MediumType]:
                                hroot.findall(f'.//{MediumType}[@{pa}]')[recordMedium[MediumType]].attrib[pa] = self._config[i][j]
                                print(f'---{pa} = {self._config[i][j]}')
                            continue
                        if ('_' in j) and (j.split('_')[0] in self._config[i]['options']):
                            ta = j.split('_')[0]
                            pa = j.split('_')[1]
                            if hroot.findall(f'.//{MediumType}')[recordMedium[MediumType]].findall(f'.//{ta}//*[@{pa}]'):
                                hroot.findall(f'.//{MediumType}')[recordMedium[MediumType]].findall(f'.//{ta}//*[@{pa}]')[0].attrib[pa] = self._config[i][j]
                                print(f'---{ta}: {pa} = {self._config[i][j]}')
                            continue
                else:
                    print(f'! Some wrong, do not find enough {MediumType} in .ski file')
                    
            if i == 'InstrumentSystem':
                print(f'\nSet {i}:')
                for j in self._config[i]:
                    if ('_' in j) and (j[:11]=='instruments'):
                        instrunum = j.split('_')[0]
                        instruType = self._config[i][instrunum]
                        if instruType in recordInstrument:
                            if instrunum not in recordInstrument[instruType]:
                                recordInstrument[instruType].append(instrunum)
                        else:
                            recordInstrument[instruType] = []
                            recordInstrument[instruType].append(instrunum)
                        ind = recordInstrument[instruType].index(instrunum)
                        thisinstru = hroot.findall(f'.//InstrumentSystem//{instruType}')[ind]
                        
                        if len(j.split('_'))==2:
                            pa = j.split('_')[1]
                            if thisinstru.findall(f'[@{pa}]'):
                                thisinstru.findall(f'[@{pa}]')[0].attrib[pa] = self._config[i][j]
                                print(f'---{instrunum}: {pa} = {self._config[i][j]}')
                            continue
                        if len(j.split('_'))==3:
                            ta = j.split('_')[1]
                            pa = j.split('_')[2]
                            if thisinstru.findall(f'.//{ta}//*[@{pa}]'):
                                thisinstru.findall(f'.//{ta}//*[@{pa}]')[0].attrib[pa] = self._config[i][j]
                                print(f'---{instrunum}: {ta}: {pa} = {self._config[i][j]}')
                            continue
                        continue
                    if ('_' in j):
                        ta = j.split('_')[0]
                        pa = j.split('_')[1]
                        if hroot.findall(f'.//InstrumentSystem//{ta}//*[@{pa}]'):
                            hroot.findall(f'.//InstrumentSystem//{ta}//*[@{pa}]')[0].attrib[pa] = self._config[i][j]
                            print(f'---{ta}: {pa} = {self._config[i][j]}')
                        continue
                
        self._skirt.write_ski(filename=file,xmlelemenmt=h) 
        return
    
    
    def _pa_MonteCarloSimulation(self,hroot,i):
        print('\nSet MonteCarloSimulation')
        if len(hroot.findall(f'.//{i}'))==1:
            for j in self._config[i]:
                if hroot.findall(f'.//MonteCarloSimulation[@{j}]'):
                    hroot.findall(f'.//MonteCarloSimulation[@{j}]')[0].attrib[j] = self._config[i][j]
                    print(f'---{j} = {self._config[i][j]}')
        elif len(hroot.findall(f'.//{i}')) > 1:
            print('! Some wrong, find more than one MonteCarloSimulation')
        else:
            print('! Some wrong, do not find MonteCarloSimulation in .ski file')
        return
    
    def _pa_Random(self,hroot,i):
        print('\nSet Random:')
        if len(hroot.findall(f'.//random'))==1:
            for j in self._config[i]:
                if hroot.findall(f'.//random//*[@{j}]'):
                    hroot.findall(f'.//random//*[@{j}]')[0].attrib[j] = self._config[i][j]
                    print(f'---{j} = {self._config[i][j]}')
        elif len(hroot.findall(f'.//{i}')) > 1:
            print('! Some wrong, find more than one random')
        else:
            print('! Some wrong, do not find random in .ski file')
        return
    
    def _pa_Units(self,hroot,i):
        print('\nSet Units:')
        if len(hroot.findall(f'.//units'))==1:
            for j in self._config[i]:
                if hroot.findall(f'.//units//*[@{j}]'):
                    hroot.findall(f'.//units//*[@{j}]')[0].attrib[j] = self._config[i][j]
                    print(f'---{j} = {self._config[i][j]}')
        elif len(hroot.findall(f'.//{i}')) > 1:
            print('! Some wrong, find more than one units')
        else:
            print('! Some wrong, do not find units in .ski file')
        return
    
    def _pa_Cosmology(self,hroot,i):
        print('\nSet Cosmology:')
        if len(hroot.findall(f'.//cosmology'))==1:
            for j in self._config[i]:
                if hroot.findall(f'.//cosmology//*[@{j}]'):
                    hroot.findall(f'.//cosmology//*[@{j}]')[0].attrib[j] = self._config[i][j]
                    print(f'---{j} = {self._config[i][j]}')
        elif len(hroot.findall(f'.//{i}')) > 1:
            print('! Some wrong, find more than one cosmology')
        else:
            print('! Some wrong, do not find cosmology in .ski file')
        return
    
    def _pa_SourceSystem(self,hroot,i):
        print('\nSet SourceSystem:')
        if len(hroot.findall(f'.//sourceSystem//SourceSystem'))==1:
            for j in self._config[i]:
                if hroot.findall(f'.//sourceSystem//SourceSystem[@{j}]'):
                    hroot.findall(f'.//sourceSystem//SourceSystem[@{j}]')[0].attrib[j] = self._config[i][j]
                    print(f'---{j} = {self._config[i][j]}')
        elif len(hroot.findall(f'.//{i}//SourceSystem')) > 1:
            print('! Some wrong, find more than one cosmology')
        else:
            print('! Some wrong, do not find cosmology in .ski file')
        return
    def _pa_MediumSystem(self,hroot,i):
        print('\nSet MediumSystem:')
        if len(hroot.findall(f'.//mediumSystem//MediumSystem'))==1:
            for j in self._config[i]:
                if hroot.findall(f'.//mediumSystem//MediumSystem//*[@{j}]'):
                    hroot.findall(f'.//mediumSystem//MediumSystem//*[@{j}]')[0].attrib[j] = self._config[i][j]
                    print(f'---{j} = {self._config[i][j]}')
                    continue
                if ('_' in j) and (j.split('_')[0] in self._config[i]['options']):
                    ta = j.split('_')[0]
                    pa = j.split('_')[1]
                    if hroot.findall(f'.//mediumSystem//MediumSystem//{ta}//*[@{pa}]'):
                        hroot.findall(f'.//mediumSystem//MediumSystem//{ta}//*[@{pa}]')[0].attrib[pa] = self._config[i][j]
                        print(f'---{ta}: {pa} = {self._config[i][j]}')
                    continue
                if ('_' in j):
                    ta = j.split('_')[0]
                    pa = j.split('_')[1]
                    if hroot.findall(f'.//mediumSystem//MediumSystem//{ta}//*[@{pa}]'):
                        hroot.findall(f'.//mediumSystem//MediumSystem//{ta}//*[@{pa}]')[0].attrib[pa] = self._config[i][j]
                        print(f'---{ta}: {pa} = {self._config[i][j]}')
                    continue
        elif len(hroot.findall(f'.//{i}//MediumSystem')) > 1:
            print('! Some wrong, find more than one MediumSystem')
        else:
            print('! Some wrong, do not find MediumSystem in .ski file')
        return
    #TODO set_probeSystem
    
    
    
    
    
        
                
        