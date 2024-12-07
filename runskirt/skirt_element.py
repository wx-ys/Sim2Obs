import xml.dom.minidom as xdm



def ski_append(main,ele):
    #tn = xdm.Document().createTextNode('\n         ')
    if main.childNodes:
        main.appendChild(ele)
     #   main.appendChild(tn)
    else:
     #   main.appendChild(tn)
        main.appendChild(ele)
    #    main.appendChild(tn)
    return main




class Element_ski:
    '''
    Generates template for each element in skirt
    '''
    _element_registry={}

    def __init__(self):
        self._element={}
        self._default={}

    
    
    def __getitem__(self, name):
        """Return the element of a given kind"""
        if name in self._element:
            return self._element[name]
        else:
            return self._get_element(name)
        
    def _get_element(self,name):
        return Element_ski._element_registry[name](self)
    
    @staticmethod
    def ski_element(fn):
        Element_ski._element_registry[fn.__name__] = fn
        return fn
    
    def get_default(self,key):
        if key in self._default:
            if self._default[key]:
                return self._default[key].pop(0)
            else:
                return 
        return 
    

    def write_ski(self,filename,xmlelemenmt = None):
        if xmlelemenmt and isinstance(filename,str):
            if hasattr(xmlelemenmt,'writexml'):
                with open(filename,'w') as f:
                    xmlelemenmt.writexml(f,newl='\n',indent='',addindent='\t',encoding="UTF-8")
                print(f'{xmlelemenmt} is successfully saved in {filename}')
                
            elif hasattr(xmlelemenmt,'write'):
                xmlelemenmt.write(filename,encoding='UTF-8', xml_declaration=True)
                with open(filename,'r+') as f:
                    cont = f.readlines()
                    cont.insert(1,'<!-- A SKIRT parameter file © Astronomical Observatory, Ghent University -->\n')
                    f.seek(0)
                    f.writelines(cont)
                print(f'{xmlelemenmt} is successfully saved in {filename}')
            else:
                raise ValueError(f'{xmlelemenmt} is not writable')
        elif isinstance(filename,str):
            hh = self['skirt']
            with open(filename,'w') as f:
                hh.writexml(f,newl='\n',indent='',addindent='\t',encoding="UTF-8")
            print(f'ski file init: {hh} is successfully saved in {filename}')
        else:
            raise ValueError(f'{filename} must be str')
        return
@Element_ski.ski_element
def skirt(ski):
    import datetime
    ele = xdm.Document()
    head = xdm.Document().createComment(' A SKIRT parameter file © Astronomical Observatory, Ghent University ')
    ele.appendChild(head)
    mai = xdm.Document().createElement('skirt-simulation-hierarchy')
    mai.setAttribute('type', 'MonteCarloSimulation')
    mai.setAttribute('format', '9')
    mai.setAttribute('producer', 'Sim2Obs')
    mai.setAttribute('time', datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S.%f')[:-3])
    mai.appendChild(ski['MonteCarloSimulation'])
    ele.appendChild(mai)
    return ele