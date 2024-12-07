import xml.dom.minidom as xdm

from .skirt_element import Element_ski, ski_append


@Element_ski.ski_element
def MonteCarloSimulation(ski):
    '''
    userLevel: the user experience level
        Basic,      for beginning users (hides many options)
        Regular,    for regular users (hides esoteric options)
        Expert,     for expert users (hides no options)
        
    simulationMode: the overall simulation mode
        OligoNoMedium,          No medium - oligochromatic regime (a few discrete wavelengths)
        OligoExtinctionOnly,    Extinction only - oligochromatic regime (a few discrete wavelengths)
        NoMedium,           No medium (primary sources only)
        ExtinctionOnly,     Extinction only (no secondary emission)
        LyaExtinctionOnly,  Extinction only with Lyman-alpha line transfer
        DustEmission,       With secondary emission from dust
        GasEmission,        With secondary emission from gas
        DustAndGasEmission, With secondary emission from dust and gas
        
    iteratePrimaryEmission: bool, iterate over primary emission for self-consistent calculation
    
    iterateSecondaryEmission: bool, iterate over secondary emission for self-consistent calculation
    
    numPackets: Double, the default number of photon packets launched per simulation segment
    '''
    keysall = {'userLevel': 'Regular', 
                'simulationMode': 'DustEmission', 
                'iteratePrimaryEmission': 'false', 
                'iterateSecondaryEmission': 'false', 
                'numPackets': '1e7'}
    ele = xdm.Document().createElement('MonteCarloSimulation')
    for i in keysall:
        ele.setAttribute(i, keysall[i])                 #TODO
    for i in ['random','units','cosmology','sourceSystem','mediumSystem','instrumentSystem','probeSystem']:
        print(f'Setting {i}')
        ele=ski_append(ele,ski[i])
        
    return ele

@Element_ski.ski_element
def random(ski):
    '''
    the random number generator
    '''
    ele = xdm.Document().createElement('random')
    ele.setAttribute('type', 'Random')
    ele = ski_append(ele,ski['Random'])
    return ele

@Element_ski.ski_element
def Random(ski):
    '''
    the default random generator
    seed: int, the seed for the random generator
    '''
    ele = xdm.Document().createElement('Random')
    ele.setAttribute('seed', '0')                   #TODO
    return ele

@Element_ski.ski_element
def units(ski):
    '''
    the units system
    '''
    ele = xdm.Document().createElement('units')
    ele.setAttribute('type', 'Units')
    ele = ski_append(ele,ski['Units'])
    return ele

@Element_ski.ski_element
def Units(ski):
    '''
    type: a units system
            ExtragalacticUnits:     extragalactic units (length in pc, distance in Mpc)
            SIUnits:                SI units
            StellarUnits:           stellar units (length in AU, distance in pc)
    '''
    subkey = ski.get_default('Units')  #TODO
    ele = ski[subkey]    #TODO 
    ele.setAttribute('wavelengthOutputStyle','Wavelength')  #TODO
    ele.setAttribute('fluxOutputStyle','Frequency')         #TODO
    return ele

@Element_ski.ski_element
def ExtragalacticUnits(ski):
    '''
    extragalactic units (length in pc, distance in Mpc)
    wavelengthOutputStyle: the output style for wavelengths
        Wavelength:             as photon wavelength: λ
        Frequency:              as photon frequency: ν
        Energy:                 as photon energy: E

    fluxOutputStyle:        the output style for flux density and surface brightness
        Neutral:                neutral: λ F_λ = ν F_ν
        Wavelength:             per unit of wavelength: F_λ
        Frequency:              per unit of frequency: F_ν
        Energy:                 counts per unit of energy: F_E
    '''
    ele = xdm.Document().createElement('ExtragalacticUnits')
    
    ele.setAttribute('wavelengthOutputStyle','Wavelength')  #TODO
    ele.setAttribute('fluxOutputStyle','Frequency')         #TODO
    return ele

@Element_ski.ski_element
def SIUnits(ski):
    '''
    SI units
    wavelengthOutputStyle: the output style for wavelengths
        Wavelength:             as photon wavelength: λ
        Frequency:              as photon frequency: ν
        Energy:                 as photon energy: E

    fluxOutputStyle:        the output style for flux density and surface brightness
        Neutral:                neutral: λ F_λ = ν F_ν
        Wavelength:             per unit of wavelength: F_λ
        Frequency:              per unit of frequency: F_ν
        Energy:                 counts per unit of energy: F_E
    '''
    ele = xdm.Document().createElement('SIUnits')
    
    ele.setAttribute('wavelengthOutputStyle','Wavelength')  #TODO
    ele.setAttribute('fluxOutputStyle','Frequency')         #TODO
    return ele
    
@Element_ski.ski_element
def StellarUnits(ski):
    '''
    stellar units (length in AU, distance in pc)
    wavelengthOutputStyle: the output style for wavelengths
        Wavelength:             as photon wavelength: λ
        Frequency:              as photon frequency: ν
        Energy:                 as photon energy: E

    fluxOutputStyle:        the output style for flux density and surface brightness
        Neutral:                neutral: λ F_λ = ν F_ν
        Wavelength:             per unit of wavelength: F_λ
        Frequency:              per unit of frequency: F_ν
        Energy:                 counts per unit of energy: F_E
    '''
    ele = xdm.Document().createElement('StellarUnits')
    
    ele.setAttribute('wavelengthOutputStyle','Wavelength')  #TODO
    ele.setAttribute('fluxOutputStyle','Frequency')         #TODO
    return ele

@Element_ski.ski_element
def cosmology(ski):
    '''
    the cosmology parameters
    '''
    ele = xdm.Document().createElement('cosmology')
    ele.setAttribute('type', 'Cosmology')
    ele = ski_append(ele,ski['Cosmology'])
    return ele

@Element_ski.ski_element
def Cosmology(ski):
    '''
    a set of cosmology parameters, including redshift
    '''
    subkey = ski.get_default('Cosmology')  #TODO
    
    ele = ski[subkey]    #TODO 
    return ele

@Element_ski.ski_element
def FlatUniverseCosmology(ski):
    '''
    the model is at a given redshift in a flat universe
    redshift, double,    the redshift z of the model coordinate frame
    reducedHubbleConstant, double,  the reduced Hubble constant h
    matterDensityFraction, double,  the cosmological matter density fraction Ω_m
    '''
    ele = xdm.Document().createElement('FlatUniverseCosmology')
    ele.setAttribute('redshift', '0.01')            #TODO
    ele.setAttribute('reducedHubbleConstant', '0.6774')            #TODO
    ele.setAttribute('matterDensityFraction', '0.3075')            #TODO
    return ele

@Element_ski.ski_element
def LocalUniverseCosmology(ski):
    '''
    the model is at redshift zero in the Local Universe
    '''
    ele = xdm.Document().createElement('LocalUniverseCosmology')
    return ele
