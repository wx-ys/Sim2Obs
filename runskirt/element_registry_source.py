
import xml.dom.minidom as xdm

from .skirt_element import Element_ski, ski_append



@Element_ski.ski_element
def sourceSystem(ski):
    '''
    the source system
    '''
    ele = xdm.Document().createElement('sourceSystem')
    ele.setAttribute('type', 'SourceSystem')
    ele = ski_append(ele,ski['SourceSystem'])
    return ele

@Element_ski.ski_element
def SourceSystem(ski):
    '''
    a primary source system
    
    minWavelength	Double	the shortest wavelength of photon packets launched from primary sources
    maxWavelength	Double	the longest wavelength of photon packets launched from primary sources
    wavelengths	    DoubleList	the discrete wavelengths of photon packets launched from primary sources
    sourceBias	    Double	the fraction of photon packets distributed uniformly across primary sources
    
    containing: sources
    '''
    ele = xdm.Document().createElement('SourceSystem')
    
    keysall = {'minWavelength': '0.08 micron', 
                   'maxWavelength': '1000 micron', 
                   'wavelengths': '0.55 micron', 
                   'sourceBias': '0.5'}
    for i in keysall:
        ele.setAttribute(i, keysall[i])             #TODO
        
    ele = ski_append(ele,ski['sources'])
    return ele

@Element_ski.ski_element
def sources(ski):
    '''
    the primary sources
    '''
    ele = xdm.Document().createElement('sources')
    ele.setAttribute('type', 'Source')
    for i in range(100):
        subkey = ski.get_default('sources')                    #TODO
        if subkey:
            ele = ski_append(ele,ski[subkey])
        else:
            break
    return ele

@Element_ski.ski_element
def AdaptiveMeshSource(ski):
    '''
    a primary source imported from data represented on an adaptive mesh (AMR grid)
    
    filename	String	the name of the file to be imported
    minX	Double	the start point of the domain in the X direction
    maxX	Double	the end point of the domain in the X direction
    minY	Double	the start point of the domain in the Y direction
    maxY	Double	the end point of the domain in the Y direction
    minZ	Double	the start point of the domain in the Z direction
    maxZ	Double	the end point of the domain in the Z direction
    importVelocity	Bool	import velocity components (3 columns)
    importVelocityDispersion	Bool	import velocity dispersion (spherically symmetric)
    importCurrentMass	Bool	import current mass
    useColumns	String	a list of names corresponding to columns in the file to be imported
    sourceWeight	Double	the weight of this source for the number of photon packets launched
    wavelengthBias	Double	the fraction of photon packet wavelengths sampled from a bias distribution
    
    containing: sedFamily, wavelengthBiasDistribution
    '''
    keysall={'filename': 'AdaptiveMeshSource.txt',
             'minX': '-1 kpc',
             'maxX': '1 kpc',
             'minY': '-1 kpc',
             'maxY': '1 kpc',
             'minZ': '-1 kpc',
             'maxZ': '1 kpc',
             'importVelocity': 'false',
             'importVelocityDispersion': 'false',
             'importCurrentMass': 'false',
             'useColumns': '',
             'sourceWeight': '1',
             'wavelengthBias': '0.5',
    }
    ele = xdm.Document().createElement('AdaptiveMeshSource')
    for i in keysall:
        ele.setAttribute(i, keysall[i])         #TODO

    con = ski.get_default('sources_options').split(',')
    for i in con:
        ele = ski_append(ele,ski[i])
    
    return ele

@Element_ski.ski_element
def CellSource(ski):
    '''
    a primary source imported from cuboidal cell data
    
    filename	String	the name of the file to be imported
    importVelocity	Bool	import velocity components (3 columns)
    importVelocityDispersion	Bool	import velocity dispersion (spherically symmetric)
    importCurrentMass	Bool	import current mass
    useColumns	String	a list of names corresponding to columns in the file to be imported
    sourceWeight	Double	the weight of this source for the number of photon packets launched
    wavelengthBias	Double	the fraction of photon packet wavelengths sampled from a bias distribution
    
    containing: sedFamily, wavelengthBiasDistribution
    '''
    keysall={'filename': 'AdaptiveMeshSource.txt',
             'importVelocity': 'false',
             'importVelocityDispersion': 'false',
             'importCurrentMass': 'false',
             'useColumns': '',
             'sourceWeight': '1',
             'wavelengthBias': '0.5',
    }
    ele = xdm.Document().createElement('CellSource')
    for i in keysall:
        ele.setAttribute(i, keysall[i])         #TODO

    con = ski.get_default('sources_options').split(',')
    for i in con:
        ele = ski_append(ele,ski[i])
    return ele

@Element_ski.ski_element
def CubicalBackgroundSource(ski):
    '''
    a cubical background source with an anisotropic inward radiation field
    
    centerX	Double	the center of the source, x component
    centerY	Double	the center of the source, y component
    centerZ	Double	the center of the source, z component
    edgeLength	Double	the edge length of the background cube
    velocityX	Double	the bulk velocity of the source, x component
    velocityY	Double	the bulk velocity of the source, y component
    velocityZ	Double	the bulk velocity of the source, z component
    useColumns	String	a list of names corresponding to columns in the file to be imported
    sourceWeight	Double	the weight of this source for the number of photon packets launched
    wavelengthBias	Double	the fraction of photon packet wavelengths sampled from a bias distribution
    
    containing: sed, normalization, wavelengthBiasDistribution
    '''
    keysall={'centerX': '0 kpc',
             'centerY': '0 kpc',
             'centerZ': '0 kpc',
             'edgeLength': '1 kpc',
             'velocityX': '0 km/s',
             'velocityY': '0 km/s',
             'velocityZ': '0 km/s',
             'useColumns': '',
             'sourceWeight': '1',
             'wavelengthBias': '0.5',
    }
    ele = xdm.Document().createElement('CubicalBackgroundSource')
    for i in keysall:
        ele.setAttribute(i, keysall[i])         #TODO

    con = ski.get_default('sources_options').split(',')
    for i in con:
        ele = ski_append(ele,ski[i])
    return ele

@Element_ski.ski_element
def FilePolarizedPointSource(ski):
    '''
    a primary point source with a polarized spectrum read from file
    
    filename	String	the name of the file to be imported
    positionX	Double	the position of the point source, x component
    positionY	Double	the position of the point source, y component
    positionZ	Double	the position of the point source, z component
    symmetryX	Double	the direction of the positive symmetry axis, x component
    symmetryY	Double	the direction of the positive symmetry axis, y component
    symmetryZ	Double	the direction of the positive symmetry axis, z component
    velocityX	Double	the bulk velocity of the point source, x component
    velocityY	Double	the bulk velocity of the point source, y component
    velocityZ	Double	the bulk velocity of the point source, z component
    sourceWeight	Double	the weight of this source for the number of photon packets launched
    wavelengthBias	Double	the fraction of photon packet wavelengths sampled from a bias distribution
    
    containing: normalization, wavelengthBiasDistribution
    '''
    keysall={'filename': 'FilePolarizedPointSource.txt',
             'positionX': '0 kpc',
             'positionY': '0 kpc',
             'positionZ': '0 kpc',
             'symmetryX': '1',
             'symmetryY': '1',
             'symmetryZ': '1',
             'velocityX': '0 km/s',
             'velocityY': '0 km/s',
             'velocityZ': '0 km/s',
             'sourceWeight': '1',
             'wavelengthBias': '0.5',
    }
    ele = xdm.Document().createElement('FilePolarizedPointSource')
    for i in keysall:
        ele.setAttribute(i, keysall[i])         #TODO
    
    con = ski.get_default('sources_options').split(',')
    for i in con:
        ele = ski_append(ele,ski[i])
    return ele

@Element_ski.ski_element
def GeometricSource(ski):
    '''
    a primary source with a built-in geometry
    
    velocityMagnitude	Double	the magnitude of the velocity (multiplier)
    sourceWeight	Double	the weight of this source for the number of photon packets launched
    wavelengthBias	Double	the fraction of photon packet wavelengths sampled from a bias distribution

    containing: geometry, velocityDistribution, sed, normalization, wavelengthBiasDistribution
    '''
    keysall={'velocityMagnitude': '1',
             'sourceWeight': '1',
             'wavelengthBias': '0.5',
    }
    ele = xdm.Document().createElement('GeometricSource')
    for i in keysall:
        ele.setAttribute(i, keysall[i])         #TODO

    con = ski.get_default('sources_options').split(',')
    for i in con:
        ele = ski_append(ele,ski[i])
    return ele

@Element_ski.ski_element
def ParticleSource(ski):
    '''
    a primary source imported from smoothed particle data
    
    filename	String	the name of the file to be imported
    importVelocity	Bool	import velocity components (3 columns)
    importVelocityDispersion	Bool	import velocity dispersion (spherically symmetric)
    importCurrentMass	Bool	import current mass
    useColumns	String	a list of names corresponding to columns in the file to be imported
    sourceWeight	Double	the weight of this source for the number of photon packets launched
    wavelengthBias	Double	the fraction of photon packet wavelengths sampled from a bias distribution

    containing: smoothingKernel, sedFamily, wavelengthBiasDistribution
    '''
    keysall={'filename': 'ParticleSource.txt',
             'importVelocity': 'false',
             'importVelocityDispersion': 'false',
             'importCurrentMass': 'false',
             'useColumns': '',
             'sourceWeight': '1',
             'wavelengthBias': '0.5',
    }
    ele = xdm.Document().createElement('ParticleSource')
    for i in keysall:
        ele.setAttribute(i, keysall[i])         #TODO

    con = ski.get_default('sources_options').split(',')
    for i in con:
        ele = ski_append(ele,ski[i])
    return ele

@Element_ski.ski_element
def SphericalBackgroundSource(ski):
    '''
    a spherical background source with an anisotropic inward radiation field
    
    centerX	Double	the center of the source, x component
    centerY	Double	the center of the source, y component
    centerZ	Double	the center of the source, z component
    backgroundRadius	Double	the radius of the background sphere
    velocityX	Double	the bulk velocity of the source, x component
    velocityY	Double	the bulk velocity of the source, y component
    velocityZ	Double	the bulk velocity of the source, z component
    sourceWeight	Double	the weight of this source for the number of photon packets launched
    wavelengthBias	Double	the fraction of photon packet wavelengths sampled from a bias distribution


    containing: sed, normalization, wavelengthBiasDistribution
    '''
    keysall={'centerX': '0 kpc',
             'centerY': '0 kpc',
             'centerZ': '0 kpc',
             'backgroundRadius': '1 kpc',
             'velocityX': '0 km/s',
             'velocityY': '0 km/s',
             'velocityZ': '0 km/s',
             'sourceWeight': '1',
             'wavelengthBias': '0.5',
    }
    ele = xdm.Document().createElement('SphericalBackgroundSource')
    for i in keysall:
        ele.setAttribute(i, keysall[i])         #TODO
        
    con = ski.get_default('sources_options').split(',')
    for i in con:
        ele = ski_append(ele,ski[i])
    return ele

@Element_ski.ski_element
def StellarSurfaceSource(ski):
    '''
    a stellar surface source with an anisotropic outward radiation field
    
    centerX	Double	the center of the source, x component
    centerY	Double	the center of the source, y component
    centerZ	Double	the center of the source, z component
    stellarRadius	Double	the stellar radius
    velocityX	Double	the bulk velocity of the source, x component
    velocityY	Double	the bulk velocity of the source, y component
    velocityZ	Double	the bulk velocity of the source, z component
    sourceWeight	Double	the weight of this source for the number of photon packets launched
    wavelengthBias	Double	the fraction of photon packet wavelengths sampled from a bias distribution


    containing: sed, normalization, wavelengthBiasDistribution
    '''
    keysall={'centerX': '0 kpc',
             'centerY': '0 kpc',
             'centerZ': '0 kpc',
             'stellarRadius': '1 AU',
             'velocityX': '0 km/s',
             'velocityY': '0 km/s',
             'velocityZ': '0 km/s',
             'sourceWeight': '1',
             'wavelengthBias': '0.5',
    }
    ele = xdm.Document().createElement('SphericalBackgroundSource')
    for i in keysall:
        ele.setAttribute(i, keysall[i])         #TODO
        
    con = ski.get_default('sources_options').split(',')
    for i in con:
        ele = ski_append(ele,ski[i])
    return ele

@Element_ski.ski_element
def VoronoiMeshSource(ski):
    '''
    a primary source imported from data represented on a Voronoi mesh
    
    filename	String	the name of the file to be imported
    minX	Double	the start point of the domain in the X direction
    maxX	Double	the end point of the domain in the X direction
    minY	Double	the start point of the domain in the Y direction
    maxY	Double	the end point of the domain in the Y direction
    minZ	Double	the start point of the domain in the Z direction
    maxZ	Double	the end point of the domain in the Z direction
    importVelocity	Bool	import velocity components (3 columns)
    importVelocityDispersion	Bool	import velocity dispersion (spherically symmetric)
    importCurrentMass	Bool	import current mass
    useColumns	String	a list of names corresponding to columns in the file to be imported
    sourceWeight	Double	the weight of this source for the number of photon packets launched
    wavelengthBias	Double	the fraction of photon packet wavelengths sampled from a bias distribution


    containing: sedFamily, wavelengthBiasDistribution
    '''
    keysall={'filename': 'VoronoiMeshSource.txt',
             
             'minX': '-1 kpc',
             'maxX': '1 kpc',
             'minY': '-1 kpc',
             'maxY': '1 kpc',
             'minZ': '-1 kpc',
             'maxZ': '1 kpc',
             'importVelocity': 'false',
             'importVelocityDispersion': 'false',
             'importCurrentMass': 'false',
             'useColumns': '',
             'sourceWeight': '1',
             'wavelengthBias': '0.5',
    }
    ele = xdm.Document().createElement('SphericalBackgroundSource')
    for i in keysall:
        ele.setAttribute(i, keysall[i])         #TODO
        
    con = ski.get_default('sources_options').split(',')
    for i in con:
        ele = ski_append(ele,ski[i])
    return ele

@Element_ski.ski_element
def smoothingKernel(ski):
    '''
    the kernel for interpolating the smoothed particles
    '''
    ele = xdm.Document().createElement('smoothingKernel')
    ele.setAttribute('type','SmoothingKernel')
    ele =ski_append(ele,ski['SmoothingKernel'])
    
    return ele

@Element_ski.ski_element
def SmoothingKernel(ski):
    '''
    a smoothing kernel
    CubicSplineSmoothingKernel, ScaledGaussianSmoothingKernel, UniformSmoothingKernel
    '''
    subkey = ski.get_default('smoothingKernel') #TODO
    ele = ski[subkey]
    return ele

@Element_ski.ski_element
def CubicSplineSmoothingKernel(ski):
    '''
    a cubic spline smoothing kernel
    '''
    ele = xdm.Document().createElement('CubicSplineSmoothingKernel')
    return ele

@Element_ski.ski_element
def ScaledGaussianSmoothingKernel(ski):
    '''
    a scaled Gaussian smoothing kernel
    '''
    ele = xdm.Document().createElement('ScaledGaussianSmoothingKernel')
    return ele

@Element_ski.ski_element
def UniformSmoothingKernel(ski):
    '''
    a uniform smoothing kernel
    '''
    ele = xdm.Document().createElement('UniformSmoothingKernel')
    return ele

@Element_ski.ski_element
def sedFamily(ski):
    '''
    the SED family for assigning spectra to the imported sources
    '''
    ele = xdm.Document().createElement('sedFamily')
    ele.setAttribute('type','SEDFamily')
    ele =ski_append(ele,ski['SEDFamily'])
    return ele

@Element_ski.ski_element
def SEDFamily(ski):
    '''
    the SED family for assigning spectra to the imported sources
    '''
  
    subkey = ski.get_default('sedFamily')
    if subkey:
        ele = ski[subkey]           #TODO
        return ele
    else:
        print(f'Some wrong! use SEDFamily: BruzualCharlotSEDFamily')
        ele = ski['BruzualCharlotSEDFamily']          
        return ele
@Element_ski.ski_element
def BlackBodySEDFamily(ski):
    '''
    a black body SED family
    '''
    ele = xdm.Document().createElement('BlackBodySEDFamily')
    return ele

@Element_ski.ski_element
def BpassSEDFamily(ski):
    '''
    a BPASS SED family for single stellar populations
    '''
    ele = xdm.Document().createElement('BpassSEDFamily')
    return ele

@Element_ski.ski_element
def BruzualCharlotSEDFamily(ski):
    '''
    a Bruzual-Charlot SED family for single stellar populations
    imf, the assumed initial mass function, (Chabrier, Salpeter)
    resolution, the wavelength resolution, (Low, High), 1221 points, 6900 points
    '''
    keysall = {'imf': 'Chabrier',
               'resolution': 'High'}
    ele = xdm.Document().createElement('BruzualCharlotSEDFamily')
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    return ele

@Element_ski.ski_element
def CastelliKuruczSEDFamily(ski):
    '''
    a Castelli-Kurucz SED family for stellar atmospheres
    '''
    ele = xdm.Document().createElement('CastelliKuruczSEDFamily')
    return ele

@Element_ski.ski_element
def FSPSSEDFamily(ski):
    '''
    an FSPS SED family for single stellar populations
    imf : the assumed initial mass function, (Chabrier, Kroupa, Salpeter)
    '''
    keysall = {'imf': 'Chabrier'}
    ele = xdm.Document().createElement('FSPSSEDFamily')
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    return ele

@Element_ski.ski_element
def FileIndexedSEDFamily(ski):
    '''
    a user-provided, indexed SED family
    filename : the name of the stored table file listing the SEDs
    '''
    keysall = {'filename': 'filename.txt'}
    ele = xdm.Document().createElement('FileIndexedSEDFamily')
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    return ele

@Element_ski.ski_element
def FileSSPSEDFamily(ski):
    '''
    a user-provided SED family for single stellar populations
    filename :String the name of the stored table file defining the SED templates
    hasIonizationParameter:Bool include ionization U as an extra parameter for the SED family
    '''
    keysall = {'filename': 'filename.txt',
               'hasIonizationParameter': 'false'}
    ele = xdm.Document().createElement('FileSSPSEDFamily')
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    return ele

@Element_ski.ski_element
def LyaDoublePeakedSEDFamily(ski):
    '''
    a family of double-peaked spectra around the central Lyman-alpha wavelength
    '''
    ele = xdm.Document().createElement('LyaDoublePeakedSEDFamily')
    return ele

@Element_ski.ski_element
def LyaGaussianSEDFamily(ski):
    '''
    a family of Gaussian spectra around the central Lyman-alpha wavelength
    '''
    ele = xdm.Document().createElement('LyaGaussianSEDFamily')
    return ele

@Element_ski.ski_element
def LyaSEDFamilyDecorator(ski):
    '''
    an SED family decorator replacing ionizing radiation with Lyman-alpha line emission
    conversionFraction	Double	the fraction of ionizing radiation replaced by Lyman-alpha emission

    containing: sedFamilyOriginal, sedLymanAlpha
    '''
    ele = xdm.Document().createElement('LyaSEDFamilyDecorator')
    ele.setAttribute('conversionFraction', '0.2')  #TODO
    ele = ski_append(ele, ski['sedFamilyOriginal'])
    ele = ski_append(ele, ski['sedLymanAlpha'])
    return ele

@Element_ski.ski_element
def MappingsSEDFamily(ski):
    '''
    a MAPPINGS III SED family for star-forming regions
    '''
    ele = xdm.Document().createElement('MappingsSEDFamily')
    return ele

@Element_ski.ski_element
def MarastonSEDFamily(ski):
    '''
    a Maraston SED family for single stellar populations
    imf: the assumed initial mass function, (Kroupa, Salpeter)
    '''
    ele = xdm.Document().createElement('MarastonSEDFamily')
    ele.setAttribute('imf', 'Kroupa')
    return ele

@Element_ski.ski_element
def SpinFlipSEDFamily(ski):
    '''
    a family of Gaussian spectra around the central spin-flip wavelength
    '''
    ele = xdm.Document().createElement('SpinFlipSEDFamily')
    return ele

@Element_ski.ski_element
def Starburst99SEDFamily(ski):
    '''
    a Starburst99 SED family for single stellar populations
    '''
    ele = xdm.Document().createElement('Starburst99SEDFamily')
    return ele

@Element_ski.ski_element
def ToddlersSEDFamily(ski):
    '''
    a Toddlers SED family for emission from star-forming regions
    pahfraction: the maximum PAH-to-dust fraction, (High,Low), High PAH-to-dust fraction (4.6%), Low PAH-to-dust fraction (1%)
    resolution: he wavelength resolution, (Low,High), (continuum and lines at R=300), (continuum at R=300 and lines at R=5e4)
    '''
    ele = xdm.Document().createElement('ToddlersSEDFamily')
    ele.setAttribute('pahfraction', 'High') # TODO
    ele.setAttribute('resolution', 'High') # TODO
    return ele

@Element_ski.ski_element
def wavelengthBiasDistribution(ski):
    '''
    the bias distribution for sampling photon packet wavelengths
    DefaultWavelengthDistribution, DiscreteWavelengthDistribution, FileWavelengthDistribution, LinWavelengthDistribution
    ListWavelengthDistribution, LogWavelengthDistribution
    '''
    ele = xdm.Document().createElement('wavelengthBiasDistribution')
    ele.setAttribute('type', 'WavelengthDistribution')
    subkey = ski.get_default('wavelengthBiasDistribution')
    if subkey:
        ele = ski[subkey]
        return ele
    else:
        print('Some wrong, use wavelengthBiasDistribution: LogWavelengthDistribution')
        ele = ski['LogWavelengthDistribution']
        return ele


@Element_ski.ski_element
def WavelengthDistribution(ski):
    '''
    a wavelength probability distribution
    '''
    subkey = ski.get_default('WavelengthDistribution') #TODO
    if subkey:
        ele = ski[subkey]
        return ele
    else:
        print('Some wrong, use WavelengthDistribution: LogWavelengthDistribution')
        ele = ski['LogWavelengthDistribution']
        return ele
@Element_ski.ski_element
def DefaultWavelengthDistribution(ski):
    '''
    the default (logarithmic) wavelength probability distribution
    '''
    ele = xdm.Document().createElement('DefaultWavelengthDistribution')
    return ele

@Element_ski.ski_element
def DiscreteWavelengthDistribution(ski):
    '''
    a discrete wavelength probability distribution derived from a wavelength grid
    containing: wavelengthGrid
    '''
    ele = xdm.Document().createElement('DiscreteWavelengthDistribution')
    ele = ski_append(ele,ski['wavelengthGrid'])
    return ele

@Element_ski.ski_element
def wavelengthGrid(ski):
    '''
    the wavelength grid for this instrument
    containing: CompositeWavelengthGrid	 ConfigurableBandWavelengthGrid FileBorderWavelengthGrid FileWavelengthGrid LinBorderWavelengthGrid LinWavelengthGrid
        ListBorderWavelengthGrid  ListWavelengthGrid LogBorderWavelengthGrid LogWavelengthGrid NestedLogWavelengthGrid PredefinedBandWavelengthGrid
        ResolutionBorderWavelengthGrid ResolutionWavelengthGrid
    '''
    ele = xdm.Document().createElement('wavelengthGrid')
    ele.setAttribute('type', 'WavelengthGrid')
    subkey = ski.get_default('wavelengthGrid')    #TODO
    if subkey:
        ele = ski_append(ele, ski[subkey])
        return ele
    else:
        print('Some wrong! use wavelengthGrid: LogWavelengthGrid')
        ele = ski_append(ele, ski['LogWavelengthGrid'])
        return ele

@Element_ski.ski_element
def FileWavelengthDistribution(ski):
    '''
    a wavelength probability distribution loaded from a text file
    filename String  the name of the file with the wavelength probability distribution
    '''
    ele = xdm.Document().createElement('FileWavelengthDistribution')
    ele.setAttribute('filename', 'FileWavelengthDistribution.txt')  #TODO
    return ele

@Element_ski.ski_element
def LinWavelengthDistribution(ski):
    '''
    a linear wavelength probability distribution
    minWavelength	Double	the shortest wavelength of the wavelength probability distribution
    maxWavelength	Double	the longest wavelength of the wavelength probability distribution
    '''
    ele = xdm.Document().createElement('LinWavelengthDistribution')
    ele.setAttribute('minWavelength', '0.001 micron')  #TODO
    ele.setAttribute('maxWavelength', '1e4 micron')  #TODO
    return ele

@Element_ski.ski_element
def ListWavelengthDistribution(ski):
    '''
    a wavelength probability distribution specified inside the configuration file
    wavelengths	    DoubleList	the wavelengths at which to specify the probability
    unitStyle	the probability unit style, (neutralmonluminosity, wavelengthmonluminosity, frequencymonluminosity, energymonluminosity)
    probabilities   DoubleList  the probabilities at each of the given wavelengths
    '''
    ele = xdm.Document().createElement('ListWavelengthDistribution')
    
    ele.setAttribute('wavelengths', '')  #TODO
    ele.setAttribute('unitStyle', 'wavelengthmonluminosity')  #TODO
    ele.setAttribute('probabilities', '')  #TODO
    return ele

@Element_ski.ski_element
def LogWavelengthDistribution(ski):
    '''
    a logarithmic wavelength probability distribution
    minWavelength	Double	the shortest wavelength of the wavelength probability distribution
    maxWavelength	Double	the longest wavelength of the wavelength probability distribution
    '''
    ele = xdm.Document().createElement('LogWavelengthDistribution')
    
    ele.setAttribute('minWavelength', '0.001 micron')  #TODO
    ele.setAttribute('maxWavelength', '1e4 micron')  #TODO
    return ele

@Element_ski.ski_element
def sed(ski):
    '''
    the spectral energy distribution for the source
    
    '''
    ele = xdm.Document().createElement('sed')
    ele.setAttribute('type', 'SED')
    ele = ski_append(ele, ski['SED'])
    return ele

@Element_ski.ski_element
def SED(ski):
    '''
    a spectral energy distribution
    BlackBodySED, BpassSED, BruzualCharlotSED, CastelliKuruczSED, FSPSSED, FileLineSED
    FileSED, ListLineSED, ListSED, LyaDoublePeakedSED, LyaGaussianSED, LyaSEDDecorator
    MappingsSED, MarastonSED, QuasarSED, SingleWavelengthSED, Starburst99SED, SunSED, ToddlersSED
    '''
    subkey = ski.get_default('SED') #TODO
    if subkey:
        ele = ski[subkey]
        return ele
    else:
        print('Some Wrong, use SED: BpassSED')
        ele = ski['BpassSED']
        return ele

@Element_ski.ski_element
def BlackBodySED(ski):
    '''
    a black-body spectral energy distribution
    temperature Double the black body temperature
    '''
    ele = xdm.Document().createElement('BlackBodySED')
    ele.setAttribute('temperature', '6000 K')       #TODO
    return ele

@Element_ski.ski_element
def BpassSED(ski):
    '''
    a BPASS single stellar population SED
    
    metallicity	Double	the metallicity of the SSP
    age	Double	the age of the SSP
    '''
    ele = xdm.Document().createElement('BpassSED')
    ele.setAttribute('metallicity', '0.02')
    ele.setAttribute('age', '5 Gyr')    #TODO
    return ele


@Element_ski.ski_element
def BruzualCharlotSED(ski):
    '''
    a Bruzual-Charlot simple stellar population SED
    imf, the assumed initial mass function, (Chabrier, Salpeter)
    metallicity	Double	the metallicity of the SSP
    age	Double	the age of the SSP
    '''
    ele = xdm.Document().createElement('BruzualCharlotSED')
    ele.setAttribute('imf', 'Chabrier')
    ele.setAttribute('metallicity', '0.02')
    ele.setAttribute('age', '5 Gyr')    #TODO
    return ele

@Element_ski.ski_element
def CastelliKuruczSED(ski):
    '''
    a Castelli-Kurucz stellar atmosphere SED
    metallicity	Double	the metallicity of the SSP
    age	Double	the age of the SSP
    gravity	Double	the surface gravity
    '''
    ele = xdm.Document().createElement('CastelliKuruczSED')
    ele.setAttribute('metallicity', '0.02')
    ele.setAttribute('age', '5 Gyr')    #TODO
    ele.setAttribute('gravity', '')
    return ele

@Element_ski.ski_element
def FSPSSED(ski):
    '''
    an FSPS simple stellar population SED
    imf, the assumed initial mass function, (Chabrier, Kroupa, Salpeter)
    metallicity	Double	the metallicity of the SSP
    age	Double	the age of the SSP
    '''
    ele = xdm.Document().createElement('FSPSSED')
    ele.setAttribute('imf', 'Chabrier')
    ele.setAttribute('metallicity', '0.02')
    ele.setAttribute('age', '5 Gyr')    #TODO
    return ele

@Element_ski.ski_element
def FileLineSED(ski):
    '''
    a discrete line SED loaded from a text file
    filename	String	the name of the file with the line definitions
    '''
    ele = xdm.Document().createElement('FileLineSED')
    ele.setAttribute('filename', 'FileLineSED.txt')    #TODO
    return ele

@Element_ski.ski_element
def FileSED(ski):
    '''
    a spectral energy distribution loaded from a text file
    filename	String	the name of the file with the spectral energy distribution
    '''
    ele = xdm.Document().createElement('FileSED')
    ele.setAttribute('filename', 'FileSED.txt')    #TODO
    return ele

@Element_ski.ski_element
def ListLineSED(ski):
    '''
    a discrete line SED specified inside the configuration file
    wavelengths 	DoubleList	the line wavelengths
    luminosities	DoubleList	the line luminosities
    '''
    ele = xdm.Document().createElement('ListLineSED')
    ele.setAttribute('wavelengths', '')    #TODO
    ele.setAttribute('luminosities', '')    #TODO
    return ele

@Element_ski.ski_element
def ListSED(ski):
    '''
    a spectral energy distribution specified inside the configuration file
    wavelengths 	DoubleList	the line wavelengths
    unitStyle (neutralmonluminosity, wavelengthmonluminosity, frequencymonluminosity, energymonluminosity)
    specificLuminosities	DoubleList	the specific luminosities at each of the given wavelengths
    '''
    ele = xdm.Document().createElement('ListSED')
    ele.setAttribute('wavelengths', '')    #TODO
    ele.setAttribute('unitStyle', 'wavelengthmonluminosity')    #TODO
    ele.setAttribute('specificLuminosities', '')    #TODO
    return ele

@Element_ski.ski_element
def LyaDoublePeakedSED(ski):
    '''
    a double-peaked spectrum around the central Lyman-alpha wavelength
    scale	Double	the velocity scale
    '''
    ele = xdm.Document().createElement('LyaDoublePeakedSED')
    ele.setAttribute('scale', '1')    #TODO
    return ele

@Element_ski.ski_element
def LyaGaussianSED(ski):
    '''
    a Gaussian spectrum around the central Lyman-alpha wavelength
    dispersion	Double	the Gaussian velocity dispersion
    '''
    ele = xdm.Document().createElement('LyaGaussianSED')
    ele.setAttribute('dispersion', '1')    #TODO
    return ele

@Element_ski.ski_element
def LyaSEDDecorator(ski):
    '''
    an SED decorator replacing ionizing radiation with Lyman-alpha line emission
    conversionFraction	Double	the fraction of ionizing radiation replaced by Lyman-alpha emission
    containing: sedOriginal, sedLymanAlpha
    '''
    ele = xdm.Document().createElement('LyaSEDDecorator')
    ele.setAttribute('conversionFraction', '0.2')    #TODO
    ele = ski_append(ele, ski['sedOriginal'])
    ele = ski_append(ele, ski['sedLymanAlpha'])
    return ele

@Element_ski.ski_element
def MappingsSED(ski):
    '''
    a star-forming region SED from the MAPPINGS III model
    metallicity	Double	the metallicity
    compactness	Double	the logarithm of the compactness parameter
    pressure	Double	the ISM pressure
    coveringFactor	Double	the PDR covering factor
    '''
    ele = xdm.Document().createElement('MappingsSED')
    ele.setAttribute('metallicity', '0.02')    #TODO
    ele.setAttribute('compactness', '5')    #TODO
    ele.setAttribute('pressure', '1.38e-12')    #TODO
    ele.setAttribute('coveringFactor', '0.2')    #TODO
    return ele

@Element_ski.ski_element
def MarastonSED(ski):
    '''
    a Maraston simple stellar population SED
    imf (Kroupa, Salpeter)
    metallicity	Double	the metallicity of the SSP
    age	Double	the age of the SSP
    '''
    ele = xdm.Document().createElement('MarastonSED')
    ele.setAttribute('metallicity', '0.02')    #TODO
    ele.setAttribute('age', '5 Gyr')    #TODO
    ele.setAttribute('imf', 'Kroupa')    #TODO
    return ele


@Element_ski.ski_element
def SingleWavelengthSED(ski):
    '''
    a single-wavelength SED in the form of a Dirac-delta function
    wavelength	Double	the single emission wavelength
    '''
    ele = xdm.Document().createElement('SingleWavelengthSED')
    ele.setAttribute('wavelength', '0.55 micron') #TODO
    return ele

@Element_ski.ski_element
def Starburst99SED(ski):
    '''
    a Starburst99 simple stellar population SED
    metallicity		Double	the metallicity of the SSP
    age	Double	the age of the SSP
    '''
    ele = xdm.Document().createElement('Starburst99SED')
    ele.setAttribute('metallicity', '0.02') #TODO
    ele.setAttribute('age', '5 Gyr') #TODO
    return ele

@Element_ski.ski_element
def SunSED(ski):
    '''
    the spectral energy distribution of the Sun
    '''
    ele = xdm.Document().createElement('SunSED')
    return ele


@Element_ski.ski_element
def ToddlersSED(ski):
    '''
    a Toddlers SED for emission from star-forming regions
    pahfraction the maximum PAH-to-dust fraction (High,Low ) High PAH-to-dust fraction (4.6%) Low	Low PAH-to-dust fraction (1%)
    resolution (Low, High), (continuum and lines at R=300), (continuum at R=300 and lines at R=5e4)
    age	Double	system age
    metallicity	Double	system metallicity
    SFE	Double	star formation efficiency
    cloudNumDensity	Double	natal cloud number density
    '''
    ele = xdm.Document().createElement('ToddlersSED')   
    ele.setAttribute('pahfraction', 'High')             #TODO
    ele.setAttribute('resolution', 'High')              #TODO
    ele.setAttribute('age', '2.5 Myr')                  #TODO
    ele.setAttribute('metallicity', '0.02')             #TODO
    ele.setAttribute('SFE', '0.025')                    #TODO
    ele.setAttribute('cloudNumDensity', '320 /cm3')     #TODO
    return ele


@Element_ski.ski_element
def QuasarSED(ski):
    '''
    the spectral energy distribution of a typical quasar
    '''
    ele = xdm.Document().createElement('QuasarSED')
    return ele
