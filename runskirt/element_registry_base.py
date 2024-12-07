import xml.dom.minidom as xdm

from .skirt_element import Element_ski, ski_append

@Element_ski.ski_element
def instrumentSystem(ski):
    '''
    the instrument system
    '''
    ele = xdm.Document().createElement('instrumentSystem')
    ele.setAttribute('type', 'InstrumentSystem')
    ele = ski_append(ele, ski['InstrumentSystem'])
    return ele

@Element_ski.ski_element
def InstrumentSystem(ski):
    '''
    an instrument system
    '''
    ele = xdm.Document().createElement('InstrumentSystem')
    ele = ski_append(ele, ski['defaultWavelengthGrid'])
    ele = ski_append(ele, ski['instruments'])
    return ele


@Element_ski.ski_element
def defaultWavelengthGrid(ski):
    '''
    the default instrument wavelength grid
    containing: CompositeWavelengthGrid	 ConfigurableBandWavelengthGrid FileBorderWavelengthGrid FileWavelengthGrid LinBorderWavelengthGrid LinWavelengthGrid
        ListBorderWavelengthGrid  ListWavelengthGrid LogBorderWavelengthGrid LogWavelengthGrid NestedLogWavelengthGrid PredefinedBandWavelengthGrid
        ResolutionBorderWavelengthGrid ResolutionWavelengthGrid
    '''
    ele = xdm.Document().createElement('defaultWavelengthGrid')
    ele.setAttribute('type', 'WavelengthGrid')
    subkey = ski.get_default('defaultWavelengthGrid')    #TODO
    ele = ski_append(ele, ski[subkey])
    return ele

@Element_ski.ski_element
def instruments(ski):
    '''
    the instruments
    '''
    ele = xdm.Document().createElement('instruments')
    ele.setAttribute('type', 'Instrument')
    ele = ski_append(ele, ski['Instrument'])
    return ele

@Element_ski.ski_element
def Instrument(ski):
    '''
    an instrument
    containing: AllSkyInstrument FrameInstrument FullInstrument HEALPixSkyInstrument PerspectiveInstrument SEDInstrument
    '''
    ele = xdm.Document().createElement('Instrument')
    for i in range(100):
        subkey = ski.get_default('Instrument') #TODO ItemList
        if subkey:
            ele = ski_append(ele, ski[subkey])
        else:
            break
    return ele


@Element_ski.ski_element
def AllSkyInstrument(ski):
    '''
    an all-sky instrument (for observing inside a model)
    instrumentName	String	the name for this instrument
    numPixelsY	Int	the number of image pixels in the vertical (shortest) direction
    radius	Double	the radius of the observer's all-sky sphere
    observerX	Double	the position of the observer, x component
    observerY	Double	the position of the observer, y component
    observerZ	Double	the position of the observer, z component
    crossX	Double	the position of the crosshair, x component
    crossY	Double	the position of the crosshair, y component
    crossZ	Double	the position of the crosshair, z component
    upX	Double	the upwards direction, x component
    upY	Double	the upwards direction, y component
    upZ	Double	the upwards direction, z component
    recordComponents	Bool	record flux components separately
    numScatteringLevels	Int	the number of individually recorded scattering levels
    recordPolarization	Bool	record polarization (Stokes vector elements)
    recordStatistics	Bool	record information for calculating statistical properties
    '''
    keysall={'instrumentName': 'instrumentName',
             'numPixelsY': '250',
             'radius': '30 kpc',
             'observerX': '10 kpc',
             'observerY': '10 kpc',
             'observerZ': '10 kpc',
             'crossX': '5 kpc',
             'crossY': '5 kpc',
             'crossZ': '5 kpc',
             'upX': '5 kpc',
             'upY': '5 kpc',
             'upZ': '5 kpc',
             'recordComponents': 'false',
             'numScatteringLevels': '0',
             'recordPolarization': 'false',
             'recordStatistics': 'false',}
    ele = xdm.Document().createElement('AllSkyInstrument')
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
        
    ele = ski_append(ele, ski['wavelengthGrid'])
    ele = ski_append(ele, ski['projection'])
    return ele

@Element_ski.ski_element
def projection(ski):
    '''
    the projection used for mapping the sky to a rectangle
    '''
    ele = xdm.Document().createElement('projection')
    ele.setAttribute('type', 'AllSkyProjection')
    ele = ski_append(ele, ski['AllSkyProjection'])
    return ele

@Element_ski.ski_element
def AllSkyProjection(ski):
    '''
    an all-sky projection
    containing: HammerAitoffProjection MollweideProjection
    '''
    ele = xdm.Document().createElement('AllSkyProjection')
    subkey = ski.xxx #TODO
    ele = ski_append(ele, ski[subkey])
    return ele

@Element_ski.ski_element
def HammerAitoffProjection(ski):
    '''
    the Hammer-Aitoff all-sky projection
    '''
    ele = xdm.Document().createElement('HammerAitoffProjection')
    return ele

@Element_ski.ski_element
def MollweideProjection(ski):
    '''
    the Mollweide all-sky projection
    '''
    ele = xdm.Document().createElement('MollweideProjection')
    return ele

@Element_ski.ski_element
def FrameInstrument(ski):
    '''
    a distant instrument that outputs the surface brightness in every pixel as a data cube
    instrumentName	String	the name for this instrument
    distance	Double	the distance to the system
    inclination	Double	the inclination angle θ of the detector
    azimuth	Double	the azimuth angle φ of the detector
    roll	Double	the roll angle ω of the detector
    fieldOfViewX	Double	the total field of view in the horizontal direction
    numPixelsX	Int	the number of pixels in the horizontal direction
    centerX	Double	the center of the frame in the horizontal direction
    fieldOfViewY	Double	the total field of view in the vertical direction
    numPixelsY	Int	the number of pixels in the vertical direction
    centerY	Double	the center of the frame in the vertical direction
    recordComponents	Bool	record flux components separately
    numScatteringLevels	Int	the number of individually recorded scattering levels
    recordPolarization	Bool	record polarization (Stokes vector elements)
    recordStatistics	Bool	record information for calculating statistical properties
    '''
    keysall={'instrumentName': 'instrumentName',
             'distance': '10 Mpc',
             'inclination': '0 deg',
             'azimuth': '0 deg',
             'roll': '0 deg',
             'fieldOfViewX': '50 kpc',
             'numPixelsX': '250',
             'centerX': '0 kpc',
             'fieldOfViewY': '50 kpc',
             'numPixelsY': '250',
             'centerY': '0 kpc',
             'recordComponents': 'false',
             'numScatteringLevels': '0',
             'recordPolarization': 'false',
             'recordStatistics': 'false',
    }
    ele = xdm.Document().createElement('FrameInstrument')
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    ele = ski_append(ele, ski['wavelengthGrid'])
    return ele


@Element_ski.ski_element
def FullInstrument(ski):
    '''
    a distant instrument that outputs both the flux density (SED) and surface brightness (data cube)
    instrumentName	String	the name for this instrument
    distance	Double	the distance to the system
    inclination	Double	the inclination angle θ of the detector
    azimuth	Double	the azimuth angle φ of the detector
    roll	Double	the roll angle ω of the detector
    fieldOfViewX	Double	the total field of view in the horizontal direction
    numPixelsX	Int	the number of pixels in the horizontal direction
    centerX	Double	the center of the frame in the horizontal direction
    fieldOfViewY	Double	the total field of view in the vertical direction
    numPixelsY	Int	the number of pixels in the vertical direction
    centerY	Double	the center of the frame in the vertical direction
    recordComponents	Bool	record flux components separately
    numScatteringLevels	Int	the number of individually recorded scattering levels
    recordPolarization	Bool	record polarization (Stokes vector elements)
    recordStatistics	Bool	record information for calculating statistical properties
    '''
    keysall={'instrumentName': 'instrumentName',
             'distance': '10 Mpc',
             'inclination': '0 deg',
             'azimuth': '0 deg',
             'roll': '0 deg',
             'fieldOfViewX': '50 kpc',
             'numPixelsX': '250',
             'centerX': '0 kpc',
             'fieldOfViewY': '50 kpc',
             'numPixelsY': '250',
             'centerY': '0 kpc',
             'recordComponents': 'false',
             'numScatteringLevels': '0',
             'recordPolarization': 'false',
             'recordStatistics': 'false',
    }
    ele = xdm.Document().createElement('FullInstrument')
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    ele = ski_append(ele, ski['wavelengthGrid'])
    return ele

@Element_ski.ski_element
def HEALPixSkyInstrument(ski):
    '''
    a HEALPix all-sky instrument (for Planck-like observations inside a model)
    instrumentName	String	the name for this instrument
    order	Int	HEALPix order
    radius	Double	the radius of the observer's all-sky sphere
    observerX	Double	the position of the observer, x component
    observerY	Double	the position of the observer, y component
    observerZ	Double	the position of the observer, z component
    crossX	Double	the position of the crosshair, x component
    crossY	Double	the position of the crosshair, y component
    crossZ	Double	the position of the crosshair, z component
    upX	Double	the upwards direction, x component
    upY	Double	the upwards direction, y component
    upZ	Double	the upwards direction, z component
    recordComponents	Bool	record flux components separately
    numScatteringLevels	Int	the number of individually recorded scattering levels
    recordPolarization	Bool	record polarization (Stokes vector elements)
    recordStatistics	Bool	record information for calculating statistical properties
    '''
    keysall={'instrumentName': 'instrumentName',
             'order': '6',
             'radius': '30 kpc',
             'observerX': '10 kpc',
             'observerY': '10 kpc',
             'observerZ': '10 kpc',
             'crossX': '5 kpc',
             'crossY': '5 kpc',
             'crossZ': '5 kpc',
             'upX': '5 kpc',
             'upY': '5 kpc',
             'upZ': '5 kpc',
             'recordComponents': 'false',
             'numScatteringLevels': '0',
             'recordPolarization': 'false',
             'recordStatistics': 'false',}
    ele = xdm.Document().createElement('HEALPixSkyInstrument')
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    ele = ski_append(ele, ski['wavelengthGrid'])
    return ele

@Element_ski.ski_element
def PerspectiveInstrument(ski):
    '''
    a perspective instrument (mostly for making movies)
    instrumentName	String	the name for this instrument
    numPixelsX	Int	the number of viewport pixels in the horizontal direction
    numPixelsY	Int	the number of viewport pixels in the vertical direction
    width	Double	the width of the viewport
    viewX	Double	the position of the viewport origin, x component
    viewY	Double	the position of the viewport origin, y component
    viewZ	Double	the position of the viewport origin, z component
    crossX	Double	the position of the crosshair, x component
    crossY	Double	the position of the crosshair, y component
    crossZ	Double	the position of the crosshair, z component
    upX	Double	the upwards direction, x component
    upY	Double	the upwards direction, y component
    upZ	Double	the upwards direction, z component
    focal	Double	the distance from the eye to the viewport origin
    recordComponents	Bool	record flux components separately
    numScatteringLevels	Int	the number of individually recorded scattering levels
    recordPolarization	Bool	record polarization (Stokes vector elements)
    recordStatistics	Bool	record information for calculating statistical properties
    '''
    keysall={'instrumentName': 'instrumentName',
             'numPixelsX': '250',
             'numPixelsY': '250',
             'width': '30 kpc',
             'viewX': '10 kpc',
             'viewY': '10 kpc',
             'viewZ': '10 kpc',
             'crossX': '5 kpc',
             'crossY': '5 kpc',
             'crossZ': '5 kpc',
             'upX': '5 kpc',
             'upY': '5 kpc',
             'upZ': '5 kpc',
             'focal': '10 kpc',
             'recordComponents': 'false',
             'numScatteringLevels': '0',
             'recordPolarization': 'false',
             'recordStatistics': 'false',}
    ele = xdm.Document().createElement('PerspectiveInstrument')
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    ele = ski_append(ele, ski['wavelengthGrid'])
    return ele

@Element_ski.ski_element
def SEDInstrument(ski):
    '''
    a distant instrument that outputs the spatially integrated flux density as an SED
    instrumentName	String	the name for this instrument
    distance	Double	the distance to the system
    inclination	Double	the inclination angle θ of the detector
    azimuth	Double	the azimuth angle φ of the detector
    roll	Double	the roll angle ω of the detector
    radius	Double	the radius of the circular aperture, or zero for no aperture
    recordComponents	Bool	record flux components separately
    numScatteringLevels	Int	the number of individually recorded scattering levels
    recordPolarization	Bool	record polarization (Stokes vector elements)
    recordStatistics	Bool	record information for calculating statistical properties
    '''
    keysall={'instrumentName': 'instrumentName',
             'distance': '10 Mpc',
             'inclination': '0 deg',
             'azimuth': '0 deg',
             'roll': '0 deg',
             'radius': '50 kpc',
             'recordComponents': 'false',
             'numScatteringLevels': '0',
             'recordPolarization': 'false',
             'recordStatistics': 'false',
    }
    ele = xdm.Document().createElement('SEDInstrument')
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    ele = ski_append(ele, ski['wavelengthGrid'])
    return ele

@Element_ski.ski_element
def probeSystem(ski):
    '''
    the probe system
    '''
    ele = xdm.Document().createElement('probeSystem')
    ele.setAttribute('type', 'ProbeSystem')
    ele = ski_append(ele, ski['ProbeSystem'])
    return ele

@Element_ski.ski_element
def ProbeSystem(ski):
    '''
    a probe system
    '''
    ele = xdm.Document().createElement('ProbeSystem')
    if ski.get_default('ProbeSystem'):
        ele = ski_append(ele, ski['probes'])
        return ele
    else:
        return ele

@Element_ski.ski_element
def probes(ski):
    '''
    the probes
    containing: ConvergenceCutsProbe ConvergenceInfoProbe CustomStateProbe DensityProbe DustAbsorptionPerCellProbe DustEmissionWavelengthGridProbe DustEmissivityProbe
    DustGrainPopulationsProbe DustGrainSizeDistributionProbe ImportedMediumDensityProbe ImportedMediumMetallicityProbe ImportedMediumTemperatureProbe ImportedMediumVelocityProbe
    ImportedSourceAgeProbe  ImportedSourceDensityProbe ImportedSourceLuminosityProbe ImportedSourceMetallicityProbe ImportedSourceVelocityProbe InstrumentWavelengthGridProbe
    LaunchedPacketsProbe LuminosityProbe MagneticFieldProbe MetallicityProbe OpacityProbe OpticalMaterialPropertiesProbe RadiationFieldProbe RadiationFieldWavelengthGridProbe
    SecondaryDustLuminosityProbe SecondaryLineLuminosityProbe SpatialCellPropertiesProbe SpatialGridPlotProbe SpatialGridSourceDensityProbe TemperatureProbe
    TreeSpatialGridTopologyProbe VelocityProbe
    '''
    ele = xdm.Document().createElement('probes')
    ele.setAttribute('type', 'Probe')
    subkey = ski.get_default('probes')    #TODO   ItemList
    ele = ski_append(ele, ski[subkey])
    return ele

@Element_ski.ski_element
def ConvergenceCutsProbe(ski):
    '''
    convergence: cuts of the medium density along the coordinate planes
    probeName	String	the name for this probe
    probeAfter	Enum	perform the probe after
    --> Setup	after setup
    --> Run	after the complete simulation run
    '''
    keysall = {'probeName': 'ConvergenceCutsProbe',
               'probeAfter': 'Setup',}
    ele = xdm.Document().createElement('ConvergenceCutsProbe')
    for i in keysall:
        ele.setAttribute(i, keysall[i])
    return ele

@Element_ski.ski_element
def ConvergenceInfoProbe(ski):
    '''
    convergence: information on the spatial grid
    probeName	String	the name for this probe
    wavelength Double the wavelength at which to determine the optical depth
    probeAfter	Enum	perform the probe after
    --> Setup	after setup
    --> Run	after the complete simulation run
    '''
    keysall = {'probeName': 'ConvergenceInfoProbe',
               'wavelength': '0.55 micron',
               'probeAfter': 'Setup',}
    ele = xdm.Document().createElement('ConvergenceInfoProbe')
    for i in keysall:
        ele.setAttribute(i, keysall[i])
    return ele

@Element_ski.ski_element
def CustomStateProbe(ski):
    '''
    internal spatial grid: custom medium state quantities
    probeName	String	the name for this probe
    indices	String	zero-based indices or index ranges of quantities to probe, or empty for all
    probeAfter	Enum	perform the probe after
    --> Setup	after setup
    --> Run	after the complete simulation run
    --> Primary	after each iteration over primary emission
    --> Secondary	after each iteration over secondary emission
    '''
    keysall = {'probeName': 'CustomStateProbe',
               'indices': '',
               'probeAfter': 'Setup',}
    ele = xdm.Document().createElement('CustomStateProbe')
    for i in keysall:
        ele.setAttribute(i, keysall[i])
    return ele

@Element_ski.ski_element
def DensityProbe(ski):
    '''
    internal spatial grid: density of the medium
    probeName	String	the name for this probe
    aggregation	Enum	how to aggregate the density
    --> Fragment	per fragment (dust grain material type and/or size bin)
    --> Component	per medium component
    --> Type	per medium type (dust, electrons, gas)
    probeAfter	Enum	perform the probe after
    --> Setup	after setup
    --> Run	after the complete simulation run
    --> Primary	after each iteration over primary emission
    --> Secondary	after each iteration over secondary emission
    '''
    keysall = {'probeName': 'DensityProbe',
               'aggregation': 'Fragment',
               'probeAfter': 'Setup',}
    ele = xdm.Document().createElement('DensityProbe')
    for i in keysall:
        ele.setAttribute(i, keysall[i])
    ele = ski_append(ele, ski['form'])
    return ele


@Element_ski.ski_element
def form(ski):
    '''
    the form describing how this quantity should be probed
    '''
    ele = xdm.Document().createElement('form')
    ele.setAttribute('type', 'Form')
    ele = ski_append(ele, ski['Form'])
    return ele

@Element_ski.ski_element
def Form(ski):
    '''
    a probe form
    containing: AllSkyProjectionForm AtPositionsForm DefaultCutsForm LinearCutForm MeridionalCutForm ParallelProjectionForm
    PerCellForm PlanarCutsForm
    '''
    ele = xdm.Document().createElement('Form')
    subkey = ski.xxx    #TODO
    ele = ski_append(ele, ski[subkey])
    return ele


@Element_ski.ski_element
def AllSkyProjectionForm(ski):
    '''
    all-sky projection at a position inside the model
    numPixelsY	Int	the number of image pixels in the vertical (shortest) direction
    observerX	Double	the position of the observer, x component
    observerY	Double	the position of the observer, y component
    observerZ	Double	the position of the observer, z component
    crossX	Double	the position of the crosshair, x component
    crossY	Double	the position of the crosshair, y component
    crossZ	Double	the position of the crosshair, z component
    upX	Double	the upwards direction, x component
    upY	Double	the upwards direction, y component
    upZ	Double	the upwards direction, z component
    '''
    keysall={'numPixelsY': '250',
             'observerX': '10 kpc',
             'observerY': '10 kpc',
             'observerZ': '10 kpc',
             'crossX': '5 kpc',
             'crossY': '5 kpc',
             'crossZ': '5 kpc',
             'upX': '5 kpc',
             'upY': '5 kpc',
             'upZ': '5 kpc',}
    ele = xdm.Document().createElement('AllSkyProjectionForm')
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    ele = ski_append(ele, ski['projection'])
    return ele

@Element_ski.ski_element
def AtPositionsForm(ski):
    '''
    a text column file with values for each imported position
    filename	String	the name of the file listing the positions
    useColumns	String	a list of names corresponding to columns in the file to be imported
    '''
    keysall={'filename': 'AtPositionsForm',
             'useColumns': ''}
    
    ele = xdm.Document().createElement('AtPositionsForm')
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    return ele

@Element_ski.ski_element
def DefaultCutsForm(ski):
    '''
    default planar cuts along the coordinate planes
    '''
    ele = xdm.Document().createElement('DefaultCutsForm')
    return ele

@Element_ski.ski_element
def LinearCutForm(ski):
    '''
    a text column file with values along a given line segment
    numSamples	Int	the number of samples along the line segment
    startX	Double	the position of the starting point, x component
    startY	Double	the position of the starting point, y component
    startZ	Double	the position of the starting point, z component
    endX	Double	the position of the ending point, x component
    endY	Double	the position of the ending point, y component
    endZ	Double	the position of the ending point, z component
    '''
    keysall={'numSamples': '250',
             'startX': '-10 kpc',
             'startY': '-10 kpc',
             'startZ': '-10 kpc',
             'endX': '10 kpc',
             'endY': '10 kpc',
             'endZ': '10 kpc',}
    
    ele = xdm.Document().createElement('LinearCutForm')
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    return ele

@Element_ski.ski_element
def MeridionalCutForm(ski):
    '''
    a text column file with values along a meridian
    numSamples	Int	the number of samples along the meridian
    radius	Double	the radius of the circle containing the meridian
    azimuth	Double	the azimuth angle φ of the meridian
    '''
    keysall={'numSamples': '250',
             'radius': '10 kpc',
             'azimuth': '0 deg',}
    
    ele = xdm.Document().createElement('MeridionalCutForm')
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    return ele

@Element_ski.ski_element
def ParallelProjectionForm(ski):
    '''
    parallel projection on a distant plane
    inclination	Double	the inclination angle θ of the projection
    azimuth	Double	the azimuth angle φ of the projection
    roll	Double	the roll angle ω of the projection
    fieldOfViewX	Double	the total field of view in the horizontal direction
    numPixelsX	Int	the number of pixels in the horizontal direction
    centerX	Double	the center of the frame in the horizontal direction
    fieldOfViewY	Double	the total field of view in the vertical direction
    numPixelsY	Int	the number of pixels in the vertical direction
    centerY	Double	the center of the frame in the vertical direction
    numSampling	Int	the oversampling rate in each direction
    '''
    keysall={'inclination': '0 deg',
             'azimuth': '0 deg',
             'roll': '0 deg',
             'fieldOfViewX': '50 kpc',
             'numPixelsX': '250',
             'centerX': '0 kpc',
             'fieldOfViewY': '50 kpc',
             'centerY': '0 kpc',
             'numSampling': '250',}
    ele = xdm.Document().createElement('ParallelProjectionForm')
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    return ele

@Element_ski.ski_element
def PerCellForm(ski):
    '''
    a text column file with values for each spatial cell
    '''
    ele = xdm.Document().createElement('PerCellForm')
    return ele

@Element_ski.ski_element
def PlanarCutsForm(ski):
    '''
    configurable planar cuts parallel to the coordinate planes
    minX	Double	the start point of the box in the X direction
    maxX	Double	the end point of the box in the X direction
    minY	Double	the start point of the box in the Y direction
    maxY	Double	the end point of the box in the Y direction
    minZ	Double	the start point of the box in the Z direction
    maxZ	Double	the end point of the box in the Z direction
    positionX	Double	the x position of the cut parallel to the yz plane
    positionY	Double	the y position of the cut parallel to the xz plane
    positionZ	Double	the z position of the cut parallel to the xy plane
    numPixelsX	Int	the number of pixels in the x direction
    numPixelsY	Int	the number of pixels in the y direction
    numPixelsZ	Int	the number of pixels in the z direction
    '''
    keysall={'minX': '-10 kpc',
             'maxX': '10 kpc',
             'minY': '-10 kpc',
             'maxY': '10 kpc',
             'minZ': '-10 kpc',
             'maxZ': '10 kpc',
             'positionX': '0 kpc',
             'positionY': '0 kpc',
             'positionZ': '0 kpc',
             'numPixelsX': '250',
             'numPixelsY': '250',
             'numPixelsZ': '250',}
    ele = xdm.Document().createElement('PlanarCutsForm')
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    return ele

@Element_ski.ski_element
def DustAbsorptionPerCellProbe(ski):
    '''
    specialty: spectral luminosity absorbed by dust for each spatial cell
    probeName	String	the name for this probe
    writeWavelengthGrid	Bool	output a text file with the radiation field wavelength grid
    '''
    keysall={'probeName': 'DustAbsorptionPerCellProbe',
             'writeWavelengthGrid': 'false',}
    
    ele = xdm.Document().createElement('DustAbsorptionPerCellProbe')
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    return ele

@Element_ski.ski_element
def DustEmissionWavelengthGridProbe(ski):
    '''
    wavelength grid: dust emission
    probeName	String	the name for this probe
    '''
    keysall={'probeName': 'DustEmissionWavelengthGridProbe'}
    
    ele = xdm.Document().createElement('DustEmissionWavelengthGridProbe')
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    return ele

@Element_ski.ski_element
def DustEmissivityProbe(ski):
    '''
    properties: emissivity spectrum for each dust mix in a range of standard fields
    probeName	String	the name for this probe
    writeWavelengthGrid	Bool	output a text file with the emission spectrum wavelength grid
    '''
    keysall={'probeName': 'DustEmissivityProbe',
             'writeWavelengthGrid': 'false',}
    
    ele = xdm.Document().createElement('DustEmissivityProbe')
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    return ele

@Element_ski.ski_element
def DustGrainPopulationsProbe(ski):
    '''
    properties: dust grain population mass and size information
    probeName	String	the name for this probe
    '''
    keysall={'probeName': 'DustGrainPopulationsProbe'}
    
    ele = xdm.Document().createElement('DustGrainPopulationsProbe')
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    return ele

@Element_ski.ski_element
def DustGrainSizeDistributionProbe(ski):
    '''
     properties: dust grain size distribution
    probeName	String	the name for this probe
    numSamples	Int	the number of samples in the size distribution table
    '''
    keysall={'probeName': 'DustGrainSizeDistributionProbe',
             'numSamples': '250',}
    
    ele = xdm.Document().createElement('DustGrainSizeDistributionProbe')
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    return ele

@Element_ski.ski_element
def ImportedMediumDensityProbe(ski):
    '''
    imported medium: mass or number density
    probeName	String	the name for this probe
    '''
    keysall={'probeName': 'ImportedMediumDensityProbe'}
    
    ele = xdm.Document().createElement('ImportedMediumDensityProbe')
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    ele = ski_append(ele, ski['form'])
    return ele



#TODO  many probes, but maybe rarely use