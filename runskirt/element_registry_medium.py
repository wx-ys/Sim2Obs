import xml.dom.minidom as xdm

from .skirt_element import Element_ski, ski_append



@Element_ski.ski_element
def mediumSystem(ski):
    '''
    the medium system
    '''
    ele = xdm.Document().createElement('mediumSystem')
    ele.setAttribute('type', 'MediumSystem')
    ele = ski_append(ele, ski['MediumSystem'])
    
    return ele

@Element_ski.ski_element
def MediumSystem(ski):
    '''
    a medium system
    containing: photonPacketOptions lyaOptions dynamicStateOptions radiationFieldOptions secondaryEmissionOptions
        iterationOptions dustEmissionOptions media samplingOptions grid
    '''
    ele = xdm.Document().createElement('MediumSystem')
    containall = ski.get_default('mediumSystem_options').split(',')
    for i in containall:
        ele = ski_append(ele, ski[i])           #TODO
    return ele

@Element_ski.ski_element
def photonPacketOptions(ski):
    '''
    the photon packet options
    '''
    ele = xdm.Document().createElement('photonPacketOptions')
    ele.setAttribute('type', 'PhotonPacketOptions')
    ele = ski_append(ele, ski['PhotonPacketOptions'])
    
    return ele

@Element_ski.ski_element
def PhotonPacketOptions(ski):
    '''
    a set of options related to the photon packet lifecycle
    explicitAbsorption	Bool	use explicit absorption to allow negative absorption (stimulated emission)
    forceScattering 	Bool	use forced scattering to reduce noise
    minWeightReduction	Double	the minimum weight reduction factor before a photon packet is terminated
    minScattEvents  	Int	the minimum number of forced scattering events before a photon packet is terminated
    pathLengthBias  	Double	the fraction of path lengths sampled from a stretched distribution
    '''
    ele = xdm.Document().createElement('PhotonPacketOptions')
    keysall = {'explicitAbsorption': 'false',
               'forceScattering': 'true',
               'minWeightReduction': '1e4',
               'minScattEvents': '0',
               'pathLengthBias': '0.5',
    }
    for i in keysall:
        ele.setAttribute(i, keysall[i])   #TODO
    return ele

@Element_ski.ski_element
def lyaOptions(ski):
    '''
    the Lyman-alpha line transfer options
    '''
    ele = xdm.Document().createElement('lyaOptions')
    ele.setAttribute('type', 'LyaOptions')
    ele = ski_append(ele, ski['LyaOptions'])
    
    return ele


@Element_ski.ski_element
def LyaOptions(ski):
    '''
    a set of options related to Lyman-alpha line transfer
    lyaAccelerationScheme the Lyman-alpha line transfer acceleration scheme,
        (None, Constant, Variable), no acceleration, acceleration scheme with a constant critical value, acceleration scheme depending on local gas temperature and density
    lyaAccelerationStrength	Double	the acceleration strength; higher is faster but less accurate
    includeHubbleFlow   	Bool	include the Doppler shift caused by the expansion of the universe   
    '''
    ele = xdm.Document().createElement('LyaOptions')
    keysall = {'lyaAccelerationScheme': 'Variable',
               'lyaAccelerationStrength': '1',
               'includeHubbleFlow': 'false',
    }
    for i in keysall:
        ele.setAttribute(i, keysall[i])   #TODO
    return ele

@Element_ski.ski_element
def dynamicStateOptions(ski):
    '''
    the dynamic medium state options
    '''
    ele = xdm.Document().createElement('dynamicStateOptions')
    ele.setAttribute('type', 'DynamicStateOptions')
    ele = ski_append(ele, ski['DynamicStateOptions'])
    
    return ele

@Element_ski.ski_element
def DynamicStateOptions(ski):
    '''
    a set of options for dynamically adjusting the medium state
    '''
    ele = xdm.Document().createElement('DynamicStateOptions')
    ele = ski_append(ele, ski['recipes'])
    
    return ele

@Element_ski.ski_element
def recipes(ski):
    '''
    the dynamic medium state recipes
    '''
    ele = xdm.Document().createElement('recipes')
    ele.setAttribute('type', 'DynamicStateRecipe')
    ele = ski_append(ele, ski['ClearDensityRecipe'])
    ele = ski_append(ele, ski['LinearDustDestructionRecipe'])
    return ele


@Element_ski.ski_element
def ClearDensityRecipe(ski):
    '''
    a recipe that clears material in cells above a given radiation field strength
    maxNotConvergedCells	Int	the number of spatial cells allowed to not converge
    fieldStrengthThreshold	Double	the field strength above which material is cleared from a cell
    '''
    ele = xdm.Document().createElement('ClearDensityRecipe')
    keysall = {'maxNotConvergedCells': '0',
               'fieldStrengthThreshold': '1',
    }
    for i in keysall:
        ele.setAttribute(i, keysall[i])   #TODO
    return ele

@Element_ski.ski_element
def LinearDustDestructionRecipe(ski):
    '''
    a dust destruction recipe using a linear temperature dependence
    maxNotConvergedCells	Int	the number of spatial cells allowed to not converge
    densityFractionTolerance	Double	the convergence tolerance on the dynamic density fraction
    minSilicateTemperature	Double	the temperature below which silicate grains are not destroyed
    maxSilicateTemperature	Double	the temperature above which all silicate grains are destroyed
    minGraphiteTemperature	Double	the temperature below which graphite grains are not destroyed
    maxGraphiteTemperature	Double	the temperature above which all graphite grains are destroyed
    '''
    ele = xdm.Document().createElement('LinearDustDestructionRecipe')
    keysall = {'maxNotConvergedCells': '0',
               'densityFractionTolerance': '0.05',
               'minSilicateTemperature': '1200 K',
               'maxSilicateTemperature': '1200 K',
               'minGraphiteTemperature': '2000 K',
               'maxGraphiteTemperature': '2000 K',
    }
    for i in keysall:
        ele.setAttribute(i, keysall[i])   #TODO
    return ele
    
@Element_ski.ski_element
def radiationFieldOptions(ski):
    '''
     the radiation field options
    '''
    ele = xdm.Document().createElement('radiationFieldOptions')
    ele.setAttribute('type', 'RadiationFieldOptions')
    ele = ski_append(ele, ski['RadiationFieldOptions'])
    return ele
@Element_ski.ski_element
def RadiationFieldOptions(ski):
    '''
    a set of options related to the radiation field
    '''
    ele = xdm.Document().createElement('RadiationFieldOptions')
    ele.setAttribute('storeRadiationField', 'false')   #TODO
    ele = ski_append(ele, ski['radiationFieldWLG'])
    return ele

@Element_ski.ski_element
def radiationFieldWLG(ski):
    '''
    the wavelength grid for storing the radiation field
    '''
    ele = xdm.Document().createElement('radiationFieldWLG')
    ele.setAttribute('type', 'DisjointWavelengthGrid')
    subkey = ski.get_default('radiationFieldWLG')
    ele = ski_append(ele, ski[subkey]) #TODO
    return ele

@Element_ski.ski_element
def CompositeWavelengthGrid(ski):
    '''
    a wavelength grid composited from a list of wavelength grids
    '''
    ele = xdm.Document().createElement('CompositeWavelengthGrid')
    ele.setAttribute('log' 'true')  #TODO
    ele = ski_append(ele, ski['wavelengthGrids'])
    return ele

@Element_ski.ski_element
def wavelengthGrids(ski):
    '''
    the wavelength grids to be composited
    '''
    ele = xdm.Document().createElement('wavelengthGrids')
    ele.setAttribute('type', 'DisjointWavelengthGrid')
    grikey = ski.xxx
    ele = ski_append(ele, ski[grikey]) #TODO
    return ele

@Element_ski.ski_element
def FileBorderWavelengthGrid(ski):
    '''
    a wavelength grid loaded from a text file listing bin borders
    filename    the name of the file with the wavelength bin borders
    characteristic  (Linear, Logarithmic, Specified)
    '''
    ele = xdm.Document().createElement('FileBorderWavelengthGrid')
    ele.setAttribute('filename', 'FileBorderWavelengthGrid.txt')
    ele.setAttribute('characteristic', 'Logarithmic')           #TODO
    return ele

@Element_ski.ski_element
def FileWavelengthGrid(ski):
    '''
    a wavelength grid loaded from a text file
    filename	String	the name of the file with the characteristic wavelengths
    relativeHalfWidth	Double	the relative half width for discrete bins, or zero for a consecutive range
    log	docs	use logarithmic scale
    '''
    keysall = {'filename':  'FileWavelengthGrid.txt',
               'relativeHalfWidth': '0',
               'log': 'true'}
    ele = xdm.Document().createElement('FileBorderWavelengthGrid')
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    return ele

@Element_ski.ski_element
def LinBorderWavelengthGrid(ski):
    '''
    a linear wavelength grid with given outer borders
    minWavelength	Double	the shortest wavelength
    maxWavelength	Double	the longest wavelength
    numWavelengthBins	Int	the number of wavelength bins
    '''
    keysall = {'minWavelength':  '1 pm',
               'maxWavelength': '1 m',
               'numWavelengthBins': '25'}
    ele = xdm.Document().createElement('LinBorderWavelengthGrid')
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    return ele

@Element_ski.ski_element
def LinWavelengthGrid(ski):
    '''
    a linear wavelength grid
    minWavelength	Double	the shortest wavelength
    maxWavelength	Double	the longest wavelength
    numWavelengths	Int	the number of wavelength grid points
    '''
    keysall = {'minWavelength':  '1 pm',
               'maxWavelength': '1 m',
               'numWavelengths': '25'}
    ele = xdm.Document().createElement('LinWavelengthGrid')
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    return ele

@Element_ski.ski_element
def ListBorderWavelengthGrid(ski):
    '''
    a wavelength grid configured as a list of bin borders
    wavelengths DoubleList	the wavelength bin borders
    characteristic (Linear, Logarithmic, Specified)
    '''
    keysall = {'wavelengths':  '',
               'characteristic': 'Logarithmic',}
    ele = xdm.Document().createElement('ListBorderWavelengthGrid')
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    return ele

@Element_ski.ski_element
def ListWavelengthGrid(ski):
    '''
    a wavelength grid configured as a list
    wavelengths	DoubleList	the characteristic wavelength for each bin
    relativeHalfWidth	Double	the relative half width for discrete bins, or zero for a consecutive range
    log	Bool	use logarithmic scale
    '''
    keysall = {'wavelengths':  '',
               'relativeHalfWidth': '0',
               'log': 'true'}
    ele = xdm.Document().createElement('ListWavelengthGrid')
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    return ele

@Element_ski.ski_element
def LogBorderWavelengthGrid(ski):
    '''
    a logarithmic wavelength grid with given outer borders
    minWavelength	Double	the shortest wavelength
    maxWavelength	Double	the longest wavelength
    numWavelengthBins	Int	the number of wavelength bins
    '''
    keysall = {'minWavelength':  '1 pm',
               'maxWavelength': '1 m',
               'numWavelengthBins': '25'}
    ele = xdm.Document().createElement('LogBorderWavelengthGrid')
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    return ele

@Element_ski.ski_element
def LogWavelengthGrid(ski):
    '''
    a logarithmic wavelength grid
    minWavelength	Double	the shortest wavelength
    maxWavelength	Double	the longest wavelength
    numWavelengths	Int	the number of wavelength grid points
    '''
    keysall = {'minWavelength':  '1 pm',
               'maxWavelength': '1 m',
               'numWavelengths': '25'}
    ele = xdm.Document().createElement('LogWavelengthGrid')
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    return ele

@Element_ski.ski_element
def NestedLogWavelengthGrid(ski):
    '''
    a nested logarithmic wavelength grid
    minWavelengthBaseGrid	Double	the shortest wavelength of the low-resolution grid
    maxWavelengthBaseGrid	Double	the longest wavelength of the low-resolution grid
    numWavelengthsBaseGrid	Int	the number of wavelength grid points in the low-resolution grid
    minWavelengthSubGrid	Double	the shortest wavelength of the high-resolution subgrid
    maxWavelengthSubGrid	Double	the longest wavelength of the high-resolution subgrid
    numWavelengthsSubGrid	Int	the number of wavelength grid points in the high-resolution subgrid
    '''
    keysall = {'minWavelengthBaseGrid':  '1 pm',
               'maxWavelengthBaseGrid': '1 m',
               'numWavelengthsBaseGrid': '25',
               'minWavelengthSubGrid':  '1 pm',
               'maxWavelengthSubGrid': '1 m',
               'numWavelengthsSubGrid': '25'}
    ele = xdm.Document().createElement('NestedLogWavelengthGrid')
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    return ele

@Element_ski.ski_element
def ResolutionBorderWavelengthGrid(ski):
    '''
    a logarithmic wavelength grid with given spectral resolution and outer borders
    minWavelength	docs	Double	the shortest wavelength
    maxWavelength	docs	Double	the longest wavelength
    resolution	docs	Double	the spectral resolution R of the grid
    '''
    keysall = {'minWavelength':  '1 pm',
               'maxWavelength': '1 m',
               'resolution': '10',}
    ele = xdm.Document().createElement('ResolutionBorderWavelengthGrid')
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    return ele

@Element_ski.ski_element
def ResolutionWavelengthGrid(ski):
    '''
    a logarithmic wavelength grid with given spectral resolution
    minWavelength	docs	Double	the shortest wavelength
    maxWavelength	docs	Double	the longest wavelength
    resolution	docs	Double	the spectral resolution R of the grid
    '''
    keysall = {'minWavelength':  '1 pm',
               'maxWavelength': '1 m',
               'resolution': '10',}
    ele = xdm.Document().createElement('ResolutionWavelengthGrid')
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    return ele

@Element_ski.ski_element
def secondaryEmissionOptions(ski):
    '''
    the secondary emission options
    '''
    ele = xdm.Document().createElement('secondaryEmissionOptions')
    ele.setAttribute('type', 'SecondaryEmissionOptions')
    ele = ski_append(ele, ski['SecondaryEmissionOptions'])
    return ele

@Element_ski.ski_element
def SecondaryEmissionOptions(ski):
    '''
    a set of options related to secondary emission    
    storeEmissionRadiationField	Bool	store the radiation field during emission so that it can be probed for output
    secondaryPacketsMultiplier	Double	the multiplier on the number of photon packets launched for secondary emission
    spatialBias	Double	the fraction of secondary photon packets distributed uniformly across spatial cells
    sourceBias	Double	the fraction of photon packets distributed uniformly across secondary sources
    '''
    keysall = {'storeEmissionRadiationField': 'false',
               'secondaryPacketsMultiplier': '1',
               'spatialBias': '0.5',
               'sourceBias': '0.5',
        
    }
    ele = xdm.Document().createElement('secondaryEmissionOptions')
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    return ele

@Element_ski.ski_element
def iterationOptions(ski):
    '''
    the primary and/or secondary emission iteration options
    '''
    ele = xdm.Document().createElement('iterationOptions')
    ele.setAttribute('type', 'IterationOptions')
    ele = ski_append(ele, ski['IterationOptions'])
    return ele
@Element_ski.ski_element
def IterationOptions(ski):
    '''
    a set of options for configuring iterations during primary and/or secondary emission
    minPrimaryIterations	Int	the minimum number of iterations during primary emission
    maxPrimaryIterations	Int	the maximum number of iterations during primary emission
    minSecondaryIterations	Int	the minimum number of iterations during secondary emission
    maxSecondaryIterations	Int	the maximum number of iterations during secondary emission
    includePrimaryEmission	Bool	include primary emission in the secondary emission iterations
    primaryIterationPacketsMultiplier	Double	the multiplier on the number of photon packets launched for each primary emission iteration
    secondaryIterationPacketsMultiplier	Double	the multiplier on the number of photon packets launched for each secondary emission iteration
    '''
    keysall = {'minPrimaryIterations': '1',
               'maxPrimaryIterations': '10',
               'minSecondaryIterations': '1',
               'maxSecondaryIterations': '10',
               'includePrimaryEmission': 'false',
               'primaryIterationPacketsMultiplier': '1',
               'secondaryIterationPacketsMultiplier': '1',
    }
    ele = xdm.Document().createElement('IterationOptions')
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    return ele

@Element_ski.ski_element
def dustEmissionOptions(ski):
    '''
    the dust emission options
    '''
    ele = xdm.Document().createElement('dustEmissionOptions')
    ele.setAttribute('type', 'DustEmissionOptions')
    ele = ski_append(ele, ski['DustEmissionOptions'])
    return ele
    
@Element_ski.ski_element
def DustEmissionOptions(ski):
    '''
    a set of options related to thermal emission from dust
    dustEmissionType    the method used for dust emission calculations (Equilibrium, Stochastic)
    includeHeatingByCMB	Bool	add the cosmic microwave background (CMB) as a dust heating source term
    maxFractionOfPrimary	Double	convergence is reached when the total absorbed dust luminosity is less than this fraction of the total absorbed primary luminosity
    maxFractionOfPrevious	Double	convergence is reached when the total absorbed dust luminosity has changed by less than this fraction compared to the previous iteration
    sourceWeight	Double	the weight of dust emission for the number of photon packets launched
    wavelengthBias	Double	the fraction of secondary photon packet wavelengths sampled from a bias distribution
    '''
    keysall = {'dustEmissionType': 'Equilibrium',
               'includeHeatingByCMB': 'false',
               'maxFractionOfPrimary': '0.01',
               'maxFractionOfPrevious': '0.03',
               'sourceWeight': '1',
               'wavelengthBias': '0.5',
    }
    ele = xdm.Document().createElement('DustEmissionOptions')
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    return ele

@Element_ski.ski_element
def cellLibrary(ski):
    '''
    the spatial cell grouping scheme for calculating dust emission
    containing AllCellsLibrary, FieldStrengthCellLibrary, TemperatureWavelengthCellLibrary
    '''
    ele = xdm.Document().createElement('cellLibrary')
    ele.setAttribute('type', 'SpatialCellLibrary')
    subkey = ski.get_default('cellLibrary') #TODO
    ele = ski_append(ele,ski[subkey])
    return ele
@Element_ski.ski_element
def AllCellsLibrary(ski):
    '''
    a library scheme that has a separate entry for every spatial cell
    '''
    ele = xdm.Document().createElement('AllCellsLibrary')
    return ele

@Element_ski.ski_element
def FieldStrengthCellLibrary(ski):
    '''
    a library scheme for grouping spatial cells based on radiation field strength
    numFieldStrengths int the number of field strength bins
    '''
    ele = xdm.Document().createElement('FieldStrengthCellLibrary')
    ele.setAttribute('elnumFieldStrengthse', '1000')        #TODO
    return ele

@Element_ski.ski_element
def TemperatureWavelengthCellLibrary(ski):
    '''
    a library scheme for grouping spatial cells based on indicative temperature and wavelength
    numTemperatures	Int	the number of temperature bins
    numWavelengths	Int	the number of wavelength bins
    '''
    ele = xdm.Document().createElement('TemperatureWavelengthCellLibrary')
    ele.setAttribute('numTemperatures', '40')
    ele.setAttribute('numWavelengths', '25')
    return ele

@Element_ski.ski_element
def dustEmissionWLG(ski):
    '''
    the wavelength grid for calculating the dust emission spectrum
    '''
    ele = xdm.Document().createElement('dustEmissionWLG')
    ele.setAttribute('type', 'DisjointWavelengthGrid')
    subkey = ski.get_default('dustEmissionWLG')
    ele = ski_append(ele,ski[subkey])
    return ele

@Element_ski.ski_element
def media(ski):
    '''
    the transfer media
    AdaptiveMeshMedium, CellMedium, GeometricMedium, ParticleMedium, VoronoiMeshMedium
    '''
    ele = xdm.Document().createElement('media')
    ele.setAttribute('type', 'Medium')
    
    for i in range(100):
        subkey = ski.get_default('mediums')    
        if subkey:    
            ele = ski_append(ele, ski[subkey])
        else:
            break
    return ele

@Element_ski.ski_element
def AdaptiveMeshMedium(ski):
    '''
    a transfer medium imported from data represented on an adaptive mesh (AMR grid)
    filename	String	the name of the file to be imported
    minX	Double	the start point of the domain in the X direction
    maxX	Double	the end point of the domain in the X direction
    minY	Double	the start point of the domain in the Y direction
    maxY	Double	the end point of the domain in the Y direction
    minZ	Double	the start point of the domain in the Z direction
    maxZ	Double	the end point of the domain in the Z direction
    massType    the type of mass quantity to be imported
        (MassDensity, Mass, MassDensityAndMass, NumberDensity, Number, NumberDensityAndNumber)
    massFraction	Double	the fraction of the mass to be included (or one to include all)
    importMetallicity	Bool	import a metallicity column
    importTemperature	Bool	import a temperature column
    maxTemperature	Double	the maximum temperature for included mass (or zero to include all)
    importVelocity	Bool	import velocity components (3 columns)
    importMagneticField	Bool	import magnetic field components (3 columns)
    importVariableMixParams	Bool	import parameter(s) to select a spatially varying material mix
    useColumns	String	a list of names corresponding to columns in the file to be imported
    '''
    keysall = {'filename':  'AdaptiveMeshMedium.txt',
               'minX':  '-10 kpc',
               'maxX':  '10 kpc',
               'minY':  '-10 kpc',
               'maxY':  '10 kpc',
               'minZ':  '-10 kpc',
               'maxZ':  '10 kpc',
               'massType':  'Mass',
               'massFraction':  '1',
               'importMetallicity':  'false',
               'importTemperature':  'false',
               'maxTemperature':  '0 K',
               'importVelocity':  'false',
               'importMagneticField':  'false',
               'importVariableMixParams':  'false',
               'useColumns':  '',
        
    }
    ele = xdm.Document().createElement('AdaptiveMeshMedium')
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    ele = ski_append(ele, ski['materialMix'])
    ele = ski_append(ele, ski['materialMixFamily'])
    
    return ele

@Element_ski.ski_element
def materialMix(ski):
    '''
    the material type and properties throughout the medium
    containing: ConfigurableDustMix DraineLiDustMix ElectronMix  FragmentDustMixDecorator LyaNeutralHydrogenGasMix MRNDustMix MeanFileDustMix MeanInterstellarDustMix
        MeanIvezicBenchmarkDustMix  MeanListDustMix MeanPascucciBenchmarkDustMix MeanPinteBenchmarkDustMix MeanTrustBenchmarkDustMix NonLTELineGasMix SpinFlipAbsorptionMix
        SpinFlipHydrogenGasMix ThemisDustMix TrivialGasMix TrustBenchmarkDustMix WeingartnerDraineDustMix XRayAtomicGasMix ZubkoDustMix
    '''
    ele = xdm.Document().createElement('materialMix')
    ele.setAttribute('type','MaterialMix')
    subkey = ski.get_default('materialMix')        #TODO
    ele = ski_append(ele, ski[subkey])
    return ele

@Element_ski.ski_element
def ConfigurableDustMix(ski):
    '''
    a configurable dust mix with one or more grain populations
    scatteringType  the type of scattering to be implemented
    (HenyeyGreenstein, MaterialPhaseFunction, SphericalPolarization, SpheroidalPolarization)
    '''
    ele = xdm.Document().createElement('ConfigurableDustMix')
    ele.setAttribute('scatteringType','HenyeyGreenstein')       #TODO
    ele = ski_append(ele, ski['populations'])
    return ele


@Element_ski.ski_element
def populations(ski):
    '''
    the grain populations
    
    '''
    ele = xdm.Document().createElement('populations')
    ele.setAttribute('type','GrainPopulation') 
    ele = ski_append(ele , ski['GrainPopulation'])
    return ele

@Element_ski.ski_element
def GrainPopulation(ski):
    '''
    a dust grain population
    numSizes	Int	the number of grain size bins
    normalizationType	Enum	the mechanism for specifying the amount of dust in the population
    --> DustMassPerHydrogenAtom	an absolute dust mass per hydrogen atom
    --> DustMassPerHydrogenMass	a ratio of dust mass per hydrogen mass
    --> FactorOnSizeDistribution	a proportionality factor on the size distribution
    dustMassPerHydrogenAtom	Double	the dust mass per hydrogen atom
    dustMassPerHydrogenMass	Double	the dust mass per hydrogen mass
    factorOnSizeDistribution	Double	the proportionality factor on the size distribution
    '''
    keysall = {'numSizes': '8',
               'normalizationType': 'DustMassPerHydrogenMass',
               'dustMassPerHydrogenAtom': '0',
               'dustMassPerHydrogenMass': '1',
               'factorOnSizeDistribution': '1',}
    ele = xdm.Document().createElement('GrainPopulation')
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    ele = ski_append(ele, ski['composition'])
    ele = ski_append(ele, ski['sizeDistribution'])
    
@Element_ski.ski_element
def composition(ski):
    '''
    the dust grain composition
    containing: BegemannPorousAluminaGrainComposition CrystalEnstatiteGrainComposition CrystalForsteriteGrainComposition DorschnerOlivineGrainComposition
    DraineGraphiteGrainComposition DraineIonizedPAHGrainComposition DraineNeutralPAHGrainComposition DraineSilicateGrainComposition DustEmGrainComposition
    HofmeisterPericlaseGrainComposition  MieSilicateGrainComposition MinSilicateGrainComposition PolarizedGraphiteGrainComposition PolarizedSilicateGrainComposition
    SpheroidalGraphiteGrainComposition SpheroidalSilicateGrainComposition TrustGraphiteGrainComposition TrustNeutralPAHGrainComposition TrustSilicateGrainComposition
    '''
    ele = xdm.Document().createElement('composition')
    ele.setAttribute('type', 'GrainComposition')
    subkey = ski.xxx #TODO
    ele = ski_append(ele, ski[subkey])
    return ele

@Element_ski.ski_element
def BegemannPorousAluminaGrainComposition(ski):
    '''
    a Begemann porous alumina dust grain composition
    '''
    ele = xdm.Document().createElement('BegemannPorousAluminaGrainComposition')
    return ele

@Element_ski.ski_element
def CrystalEnstatiteGrainComposition(ski):
    '''
    a crystalline Enstatite dust grain composition
    '''
    ele = xdm.Document().createElement('CrystalEnstatiteGrainComposition')
    return ele

@Element_ski.ski_element
def CrystalForsteriteGrainComposition(ski):
    '''
    a crystalline Forsterite dust grain composition
    '''
    ele = xdm.Document().createElement('CrystalForsteriteGrainComposition')
    return ele

@Element_ski.ski_element
def DorschnerOlivineGrainComposition(ski):
    '''
    a Dorschner olivine dust grain composition
    '''
    ele = xdm.Document().createElement('DorschnerOlivineGrainComposition')
    return ele

@Element_ski.ski_element
def DraineGraphiteGrainComposition(ski):
    '''
    a Draine graphite dust grain composition
    '''
    ele = xdm.Document().createElement('DraineGraphiteGrainComposition')
    return ele

@Element_ski.ski_element
def DraineIonizedPAHGrainComposition(ski):
    '''
    a Draine ionized PAH dust grain composition
    '''
    ele = xdm.Document().createElement('DraineIonizedPAHGrainComposition')
    return ele

@Element_ski.ski_element
def DraineNeutralPAHGrainComposition(ski):
    '''
    a Draine neutral PAH dust grain composition
    '''
    ele = xdm.Document().createElement('DraineNeutralPAHGrainComposition')
    return ele

@Element_ski.ski_element
def DraineSilicateGrainComposition(ski):
    '''
    a Draine silicate dust grain composition
    '''
    ele = xdm.Document().createElement('DraineSilicateGrainComposition')
    return ele

@Element_ski.ski_element
def DustEmGrainComposition(ski):
    '''
    a dust grain composition based on DustEM data
    grainType	Enum	the DustEM grain type
    --> aSil	astronomical silicate grains (Draine & Li 2007)
    --> Gra	graphite grains (Draine & Li 2001; Li & Draine 2001)
    --> PAH0DL07	neutral PAHs (Draine & Li 2007)
    --> PAH1DL07	ionized PAHs (Draine & Li 2007)
    --> PAH0MC10	neutral PAHs (Compiègne et al. 2011)
    --> PAH1MC10	ionized PAHs (Compiègne et al. 2011)
    --> CM20	amorphous hydro-carbon grains (Jones et al. 2013)
    --> aOlM5	amorphous olivine/forsterite grains (Koehler et al. 2014)
    --> aPyM5	amorphous pyroxene/enstatite grains (Koehler et al. 2014)
    bulkMassDensity	Double	the bulk mass density for this grain material
    '''
    keysall = {'grainType': 'aSil',
               'bulkMassDensity': '3500 kg/m3'
    }
    ele = xdm.Document().createElement('DustEmGrainComposition')
    for i in keysall:
        ele.setAttribute(i, keysall[i])     #TODO
    return ele

@Element_ski.ski_element
def HofmeisterPericlaseGrainComposition(ski):
    '''
     a Hofmeister periclase dust grain composition
    '''
    ele = xdm.Document().createElement('HofmeisterPericlaseGrainComposition')
    return ele

@Element_ski.ski_element
def MieSilicateGrainComposition(ski):
    '''
     a MieX-based Draine silicate dust grain composition
    '''
    ele = xdm.Document().createElement('MieSilicateGrainComposition')
    return ele

@Element_ski.ski_element
def MinSilicateGrainComposition(ski):
    '''
     a Min 2007 amorphous silicate dust grain composition
    '''
    ele = xdm.Document().createElement('MinSilicateGrainComposition')
    return ele

@Element_ski.ski_element
def PolarizedGraphiteGrainComposition(ski):
    '''
     a graphite dust grain composition with support for polarization
    '''
    ele = xdm.Document().createElement('PolarizedGraphiteGrainComposition')
    return ele

@Element_ski.ski_element
def PolarizedSilicateGrainComposition(ski):
    '''
     a silicate dust grain composition with support for polarization
    '''
    ele = xdm.Document().createElement('PolarizedSilicateGrainComposition')
    return ele

@Element_ski.ski_element
def SpheroidalGraphiteGrainComposition(ski):
    '''
    a spheroidal graphite dust grain composition with support for polarization
    tableType	Enum	the type of emission tables to use
    --> Builtin	builtin resources
    --> OneTable	single custom table
    --> TwoTables	two custom tables with interpolation
    emissionTable	String	the name of the file tabulating properties for polarized emission by arbitrarily aligned spheroidal grains
    alignedEmissionTable	String	the name of the file tabulating properties for polarized emission by perfectly aligned spheroidal grains
    nonAlignedEmissionTable	String	the name of the file tabulating properties for polarized emission by non-aligned spheroidal grains
    alignmentFraction	Double	the alignment fraction of the spheroidal grains with the local magnetic field

    '''
    keysall = {'tableType': 'Builtin',
               'emissionTable': 'emissionTable.txt',
               'alignedEmissionTable': 'alignedEmissionTable.txt',
               'nonAlignedEmissionTable': 'nonAlignedEmissionTable.txt',
               'alignmentFraction': '1.',
    }
    ele = xdm.Document().createElement('SpheroidalGraphiteGrainComposition')
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    return ele

@Element_ski.ski_element
def SpheroidalSilicateGrainComposition(ski):
    '''
     a spheroidal silicate dust grain composition with support for polarization
     tableType	Enum	the type of emission tables to use
    --> Builtin	builtin resources
    --> OneTable	single custom table
    --> TwoTables	two custom tables with interpolation
    emissionTable	String	the name of the file tabulating properties for polarized emission by arbitrarily aligned spheroidal grains
    alignedEmissionTable	String	the name of the file tabulating properties for polarized emission by perfectly aligned spheroidal grains
    nonAlignedEmissionTable	String	the name of the file tabulating properties for polarized emission by non-aligned spheroidal grains
    alignmentFraction	Double	the alignment fraction of the spheroidal grains with the local magnetic field
    '''
    keysall = {'tableType': 'Builtin',
               'emissionTable': 'emissionTable.txt',
               'alignedEmissionTable': 'alignedEmissionTable.txt',
               'nonAlignedEmissionTable': 'nonAlignedEmissionTable.txt',
               'alignmentFraction': '1.',}
    ele = xdm.Document().createElement('SpheroidalSilicateGrainComposition')
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    return ele

@Element_ski.ski_element
def TrustGraphiteGrainComposition(ski):
    '''
     a TRUST benchmark graphite dust grain composition
    '''
    ele = xdm.Document().createElement('TrustGraphiteGrainComposition')
    return ele

@Element_ski.ski_element
def TrustNeutralPAHGrainComposition(ski):
    '''
    a TRUST benchmark neutral PAH dust grain composition
    '''
    ele = xdm.Document().createElement('TrustNeutralPAHGrainComposition')
    return ele

@Element_ski.ski_element
def TrustSilicateGrainComposition(ski):
    '''
    a TRUST benchmark silicate dust grain composition
    '''
    ele = xdm.Document().createElement('TrustSilicateGrainComposition')
    return ele

@Element_ski.ski_element
def sizeDistribution(ski):
    '''
    the dust grain size distribution
    containing:  FileGrainSizeDistribution HirashitaLogNormalGrainSizeDistribution ListGrainSizeDistribution LogNormalGrainSizeDistribution ModifiedLogNormalGrainSizeDistribution
        ModifiedPowerLawGrainSizeDistribution PowerLawGrainSizeDistribution SingleGrainSizeDistribution ZubkoGraphiteGrainSizeDistribution ZubkoPAHGrainSizeDistribution
        ZubkoSilicateGrainSizeDistribution
    '''
    ele = xdm.Document().createElement('sizeDistribution')
    ele.setAttribute('type', 'GrainSizeDistribution')
    subkey =ski.xx #TODO
    ele = ski_append(ele, ski[subkey])
    return ele

@Element_ski.ski_element
def FileGrainSizeDistribution(ski):
    '''
    a dust grain size distribution loaded from a text file
    filename	String	the name of the file with the dust grain size distribution
    '''
    ele = xdm.Document().createElement('FileGrainSizeDistribution')
    ele.setAttribute('filename', 'FileGrainSizeDistribution.txt' ) #TODO
    return ele

@Element_ski.ski_element
def HirashitaLogNormalGrainSizeDistribution(ski):
    '''
    a Hirashita (2015) log-normal dust grain size distribution
    minSize	Double	the minimum grain size for this distribution
    maxSize	Double	the maximum grain size for this distribution
    centroid	Double	the centroid a0 of the log-normal law
    width	Double	the width σ of the log-normal law
    '''
    keysall = {'minSize': '1 mm',
               'maxSize': '1 mm',
               'centroid': '1 nm',
               'width': '0.4'}
    ele = xdm.Document().createElement('HirashitaLogNormalGrainSizeDistribution')
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    return ele

@Element_ski.ski_element
def ListGrainSizeDistribution(ski):
    '''
    a dust grain size distribution specified inside the configuration file
    sizes	DoubleList	the grain sizes at which to specify the size distribution
    sizeDistributionValues	DoubleList	the size distribution values at each of the given grain sizes
    '''
    keysall = {'sizes': '',
               'sizeDistributionValues': '',}
    ele = xdm.Document().createElement('ListGrainSizeDistribution')
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    return ele

@Element_ski.ski_element
def LogNormalGrainSizeDistribution(ski):
    '''
    a log-normal dust grain size distribution
    minSize	Double	the minimum grain size for this distribution
    maxSize	Double	the maximum grain size for this distribution
    centroid	Double	the centroid a0 of the log-normal law
    width	Double	the width σ of the log-normal law
    '''
    keysall = {'minSize': '1 mm',
               'maxSize': '1 mm',
               'centroid': '1 nm',
               'width': '0.4'}
    ele = xdm.Document().createElement('LogNormalGrainSizeDistribution')
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    return ele

@Element_ski.ski_element
def ModifiedLogNormalGrainSizeDistribution(ski):
    '''
    a modified log-normal dust grain size distribution
    minSize	Double	the minimum grain size for this distribution
    maxSize	Double	the maximum grain size for this distribution
    centroid	Double	the centroid a0 of the log-normal law
    width	Double	the width σ of the log-normal law
    firstMixingParameter	Double	the first mixing parameter y0
    secondMixingParameter	Double	the second mixing parameter y1
    '''
    keysall = {'minSize': '1 mm',
               'maxSize': '1 mm',
               'centroid': '1 nm',
               'width': '0.4',
               'firstMixingParameter': '1',
               'secondMixingParameter': '1',}
    ele = xdm.Document().createElement('ModifiedLogNormalGrainSizeDistribution')
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    return ele

@Element_ski.ski_element
def ModifiedPowerLawGrainSizeDistribution(ski):
    '''
    minSize	Double	the minimum grain size for this distribution
    maxSize	Double	the maximum grain size for this distribution
    powerLawIndex	Double	the index α of the power law
    turnOffPoint	Double	the turn-off point a_t in the exponential decay term
    scaleExponentialDecay	Double	the scale a_c in the exponential decay term
    exponentExponentialDecay	Double	the exponent γ in the exponential decay term
    scaleCurvature	Double	the scale a_u in the curvature term
    strengthCurvature	Double	the strength ζ in the curvature term
    exponentCurvature	Double	the exponent η in the curvature term
    '''
    keysall = {'minSize': '1 mm',
               'maxSize': '1 mm',
               'powerLawIndex': '-3.5',
               'turnOffPoint': '0.1 micron',
               'scaleExponentialDecay': '0.1 micron',
               'exponentExponentialDecay': '3',
               'scaleCurvature': '0.1 micron',
               'strengthCurvature': '0.3',
               'exponentCurvature': '1',}
    ele = xdm.Document().createElement('ModifiedPowerLawGrainSizeDistribution')
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    return ele

@Element_ski.ski_element
def PowerLawGrainSizeDistribution(ski):
    '''
    a power-law dust grain size distribution
    minSize	Double	the minimum grain size for this distribution
    maxSize	Double	the maximum grain size for this distribution
    exponent	Double	the (absolute value of the) exponent in the power-law distribution function

    '''
    keysall = {'minSize': '1 mm',
               'maxSize': '1 mm',
               'exponent': '3.5',}
    ele = xdm.Document().createElement('PowerLawGrainSizeDistribution')
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    return ele

@Element_ski.ski_element
def SingleGrainSizeDistribution(ski):
    '''
    a single-size dust grain size distribution
    size	Double	the single grain size for this distribution
    '''
    keysall = {'size': '1 mm',}
    ele = xdm.Document().createElement('SingleGrainSizeDistribution')
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    return ele

@Element_ski.ski_element
def ZubkoGraphiteGrainSizeDistribution(ski):
    '''
    a Zubko, Dwek & Arendt size distribution for graphite dust grains
    '''
    ele = xdm.Document().createElement('ZubkoGraphiteGrainSizeDistribution')
    return ele

@Element_ski.ski_element
def ZubkoPAHGrainSizeDistribution(ski):
    '''
    a Zubko, Dwek & Arendt size distribution for PAH molecules
    '''
    ele = xdm.Document().createElement('ZubkoPAHGrainSizeDistribution')
    return ele

@Element_ski.ski_element
def ZubkoSilicateGrainSizeDistribution(ski):
    '''
    a Zubko, Dwek & Arendt size distribution for silicate dust grains
    '''
    ele = xdm.Document().createElement('ZubkoSilicateGrainSizeDistribution')
    return ele

@Element_ski.ski_element
def DraineLiDustMix(ski):
    '''
    a Draine and Li (2007) dust mix
    numSilicateSizes	Int	the number of silicate grain size bins
    numGraphiteSizes	Int	the number of graphite grain size bins (for each of two populations)
    numPAHSizes	Int	the number of neutral and ionized PAH size bins (each)
    '''
    keysall = {'numSilicateSizes': '5',
               'numGraphiteSizes': '5',
               'expnumPAHSizesonent': '5',}
    ele = xdm.Document().createElement('DraineLiDustMix')       
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    
    return ele

@Element_ski.ski_element
def ElectronMix(ski):
    '''
    a population of electrons
    includePolarization	Bool	include support for polarization
    includeThermalDispersion	Bool	include thermal velocity dispersion
    defaultTemperature	Double	the default temperature of the electron population

    '''
    keysall = {'includePolarization': 'false',
               'includeThermalDispersion': 'false',
               'defaultTemperature': '1e4',}
    ele = xdm.Document().createElement('ElectronMix')       
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    return ele
    
@Element_ski.ski_element
def FragmentDustMixDecorator(ski):
    '''
    a dust mix decorator that manages separate densities for fragments of another dust mix
    fragmentSizeBins	Bool	fragment each dust population into its grain size bins
    hasDynamicDensities	Bool	allow the fragment densities to be adjusted dynamically
    initialDensityFraction	Double	the initial value of the dynamic density fraction
    '''
    keysall = {'fragmentSizeBins': 'false',
               'hasDynamicDensities': 'false',
               'initialDensityFraction': '0',}
    ele = xdm.Document().createElement('FragmentDustMixDecorator')       
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    ele = ski_append(ele, ski['dustMix'])
    return ele

@Element_ski.ski_element
def dustMix(ski):
    '''
    a dust mix with one or more grain populations
    containing: ConfigurableDustMix DraineLiDustMix DraineLiDustMix MRNDustMix ThemisDustMix TrustBenchmarkDustMix WeingartnerDraineDustMix ZubkoDustMix
    '''
    ele = xdm.Document().createElement('dustMix')       
    ele.setAttribute('type', 'MultiGrainDustMix')
    subkey =ski.xxx #TODO
    ele = ski_append(ele, ski[subkey])
    return ele

@Element_ski.ski_element
def MRNDustMix(ski):
    '''
    an MRN (1997) dust mix
    numSilicateSizes	Int	the number of silicate grain size bins
    numGraphiteSizes	Int	the number of graphite grain size bins
    '''
    keysall = {'numSilicateSizes': '5',
               'numGraphiteSizes': '5',}
    ele = xdm.Document().createElement('MRNDustMix')       
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    return ele

@Element_ski.ski_element
def ThemisDustMix(ski):
    '''
    a THEMIS (Jones et al. 2017) dust mix
    numSilicateSizes	Int	the number of grain size bins for each of the silicate populations
    numHydrocarbonSizes	Int	the number of grain size bins for each of the hydrocarbon populations
    '''
    keysall = {'numSilicateSizes': '5',
               'numHydrocarbonSizes': '5',}
    ele = xdm.Document().createElement('ThemisDustMix')       
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    return ele

@Element_ski.ski_element
def TrustBenchmarkDustMix(ski):
    '''
    a TRUST benchmark dust mix
    numSilicateSizes	Int	the number of silicate grain size bins
    numGraphiteSizes	Int	the number of graphite grain size bins
    numPAHSizes	Int	the number of PAH size bins
    '''
    keysall = {'numSilicateSizes': '5',
               'numGraphiteSizes': '5',
               'numPAHSizes': '5',}
    ele = xdm.Document().createElement('TrustBenchmarkDustMix')       
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    return ele

@Element_ski.ski_element
def WeingartnerDraineDustMix(ski):
    '''
    a Weingartner and Draine (2001) dust mix
    environment	Enum	the environment determining the dust model
    --> MilkyWay	the Milky Way
    --> LMC	the Large Magellanic Cloud
    --> SMC	the Small Magellanic Cloud
    numSilicateSizes	Int	the number of silicate grain size bins
    numGraphiteSizes	Int	the number of graphite grain size bins
    numPAHSizes	Int	the number of neutral and ionized PAH size bins (each)
    '''
    keysall = {'environment': 'MilkyWay',
                'numSilicateSizes': '5',
               'numGraphiteSizes': '5',
               'numPAHSizes': '5',}
    ele = xdm.Document().createElement('WeingartnerDraineDustMix')       
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    return ele

@Element_ski.ski_element
def ZubkoDustMix(ski):
    '''
    a Zubko et al. (2004) dust mix
    numSilicateSizes	Int	the number of silicate grain size bins
    numGraphiteSizes	Int	the number of graphite grain size bins
    numPAHSizes	Int	the number of neutral and ionized PAH size bins (each)
    '''
    keysall = {'numSilicateSizes': '5',
               'numGraphiteSizes': '5',
               'numPAHSizes': '5',}
    ele = xdm.Document().createElement('ZubkoDustMix')       
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    return ele

@Element_ski.ski_element
def LyaNeutralHydrogenGasMix(ski):
    '''
    neutral hydrogen for Lyman-alpha line transfer
    defaultTemperature	Double	the default temperature of the neutral hydrogen gas
    includePolarization	Bool	include support for polarization
    '''
    keysall = {'defaultTemperature': '1e4',
               'includePolarization': 'false',}
    ele = xdm.Document().createElement('LyaNeutralHydrogenGasMix')       
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    return ele

@Element_ski.ski_element
def MeanFileDustMix(ski):
    '''
    a dust mix with mean properties loaded from a text file
    filename	String	the name of the file with the optical properties for the dust mix
    '''
    ele = xdm.Document().createElement('MeanFileDustMix')
    ele.setAttribute('filename', 'MeanFileDustMix.txt') #TODO
    return ele

@Element_ski.ski_element
def MeanInterstellarDustMix(ski):
    '''
    a typical interstellar dust mix (mean properties)
    '''
    ele = xdm.Document().createElement('MeanInterstellarDustMix')
    return ele

@Element_ski.ski_element
def MeanIvezicBenchmarkDustMix(ski):
    '''
    an Ivezic 1D benchmark dust mix (mean properties)
    '''
    ele = xdm.Document().createElement('MeanIvezicBenchmarkDustMix')
    return ele

@Element_ski.ski_element
def MeanListDustMix(ski):
    '''
    a dust mix with mean properties specified inside the configuration file
    wavelengths	DoubleList	the wavelengths at which to specify the optical properties
    extinctionCoefficients	DoubleList	the extinction mass coefficients at each of the given wavelengths
    albedos	DoubleList	the scattering albedos at each of the given wavelengths
    asymmetryParameters	DoubleList	the scattering asymmetry parameters at each of the given wavelengths
    '''
    keysall = {'wavelengths': '',
               'extinctionCoefficients': '',
               'albedos': '',
               'asymmetryParameters': '',}
    ele = xdm.Document().createElement('MeanListDustMix')       
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    return ele

@Element_ski.ski_element
def MeanPascucciBenchmarkDustMix(ski):
    '''
    a Pascucci 2D benchmark dust mix (mean properties)
    '''
    ele = xdm.Document().createElement('MeanPascucciBenchmarkDustMix')
    return ele

@Element_ski.ski_element
def MeanPinteBenchmarkDustMix(ski):
    '''
    a Pinte 2D benchmark dust mix (mean properties, optionally with polarization at 1 micron)
    scatteringType	Enum	the type of scattering to be implemented
    --> HenyeyGreenstein	use the Henyey-Greenstein phase function (unpolarized)
    --> MaterialPhaseFunction	use the phase function derived from actual material properties (unpolarized)
    --> SphericalPolarization	support polarization through scattering by spherical grains
    '''
    keysall = {'scatteringType': 'HenyeyGreenstein',}
    ele = xdm.Document().createElement('MeanPinteBenchmarkDustMix')       
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    return ele

@Element_ski.ski_element
def MeanTrustBenchmarkDustMix(ski):
    '''
    a TRUST benchmark dust mix (mean properties, optionally with polarization)
    scatteringType	Enum	the type of scattering to be implemented
    --> HenyeyGreenstein	use the Henyey-Greenstein phase function (unpolarized)
    --> MaterialPhaseFunction	use the phase function derived from actual material properties (unpolarized)
    --> SphericalPolarization	support polarization through scattering by spherical grains
    '''
    keysall = {'scatteringType': 'HenyeyGreenstein',}
    ele = xdm.Document().createElement('MeanTrustBenchmarkDustMix')       
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    return ele

@Element_ski.ski_element
def NonLTELineGasMix(ski):
    '''
    A gas mix supporting rotational transitions in specific molecules and atoms
    sourceWeight	Double	the weight of this secondary source for the number of photon packets launched
    wavelengthBias	Double	the fraction of secondary photon packet wavelengths sampled from a bias distribution
    species	Enum	the molecular or atomic species being represented
    --> Test	Fictive two-level test molecule (TT)
    --> Hydroxyl	Hydroxyl radical (OH)
    --> HydroxylHFS	Hydroxyl radical (OH) with hyperfine structure
    --> Formyl	Formyl cation (HCO+)
    --> CarbonMonoxide	Carbon monoxide (CO)
    --> AtomicCarbon	Atomic carbon (C)
    --> IonizedCarbon	Ionized carbon (C+)
    --> MolecularHydrogen	Molecular hydrogen (H2)
    numEnergyLevels	Int	the number of energy levels used (or 999 for all supported)
    defaultTemperature	Double	the default temperature of the gas
    defaultCollisionPartnerRatios	DoubleList	the default relative abundances of the collisional partners
    defaultTurbulenceVelocity	Double	the default (non-thermal) turbulence velocity
    maxChangeInLevelPopulations	Double	the maximum relative change for the level populations in a cell to be considered converged
    maxFractionNotConvergedCells	Double	the maximum fraction of not-converged cells for all cells to be considered converged
    maxChangeInGlobalLevelPopulations	Double	the maximum relative change for the global level populations to be considered converged
    lowestOpticalDepth	Double	Lower limit of (negative) optical depth along a cell diagonal
    storeMeanIntensities	Bool	store the mean radiation field intensity at each transition line
    initialLevelPopsFilename	String	the name of the file with initial level populations
    '''
    keysall = {'sourceWeight': '1',
               'wavelengthBias': '0.5',
               'species': 'CarbonMonoxide',
               'numEnergyLevels': '999',
               'defaultTemperature': '1000',
               'defaultCollisionPartnerRatios': '1e4',
               'defaultTurbulenceVelocity': '0 km/s',
               'maxChangeInLevelPopulations': '0.05',
               'maxFractionNotConvergedCells': '0.01',
               'lowestOpticalDepth': '-2',
               'storeMeanIntensities': 'false',
               'initialLevelPopsFilename': 'false',}
    ele = xdm.Document().createElement('NonLTELineGasMix')       
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    ele = ski_append(ele, ski['wavelengthBiasDistribution']) #TODO
    return ele

@Element_ski.ski_element
def SpinFlipAbsorptionMix(ski):
    '''
    A gas mix supporting the spin-flip 21 cm hydrogen absorption
    defaultTemperature	Double	the default temperature of the gas
    '''
    keysall = {'defaultTemperature': '1e4',}
    ele = xdm.Document().createElement('SpinFlipAbsorptionMix')       
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    return ele


@Element_ski.ski_element
def SpinFlipHydrogenGasMix(ski):
    '''
    A gas mix supporting the spin-flip 21 cm hydrogen transition
    sourceWeight	Double	the weight of this secondary source for the number of photon packets launched
    wavelengthBias	Double	the fraction of secondary photon packet wavelengths sampled from a bias distribution
    defaultMetallicity	Double	the default metallicity of the gas
    defaultTemperature	Double	the default temperature of the gas
    defaultNeutralSurfaceDensity	Double	the default neutral hydrogen surface density
    '''
    keysall = {'sourceWeight': '1',
               'wavelengthBias': '0.5',
               'defaultMetallicity': '0.02',
               'defaultNeutralSurfaceDensity': '10 Msun/pc2',}
    ele = xdm.Document().createElement('SpinFlipHydrogenGasMix')       
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    ele = ski_append(ele, ski['wavelengthBiasDistribution']) #TODO
    return ele

@Element_ski.ski_element
def TrivialGasMix(ski):
    '''
    A trivial gas mix for testing purposes
    absorptionCrossSection	Double	the absorption cross section per hydrogen atom
    scatteringCrossSection	Double	the scattering cross section per hydrogen atom
    asymmetryParameter	Double	the scattering asymmetry parameter

    '''
    keysall = {'absorptionCrossSection': '',
               'scatteringCrossSection': '0',
               'asymmetryParameter': '0',}
    ele = xdm.Document().createElement('TrivialGasMix')       
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    return ele

@Element_ski.ski_element
def XRayAtomicGasMix(ski):
    '''
    A gas mix supporting photo-absorption and fluorescence for X-ray wavelengths
    abundancies	DoubleList	the abundancies for the elements with atomic number Z = 1,...,30
    temperature	Double	the temperature of the gas or zero to disable thermal dispersion
    scatterBoundElectrons	Enum	implementation of scattering by bound electrons
    --> None	ignore bound electrons
    --> Free	use free-electron Compton scattering
    --> FreeWithPolarization	use free-electron Compton scattering with support for polarization
    --> Good	use smooth Rayleigh scattering and exact bound-Compton scattering
    --> Exact	use anomalous Rayleigh scattering and exact bound-Compton scattering
    '''
    keysall = {'abundancies': '1',
               'temperature': '1e4',
               'scatterBoundElectrons': 'Good',}
    ele = xdm.Document().createElement('XRayAtomicGasMix')       
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    return ele

@Element_ski.ski_element
def materialMixFamily(ski):
    '''
    a family of material mixes
    '''
    ele = xdm.Document().createElement('materialMixFamily')   
    ele.setAttribute('type', 'MaterialMixFamily')
    ele = ski_append(ele, ski['SelectDustMixFamily'])
    return ele

@Element_ski.ski_element
def SelectDustMixFamily(ski):
    '''
     a family of dust mixes specified in the configuration
    '''
    ele = xdm.Document().createElement('SelectDustMixFamily')   
    ele = ski_append(ele, ski['dustMixes'])
    return ele


@Element_ski.ski_element
def dustMixes(ski):
    '''
     a family of dust mixes specified in the configuration
     containing: ConfigurableDustMix DraineLiDustMix MRNDustMix MeanFileDustMix MeanInterstellarDustMix MeanIvezicBenchmarkDustMix MeanListDustMix
     MeanPascucciBenchmarkDustMix MeanPinteBenchmarkDustMix MeanTrustBenchmarkDustMix ThemisDustMix TrustBenchmarkDustMix WeingartnerDraineDustMix ZubkoDustMix
    '''
    ele = xdm.Document().createElement('dustMixes')   
    ele.setAttribute('type', 'DustMix')
    subkey = ski.get_default('dustMixes') #TODO
    ele = ski_append(ele, ski[subkey])
    return ele

@Element_ski.ski_element
def CellMedium(ski):
    '''
     a transfer medium imported from cuboidal cell data
     
    Scalar Property		Type	Description
    filename	String	the name of the file to be imported
    massType	Enum	the type of mass quantity to be imported
    --> MassDensity	mass density
    --> Mass	mass (volume-integrated)
    --> NumberDensity	number density
    --> Number	number (volume-integrated)
    massFraction	Double	the fraction of the mass to be included (or one to include all)
    importMetallicity	Bool	import a metallicity column
    importTemperature	Bool	import a temperature column
    maxTemperature	Double	the maximum temperature for included mass (or zero to include all)
    importVelocity	Bool	import velocity components (3 columns)
    importMagneticField	Bool	import magnetic field components (3 columns)
    importVariableMixParams	Bool	import parameter(s) to select a spatially varying material mix
    useColumns	String	a list of names corresponding to columns in the file to be imported
    # materialMix materialMixFamily
    '''
    keysall = {'filename':  'AdaptiveMeshMedium.txt',
               'massType':  'Mass',
               'massFraction':  '1',
               'importMetallicity':  'false',
               'importTemperature':  'false',
               'maxTemperature':  '0 K',
               'importVelocity':  'false',
               'importMagneticField':  'false',
               'importVariableMixParams':  'false',
               'useColumns':  '',
        
    }
    ele = xdm.Document().createElement('CellMedium')   
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    for i in ski.get_default('mediums_options').split(','):
        ele = ski_append(ele, ski[i])
    return ele

@Element_ski.ski_element
def GeometricMedium(ski):
    '''
    a transfer medium with a built-in geometry
    velocityMagnitude	docs	Double	the magnitude of the velocity (multiplier)
    magneticFieldStrength	docs	Double	the strength of the magnetic field (multiplier)
    # geometry, materialMix normalization velocityDistribution magneticFieldDistribution
    '''
    keysall = {'velocityMagnitude':  '0 km/s.txt',
               'magneticFieldStrength':  '0 T',}
    ele = xdm.Document().createElement('GeometricMedium')
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    for i in ski.get_default('mediums_options').split(','):
        ele = ski_append(ele, ski[i])
    return ele

@Element_ski.ski_element
def ParticleMedium(ski):
    '''
    a transfer medium imported from smoothed particle data
    
    Scalar Property		Type	Description
    filename	String	the name of the file to be imported
    massType	Enum	the type of mass quantity to be imported
    --> MassDensity	mass density
    --> Mass	mass (volume-integrated)
    --> NumberDensity	number density
    --> Number	number (volume-integrated)
    massFraction	Double	the fraction of the mass to be included (or one to include all)
    importMetallicity	Bool	import a metallicity column
    importTemperature	Bool	import a temperature column
    maxTemperature	Double	the maximum temperature for included mass (or zero to include all)
    importVelocity	Bool	import velocity components (3 columns)
    importMagneticField	Bool	import magnetic field components (3 columns)
    importVariableMixParams	Bool	import parameter(s) to select a spatially varying material mix
    useColumns	String	a list of names corresponding to columns in the file to be imported
    '''
    keysall = {'filename':  'ParticleMedium.txt',
               'massType':  'Mass',
               'massFraction':  '1',
               'importMetallicity':  'false',
               'importTemperature':  'false',
               'maxTemperature':  '0 K',
               'importVelocity':  'false',
               'importMagneticField':  'false',
               'importVariableMixParams':  'false',
               'useColumns':  '',
        
    }
    ele = xdm.Document().createElement('ParticleMedium')   
    for i in keysall:
        ele.setAttribute(i, keysall[i]) #TODO
    for i in ski.get_default('mediums_options').split(','):
        ele = ski_append(ele, ski[i])
    return ele

@Element_ski.ski_element
def VoronoiMeshMedium(ski):
    '''
    a transfer medium imported from data represented on an adaptive mesh (AMR grid)
    filename	String	the name of the file to be imported
    minX	Double	the start point of the domain in the X direction
    maxX	Double	the end point of the domain in the X direction
    minY	Double	the start point of the domain in the Y direction
    maxY	Double	the end point of the domain in the Y direction
    minZ	Double	the start point of the domain in the Z direction
    maxZ	Double	the end point of the domain in the Z direction
    massType    the type of mass quantity to be imported
        (MassDensity, Mass, MassDensityAndMass, NumberDensity, Number, NumberDensityAndNumber)
    massFraction	Double	the fraction of the mass to be included (or one to include all)
    importMetallicity	Bool	import a metallicity column
    importTemperature	Bool	import a temperature column
    maxTemperature	Double	the maximum temperature for included mass (or zero to include all)
    importVelocity	Bool	import velocity components (3 columns)
    importMagneticField	Bool	import magnetic field components (3 columns)
    importVariableMixParams	Bool	import parameter(s) to select a spatially varying material mix
    useColumns	String	a list of names corresponding to columns in the file to be imported
    '''
    keysall = {'filename':  'VoronoiMeshMedium.txt',
               'minX':  '-10 kpc',
               'maxX':  '10 kpc',
               'minY':  '-10 kpc',
               'maxY':  '10 kpc',
               'minZ':  '-10 kpc',
               'maxZ':  '10 kpc',
               'massType':  'Mass',
               'massFraction':  '1',
               'importMetallicity':  'false',
               'importTemperature':  'false',
               'maxTemperature':  '0 K',
               'importVelocity':  'false',
               'importMagneticField':  'false',
               'importVariableMixParams':  'false',
               'useColumns':  '',
        
    }
    ele = xdm.Document().createElement('VoronoiMeshMedium')
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    for i in ski.get_default('mediums_options').split(','):
        ele = ski_append(ele, ski[i])
    
    return ele

@Element_ski.ski_element
def samplingOptions(ski):
    '''
    a set of options related to media sampling for the spatial grid
    '''
    ele = xdm.Document().createElement('samplingOptions')
    ele.setAttribute('type', 'SamplingOptions')
    ele = ski_append(ele, ski['SamplingOptions'])
    return ele

@Element_ski.ski_element
def SamplingOptions(ski):
    '''
    a set of options related to media sampling for the spatial grid
    numDensitySamples	Int	the number of random density samples for determining spatial cell mass
    numPropertySamples	Int	the number of random samples for determining other medium properties
    aggregateVelocity	Enum	aggregating the bulk velocity from multiple medium components
    --> Average	Use the density-weighted average; missing values are taken to be zero
    --> Maximum	Use the vector with largest magnitude; missing values are taken to be zero
    --> First	Use the vector of the first medium component for which one is available
    '''
    keysall = {'numDensitySamples': '100',
               'numPropertySamples': '1',
               'aggregateVelocity': 'Average',}
    ele = xdm.Document().createElement('SamplingOptions')
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    return ele

@Element_ski.ski_element
def grid(ski):
    '''
    the spatial grid
    containing: AdaptiveMeshSpatialGrid CartesianSpatialGrid Cylinder2DSpatialGrid FileTreeSpatialGrid PolicyTreeSpatialGrid Sphere1DSpatialGrid Sphere2DSpatialGrid
        VoronoiMeshSpatialGrid
    '''
    ele = xdm.Document().createElement('grid')
    ele.setAttribute('type', 'SpatialGrid')
    subkey = ski.get_default('grid') #TODO
    ele = ski_append(ele, ski[subkey])
    return ele
    
@Element_ski.ski_element
def AdaptiveMeshSpatialGrid(ski):
    '''
    a spatial grid taken from an imported adaptive mesh snapshot
    '''
    ele = xdm.Document().createElement('AdaptiveMeshSpatialGrid')
    return ele

@Element_ski.ski_element
def CartesianSpatialGrid(ski):
    '''
    a Cartesian spatial grid
    minX	Double	the start point of the box in the X direction
    maxX	Double	the end point of the box in the X direction
    minY	Double	the start point of the box in the Y direction
    maxY	Double	the end point of the box in the Y direction
    minZ	Double	the start point of the box in the Z direction
    maxZ	Double	the end point of the box in the Z direction
    '''
    keysall = {'minX': '-10 kpc',
               'maxX': '10 kpc',
               'minY': '-10 kpc',
               'maxY': '10 kpc',
               'minZ': '-10 kpc',
               'maxZ': '10 kpc',}
    ele = xdm.Document().createElement('CartesianSpatialGrid')
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    ele = ski_append(ele,ski['meshX'])
    ele = ski_append(ele,ski['meshY'])
    ele = ski_append(ele,ski['meshZ'])
    return ele

@Element_ski.ski_element
def meshX(ski):
    '''
    the bin distribution in the X direction
    '''
    ele = xdm.Document().createElement('meshX')
    ele.setAttribute('type', 'MoveableMesh')
    subkey = ski.xxx #TODO
    ele = ski_append(ele, ski[subkey])
    return ele

@Element_ski.ski_element
def meshY(ski):
    '''
    the bin distribution in the Y direction
    '''
    ele = xdm.Document().createElement('meshY')
    ele.setAttribute('type', 'MoveableMesh')
    subkey = ski.xxx #TODO
    ele = ski_append(ele, ski[subkey])
    return ele

@Element_ski.ski_element
def meshZ(ski):
    '''
    the bin distribution in the Z direction
    '''
    ele = xdm.Document().createElement('meshZ')
    ele.setAttribute('type', 'MoveableMesh')
    subkey = ski.xxx #TODO
    ele = ski_append(ele, ski[subkey])
    return ele

@Element_ski.ski_element
def FileMesh(ski):
    '''
    a linear mesh
    numBins	Int	the number of bins in the mesh
    filename	String	the name of the file with the mesh border points
    '''
    keysall={'numBins': '100',
             'filename': 'FileMesh.txt'}
    ele = xdm.Document().createElement('FileMesh')
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    return ele

@Element_ski.ski_element
def LinMesh(ski):
    '''
    a linear mesh
    numBins	Int	the number of bins in the mesh
    '''
    keysall={'numBins': '100',}
    ele = xdm.Document().createElement('LinMesh')
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    return ele

@Element_ski.ski_element
def ListMesh(ski):
    '''
    a mesh specified inside the configuration file
    numBins	docs	Int	the number of bins in the mesh
    points	docs	DoubleList	the mesh border points
    '''
    ele = xdm.Document().createElement('ListMesh')
    keysall={'numBins': '100',
             'points': '0'}
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    return ele

@Element_ski.ski_element
def PowMesh(ski):
    '''
    a power-law mesh
    numBins	Int	the number of bins in the mesh
    ratio	Double	the bin width ratio between the last and the first bin
    '''
    ele = xdm.Document().createElement('PowMesh')
    keysall={'numBins': '100',
             'ratio': '1'}
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    return ele

@Element_ski.ski_element
def SymPowMesh(ski):
    '''
    a symmetric power-law mesh
    numBins	Int	the number of bins in the mesh
    ratio	Double	the bin width ratio between the last and the first bin
    '''
    ele = xdm.Document().createElement('SymPowMesh')
    keysall={'numBins': '100',
             'ratio': '1'}
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    return ele

@Element_ski.ski_element
def Cylinder2DSpatialGrid(ski):
    '''
    an axisymmetric spatial grid in cylindrical coordinates
    maxRadius	Double	the cylindrical radius of the grid
    minZ	Double	the start point of the cylinder in the Z direction
    maxZ	Double	the end point of the cylinder in the Z direction
    '''
    ele = xdm.Document().createElement('Cylinder2DSpatialGrid')
    keysall={'maxRadius': '20 kpc',
             'minZ': '-10 kpc',
             'maxZ': '10 kpc',}
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    ele = ski_append(ele, ski['meshRadial'])
    ele = ski_append(ele, ski['meshZ'])
    return ele

@Element_ski.ski_element
def meshRadial(ski):
    '''
    the bin distribution in the radial direction
    '''
    ele = xdm.Document().createElement('meshRadial')
    ele.setAttribute('type', 'Mesh')
    subkey = ski.xxx #TODO
    ele = ski_append(ele, ski[subkey])
    return ele

@Element_ski.ski_element
def FileTreeSpatialGrid(ski):
    '''
    a tree-based spatial grid loaded from a topology data file
    minX	Double	the start point of the box in the X direction
    maxX	Double	the end point of the box in the X direction
    minY	Double	the start point of the box in the Y direction
    maxY	Double	the end point of the box in the Y direction
    minZ	Double	the start point of the box in the Z direction
    maxZ	Double	the end point of the box in the Z direction
    filename	String	the name of the file with the tree topology data
    '''
    ele = xdm.Document().createElement('FileTreeSpatialGrid')
    keysall={'minX': '-10 kpc',
             'maxX': '10 kpc',
             'minY': '-10 kpc',
             'maxY': '10 kpc',
             'minZ': '-10 kpc',
             'maxZ': '10 kpc',
             'filename': 'FileTreeSpatialGrid.txt'}
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    return ele
    
@Element_ski.ski_element
def PolicyTreeSpatialGrid(ski):
    '''
    a tree-based spatial grid loaded from a topology data file
    minX	Double	the start point of the box in the X direction
    maxX	Double	the end point of the box in the X direction
    minY	Double	the start point of the box in the Y direction
    maxY	Double	the end point of the box in the Y direction
    minZ	Double	the start point of the box in the Z direction
    maxZ	Double	the end point of the box in the Z direction
    treeType	Enum	the type of tree
    --> OctTree	an octtree (8 children per node)
    --> BinTree	a binary tree (2 children per node)
    '''
    ele = xdm.Document().createElement('PolicyTreeSpatialGrid')
    keysall={'minX': '-10 kpc',
             'maxX': '10 kpc',
             'minY': '-10 kpc',
             'maxY': '10 kpc',
             'minZ': '-10 kpc',
             'maxZ': '10 kpc',
             'treeType': 'OctTree'}
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    ele = ski_append(ele, ski['policy'])
    return ele
@Element_ski.ski_element
def policy(ski):
    '''
     the tree construction policy (configuration options)
    '''
    ele = xdm.Document().createElement('policy')
    ele.setAttribute('type', 'TreePolicy')
    subkey = ski.get_default('policy')    #TODO
    ele = ski_append(ele, ski[subkey])
    return ele

@Element_ski.ski_element
def DensityTreePolicy(ski):
    '''
    a tree grid construction policy using the medium density distribution
    minLevel	Int	the minimum level of grid refinement
    maxLevel	Int	the maximum level of grid refinement
    maxDustFraction	Double	the maximum fraction of dust contained in each cell
    maxDustOpticalDepth	Double	the maximum diagonal dust optical depth for each cell
    wavelength	Double	the wavelength at which to evaluate the optical depth
    maxDustDensityDispersion	Double	the maximum dust density dispersion in each cell
    maxElectronFraction	Double	the maximum fraction of electrons contained in each cell
    maxGasFraction	Double	the maximum fraction of gas contained in each cell
    '''
    keysall={'minLevel': '7',
             'maxLevel': '11',
             'maxDustFraction': '1e-6',
             'maxDustOpticalDepth': '0',
             'wavelength': '0.55 micron',
             'maxDustDensityDispersion': '0',
             'maxElectronFraction': '0',
             'maxGasFraction': '1e-6',}
    ele = xdm.Document().createElement('DensityTreePolicy')
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    return ele

@Element_ski.ski_element
def NestedDensityTreePolicy(ski):
    '''
    a tree grid construction policy using a nested density tree policy
    minLevel	Int	the minimum level of grid refinement
    maxLevel	Int	the maximum level of grid refinement
    maxDustFraction	Double	the maximum fraction of dust contained in each cell
    maxDustOpticalDepth	Double	the maximum diagonal dust optical depth for each cell
    wavelength	Double	the wavelength at which to evaluate the optical depth
    maxDustDensityDispersion	Double	the maximum dust density dispersion in each cell
    maxElectronFraction	Double	the maximum fraction of electrons contained in each cell
    maxGasFraction	Double	the maximum fraction of gas contained in each cell
    innerMinX	Double	the start point of the inner box in the X direction
    innerMaxX	Double	the end point of the inner box in the X direction
    innerMinY	Double	the start point of the inner box in the Y direction
    innerMaxY	Double	the end point of the inner box in the Y direction
    innerMinZ	Double	the start point of the inner box in the Z direction
    innerMaxZ	Double	the end point of the inner box in the Z direction
    '''
    ele = xdm.Document().createElement('NestedDensityTreePolicy')
    keysall={'minLevel': '7',
             'maxLevel': '11',
             'maxDustFraction': '1e-6',
             'maxDustOpticalDepth': '0',
             'wavelength': '0.55 micron',
             'maxDustDensityDispersion': '0',
             'maxElectronFraction': '0',
             'maxGasFraction': '1e-6',
             'innerMinX': '-10 kpc',
             'innerMaxX': '10 kpc',
             'innerMinY': '-10 kpc',
             'innerMaxY': '10 kpc',
             'innerMinZ': '-10 kpc',
             'innerMaxZ': '10 kpc',}
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    ele = ski_append(ele, ski['innerPolicy'])
    return ele

@Element_ski.ski_element
def innerPolicy(ski):
    '''
    the density tree policy for the inner region
    DensityTreePolicy   NestedDensityTreePolicy
    '''
    ele = xdm.Document().createElement('innerPolicy')
    ele.setAttribute('type', 'DensityTreePolicy')
    subkey = ski.xxx    #TODO
    ele = ski_append(ele, ski[subkey])
    return ele

@Element_ski.ski_element
def SiteListTreePolicy(ski):
    '''
    a tree grid construction policy using positions defined by an imported medium
    minLevel	Int	the minimum level of grid refinement
    maxLevel	Int	the maximum level of grid refinement
    numExtraLevels	Int	the number of additional subdivision levels
    '''
    keysall={'minLevel': '7',
             'maxLevel': '11',
             'numExtraLevels': '0',}
    ele = xdm.Document().createElement('SiteListTreePolicy')
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    return ele

@Element_ski.ski_element
def Sphere1DSpatialGrid(ski):
    '''
    a spherically symmetric spatial grid
    maxRadius	Double	the outer radius of the grid
    minRadius	Double	the inner radius of the grid
    '''
    keysall={'maxRadius': '10 kpc',
             'minRadius': '0 kpc',}
    ele = xdm.Document().createElement('Sphere1DSpatialGrid')
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    ele = ski_append(ele, ski['meshRadial'])
    return ele


@Element_ski.ski_element
def Sphere2DSpatialGrid(ski):
    '''
    an axisymmetric spatial grid in spherical coordinates
    maxRadius	Double	the outer radius of the grid
    '''
    keysall={'maxRadius': '10 kpc'}
    ele = xdm.Document().createElement('Sphere2DSpatialGrid')
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    ele = ski_append(ele, ski['meshRadial'])
    ele = ski_append(ele, ski['meshPolar'])
    return ele

@Element_ski.ski_element
def meshPolar(ski):
    '''
    the bin distribution in the polar direction
    '''
    ele = xdm.Document().createElement('meshPolar')
    ele.setAttribute('type', 'Mesh')
    subkey = ski.xxx    #TODO
    
    ele = ski_append(ele, ele[subkey])
    return ele

@Element_ski.ski_element
def VoronoiMeshSpatialGrid(ski):
    '''
    a Voronoi tessellation-based spatial grid
    minX	Double	the start point of the box in the X direction
    maxX	Double	the end point of the box in the X direction
    minY	Double	the start point of the box in the Y direction
    maxY	Double	the end point of the box in the Y direction
    minZ	Double	the start point of the box in the Z direction
    maxZ	Double	the end point of the box in the Z direction
    policy	Enum	the policy for determining the positions of the sites
    --> Uniform	random from uniform distribution
    --> CentralPeak	random from distribution with a steep central peak
    --> DustDensity	random from dust density distribution
    --> ElectronDensity	random from electron density distribution
    --> GasDensity	random from gas density distribution
    --> File	loaded from text column data file
    --> ImportedSites	positions of particles, sites or cells in imported distribution
    --> ImportedMesh	employ imported Voronoi mesh in medium system
    numSites	Int	the number of random sites (or cells in the grid)
    filename	String	the name of the file containing the site positions
    relaxSites	Bool	perform site relaxation to avoid overly elongated cells
    '''
    keysall={'minX': '-10 kpc',
             'maxX': '10 kpc',
             'minY': '-10 kpc',
             'maxY': '10 kpc',
             'minZ': '-10 kpc',
             'maxZ': '10 kpc',
             'policy': 'DustDensity',
             'numSites': '500',
             'filename': 'VoronoiMeshSpatialGrid.txt',
             'relaxSites': 'false',}
    ele = xdm.Document().createElement('VoronoiMeshSpatialGrid')
    for i in keysall:
        ele.setAttribute(i, keysall[i])  #TODO
    return ele
