### Copied from pyabf (see below), with a little fix for Python 2

"""Code to interact with ABF files. https://github.com/swharden/pyABF/ """

import warnings
import numpy as np

np.set_printoptions(suppress=True)  # don't use scientific notation
import matplotlib.pyplot as plt

"""
ABFheader - a python module to provide simple access to header content of ABF (Axon Binary Format) files.

GOAL:
    This standalone module is minimal and simplistic to facilitate learning and porting to other languages.
    Although other scripts in this project have complex dependencies (like numpy and matplotlib) this script
    is entirely standalone and uses only standard python libraries.

NOTE ABOUT NUMPY:
    When it comes time to read values out of ABF files and scale them, numpy-optimized vector math offers
    a massive performance improvement! I'm not making numpy required to use this code. Everything will
    work without numpy. BUT, if numpy is found, it will use it when doing vector operations on signals.

RESOURCES:
    Code here is a blend of original ideas and ideas compiled from reading others' code. For a full list of
    resources, see my Unofficial ABF File Format Guide (among other things) on the project homepage:
    https://github.com/swharden/pyABF/
"""

import os
import time
import datetime
import struct
import collections
import warnings

try:
    import numpy as np
except:
    np = False


class ABFheader:
    def __init__(self, abfFileName, loadDataIntoMemory=True):
        """
        The ABFheader class provides low-level access to ABF file contents (header values and signal data).
        You can pull all the information you need from abf.header and abf.data. Additional functions help
        display header data in various formats.

        See what information is available about the ABF
            >>> abfHeader=ABFheader("filename.abf")
            >>> abfHeader.show()
            ### Header ###
            fFileSignature = b'ABF2'
            fFileVersionNumber = (0, 0, 0, 2)
            uFileInfoSize = 512
            lActualEpisodes = 187
            uFileStartDate = 20161205
            ...

            * This will list values which can be pulled from the abf.header dictionary by their names

        Get a value from the header by its name:
            >>> abf.header['lActualEpisodes'] # number of sweeps (episodes)
            187

            * Each variable in the header has a special name, but it's not always obvious what it means.
            * Check the ABF Format PDF http://mdc.custhelp.com/euf/assets/content/ABFHelp.pdf for info

        Get recorded signal data:
            >>> abf.data
            [-70.123, -70.321, ..., -70.213, -70.031]

            * If numpy is available, abf.data will be a numpy.ndarray with a 32-bit float datatype.
            * If numpy is not available, abf.data will be a list of python integers
            * Note that this doesn't work sweep by sweep, it's ALL the data in one array!
            * Dividing data into sweeps is your job.

        Separate multichannel data (it's just interleaved):
            >>> abf=ABFheader(abfFileName)
            >>> for i in range(abf.header['dataChannels']):
            >>>     print(abf.data[i::abf.header['dataChannels']])
            [-146.34 -144.96 -145.65 ..., -117.11 -116.55 -117.37]
            [ 624.13  624.06  624.19 ...,  624.32  624.45 624.45]

        Save a formatted header to file:
            >>> abf.saveHTML("someFile.html")
            >>> abf.saveMarkdown("someFile.md")

        """

        t1 = time.time()
        self._fb = open(abfFileName, 'rb')
        self._fileSignature = self._fb.read(4)
        if self._fileSignature == b'ABF ':
            self._readHeaderABF1()
        elif self._fileSignature == b'ABF2':
            self._readHeaderABF2()
        else:
            raise ValueError('invalid file (does not appear to be an ABF at all!)')

        # read the signal data into memory and scale it
        self.data = None
        if loadDataIntoMemory:
            self._fileReadData()

        # we are now done reading contents of the ABF file
        self._fb.close()

        # stop the performance counter and calculate the time it took to load/process the ABF file
        self.abfLoadTime = (time.time() - t1)

    ### HEADER READING AND VALUE TOUCHUP

    def _readHeaderABF1(self):
        """populate self.header with values read the ABF1 header. Not many extra keys are added."""
        self.header = collections.OrderedDict()
        self.header["### ABF1 Header ###"] = [None]
        for key, offset, varFormat in HEADERV1:
            self._fb.seek(offset)
            varVal = list(struct.unpack(varFormat, self._fb.read(struct.calcsize(varFormat))))
            for i, item in enumerate(varVal):
                if type(item) == bytes:
                    varVal[i] = item.decode().strip()  # TODO: this for ABF2 comments
            self.header[key] = varVal
        for key in [key for key in self.header.keys() if len(self.header[key]) == 1]:
            self.header[key] = self.header[key][0]  # flatten lists with just 1 element

        # add a few extra things I think are useful. This list isn't as extensive as ABF2
        sz = struct.calcsize("2i") if self.header['nDataFormat'] == 0 else struct.calcsize("4f")
        self.header['dataByteStart'] = self.header['lDataSectionPtr'] * 512 + self.header['nNumPointsIgnored'] * sz
        self.header['dataPointCount'] = self.header['lActualAcqLength']
        self.header['dataChannels'] = self.header['nADCNumChannels']
        self.header['timeSecPerPoint'] = self.header['fADCSampleInterval'] / 1e6
        self.header['timePointPerSec'] = 1e6 / self.header['fADCSampleInterval']
        self.header['abfFilename'] = os.path.abspath(self._fb.name)
        self.header['abfID'] = os.path.basename(self._fb.name)[:-4]
        self.header['abfDatetime'] = "ABF1 not sure"
        self.header['sweepPointCount'] = int(self.header['lNumSamplesPerEpisode'] / self.header['dataChannels'])
        self.header['rate'] = 1e6 / self.header['fADCSampleInterval']
        self.header['sweepCount'] = self.header['lActualEpisodes']
        self.header['sweepLengthSec'] = self.header['sweepPointCount'] * self.header['timeSecPerPoint']
        self.header['mode'] = "IC" if self.header['sADCUnits'][0] == "mV" else "VC"
        self.header['units'] = "mV" if self.header['mode'] == "IC" else "pA"
        self.header['unitsCommand'] = "pA" if self.header['mode'] == "IC" else "mV"
        self.header['commandHoldingByDAC'] = self.header['fEpochInitLevel']
        self.header['lEpochPulsePeriod'] = None  # pulses unsupported in ABF1
        self.header['lEpochPulseWidth'] = None  # pulses unsupported in ABF1
        self.header['nEpochDigitalOutput'] = self.header['nDigitalValue']
        self.header['dataScale'] = self.header['lADCResolution'] / 1e6
        self.header['protocolPath'] = self._readHeaderProtocol()
        self.header['protocol'] = os.path.basename(self.header['protocolPath'])
        self._calculateScaleFactor()
        # TODO: make ABF1 header items the same as ABF2 headers. There are so many commonalities.
        return

    def _readHeaderABF2(self):
        """populate self.header with values read the ABF2 header. Extra helpful keys are added too."""

        # pull values out of the header
        self.header = collections.OrderedDict()
        self._byteMap = collections.OrderedDict()
        self._fileReadStructMap(HEADER, sectionName="Header")
        self._fileReadStructMap(SECTIONS, 76, 16, sectionName="Section Map")
        self._fileReadSection('ProtocolSection', PROTO)
        self._fileReadSection('ADCSection', ADC)
        self._fileReadSection('DACSection', DAC)
        self._fileReadSection('EpochPerDACSection', EPPERDAC)
        self._fileReadSection('EpochSection', EPSEC)
        self._fileReadSection('TagSection', TAGS)

        # touch-up comments
        if 'sComment' in self.header.keys():
            self.header['sComment'] = [x.decode().strip() for x in self.header['sComment']]

        # make values that are a list with just 1 element just the element (no list required)
        for key in [key for key in self.header.keys() if len(self.header[key]) == 1]:
            self.header[key] = self.header[key][0]

        # add a few extra things I think are useful.
        self.header["### Extras ###"] = None
        self.header['abfFilename'] = os.path.abspath(self._fb.name)
        self.header['abfID'] = os.path.basename(self._fb.name)[:-4]
        dt = datetime.datetime.strptime(str(self.header['uFileStartDate']), "%Y%M%d")
        self.header['abfDatetime'] = dt + datetime.timedelta(seconds=self.header['uFileStartTimeMS'] / 1000)
        self.header['dataByteStart'] = self.header['DataSection'][0] * 512
        self.header['dataPointCount'] = self.header['DataSection'][2]
        self.header['dataChannels'] = self.header['ADCSection'][2]
        self.header['timeSecPerPoint'] = self.header['fADCSequenceInterval'] / 1e6
        self.header['timePointPerSec'] = 1e6 / self.header['fADCSequenceInterval']
        self.header['rate'] = 1e6 / self.header['fADCSequenceInterval']
        self.header['sweepCount'] = self.header['lActualEpisodes']
        self.header['sweepPointCount'] = int(self.header['lNumSamplesPerEpisode'] / self.header['dataChannels'])
        if self.header['sweepCount'] == 0:  # gap free mode
            self.header['sweepPointCount'] = int(self.header['dataPointCount'] / self.header['dataChannels'])
        self.header['sweepLengthSec'] = self.header['sweepPointCount'] * self.header['timeSecPerPoint']
        self.header['gain'] = self.header['fTelegraphAdditGain']
        self.header['mode'] = "IC" if self.header['nTelegraphMode'] else "VC"
        self.header['units'] = "mV" if self.header['mode'] == "IC" else "pA"
        self.header['unitsCommand'] = "pA" if self.header['mode'] == "IC" else "mV"
        self.header['commandHoldingByDAC'] = self.header['fDACHoldingLevel']
        self.header['protocolPath'] = self._readHeaderProtocol()
        self.header['protocol'] = os.path.basename(self.header['protocolPath'])

        self._calculateScaleFactor()

    def _readHeaderProtocol(self):
        """read the the protocol filename out of the ABF header"""
        self._fb.seek(0)
        raw = self._fb.read(self.header['dataByteStart'])
        match = b".pro"
        matchI = raw.find(match)
        if matchI:
            chunk = raw[matchI - 256:matchI + len(match)]
            proto = chunk.split(b"\x00")[-1]
            return proto.decode()
        else:
            return None

    def _calculateScaleFactor(self):
        """
        Populates header['dataScale'] with a data scale multiplier. Note this only reports for channel 0,
        and multi-channel recordings may need to calculate this for each channel individually.
        """
        dataScale = 1
        dataScale /= self._first(self.header['fInstrumentScaleFactor'])
        dataScale /= self._first(self.header['fSignalGain'])
        dataScale /= self._first(self.header['fADCProgrammableGain'])
        if self.header['nTelegraphEnable']:
            dataScale /= self._first(self.header['fTelegraphAdditGain'])
        dataScale *= self._first(self.header['fADCRange'])
        dataScale /= self._first(self.header['lADCResolution'])
        dataScale += self._first(self.header['fInstrumentOffset'])
        dataScale -= self._first(self.header['fSignalOffset'])
        self.header['dataScale'] = dataScale

    def _first(self, thing):
        """If a thing is just a thing, return it. If a thing is a list, return the first thing."""
        if type(thing) in [int, float]:
            return thing
        if type(thing) is list:
            return thing[0]
        if len(thing):
            return thing[0]
        return thing

    ### FILE READING

    def _fileReadStructMap(self, structMap, startByte=0, fixedOffset=None, sectionName=None):
        """Given a string of varName_varFormat structs, get the objects from the file."""
        if sectionName:
            self.header["### %s ###" % sectionName] = [None]
            self._byteMap["### %s (fixed byte positions) ###" % sectionName] = [None]
        self._fb.seek(startByte)
        for structCode in structMap.replace("\n", "").split(","):
            varName, varFormat = structCode.strip().split("_")
            if sectionName:
                self._byteMap.setdefault(varName, []).append("%d" % (self._fb.tell()))
            else:
                self._byteMap.setdefault(varName, []).append("+%d" % (self._fb.tell() - startByte))
            varVal = struct.unpack(varFormat, self._fb.read(struct.calcsize(varFormat)))
            varVal = varVal if len(varVal) > 1 else varVal[0]
            self.header.setdefault(varName, []).append(varVal)
            if fixedOffset:
                self._fb.read(fixedOffset - struct.calcsize(varFormat))

    def _fileReadSection(self, sectionName, structMap):
        """Read a structure map repeatedly according to its name in the section map."""
        self.header["### %s ###" % sectionName] = [None]
        self._byteMap["### %s (section byte offsets) ###" % sectionName] = [None]
        entryStartBlock, entryBytes, entryCount = self.header[sectionName][0]
        for entryNumber in range(entryCount):
            self._fileReadStructMap(structMap, entryStartBlock * 512 + entryNumber * entryBytes)

    def _fileReadData(self):
        """
        Read the full file data into memory. Scale it too. Uses numpy if available.
        If the signal is multiple channels (dataChannels) it's your responsability to reshape the ouptut.
        """
        self._fb.seek(self.header['dataByteStart'])
        pointCount = self.header['dataPointCount']
        scaleFactor = self.header['dataScale']
        if True:
            self.data = np.fromfile(self._fb, dtype=np.int16, count=pointCount)
            self.data = np.multiply(self.data, scaleFactor, dtype='float32')
        else:
            warnings.warn("Numpy is not installed. We can go without it, but performance will suffer.")
            self.data = struct.unpack("%dh" % (pointCount), self._fb.read(pointCount * 2))
            self.data = [point * scaleFactor for point in self.data]

    ### HEADER DISPLAY

    def show(self):
        """Display the contents of the header to the console in an easy to read format."""
        for key in self.header.keys():
            if key.startswith("###"):
                print("\n%s" % key)
            else:
                print("%s = %s" % (key, self.header[key]))
        print()

    def saveHTML(self, fname):
        """Generate a HTML-formatted document with all header information."""
        html = "<html><body><code>"
        for key in self.header.keys():
            if key.startswith("###"):
                html += "<br><b style='font-size: 200%%;'>%s</b><br>" % (key.replace("#", "").strip())
            else:
                html += "%s = %s<br>" % (key, self.header[key])
        html += "</code></html></body>"
        with open(fname, 'w') as f:
            f.write(html)
        print("wrote", os.path.abspath(fname))

    def saveMarkdown(self, fname):
        """Generate a markdown-formatted document with all header information."""
        out = "# ABF Header Contents\n"
        for key in self.header.keys():
            if key.startswith("###"):
                out += "\n## %s\n" % (key.replace("#", "").strip())
            else:
                out += "* %s = `%s`\n" % (key, self.header[key])
        with open(fname, 'w') as f:
            f.write(out)
        print("wrote", os.path.abspath(fname))


        ### Data structures for ABF2 files:


HEADER = """fFileSignature_4s,fFileVersionNumber_4b,uFileInfoSize_I,lActualEpisodes_I,uFileStartDate_I,
uFileStartTimeMS_I,uStopwatchTime_I,nFileType_H,nDataFormat_H,nSimultaneousScan_H,nCRCEnable_H,uFileCRC_I,
FileGUID_I,unknown1_I,unknown2_I,unknown3_I,uCreatorVersion_I,uCreatorNameIndex_I,uModifierVersion_I,
uModifierNameIndex_I,uProtocolPathIndex_I"""
SECTIONS = """ProtocolSection_IIl,ADCSection_IIl,DACSection_IIl,EpochSection_IIl,ADCPerDACSection_IIl,
EpochPerDACSection_IIl,UserListSection_IIl,StatsRegionSection_IIl,MathSection_IIl,StringsSection_IIl,
DataSection_IIl,TagSection_IIl,ScopeSection_IIl,DeltaSection_IIl,VoiceTagSection_IIl,SynchArraySection_IIl,
AnnotationSection_IIl,StatsSection_IIl"""
PROTO = """nOperationMode_h,fADCSequenceInterval_f,bEnableFileCompression_b,sUnused_3s,
uFileCompressionRatio_I,fSynchTimeUnit_f,fSecondsPerRun_f,lNumSamplesPerEpisode_i,lPreTriggerSamples_i,
lEpisodesPerRun_i,lRunsPerTrial_i,lNumberOfTrials_i,nAveragingMode_h,nUndoRunCount_h,nFirstEpisodeInRun_h,
fTriggerThreshold_f,nTriggerSource_h,nTriggerAction_h,nTriggerPolarity_h,fScopeOutputInterval_f,
fEpisodeStartToStart_f,fRunStartToStart_f,lAverageCount_i,fTrialStartToStart_f,nAutoTriggerStrategy_h,
fFirstRunDelayS_f,nChannelStatsStrategy_h,lSamplesPerTrace_i,lStartDisplayNum_i,lFinishDisplayNum_i,
nShowPNRawData_h,fStatisticsPeriod_f,lStatisticsMeasurements_i,nStatisticsSaveStrategy_h,fADCRange_f,
fDACRange_f,lADCResolution_i,lDACResolution_i,nExperimentType_h,nManualInfoStrategy_h,nCommentsEnable_h,
lFileCommentIndex_i,nAutoAnalyseEnable_h,nSignalType_h,nDigitalEnable_h,nActiveDACChannel_h,
nDigitalHolding_h,nDigitalInterEpisode_h,nDigitalDACChannel_h,nDigitalTrainActiveLogic_h,nStatsEnable_h,
nStatisticsClearStrategy_h,nLevelHysteresis_h,lTimeHysteresis_i,nAllowExternalTags_h,nAverageAlgorithm_h,
fAverageWeighting_f,nUndoPromptStrategy_h,nTrialTriggerSource_h,nStatisticsDisplayStrategy_h,
nExternalTagType_h,nScopeTriggerOut_h,nLTPType_h,nAlternateDACOutputState_h,nAlternateDigitalOutputState_h,
fCellID_3f,nDigitizerADCs_h,nDigitizerDACs_h,nDigitizerTotalDigitalOuts_h,nDigitizerSynchDigitalOuts_h,
nDigitizerType_h"""
ADC = """nADCNum_h,nTelegraphEnable_h,nTelegraphInstrument_h,fTelegraphAdditGain_f,
fTelegraphFilter_f,fTelegraphMembraneCap_f,nTelegraphMode_h,fTelegraphAccessResistance_f,nADCPtoLChannelMap_h,
nADCSamplingSeq_h,fADCProgrammableGain_f,fADCDisplayAmplification_f,fADCDisplayOffset_f,
fInstrumentScaleFactor_f,fInstrumentOffset_f,fSignalGain_f,fSignalOffset_f,fSignalLowpassFilter_f,
fSignalHighpassFilter_f,nLowpassFilterType_b,nHighpassFilterType_b,fPostProcessLowpassFilter_f,
nPostProcessLowpassFilterType_c,bEnabledDuringPN_b,nStatsChannelPolarity_h,lADCChannelNameIndex_i,
lADCUnitsIndex_i"""
DAC = """nDACNum_h,nTelegraphDACScaleFactorEnable_h,fInstrumentHoldingLevel_f,fDACScaleFactor_f,
fDACHoldingLevel_f,fDACCalibrationFactor_f,fDACCalibrationOffset_f,lDACChannelNameIndex_i,
lDACChannelUnitsIndex_i,lDACFilePtr_i,lDACFileNumEpisodes_i,nWaveformEnable_h,nWaveformSource_h,
nInterEpisodeLevel_h,fDACFileScale_f,fDACFileOffset_f,lDACFileEpisodeNum_i,nDACFileADCNum_h,nConditEnable_h,
lConditNumPulses_i,fBaselineDuration_f,fBaselineLevel_f,fStepDuration_f,fStepLevel_f,fPostTrainPeriod_f,
fPostTrainLevel_f,nMembTestEnable_h,nLeakSubtractType_h,nPNPolarity_h,fPNHoldingLevel_f,nPNNumADCChannels_h,
nPNPosition_h,nPNNumPulses_h,fPNSettlingTime_f,fPNInterpulse_f,nLTPUsageOfDAC_h,nLTPPresynapticPulses_h,
lDACFilePathIndex_i,fMembTestPreSettlingTimeMS_f,fMembTestPostSettlingTimeMS_f,nLeakSubtractADCIndex_h"""
EPPERDAC = """nEpochNum_h,nDACNum_h,nEpochType_h,fEpochInitLevel_f,fEpochLevelInc_f,
lEpochInitDuration_i,lEpochDurationInc_i,lEpochPulsePeriod_i,lEpochPulseWidth_i"""
TAGS = """lTagTime_i,sComment_56s,nTagType_h,nVoiceTagNumberorAnnotationIndex_h"""
EPSEC = """nEpochNum_h,nEpochDigitalOutput_h"""

### Data structures for ABF1 files:
HEADERV1 = [('fFileSignature', 0, '4s'), ('fFileVersionNumber', 4, 'f'), ('nOperationMode', 8, 'h'),
            ('lActualAcqLength', 10, 'i'), ('nNumPointsIgnored', 14, 'h'), ('lActualEpisodes', 16, 'i'),
            ('lFileStartTime', 24, 'i'),
            ('lDataSectionPtr', 40, 'i'), ('lTagSectionPtr', 44, 'i'), ('lNumTagEntries', 48, 'i'),
            ('lSynchArrayPtr', 92, 'i'),
            ('lSynchArraySize', 96, 'i'), ('nDataFormat', 100, 'h'), ('nADCNumChannels', 120, 'h'),
            ('fADCSampleInterval', 122, 'f'),
            ('fSynchTimeUnit', 130, 'f'), ('lNumSamplesPerEpisode', 138, 'i'), ('lPreTriggerSamples', 142, 'i'),
            ('lEpisodesPerRun', 146, 'i'), ('fADCRange', 244, 'f'), ('lADCResolution', 252, 'i'),
            ('nFileStartMillisecs', 366, 'h'),
            ('nADCPtoLChannelMap', 378, '16h'), ('nADCSamplingSeq', 410, '16h'), ('sADCChannelName', 442, '10s' * 16),
            ('sADCUnits', 602, '8s' * 16), ('fADCProgrammableGain', 730, '16f'), ('fInstrumentScaleFactor', 922, '16f'),
            ('fInstrumentOffset', 986, '16f'), ('fSignalGain', 1050, '16f'), ('fSignalOffset', 1114, '16f'),
            ('nDigitalEnable', 1436, 'h'), ('nActiveDACChannel', 1440, 'h'), ('nDigitalHolding', 1584, 'h'),
            ('nDigitalInterEpisode', 1586, 'h'), ('nDigitalValue', 2588, '10h'), ('lDACFilePtr', 2048, '2i'),
            ('lDACFileNumEpisodes', 2056, '2i'), ('fDACCalibrationFactor', 2074, '4f'),
            ('fDACCalibrationOffset', 2090, '4f'),
            ('nWaveformEnable', 2296, '2h'), ('nWaveformSource', 2300, '2h'), ('nInterEpisodeLevel', 2304, '2h'),
            ('nEpochType', 2308, '20h'), ('fEpochInitLevel', 2348, '20f'), ('fEpochLevelInc', 2428, '20f'),
            ('lEpochInitDuration', 2508, '20i'), ('lEpochDurationInc', 2588, '20i'), ('nTelegraphEnable', 4512, '16h'),
            ('fTelegraphAdditGain', 4576, '16f'), ('sProtocolPath', 4898, '384s')]


def _compareHeaders(abfFile1, abfFile2):
    """Given two ABF filenames, show how their headers are different."""
    header1 = ABFheader(abfFile1).header
    header2 = ABFheader(abfFile2).header
    for key in header1.keys():
        if not key in header2.keys():
            continue
        if key.startswith("#"):
            print("\n" + key)
        if header1[key] == header2[key]:
            continue
        if type(header1[key]) in [list, tuple]:
            continue
        else:
            print(key, header1[key], header2[key])


def _graphSomeData(abfFileName):
    """Graph data from a file. Exclusively used for testing."""
    import matplotlib.pyplot as plt
    abf = ABFheader(abfFileName)
    plt.figure(figsize=(10, 2))
    for i in range(abf.header['dataChannels']):
        Ys = abf.data[i::abf.header['dataChannels']]
        Ys = Ys[20000 * 13:20000 * 19]
        Xs = np.arange(len(Ys)) * abf.header['timeSecPerPoint']
        plt.plot(Xs, Ys, label="channel %d" % (i + 1))
    plt.legend(fontsize=8)
    plt.margins(0, .1)
    plt.tight_layout()
    plt.show()

class ABF:
    def __init__(self, abf):
        """The ABF class provides easy pythonic access to header and signal
        data in ABF files. Although it is typically instantiated with a path
        (string), you can also use an ABF or ABFheader. This will reset that
        ABF to default values without re-loading original data from the file.

        Quick start:
            >>> abf = ABF("/path/to/file.abf")
            >>> abf.setSweep(0) # load data from the first sweep
            >>> print(abf.dataY) # signal data
            >>> print(abf.dataX) # timestamps
            >>> print(abf.dataC) # command waveform

        See all the properties available to you:
            >>> abf.help()

        Developers can access the ABFheader class features:
            >>> abf._abfHeader.saveHTML()

        """

        # get our abfHeader in order depending on what type of object we were given
        if type(abf) is str:
            self._abfHeader = ABFheader(abf)
        elif str(type(abf)).endswith(".ABF'>"):
            self._abfHeader = abf._abfHeader
        elif str(type(abf)).endswith(".ABFheader'>"):
            self._abfHeader = abf
        else:
            raise ValueError('abf must be a file path (str), ABF object, or ABFheader object.')

        ### Populate meaningful ABF attributes. Think about how you will use them: abf.something
        self.ID = self._abfHeader.header['abfID']
        self.filename = self._abfHeader.header['abfFilename']
        self.datetime = self._abfHeader.header['abfDatetime']
        self.pointDurSec = self._abfHeader.header['timeSecPerPoint']
        self.pointDurMS = self._abfHeader.header['timeSecPerPoint'] * 1000.0
        self.pointsPerSweep = self._abfHeader.header['sweepPointCount']
        self.pointsPerSec = self._abfHeader.header['rate']
        self.pointsPerMS = self.pointsPerSec / 1000
        self.dataChannels = self._abfHeader.header['dataChannels']
        self.sweepCount = self._abfHeader.header['sweepCount']
        self.sweepLengthSec = self._abfHeader.header['sweepLengthSec']
        self.sweepPointCount = self._abfHeader.header['sweepPointCount']
        self.sweepList = np.arange(max(1, self.sweepCount))  # zero for gap-free files
        self.sweepTimesSec = self.sweepList * self.sweepLengthSec
        self.sweepTimesMin = self.sweepTimesSec / 60
        self.mode = self._abfHeader.header['mode']
        self.units = self._abfHeader.header['units']
        self.unitsLong = "Membrane Potential (mV)" if self.units is 'mV' else "Membrane Current (pA)"
        self.unitsCommand = self._abfHeader.header['unitsCommand']
        self.unitsCommandLong = "Command Potential (mV)" if self.unitsCommand is 'mV' else "Command Current (pA)"
        self.commandHoldingByDAC = self._abfHeader.header['commandHoldingByDAC']
        self.commandHold = self.commandHoldingByDAC[0]
        self.experimentLengthSec = self.sweepLengthSec * self.sweepCount
        self.unitsTime = "seconds"
        self.unitsTimeLong = "Signal Time (seconds)"
        self.protocolPath = self._abfHeader.header['protocolPath']
        self.protocol = self._abfHeader.header['protocol']

        ### Load Comments
        self.commentsExist, self.commentTags, self.commentTimes = False, [], []
        if "sComment" in self._abfHeader.header.keys():
            self.commentsExist = True
            self.commentTags = self._abfHeader.header["sComment"]
            self.commentTimes = self._abfHeader.header["lTagTime"]
            self.commentTimesSec = [x * self._abfHeader.header["fSynchTimeUnit"] / 1e6 for x in self.commentTimes]
            self.commentTimesMin = [x * 60 for x in self.commentTimesSec]
            self.commentSweeps = [int(x / self.sweepLengthSec) for x in self.commentTimesSec]

        ### Preload signal and time data (totalling ~10MB of memory per minute of 20kHz recording)
        self.signalData = self._abfHeader.data / self.dataChannels
        self.signalTimes = np.arange(len(self.signalData), dtype='float32') * self.pointDurSec

        ### Add information about the epochs / command waveform - we will always expect epochs to be lists.
        if "nEpochType" in self._abfHeader.header.keys():
            # ensure epochs which are just a single epoch still come out as lists
            for key in ["nEpochType", "fEpochInitLevel", "fEpochLevelInc", "lEpochInitDuration",
                        "lEpochDurationInc", "lEpochPulsePeriod", "lEpochPulseWidth", "nEpochDigitalOutput"]:
                if not type(self._abfHeader.header[key]) == list:
                    self._abfHeader.header[key] = [self._abfHeader.header[key]]
            self.epochCount = len(self._abfHeader.header['nEpochType'])
            self.epochType = self._abfHeader.header['nEpochType']
            self.epochCommand = self._abfHeader.header['fEpochInitLevel']
            self.epochCommandDelta = self._abfHeader.header['fEpochLevelInc']
            self.epochDuration = self._abfHeader.header['lEpochInitDuration']  # in points
            self.epochDurationDelta = self._abfHeader.header['lEpochDurationInc']
            self.epochPulsePeriod = self._abfHeader.header['lEpochPulsePeriod']
            self.epochPulseWidth = self._abfHeader.header['lEpochPulseWidth']
            self.epochDigOut = self._abfHeader.header['nEpochDigitalOutput']
            self.epochStartPoint = [int(self.pointsPerSweep / 64)]
            for i, duration in enumerate(self.epochDuration):
                self.epochStartPoint.append(self.epochStartPoint[-1] + duration + self.epochDurationDelta[i] * i)
            self.epochStartSec = [self.signalTimes[int(x)] for x in self.epochStartPoint]
        else:
            # this ABF has no epochs at all, so make all epoch stuff empty lists
            self.epochCount = 0
            self.epochType = []
            self.epochCommand = []
            self.epochCommandDelta = []
            self.epochDuration = []
            self.epochDurationDelta = []
            self.epochPulsePeriod = []
            self.epochPulseWidth = []
            self.epochDigOut = []
            self.epochStartSec = []
            self.epochStartPoint = []

        ### Prepare memtest placeholder
        self.memtest = self._MemtestResults(self.sweepCount)

        ### Go ahead and set sweep zero to populate command signal trace
        self.setSweep(0)

    def help(self):
        """Launch the pyABF website in a web browser."""
        import webbrowser
        webbrowser.open('https://github.com/swharden/pyABF')

    def info(self, silent=False):
        """
        Display (and return) a long message indicating what you can access
        and do with the ABF class.
        """
        functions, attributes, lists, data = [], [], [], []
        for itemName in dir(self):
            if itemName.startswith("_"):
                continue
            itemType = str(type(getattr(self, itemName))).split("'")[1]
            if itemType in ['str', 'float', 'int', 'bool']:
                attributes.append(itemName)
            elif itemType == 'list':
                lists.append(itemName)
            elif itemType == 'numpy.ndarray':
                data.append(itemName)
            elif itemType == 'method':
                functions.append(itemName)
            elif itemType in ['datetime', 'datetime.datetime']:
                continue
            else:
                print(itemType, itemName)

        msg = ""
        msg += "\n### INSTANTIATION ###\n"
        msg += "abf=pyabf.ABF(R'%s')\n" % self.filename

        msg += "\n### VALUES ###\n"
        for itemName in sorted(attributes):
            itemValue = str(getattr(self, itemName))
            msg += "* abf.%s = %s\n" % (itemName, itemValue)

        msg += "\n### LISTS ###\n"
        for itemName in sorted(lists):
            itemValue = str(getattr(self, itemName))
            msg += "* abf.%s = %s\n" % (itemName, itemValue)

        msg += "\n### SIGNAL STUFF###\n"
        for itemName in sorted(data):
            itemValue = getattr(self, itemName)
            if 'float' in str(itemValue.dtype):
                itemValue = np.array(getattr(self, itemName), dtype=np.float)
                itemValue = np.round(itemValue, decimals=5)
            msg += "* abf.%s = %s\n" % (itemName, itemValue)

        msg += "\n### FUNCTIONS ###\n"
        for itemName in sorted(functions):
            msg += "* abf.%s()\n" % (itemName)
        if not silent:
            print(msg)
        return msg

    def setSweep(self, sweepNumber=0, absoluteTime=False, channel=0):
        """
        Load a particular sweep (and times and command trace) into memory.
        This populates self.dataX, self.dataY, and self.dataC.

        If absoluteTime is used, self.dataX will reflect the time points
        (in seconds) of the sweep relative to the total experiment. If it is
        not used, self.dataX will always start at 0 (and display "sweep time").
        """
        # TODO: make function to get sweep-offset time
        # TODO: command signal not supported if using multi-channel
        if sweepNumber is None:
            return
        if sweepNumber < 0:
            sweepNumber = self.sweepList[sweepNumber]
        if not sweepNumber in self.sweepList:
            raise ValueError("Sweep %d not found (last sweep is %d)" % (sweepNumber, self.sweepList[-1]))
        self.dataSweepSelected = sweepNumber
        self.sweepSelected = sweepNumber
        pointStart = sweepNumber * (self.pointsPerSweep * self.dataChannels)
        pointEnd = pointStart + (self.pointsPerSweep * self.dataChannels)
        self.dataY = self.signalData[int(pointStart):int(pointEnd)]
        if absoluteTime:
            self.dataX = self.signalTimes[int(pointStart):int(pointEnd)]
        else:
            self.dataX = self.signalTimes[0:int(self.pointsPerSweep)]
        if self.dataChannels > 1:
            self.dataY = self.dataY[channel::self.dataChannels] * self.dataChannels
        self._updateCommandWaveform()
        if self.gaussianSigma:
            self.dataY = self._filterGaussian(self.dataY)

    def _updateCommandWaveform(self):
        """Read the epochs and populate self.dataC with the command signal."""
        # TODO: don't update if the command doesn't change from sweep to sweep
        self.dataC = np.empty(self.dataX.size)  # start as random data
        position = 0  # start at zero here for clarity
        position += int(self.pointsPerSweep / 64)  # the first 1/64th is pre-epoch (why???)
        self.dataC[:position] = self.commandHold  # fill the pre-epoch with the command holding
        for epochNumber in range(self.epochCount):
            if self.epochType[epochNumber] == 0:
                continue  # it's a disabled epoch
            if epochNumber >= len(self.epochDuration):
                print("ran out of epoch")
                break  # ran out?

            pointCount = self.epochDuration[epochNumber]
            deltaCommand = self.epochCommandDelta[epochNumber] * self.sweepSelected

            if self.epochType[epochNumber] == 1:  # STEP
                self.dataC[position:position + pointCount] = self.epochCommand[epochNumber] + deltaCommand
            elif self.epochType[epochNumber] == 2:  # RAMP
                ramp = np.arange(pointCount)*1. / pointCount
                rampStart = self.dataC[position - 1]
                rampEnd = self.epochCommand[epochNumber] + deltaCommand
                rampDiff = rampEnd - rampStart
                print (ramp,rampDiff)
                ramp *= rampDiff
                ramp += rampStart
                self.dataC[position:position + pointCount] = ramp
            else:  # UNKNOWN (so treat it like a step)
                warnings.warn("I don't know how to analyze an epoch of type %d" % self.epochType[epochNumber])
                self.dataC[position:position + pointCount] = self.epochCommand[epochNumber] + deltaCommand
            position += pointCount
        self.dataC[position:] = self.commandHold  # set the post-epoch to the command holding

    ### FILTERING

    gaussianSigma = 0
    gaussianLeft = False
    gaussianRight = True

    def _filterGaussian(self, signal, sigma=None):
        """Perform gaussian smoothing on a 1d array."""
        if self.gaussianLeft and self.gaussianRight:
            raise ValueError("can't set both gaussianLeft and gaussianRight")
        if sigma is None:
            sigma = self.gaussianSigma
        if self.gaussianLeft or self.gaussianRight:
            sigma *= 2
        size = sigma * 10
        points = np.exp(-np.power(np.arange(size) - size / 2, 2) / (2 * np.power(sigma, 2)))
        if self.gaussianLeft:
            points[-int(len(points) / 2):] = 0
            points = points[:-int(sigma)]  # offset a sigma
        if self.gaussianRight:
            points[0:int(len(points) / 2)] = 0
            points = points[int(sigma):]  # offset a sigma
        kernel = points / sum(points)
        return np.convolve(signal, kernel, mode='same')

    def _filterSimple(self, signal, sigma):
        """return a Gaussian filtered signal."""
        size = sigma * 10
        points = np.exp(-np.power(np.arange(size) - size / 2, 2) / (2 * np.power(sigma, 2)))
        kernel = points / sum(points)
        return np.convolve(signal, kernel, mode='same')

    def epochData(self, epochNumber=0):
        """Return [dataX, dataY, dataC] of a given epoch of the current sweep."""
        I1 = self.epochStartPoint[epochNumber]
        if len(self.epochStartPoint) - 1 <= epochNumber:
            I2 = -1
        else:
            I2 = self.epochStartPoint[epochNumber + 1]
        Xs = self.dataX[I1:I2]
        Ys = self.dataY[I1:I2]
        Cs = self.dataC[I1:I2]
        return Xs, Ys, Cs

    ### ANALYSIS

    def _monoExpTau(self, data, sample_rate_hz=20000, tau=.1, step=.1):
        """Given some data which decays to zero, return its time constant."""
        if len(data) == 0:
            return np.nan
        errs = [np.inf]
        normed = data / data[0]
        Xs = np.arange(len(normed)) / sample_rate_hz
        while (len(errs)) < 50:
            assert len(Xs) == len(data)
            tau = np.max((0.000001, tau))
            errs.append(np.sum(np.exp(-Xs / tau) - normed))
            if np.abs(errs[-1]) < 0.01:
                return tau
            if (errs[-1] > 0 and errs[-2] < 0) or (errs[-1] < 0 and errs[-2] > 0):
                step *= .6
            if errs[-1] < 0:
                tau += step
            elif errs[-1] > 0:
                tau -= step
        return tau

    def rms(self, chunkMS=10, quietestMS=100):
        """
        Return the RMS value of the noise floor. RMS = stdev when mean is 0.
        The noise floor is defined as the quietest parts of the signal.
        """
        chunkSize = chunkMS * self.pointsPerMS
        chunkCount = int(len(self.dataY) / chunkSize)
        stdev = np.empty(chunkCount)
        for chunkNumber in range(chunkCount):
            i1 = int(chunkSize * chunkNumber)
            i2 = int(chunkSize * (chunkNumber + 1))
            stdev[chunkNumber] = np.std(self.dataY[i1:i2])
        countToAverage = int(quietestMS / chunkMS)
        rms = np.mean(sorted(stdev)[:countToAverage])
        return rms

    def _tonic(self, data, binSize=.1, fitAboveFrac=.25):
        """
        Return a polynomial-fitted peak of the histogram of the given data.
        Only histogram points above "fitAboveFrac" are fitted.
        """
        padSize = int(200 / binSize) * 2
        pad = np.arange(padSize) * binSize
        bins = np.concatenate((data[0] - pad[::-1], data[0] + pad))
        histCount, histBins = np.histogram(data, bins=bins)
        histBins = histBins[:-1]
        validIs = histCount > np.max(histCount) * fitAboveFrac
        histBins = histBins[validIs]
        histCount = histCount[validIs]
        fit = np.poly1d(np.polyfit(histBins, histCount, 6))
        histFitVals = fit(histBins)
        histFitPeakVal = np.max(histFitVals)
        histFitPeakI = np.where(histFitVals == histFitPeakVal)[0]
        tonicValue = histBins[histFitPeakI]
        return float(tonicValue)

    def tonicPhasic(self, t1=0, t2=None):
        """
        Return [tonic, phasicNeg, and phasicPos] of the selected sweep.
        Tonic is the peak value of the polynomial-fitted histogram.
        Phasic is the average value of all points below or above the tonic
        value. All 3 are in abf.units units.
        """
        if not t2:
            t2 = self.sweepLengthSec
        i1, i2 = int(t1 * self.pointsPerSec), int(t2 * self.pointsPerSec)
        data = self.dataY[i1:i2]
        tonicValue = self._tonic(data)
        phasicNeg = tonicValue - np.average(data[data < tonicValue])
        phasicPos = np.average(data[data > tonicValue]) - tonicValue
        return [tonicValue, phasicNeg, phasicPos]

    def average(self, t1=0, t2=None, setSweep=False, tonic=False):
        """
        Return the average value of the selected sweep. If t1 and/or t2 are
        given, data between those times (seconds) will be measured. If tonic
        is True, the tonic (histogram fit) method will be used instead of
        a true average.
        """
        if not setSweep is False:
            self.setSweep(setSweep)
        if not t2:
            t2 = self.sweepLengthSec
        i1 = int(t1 * self.pointsPerSec)
        i2 = int(t2 * self.pointsPerSec)
        if tonic:
            return self._tonic(self.dataY[i1:i2])
        else:
            return np.nanmean(self.dataY[i1:i2])

    def averageEpoch(self, epoch, firstFrac=False, lastFrac=False, setSweep=False,
                     tonic=False):
        """
        Return the average of some fraction of an epoch. Similar to average(),
        but uses an epoch number and fraction rather than time points.

        Example: return the average of the last 25% of epoch B
            >>> abf.averageEpoch(epoch=1,lastFrac=.25)

        """
        if not setSweep is False:
            self.setSweep(setSweep)
        if firstFrac and lastFrac:
            raise ValueError("can't set both a first and last fraction")
        epochs = self.epochStartPoint + [self.pointsPerSweep]
        if epoch < (len(self.epochStartPoint) - 1):
            i2 = epochs[epoch + 1]
            i1 = i2 - self.epochDuration[epoch]
        else:
            i2 = self.pointsPerSweep
            i1 = epochs[epoch]
        dur = i2 - i1
        if firstFrac:
            i2 = i1 + dur * firstFrac
        if lastFrac:
            i1 = i2 - dur * lastFrac
        if tonic:
            return self._tonic(self.dataY[int(i1):int(i2)])
        else:
            return np.nanmean(self.dataY[int(i1):int(i2)])

    def stdev(self, t1=0, t2=None):
        """
        Return the standard deviation of current sweep between two times (sec)
        """
        if not t2:
            t2 = self.sweepLengthSec
        i1 = int(t1 * self.pointsPerSec)
        i2 = int(t2 * self.pointsPerSec)
        return np.nanstd(self.dataY[i1:i2])

    def stderr(self, t1=0, t2=None):
        """
        Return the standard error of current sweep between two times (sec)
        """
        if not t2:
            t2 = self.sweepLengthSec
        i1 = int(t1 * self.pointsPerSec)
        i2 = int(t2 * self.pointsPerSec)
        return np.nanstd(self.dataY[i1:i2]) / np.math.sqrt(i2 - i1)

    #    def sweepSpan(self,t1=0,t2=None):
    #        """
    #        Return just the dataY between two time points (sec)
    #        """
    #        if not t2:
    #            t2=self.sweepLengthSec
    #        i1=int(t1*self.pointsPerSec)
    #        i2=int(t2*self.pointsPerSec)
    #        return self.dataY[i1:i2]

    ### EVENT DETECTION


    def eventsDeriv(self, setSweep=None, t1=0, t2=None, dTms=1, threshold=-10,
                    alignToDerivPeak=True, alignToRawPeak=False, plot=False,
                    mustRevertMS=False):
        """
        Derivative-threshold-based event detection. Return a list of times for this
        sweep where the first derivative threshold was exceeded. This event
        detection routine is minimal and simplistic, and can be used for APs,
        EPSCs, and IPSCs.

        You may want to enable gaussian filtering before calling this function.

        setSweep:
            which sweep to use

        t1 and t2:
            time range (in seconds) to perform event detection

        dTms and threshold:
            Events are points where the threshold is exceeded over the change in
            time in milliseconds (dTms). Threshold can be positive or negative.

        alignToDerivPeak and alignToRawPeak:
            If disabled, the event times will be the first point where the
            derivative was first crossed. If enabled, the event times will be
            aligned to the peak derivative (rather than the threshold) or the
            peak of the raw trace.

        Common settings:
            Action potential detection: dTms = 1, threshold = 10, mustRevertMS = 5

        """

        # determine detection details
        if not setSweep is None:
            self.setSweep(setSweep)
        if not t2:
            t2 = self.sweepLengthSec
        i1, i2 = int(t1 * self.pointsPerSec), int(t2 * self.pointsPerSec)

        # load the data and calculate its derivative
        strip = self.dataY[i1:i2]
        Xs = self.dataX[i1:i2]
        dT = int(self.pointsPerMS * dTms)
        deriv = (strip[dT:] - strip[:-dT]) / dTms
        deriv = np.concatenate((deriv, [deriv[-1]] * dT))

        # find first-crossings of points where the derivative was crossed
        if threshold > 0:
            crossed = deriv > threshold
        else:
            crossed = deriv < threshold
        crossed[1:][crossed[:-1] & crossed[1:]] = False
        eventIs = np.where(crossed)[0]  # +dT

        # remove events which are less than one dT together
        for i in range(1, len(eventIs)):
            if eventIs[i] - eventIs[i - 1] <= dT:
                eventIs[i] = False
        eventIs = eventIs[np.where(eventIs)[0]]

        # optionally align to the peak of the first derivative
        if alignToDerivPeak:
            for n, i in enumerate(eventIs):
                if threshold > 0:
                    while deriv[i] > deriv[i - 1]:
                        i += 1
                    eventIs[n] = i - 1
                else:
                    while deriv[i] < deriv[i - 1]:
                        i += 1
                    eventIs[n] = i - 1

        # optionally align to the peak of the raw trace
        if alignToRawPeak:
            for n, i in enumerate(eventIs):
                if threshold > 0:
                    while strip[i] > strip[i - 1]:
                        i += 1
                    eventIs[n] = i
                else:
                    while strip[i] < strip[i - 1]:
                        i += 1
                    eventIs[n] = i

        if mustRevertMS:
            revertPoints = int(self.pointsPerMS * mustRevertMS)
            for n, i in enumerate(eventIs):
                if threshold > 0:
                    if not np.min(deriv[i:i + revertPoints]) < 0:
                        eventIs[n] = False
                else:
                    if not np.max(deriv[i:i + revertPoints]) > 0:
                        eventIs[n] = False
            eventIs = eventIs[np.where(eventIs)[0]]

        eventIs = np.unique(eventIs)
        eventTimes = Xs[eventIs]

        if plot:
            plt.figure()
            ax1 = plt.subplot(211)
            plt.title("sweep %d (raw signal, %s)" % (self.sweepSelected, self.units))
            plt.plot(Xs, strip)
            for eventI in [int(self.pointsPerSec * x) for x in eventTimes]:
                plt.axvline(self.dataX[eventI], color='r', alpha=.5)
            plt.subplot(212, sharex=ax1)
            plt.title("first derivative (%s / ms)" % self.units)
            plt.plot(Xs, deriv)
            plt.axhline(threshold, color='r', alpha=.5)
            for eventI in [int(self.pointsPerSec * x) for x in eventTimes]:
                plt.axvline(self.dataX[eventI], color='r', alpha=.5)
            plt.margins(0, .1)
            plt.tight_layout()
            plt.show()

        return eventTimes

    ### PLOTTING

    def plotEpochs(self):
        """
        Display a matplotlib figure of the current sweep highlighting epochs.
        """
        epochBoundsSec = self.epochStartSec + [self.sweepLengthSec]
        command = self.epochCommand + [self.commandHold]
        plt.plot(self.dataX, self.dataY, 'k-')
        for i in range(len(self.epochStartSec)):
            plt.axvspan(epochBoundsSec[i], epochBoundsSec[i + 1], alpha=.3, color=self._sweepColor(i),
                        lw=0, label=("%s: (%s %s)" % (chr(65 + i), command[i], self.unitsCommand)))
        plt.legend(fontsize=8)
        self.plotDecorate(zoomYstdev=True)

    colormap = "jet_r"

    def _sweepColor(self, sweep=None):
        """
        Return a colorcode of the current sweep using abf.colormap
        """
        if sweep is None:
            sweep = self.dataSweepSelected
        frac = sweep / self.sweepCount
        return plt.cm.get_cmap(self.colormap)(frac)

    def plotSweeps(self, sweeps=None, offsetX=0, offsetY=0, useColormap=False,
                   color='b'):
        """
        Plot signal data using matplotlib.If sweeps is a list of numbers,
        only plot those sweeps. Also accepts an integer to plot one sweep.
        """
        if type(sweeps) == list:
            pass
        elif sweeps == None or sweeps == False:
            sweeps = self.sweepList
        else:
            sweeps = [sweeps]

        for sweepNumber in sweeps:
            self.setSweep(sweepNumber)
            if useColormap:
                color = self._sweepColor()
            plt.plot(self.dataX + offsetX * sweepNumber,
                     self.dataY + offsetY * sweepNumber,
                     color=color)

        plt.margins(0, .1)

    def plotDecorate(self, command=False, title=True, xlabel=True, ylabel=True,
                     zoomYstdev=False, legend=False, axis=None):
        """
        Add axis labels and a title to the already-drawn matplotlib figure.

        Arguments:
            * zoomYstdev - optionally zoom in vertically to a certain number
                           of standard deviations from the mean.
            * axis - if provided, resize to this axis [X1,X2,Y1,Y2]
        """

        # title
        if title is True:
            plt.title(self.ID, fontsize=16)
        elif title:
            plt.title(str(title), fontsize=16)

        # x label
        if xlabel is True:
            plt.xlabel(self.unitsTimeLong)
        elif xlabel:
            plt.xlabel(str(xlabel))

        # y label
        if ylabel is True:
            if command:
                plt.ylabel(self.unitsCommandLong)
            else:
                plt.ylabel(self.unitsLong)
        elif ylabel:
            plt.ylabel(str(ylabel))

        plt.margins(0, .1)

        if zoomYstdev:
            if zoomYstdev is True:
                zoomYstdev = 3
            else:
                zoomYstdev = int(zoomYstdev)
            av = np.nanmean(self.dataY)
            stdev = np.nanstd(self.dataY)
            plt.axis([None, None, av - stdev * zoomYstdev, av + stdev * zoomYstdev])

        if legend:
            if type(legend) is int:
                plt.legend(legend)
            else:
                plt.legend()

        if axis:
            plt.axis(axis)

        plt.tight_layout()

    ### MEMBRANE TEST

    class _MemtestResults:

        class _MemtestItem:
            def __init__(self, name="Ih", units="pA", desc="Membrane Current", sweepCount=123):
                """A single membrane test item, such as holding current (numpy ndarray)"""
                self.name = name
                self.units = units
                self.desc = desc
                self.values = np.empty(sweepCount) * np.nan
                self.analyzed = False

            @property
            def average(self):
                return np.nanmean(self.values)

            @property
            def longLabel(self):
                return "%s (%s)" % (self.desc, self.units)

            @property
            def label(self):
                return "%s (%s)" % (self.name, self.units)

            def __setitem__(self, key, value):
                self.analyzed = True
                self.values[key] = value

            def __getitem__(self, key):
                if not self.analyzed:
                    warnings.warn("memtest needs to run prior to its values being accessed!")
                return self.values[key]

            def __len__(self):
                return len(self.values)

            def __repr__(self):
                return repr(self.values)

        def __init__(self, sweepCount):
            """this object stores membrane test data and common methods."""
            self.sweepCount = sweepCount
            self.Ih = self._MemtestItem("Ih", "pA", "Membrane Current", sweepCount)
            self.Rm = self._MemtestItem("Rm", "MOhm", "Membrane Resistance", sweepCount)
            self.Ra = self._MemtestItem("Ra", "MOhm", "Access Resistance", sweepCount)
            self.Cm = self._MemtestItem("Cm", "pF", "Membrane Capacitance", sweepCount)
            self.Tau = self._MemtestItem("Tau", "ms", "VC Time Constant", sweepCount)

    def _memtestThisSweep(self, tonicAnalysis=True, avgLastFrac=.75):
        """
        When called on a sweep with a voltage step in it, return a dict with membrane properties.
        Keys will have names like "Rm", "Ra", "Cm", "Ih", and "tau".
        See the cookbook to gain insight into the calculations used here.

        Ideal configuration 1 (holding: -70 mV):
            epoch A: -80 mV
            epoch B: anything (or absent)

        Ideal configuration 2 (holding: -70 mV):
            epoch A: -80 mV
            epoch B: -70 mV
            epoch C: anything (or absent)

        Ideal configuration 3 (holding: -70 mV):
            epoch A: -70 mV
            epoch B: -80 mV
            epoch C: -70 mV
            epoch D: anything (or absent)

        This method assumes a few things about your file:
            * no command deltas or time deltas are used
            * the first or second epoch is a voltage step
            * if an epoch before the step exists, its value is the same as the command current
            * the epoch after the step returns back to the command current
        """
        if not self.units == "pA":
            raise ValueError("memtest should only be run on VC traces")

        # we will calculate memtest based on two traces, so figure them out based on command steps
        if self.epochCommand[0] != self.commandHold:
            # Epoch A is the step
            dV = np.abs(self.epochCommand[0] - self.commandHold) * 1e-3
            if len(self.epochCommand) > 1:
                # Epoch A and B will be analyzed
                trace1 = self.dataY[self.epochStartPoint[0]:self.epochStartPoint[0] + self.epochDuration[0]]
                trace2 = self.dataY[self.epochStartPoint[1]:self.epochStartPoint[1] + self.epochDuration[1]]
            else:
                # Epoch A will be analyzed, B is whatever is left
                trace1 = self.dataY[self.epochStartPoint[0]:self.epochStartPoint[0] + self.epochDuration[0]]
                trace2 = self.dataY[self.epochStartPoint[0] + self.epochDuration[0]:]
        elif self.epochCommand[0] != self.epochCommand[1]:
            # Epoch A is the step
            dV = np.abs(self.epochCommand[1] - self.epochCommand[0]) * 1e-3
            if len(self.epochCommand) > 2:
                # Epoch B and C will be analyzed
                trace1 = self.dataY[self.epochStartPoint[1]:self.epochStartPoint[1] + self.epochDuration[1]]
                trace2 = self.dataY[self.epochStartPoint[2]:self.epochStartPoint[2] + self.epochDuration[2]]
            else:
                # Epoch B will be analyzed, C is whatever is left
                trace1 = self.dataY[self.epochStartPoint[1]:self.epochStartPoint[1] + self.epochDuration[1]]
                trace2 = self.dataY[self.epochStartPoint[1] + self.epochDuration[1]:]
        else:
            raise ValueError("A step memtest cannot be used on this type of ABF.")

        # this memtest dictionary is what gets returned
        mt_dict = {"dV": dV * 1e3}

        # subtract-out the steady state current so signals are centered at 0
        if tonicAnalysis:
            Ih1 = self._tonic(trace1[len(trace1) - int(avgLastFrac * len(trace1)):])
            Ih2 = self._tonic(trace2[len(trace2) - int(avgLastFrac * len(trace2)):])
        else:
            Ih1 = np.average(trace1[len(trace1) - int(avgLastFrac * len(trace1)):])
            Ih2 = np.average(trace2[len(trace2) - int(avgLastFrac * len(trace2)):])
        data1 = trace1 - Ih1
        data2 = trace2 - Ih2
        mt_dict["Ih"] = Ih2
        mt_dict["Vm"] = self.commandHold

        # Rm - compare the steady state currents to calculate membrane resistance
        dI = (np.abs(Ih2 - Ih1) * 1e-12)
        Rm = dV / dI  # Rm = dV/dI
        mt_dict["Rm"] = Rm * 1e-6

        # let's improve out data by averaging the two curves together
        point_count = np.min((len(trace1), len(trace2)))
        data = np.average((-data1[:point_count], data2[:point_count]), axis=0)

        # Find the points of the trace we intend to fit
        peakI = np.where(data == np.max(data))[0][0]
        zeroI = np.where(data[peakI:] <= 0)[0]
        if len(zeroI) == 0:
            zeroI = peakI
        else:
            zeroI = zeroI[0] + peakI

        # Fit the curve to a monoexponential equation and record tau
        tau = self._monoExpTau(data[peakI:zeroI])
        mt_dict["Tau"] = tau * 1e3

        # use tau to guess what I0 probably was at the first point after the step
        I0 = np.exp((peakI / self.pointsPerSec) / tau) * data[peakI] * 1e-12
        mt_dict["I0"] = I0 * 1e12

        # calculate Ra=dV/I0
        Ra = dV / I0
        mt_dict["Ra"] = Ra * 1e-6

        # calculate Cm=tau/Ra
        Cm = tau / Ra
        mt_dict["Cm"] = Cm * 1e12

        # populate the memtest object
        self.memtest.Ih[self.sweepSelected] = mt_dict["Ih"]
        self.memtest.Rm[self.sweepSelected] = mt_dict["Rm"]
        self.memtest.Ra[self.sweepSelected] = mt_dict["Ra"]
        self.memtest.Cm[self.sweepSelected] = mt_dict["Cm"]
        self.memtest.Tau[self.sweepSelected] = mt_dict["Tau"]

        # optionally return memtest dictionary with full details
        return mt_dict

    def memtestAnalyzeAll(self):
        """Determine membrane test properties for every sweep. Assign results to self.memtestResults"""
        for sweep in self.sweepList:
            self.setSweep(sweep)
            self._memtestThisSweep()


def _listDemoFiles(silent=False):
    """List all ABF files in the ../../data/ folder."""
    import glob
    fnames = []
    for fname in sorted(glob.glob("../../data/*.abf")):
        fnames.append(fname)
    if not silent:
        print("\n".join(fnames))
    return fnames


def _checkFirstPoint():
    """Display the first value for each ABF. A good way to ensure scaleFactor is working."""
    for fname in _listDemoFiles(silent=True):
        abf = ABF(fname)
        print(abf.ID, abf.dataY[0])


if __name__ == "__main__":
    print("do not run this script directly.")
    # _listDemoFiles()
    # _checkFirstPoint()

    # abf=ABF(R"../../data/171116sh_0014.abf") # V memtest
    # abf=ABF(R"../../data/171116sh_0019.abf") # IC steps
    # abf=ABF(R"../../data/171116sh_0011.abf") # step memtest
    # abf=ABF(R"../../data/16d05007_vc_tags.abf") # time course experiment
    # abf=ABF(R"../../data/16d22006_kim_gapfree.abf") # gap-free dual-channel file


    print("DONE")