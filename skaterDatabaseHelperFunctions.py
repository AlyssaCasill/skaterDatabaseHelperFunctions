"""
Date: Jan 7, 2019
Written for the Gamble Lab @ Albert Einstein College of Medicine

Purpose: Contains all functionality to interact with analysisID and skater databases

Repository creation: August 2, 2020
"""
__author__ =  'Alyssa D. Casill'
__version__ = '0.0.3'
__email__ = 'alyssa.casill@phd.einstein.yu.edu'

import math

def updateLine(toUpdate,index,newValue):
    """
    Function used to update specified index of 'toUpdate' with 'newValue,' used in analysis and skater databases to set attributes
    """
    if type(newValue) == list:
        templist = toUpdate.split("\t")
        templist[index] = ",".join(newValue)
    else:
        templist = toUpdate.split("\t")
        templist[index] = newValue
    newline = "\t".join([str(i) for i in templist])
    return(newline)

class AnalysisIDdata:
    """
    A class with functions to read and write lines of the analysisID database files created during SKaTERseq optimization. 
    
    The analysisID database files contain metadata describing the parameters of the analysis run. This file is most useful when testing/optimizing various parameters. All fields are typically included. Could eventually be phased out or reformatted. 
    Analysis ID fields:
        analysisID: unique identifer of this particular analysis. Typically includes the name of the experiment, the condition, and the replicates
        softwareVersion: Version of SKaTER optimize software used
        experimentName: Condition being tested (cell line, treatment)
        replicate: rep1, rep2, etc
        pathToData: path to location of database **THIS DOES NOT ACCOUNT FOR MOVING OF DATABASES TO NEW LOCATIONS**
        settingsName: name of settings file used to optimize data
        globalSettings_settingsFileInfo: comma separated list of all arguments included in settings file (except BAM files)
        BAMfiles: path to BAM files from settings file location
    """
    def __init__(self,line=""):
        """
        A class used to construct an AnalysisIDdata instance

        Parameters
		----------
		line : str
			A line of text used to construct the AnalysisIDdata. If empty, a string with the correct number of tabs will be created.
		
		Returns
		----------
		AnalysisIDdata instance
        """
        if line == "":
            self.line = "\t\t\t\t\t\t\t\t"
        else:
            self.line = line.strip()
            if len(self.line.split("\t")) != 9:
                print("\n".join(line.split("\t")))
                raise Exception("Analysis ID line does not match standard format")
    
    @property
    def analysisID(self):
        return(self.line.split("\t")[0])
    @analysisID.setter
    def analysisID(self,newValue):
        self.line = updateLine(self.line,0,newValue)
        return(self.line)

    @property
    def softwareVersion(self):
        return(self.line.split("\t")[1])
    @softwareVersion.setter
    def softwareVersion(self,newValue):
        self.line = updateLine(self.line,1,newValue)
        return(self.line)

    @property
    def experimentName(self):
        return(self.line.split("\t")[2])
    @experimentName.setter
    def experimentNawme(self,newValue):
        self.line = updateLine(self.line,2,newValue)
        return(self.line)

    @property
    def replicate(self):
        return(self.line.split("\t")[3])
    @replicate.setter
    def replicate(self,newValue):
        self.line = updateLine(self.line,3,newValue)
        return(self.line)

    @property
    def pathToData(self):
        return(self.line.split("\t")[4])
    @pathToData.setter
    def pathToData(self,newValue):
        self.line = updateLine(self.line,4,newValue)
        return(self.line)

    @property
    def settingsPath(self):
        return(self.line.split("\t")[5])
    @settingsPath.setter
    def settingsPath(self,newValue):
        self.line = updateLine(self.line,5,newValue)
        return(self.line)

    @property
    def globalSettings(self):
        return(self.line.split("\t")[6])
    @globalSettings.setter
    def globalSettings(self,newValue):
        self.line = updateLine(self.line,6,newValue)
        return(self.line)

    @property
    def BAMfile(self):
        return(self.line.split("\t")[7])
    @BAMfile.setter
    def BAMfile(self,newValue):
        self.line = updateLine(self.line,7,newValue)
        return(self.line)

    @property
    def geneDatabasePath(self):
        return(self.line.split("\t")[8])
    @geneDatabasePath.setter
    def geneDatabasePath(self,newValue):
        self.line = updateLine(self.line,8,newValue)
        return(self.line)

    def __str__(self):
        """
        Returns a formated line for a skater database file.
		"""
        return(self.line)

def addNewAnalysis(settingsFilePath,version,experiment,replicate,analysisID,geneDatabasePath):
    """ 
    Function to add new AnalysisIDdata to analysisID database.

    Parameters
    ----------
    settingsFilePath: path to settings file
    version: version of skaterOptimize software used
    experiment: Condition being tested (cell line, treatment)
    replicate: rep1, rep2, etc
    analysisID: unique identifer of this particular analysis. Typically includes the name of the experiment, the condition, and the replicates
    geneDatabasePath: path to SKaTER database file
    
    Returns
    ----------
    Line (str) of analysisID database file (can be used to write/append to analysisID database)
    """
    newAnalysis = AnalysisIDdata()
    newAnalysis.analysisID = analysisID
    newAnalysis.softwareVersion = version
    newAnalysis.experimentName = experiment
    newAnalysis.replicate = replicate
    newAnalysis.pathToData = "/".join(settingsFilePath.split("/")[:-1])
    newAnalysis.settingsPath = settingsFilePath
    bamfiles = []
    globalsettings = []
    settings = open(settingsFilePath)
    for line in settings:
        temp = line.strip().split("=") 
        if temp[0] == "gene":
            pass
        elif temp[0] == "--bamFileNames":
            bamfiles.append(temp[1])
        else:
            try:
                globalsettings.append(":".join([temp[0],temp[1]]))
            except IndexError:
                globalsettings.append(temp[0])
    newAnalysis.BAMfile = bamfiles
    newAnalysis.globalSettings = globalsettings
    newAnalysis.geneDatabasePath = geneDatabasePath
    return(str(newAnalysis))

class RunLevelData: 
    """
    A class with functions to read and write lines of a *sorted* skater database file, and to perform functions on data output for a single optimization run. 
    
    The skater database file contains data from each optimization run. One line per run. 
    skater database fields:
        geneID: gene name (database is sorted by gene name)
        analysisID: matches analysis ID in analysisID database (experiment/condition/rep)
        runID: unique identifier for this run (typically, job number from HPC)
        chrom: chromosome ("chr1")
        strand: plus or minus strand of gene
        RSS: residiual sum of squares from best solution
        DRBreleaseRate: DRB release rate
        tssPos: genomic coordinate of transcription start site, based on refflat data (will be smaller than transcription end site if gene is positive strand, will be larger than transcription end site if gene is negative strand)
        spawnRates: comma separated list of spawn rates
        speedZonePos: comma separated list of "positions" of speed zones where elongation rates are measured. If the "speedZoneType" is time, these will be the time point windows (eg 5min-10min, 10min-15min)
        speedZoneRates: comma separate list of elongation rates in each speed zone (in same order of speed zones)
        introns: comma separated, str(tuple) positions of each intron
        splicingRates: comma separated list of splicing rates (in same order as introns)
        tesPos: genomic coordinate of transcription end site, based on refflat data
        cleavageRates: comma seprated list of cleavage rates
        contamination: contamination level
        enOverIn: nascent exon over nascent intron levels 
        generation: generation at which final solution was produced
    """
    def __init__(self,line=""):
        """
        A class used to construct a RunLevelData instance, and to perform functions on an instance of RunLevelData

        Parameters
		----------
		line : str
			A line of text used to construct the RunLevelData. If empty, a string with the correct number of tabs will be created.
		
		Returns
		----------
		RunLevelData instance
        """
        if line == "":
            self.line = "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t"
        else:
            self.line = line.strip()
            if len(self.line.split("\t")) != 18:
                print(self.line)
                raise Exception("Database line does not match standard format")

    @property
    def geneID(self):
        """
		Returns a string of the gene name for this run.
		"""
        return(self.line.split("\t")[0])
    @geneID.setter
    def geneID(self,newValue):
        self.line = updateLine(self.line,0,newValue)
        return(self.line)

    @property
    def analysisID(self):
        """
		Returns a string of the analysis ID for this run.
		"""
        return(self.line.split("\t")[1])
    @analysisID.setter
    def analysisID(self,newValue):
        self.line = updateLine(self.line,1,newValue)
        return(self.line)
    
    @property
    def runID(self):
        """
		Returns a string of the run ID.
		"""
        return(self.line.split("\t")[2])
    @runID.setter
    def runID(self,newValue):
        self.line = updateLine(self.line,2,newValue)
        return(self.line)
    
    @property
    def chrom(self):
        """
		Returns a string of the chromosome (chr1, chr2, etc)
		"""
        return(self.line.split("\t")[3])
    @chrom.setter
    def chrom(self,newValue):
        self.line = updateLine(self.line,3,newValue)
        return(self.line)
    
    @property
    def strand(self):
        """
		Returns a string of the strand of this gene (+ or -)
		"""
        return(self.line.split("\t")[4])
    @strand.setter
    def strand(self,newValue):
        self.line = updateLine(self.line,4,newValue)
        return(self.line)
    
    @property
    def RSS(self):
        """
		Returns a string of the RSS for this run. 
        TO UPDATE: This should be a float?
		"""
        return(self.line.split("\t")[5])
    @RSS.setter
    def RSS(self,newValue):
        self.line = updateLine(self.line,5,newValue)
        return(self.line)

    @property
    def DRBreleaseRate(self):
        """
		Returns a float of the DRB release rate. 
		"""
        return(float(self.line.split("\t")[6]))
    @DRBreleaseRate.setter
    def DRBreleaseRate(self,newValue):
        self.line = updateLine(self.line,6,newValue)
        return(self.line)

    @property
    def tssPos(self):
        """
		Returns a list of integers representing the genomic positions of the transcription start sites.
        Currently, this will be a list of one TSS, but is a list in anticipation of solving genes with multiple TSS.
		"""
        return([int(i) for i in self.line.split("\t")[7].split(",")])
    @tssPos.setter
    def tssPos(self,newValue):
        self.line = updateLine(self.line,7,newValue)
        return(self.line)

    @property
    def spawnRates(self):
        """
		Returns a list of floats representing the spawn rates.
        Currently, this will be a list of one spawn rate, but is a list in anticipation of solving genes with multiple TSS.
		"""
        return([float(i) for i in self.line.split("\t")[8].split(",")])
    @spawnRates.setter
    def spawnRates(self,newValue):
        self.line = updateLine(self.line,8,newValue)
        return(self.line)

    @property
    def speedZonePos(self):
        """
		Returns a list of strings of speed zone positions: '(0 - 5)', '(5 - 10)'
        See speedZoneTuples() to convert string to named tuple of integers representing the positions
		"""
        return(self.line.split("\t")[9].split(","))
    @speedZonePos.setter
    def speedZonePos(self,newValue):
        self.line = updateLine(self.line,9,newValue)
        return(self.line)

    @property
    def speedZoneRates(self):
        """
		Returns a list of floats representing the elongation rates in each speed zone, in the same order as the speed zones. 
		"""
        return([float(i) for i in self.line.split("\t")[10].split(",")])
    @speedZoneRates.setter
    def speedZoneRates(self,newValue):
        self.line = updateLine(self.line,10,newValue)
        return(self.line)

    @property
    def introns(self):
        """
		Returns a list of strings of genomic positions for each intron: '(52582305 - 52587890)'
        See intronTuples() to convert string to named tuple of integers representing the positions
		"""
        return(self.line.split("\t")[11].split(","))
    @introns.setter
    def introns(self,newValue):
        self.line = updateLine(self.line,11,newValue)
        return(self.line)

    @property
    def splicingRates(self):
        """
		Returns a list of floats representing the splicing rates of each intron, in the same order as the intron positions. If the gene has no introns, the function will return NAN.
        TO UPDATE: if the gene doesn't have introns, should this just return False instead of NAN? 
		"""
        import math
        try:
            return([float(i) for i in self.line.split("\t")[12].split(",")])
        except ValueError:
            return([math.nan])
    @splicingRates.setter
    def splicingRates(self,newValue):
        self.line = updateLine(self.line,12,newValue)
        return(self.line)

    @property
    def tesPos(self):
        """
		Returns a list of integers representing the genomic positions of the transcription end sites.
        Currently, this will be a list of one TES, but is a list in anticipation of solving genes with multiple TES.
		"""
        return([int(i) for i in self.line.split("\t")[13].split(",")])
    @tesPos.setter
    def tesPos(self,newValue):
        self.line = updateLine(self.line,13,newValue)
        return(self.line)

    @property
    def cleavageRates(self):
        """
		Returns a list of floats representing the cleavage rates.
        Currently, this will be a list of one cleavage rate, but is a list in anticipation of solving genes with multiple TES.
		"""
        return([float(i) for i in self.line.split("\t")[14].split(",")])
    @cleavageRates.setter
    def cleavageRates(self,newValue):
        self.line = updateLine(self.line,14,newValue)
        return(self.line)

    @property
    def contamination(self):
        """
		Returns a float representing the contamination value.
		"""
        return(float(self.line.split("\t")[15]))
    @contamination.setter
    def contamination(self,newValue):
        self.line = updateLine(self.line,15,newValue)
        return(self.line)

    @property
    def enOverIn(self):
        """
		Returns a float representing the exon/intron value.
		"""
        try:
            return(float(self.line.split("\t")[16]))
        except TypeError:
            return(self.line.split("\t")[16])
    @enOverIn.setter
    def enOverIn(self,newValue):
        self.line = updateLine(self.line,16,newValue)
        return(self.line)

    @property
    def generation(self):
        """
		Returns an integer representing the generation where the final solution was determined.
		"""
        return(int(self.line.split("\t")[17]))
    @generation.setter
    def generation(self,newValue):
        self.line = updateLine(self.line,17,newValue)
        return(self.line)

    def __str__(self):
        """
        Returns a formated line for a skater database file.
		"""
        return(self.line)

    def geneCoverage(self):
        """
        Returns an integer of the length of a gene based on the TSS and TES of the gene. This is based on the first TSS and TES in the respective lists, which will need to be reexamined if genes with multiple TSS or TES are solved. 
        """
        if self.strand == "+":
            return(self.tesPos[0] - self.tssPos[0])
        else:
            return(self.tssPos[0] - self.tesPos[0])
        
    def speedZoneTuples(self):
        """
        Returns list of tuples of integers of speed zone positions [(0, 5),(5, 10)...]
        """
        speedZonePosTupleList = []
        for i in self.speedZonePos:
            speedZonePosTupleList.append((int(i.split("(")[1].split(" - ")[0]), int(i.split(")")[0].split(" - ")[1])))
        return(speedZonePosTupleList)
    
    def intronTuples(self):
        """
        Returns list of tuples of integers of intron genomic positions [(52582305, 52587890)]
        """
        intronTupleList = []
        for i in self.introns:
            try:
                intronTupleList.append((int(i.split("(")[1].split(" - ")[0]), int(i.split(")")[0].split(" - ")[1])))
            except IndexError:
                intronTupleList = [(0,0)] #gene with no introns
        return(intronTupleList)

    def timeTau(self, position): #position must be relative to tss
        """ 
        Function to calculate the time needed to transcribe to a specified position in a gene based on the optimize elongation rates.

        Parameters
        ----------
        position: position relative to TSS in base pairs.
        
        Returns
        ----------
        Time (in minutes) to transcribe from tss to position.
        """
        # tss = self.tssPos[0] #most upstream tss - does this need to change for multiple tss? tss are printed out in order of most upstream (so for negative strand genes the larger number is printed first)
        # tes = self.tesPos[-1]
        rates = self.speedZoneRates
        x = int(position)
        # if self.strand == "+":
        #     if x > tes:
        #         return(None)
        # else:
        #     if x < tes:
        #         return(None)
        if x > self.geneCoverage():
            # print("x error")
            return(None)
        posList = [0]
        tauList = [0]
        time = 0
        currentDist = 0
        for i in rates:
            dist = currentDist + (i*5) #this could be an argument for timepoint length
            time+=5
            posList.append(dist)
            tauList.append(time)
            currentDist = dist
        if x > posList[-1]:
            start = -2
            stop = -1
            m = (tauList[stop] - tauList[start]) / (posList[stop] - posList[start])
            b = tauList[start] - m * posList[start]
            return(m * x + b)
        for i in range(len(posList)):
            if posList[i] == x:
                return(tauList[i])
            elif posList[i] > x:
                start = i-1
                stop = i
                m = (tauList[stop] - tauList[start]) / (posList[stop] - posList[start])
                b = tauList[start] - m * posList[start]
                return(m * x + b)

def addNewRun(outfile,geneInfoDict,analysisID,errorFile):
    """ 
    Function to add new run to skater database.

    Parameters
    ----------
    outfile: path to skater database file
    geneInfoDict: key[geneName] = value(chrom,strand)
    analysisID: unique identifer of this particular analysis. Typically includes the name of the experiment, the condition, and the replicates
    errorFile: path to error file to collect information from runs that failed
    
    Returns
    ----------
    Line (str) of skater database file (can be used to write/append to skater database)
    """
    newRun = RunLevelData()
    newRun.geneID = outfile.split(".")[2].split("_")[1]
    try:
        newRun.chrom = geneInfoDict[outfile.split(".")[2].split("_")[1]][0]
        newRun.strand = geneInfoDict[outfile.split(".")[2].split("_")[1]][1]
    except KeyError:
        errordata = open(outfile).readlines()
        errorFileOutfile = open(errorFile,"a+")
        errorFileOutfile.write(outfile+"\n"+"\n".join(errordata)+"\n\n")
        errorFileOutfile.close()
        return("Error logged")
    newRun.runID = outfile.split(".")[-1].split("o")[1]
    newRun.analysisID = analysisID
    
    data = open(outfile)
    for line in data:
        if line.split(":")[0] == "generation\tstdev\tobjective\tG\tDRBrelease\tinitpos":
            try:
                names = line.split("\n")[0].split("\t")
                rates = data.readline().split("\n")[0].split("\t")
                break
            except IndexError:
                continue
    data.close()
    try:
        tssPosList = []
        spawnRatesList = []
        speedZonePosList = []
        speedZoneRatesList = []
        intronsList = []
        splicingRatesList = []
        tesPosList = []
        cleavageRatesList = []
        for i,j in zip(names,rates):
            try:
                parameter = i.split(": ")[0]
            except IndexError:
                parameter = i
                
            if parameter == "objective":
                newRun.RSS = j
            elif parameter == "generation":
                newRun.generation = j
            elif parameter == "DRBrelease":
                newRun.DRBreleaseRate = j
            elif parameter == "initpos":
                tssPosList.append(i.split(": ")[1])
                spawnRatesList.append(j)
            elif parameter == "speedZone":
                speedZonePosList.append(i.split(": ")[1])
                speedZoneRatesList.append(j)
            elif parameter == "intron":
                intronsList.append(i.split(": ")[1])
                splicingRatesList.append(j)
            elif parameter == "CleavagePos":
                tesPosList.append(i.split(": ")[1])
                cleavageRatesList.append(j)
            elif parameter == "contamination":
                newRun.contamination = j
            elif parameter == "EnOverIn":
                newRun.enOverIn = j
        newRun.tssPos = tssPosList
        newRun.spawnRates = spawnRatesList
        newRun.speedZonePos = speedZonePosList
        newRun.speedZoneRates = speedZoneRatesList 
        newRun.introns = intronsList 
        newRun.splicingRates = splicingRatesList 
        newRun.tesPos = tesPosList 
        newRun.cleavageRates = cleavageRatesList 
    except NameError:
        errordata = open(outfile).readlines()
        errorFileOutfile = open(errorFile,"a+")
        errorFileOutfile.write(outfile+"\n"+"\n".join(errordata)+"\n\n")
        errorFileOutfile.close()
        return("Error logged")
    return(str(newRun))

def varianceCalc(rateList):
    """ 
    Function to calculate variance of list of list of rates, convenience functions used in GeneLevelData

    Parameters
    ----------
    rateList: list of list of rates [[string1,string2,string3],[string1,string2,string3]]

    Returns
    ----------
    List of variances, one for each rate
    
    """
    import numpy as np
    length = len(rateList[0])
    if any(len(lst) != length for lst in [*rateList]):
        print(rateList)
        raise Exception("Rate lists are not the same length")
    var = []
    for i in zip(*rateList): #will zip for shortest list, all lists should be the same length
        data = [np.log2(float(x)) for x in i]
        var.append(np.var(data))
    return(var) #returns list of variances, one for each rate 

def bestDRBreleaseRateGenerator(gene,varLimit=math.inf,rateLB=0,rateUB=math.inf):
    """
    Generator to calculate the DRB release rate of the best optimizer run for a gene, one at a time
    To be used in a loop walking through GeneLevelData instances 

    Parameters
    ----------
    gene: instance of GeneLevelData class
    varLimit: maximum variance tolerated (variance between optimizer runs)
    rateLB: lower bound for rate
    rateUB: upper bound for rate
    Default bounds will include all rates, bounds should be set to the bounds used in the opimization argument file

    Returns
    ----------
    Tuple describing best rate: (gene, position, best rate)

    """
    rate = gene.best().DRBreleaseRate
    var = gene.drbVariance()[2]
    if var < float(varLimit) and float(rateUB)*0.95 >= rate >= float(rateLB)*1.05:
        yield((gene.geneName,int(1),float(rate))) #tuple (gene,position,best rate)

def bestSpawnRatesGenerator(gene,varLimit=math.inf,rateLB=0,rateUB=math.inf):
    """
    Generator to calculate the spawn rate of the best optimizer run for a gene, one at a time
    To be used in a loop walking through GeneLevelData instances

    Parameters
    ----------
    gene: instance of GeneLevelData class
    varLimit: maximum variance tolerated (variance between optimizer runs)
    rateLB: lower bound for rate
    rateUB: upper bound for rate
    Default bounds will include all rates, bounds should be set to the bounds used in the opimization argument file

    Returns
    ----------
    Tuple describing best rate: (gene, position, best rate)

    """ 
    pos = gene.best().tssPos
    rates = gene.best().spawnRates
    var = gene.spawnVariance()
    for i,j,k in zip(pos,rates,var):
        if k[2] < float(varLimit) and float(rateUB)*0.95 >= j >= float(rateLB)*1.05:
            yield((gene.geneName,int(i),float(j))) #tuple (gene,position,best rate)

def bestElongationRateGenerator(gene,varLimit=math.inf,rateLB=0,rateUB=math.inf): 
    """
    Generator to calculate the elongation rate of the best optimizer run for a gene, one at a time
    To be used in a loop walking through GeneLevelData instances

    Parameters
    ----------
    gene: instance of GeneLevelData class
    varLimit: maximum variance tolerated (variance between optimizer runs)
    rateLB: lower bound for rate
    rateUB: upper bound for rate
    Default bounds will include all rates, bounds should be set to the bounds used in the opimization argument file

    Returns
    ----------
    Tuple describing best rate: (gene, position, best rate)
    Will return this information for one speed zone at a time

    """
    pos = gene.best().speedZoneTuples()
    rates = gene.best().speedZoneRates
    var = gene.elongationVariance()
    for i,j,k in zip(pos,rates,var):
        if k[2] < float(varLimit) and float(rateUB)*0.95 >= j >= float(rateLB)*1.05:
            yield((gene.geneName,pos.index(i),float(j))) #tuple (gene,position,best rate)

def bestTauBasedElongationRateGenerator(gene,positionA,positionB,varLimit=math.inf,rateLB=0,rateUB=math.inf):
    """
    Generator to calculate the elongation rate within a given set of positions in a gene of the best optimizer run for a gene, one at a time
    To be used in a loop walking through GeneLevelData instances

    Parameters
    ----------
    gene: instance of GeneLevelData class
    positionA: start of zone to calculate rate within
    positionB: stop of zone to calculate rate within
    varLimit: maximum variance tolerated (variance between optimizer runs)
    rateLB: lower bound for rate
    rateUB: upper bound for rate
    Default bounds will include all rates, bounds should be set to the bounds used in the opimization argument file

    Returns
    ----------
    Tuple describing best rate: (gene, position, best rate)

    """
    var = gene.tauRateVariance(positionA,positionB)
    rate = (positionB - positionA) / (gene.best().timeTau(positionB)-gene.best().timeTau(positionA))
    if var[2] < float(varLimit) and float(rateUB)*0.95 >= rate >= float(rateLB)*1.05:
        yield((gene.geneName,(positionA,positionB),rate))

def bestSplicingRateGenerator(gene,varLimit=math.inf,rateLB=0,rateUB=math.inf):
    """
    Generator to calculate the splicing rate of the best optimizer run for a gene, one at a time
    To be used in a loop walking through GeneLevelData instances

    Parameters
    ----------
    gene: instance of GeneLevelData class
    varLimit: maximum variance tolerated (variance between optimizer runs)
    rateLB: lower bound for rate
    rateUB: upper bound for rate
    Default bounds will include all rates, bounds should be set to the bounds used in the opimization argument file

    Returns
    ----------
    Tuple describing best rate: (gene, position, best rate)
    Will return this information for one intron at a time

    """
    pos = gene.best().intronTuples()
    rates = gene.best().splicingRates
    var = gene.splicingVariance()
    if gene.best().strand == "-":
        pos.reverse()
        rates.reverse()
        var.reverse()
    for i,j,k in zip(pos,rates,var):
        if k[2] < float(varLimit) and float(rateUB)*0.95 >= j >= float(rateLB)*1.05:
            yield((gene.geneName,i,float(j))) #tuple (gene,position,best rate)
                
def bestCleavageGenerator(gene,varLimit=math.inf,rateLB=0,rateUB=math.inf):
    """
    Generator to calculate the cleavage rate of the best optimizer run for a gene, one at a time
    To be used in a loop walking through GeneLevelData instances

    Parameters
    ----------
    gene: instance of GeneLevelData class
    varLimit: maximum variance tolerated (variance between optimizer runs)
    rateLB: lower bound for rate
    rateUB: upper bound for rate
    Default bounds will include all rates, bounds should be set to the bounds used in the opimization argument file

    Returns
    ----------
    Tuple describing best rate: (gene, position, best rate)

    """
    pos = gene.best().tesPos
    rates = gene.best().cleavageRates
    var = gene.cleavageVariance()
    # elongationPos = gene.best().speedZoneTuples()
    # elongationRates = gene.best().speedZoneRates
    for i,j,k in zip(pos,rates,var):
        # overlappingSZ = overlappingSpeedZones(i,elongationPos)
        # for s in overlappingSZ:
        #     if elongationRates[s] >= float(elongationRateLimit) or s == windowIndexMax:
        #         break
        if k[2] < float(varLimit) and float(rateUB)*0.95 >= j >= float(rateLB)*1.05:
            yield((gene.geneName,int(i),float(j))) #tuple (gene,position,best rate)

def bestContaminationGenerator(gene,varLimit=math.inf):
    """
    Generator to calculate the contamination value of the best optimizer run for a gene, one at a time
    To be used in a loop walking through GeneLevelData instances

    Parameters
    ----------
    gene: instance of GeneLevelData class
    varLimit: maximum variance tolerated (variance between optimizer runs)
    rateLB: lower bound for rate
    rateUB: upper bound for rate
    Default bounds will include all rates, bounds should be set to the bounds used in the opimization argument file

    Returns
    ----------
    Tuple describing best rate: (gene, position, best rate)

    """
    rate = gene.best().contamination
    var = gene.contaminationVariance()[2]
    if var < float(varLimit):
        yield((gene.geneName,int(1),float(rate))) #tuple (gene,position,best rate)

class GeneLevelData:
    """
    A class with functions for a single gene in a SKaTERseq experiment, incorporating data from all runs (RunLevelData instances) of the optimizer for the gene. 
    """

    def __init__(self,runList): #run list is a list of RunLevelData instances that all match the same gene
        """
        A class used to construct a GeneLevelData instance, and to perform functions on an instance of GeneLevelData

        Parameters
		----------
		runList : a list of RunLevelData instances

		Returns
		----------
		GeneLevelData instance
        """
        self.runList = runList
#        self.best = bestRecord()
        self.geneName = runList[0].geneID
        
    def best(self):
        """
        Returns the RunLevelData with the lowest RSS (the "best") for this gene 
        """
        rss = []
        for record in self.runList:
            rss.append(float(record.RSS))
        bestRSS = min(rss)
        return(self.runList[rss.index(bestRSS)])
#        self.best = self.runList[rss.index(bestRSS)] #sets best record a RunLevelData instance
    
    def drbVariance(self):
        """
        Returns a tuple describing the variance of the DRB release rates for this gene (gene name, 1, float(variance)). 1 is a placeholder for the "position" since the DRB release rate does not have a position. 
        """
        rates = []
        for run in self.runList:
            rates.append([run.DRBreleaseRate])
        variance = varianceCalc(rates)
        return((self.geneName,int(1),float(variance[0])))

    def contaminationVariance(self):
        """
        Returns a tuple describing the variance of the contamination values for this gene (gene name, 1, float(variance)). 1 is a placeholder for the "position" since the contamination does not have a position. 
        """
        rates = []
        for run in self.runList:
            rates.append([run.contamination])
        variance = varianceCalc(rates)
        return((self.geneName,int(1),float(variance[0])))
        
    def spawnVariance(self):
        """
        Returns a list of tuples describing the variance of the spawn rates for this gene [(gene name, tss, float(variance)),...]
        """
        pos = self.runList[0].tssPos #pos = list positions (one position for each rate) [pos1,pos2]
        rates = []
        variance = []
        for run in self.runList:
            rates.append(run.spawnRates) #rates = list of list of rates (one list for each run)
        variance = varianceCalc(rates) #returns list of variances, one for each rate
        results = []
        for i,j in zip(pos,variance):
            results.append((self.geneName,int(i),float(j))) #results = list of tuples, (gene, position, and variance for the spawn rate at that position)
        return(results) 
    
    def elongationVariance(self):
        """
        Returns a list of tuples describing the variance of the elongation rates for this gene [(gene name, tss, float(variance)),...]
        """
        pos = self.runList[0].speedZoneTuples()
        rates = []
        variance = []
        for run in self.runList:
            rates.append(run.speedZoneRates)
        variance = varianceCalc(rates)
        results = []
        for i,j in zip(pos,variance):
            results.append((self.geneName,i,float(j)))
        return(results) 

    def tauVariance(self,position):
        """
        Returns a tuple describing the variance of the time to transcribe the gene [(gene name, position, float(variance)),...]
        """
        taus = []
        for run in self.runList:
            taus.append([run.timeTau(position)])
        variance = varianceCalc(taus)
        return((self.geneName,position,float(variance[0])))

    def tauRateVariance(self,positionA,positionB):
        """
        Returns a tuple describing the variance of the tau-based average elongation rate of the given position start stop on the gene [(gene name, position, float(variance)),...]
        """
        import numpy as np
        try:
            rates = []
            for run in self.runList:
                tau1 = run.timeTau(positionA)
                tau2 = run.timeTau(positionB)
                rate = (positionB-positionA)/(tau2-tau1)
                rates.append([np.log2(rate)])
            variance = varianceCalc(rates)
            return((self.geneName,(positionA,positionB),float(variance[0])))
        except TypeError:
            return(None)

    def bestTauBasedElongationRates(self,positionA,positionB,varLimit=math.inf,rateLB=0,rateUB=math.inf):
        """
        Function to calculate the average rate within a given region based on the best optimizer run for a gene. 

        Parameters
        ----------
        gene: instance of GeneLevelData class
        positionA: start of region (relative to TSS)
        positionB: stop of region (rel)
        varLimit: maximum variance tolerated (variance between optimizer runs)
        rateLB: lower bound for rate
        rateUB: upper bound for rate
        Default bounds will include all rates, bounds should be set to the bounds used in the opimization argument file

        Returns
        ----------
        Tuple describing tau: (gene, (positionA, positionB), rate within a region)

        """
        try:
            var = self.tauRateVariance(positionA,positionB)
            # print("var",var)
            rate = (positionB - positionA) / (self.best().timeTau(positionB)-self.best().timeTau(positionA))
            if var[2] < float(varLimit) and float(rateUB)*0.95 >= rate >= float(rateLB)*1.05:
                return((self.geneName,(positionA,positionB),rate))    
        except TypeError:
            return(None)

    def bestDeltaTimeTau(self,positionA,positionB,varLimit=math.inf):
        """
        Function to calculate the time to transcribe a given region based on the best optimizer run for a gene. 

        Parameters
        ----------
        gene: instance of GeneLevelData class
        positionA: start of region (relative to TSS)
        positionB: stop of region (rel)
        varLimit: maximum variance tolerated (variance between optimizer runs)
        rateLB: lower bound for rate
        rateUB: upper bound for rate
        Default bounds will include all rates, bounds should be set to the bounds used in the opimization argument file

        Returns
        ----------
        Tuple describing tau: (gene, (positionA, positionB), time to transcribe region)

        """
        try:
            varA = self.tauVariance(positionA)
            varB = self.tauVariance(positionB)
            # print("var",varA)
            # print("var",varB)
            deltaTau = self.best().timeTau(positionB) - self.best().timeTau(positionA)
            if varA[2] < float(varLimit) and varB[2] < float(varLimit):
                return((self.geneName,(positionA,positionB),deltaTau))
        except TypeError:
            return(None)

    def splicingVariance(self):
        """
        Returns a list of tuples describing the variance of the splicing rates for this gene [(gene name, tss, float(variance)),...]
        """
        pos = self.runList[0].intronTuples()
        rates = []
        variance = []
        for run in self.runList:
            rates.append(run.splicingRates)
        variance = varianceCalc(rates)
        results = []
        for i,j in zip(pos,variance):
            results.append((self.geneName,i,float(j)))
        return(results) 
    
    def cleavageVariance(self):
        """
        Returns a list of tuples describing the variance of the cleavage rates for this gene [(gene name, tss, float(variance)),...]
        """
        pos = self.runList[0].tesPos
        rates = []
        variance = []
        for run in self.runList:
            rates.append(run.cleavageRates)
        variance = varianceCalc(rates)
        results = []
        for i,j in zip(pos,variance):
           results.append((self.geneName,int(i),float(j)))
        return(results) 

    """
    geneLevelBestRate generators
    Generators to yield the best rate for each gene, looping through GeneLevelData instances for a given SKaTERseq experiment.
    """
    
    def geneLevelBestDRBreleaseRates(self,varLimit=math.inf,rateLB=0,rateUB=math.inf):
        for rate in bestDRBreleaseRateGenerator(self,varLimit,rateLB,rateUB):
            yield(rate)

    def geneLevelBestContamination(self,varLimit=math.inf):
        for rate in bestContaminationGenerator(self,varLimit):
            yield(rate)
    
    def geneLevelBestSpawnRates(self,varLimit=math.inf,rateLB=0,rateUB=math.inf):
        for rate in bestSpawnRatesGenerator(self,varLimit,rateLB,rateUB):
            yield(rate)
    
    def geneLevelBestElongationRates(self,varLimit=math.inf,rateLB=0,rateUB=math.inf):
        for rate in bestElongationRateGenerator(self,varLimit,rateLB,rateUB):
            yield(rate)

    def geneLevelBestSplicingRates(self,varLimit=math.inf,rateLB=0,rateUB=math.inf):
        for rate in bestSplicingRateGenerator(self,varLimit,rateLB,rateUB):
            yield(rate)
            
    def geneLevelBestCleavageRates(self,varLimit=math.inf,rateLB=0,rateUB=math.inf):
        for rate in bestCleavageGenerator(self,varLimit,rateLB,rateUB):
            yield(rate)

    def tauC(self,varLimit=math.inf,rateLB=0,rateUB=math.inf,secWithinTPfilter=0):
        rateInfo = self.bestTauBasedElongationRates(0,self.best().geneCoverage(),varLimit=varLimit,rateLB=rateLB,rateUB=rateUB)
        x = secWithinTPfilter/60
        timepoints = [5,10,15,20,25,30,35] #should be moved to an argument...
        try:
            rate = rateInfo[2]
            tau = self.best().geneCoverage()/rate
            if tau > 5 and tau < 35:
                for i in timepoints:
                    if i-x <= tau <= i or i <= tau <= i+x: 
                        # print(tau)
                        # return(None)
                        pass
                return(tau)
            else:
                return(None)
        except TypeError:
            return(None)

    def tauCrate(self,varLimit=math.inf,rateLB=0,rateUB=math.inf,secWithinTPfilter=0):
        x = secWithinTPfilter/60
        timepoints = [5,10,15,20,25,30,35] #should be moved to an argument...
        try:
            rate = self.bestTauBasedElongationRates(0,self.best().geneCoverage(),varLimit=varLimit,rateLB=rateLB,rateUB=rateUB)[2]
            tau = self.best().geneCoverage()/rate
            if tau > 5 and tau < 35:
                for i in timepoints:
                    if i-x <= tau <= i or i <= tau <= i+x: 
                        # print(tau)
                        return(None)
                return(rate)
        except TypeError:
            return(None)

    def splicingGeneAvg(self,varLimit=math.inf,rateLB=0,rateUB=math.inf,n=0):
        import numpy as np
        rates = []
        for i in self.geneLevelBestSplicingRates(varLimit=varLimit,rateLB=rateLB,rateUB=rateUB):
            rates.append(i[2])
        if len(rates) > n:
            return(np.mean(rates))
        else:
            return(None)

class dbReader:
    """
    A class with functions to open, read, and retrieve data from a SKaTERseq database file, gene by gene. SKaTERseq database file must be sorted by gene name. Records will be grouped by gene name and instances of the GeneLevelData class will be created in order to analyze the data in the database. 
    """
    def __init__(self,fileName):
        """
        A class used to construct a dbReader instance, and to perform functions on the SKaTERseq database.

        Parameters
		----------
		fileName : Path to a SKaTERseq database

		Returns
		----------
		dbReader instance
        """
        self.fileName = fileName
        
    def genes(self):
        """
        Function to open and walk through SKaTERseq database and yield one GeneLevelData instance at a time, containing all of the RunLevelData instances for that gene. SKaTERseq database must be sorted by gene name. 
        """
        file = open(self.fileName)
        firstRecord = RunLevelData(file.readline())
        recordList = [firstRecord]
        currentGene = firstRecord.geneID
        while True:
            line = file.readline()
            if line:
                try:
                    info = RunLevelData(line)
                except ValueError:
                    continue 
                gene = info.geneID
                if gene == currentGene:
                    recordList.append(info)
                    #if this is the last line of the file, it will "continue" to nothing and therefore never yield the last gene
                    continue
                else:
                    yield GeneLevelData(recordList)
                    recordList = [info]
                    currentGene = gene
            else:
                yield GeneLevelData(recordList)
                break

    """
    bestRate generators
    Generators to yield the best rate for each gene, looping through a SKaTERseq database.
    """

    def bestDRBreleaseRates(self,varLimit=math.inf,rateLB=0,rateUB=math.inf):
        for gene in self.genes():
            for rate in bestDRBreleaseRateGenerator(gene,varLimit,rateLB,rateUB):
                if rate is not None:
                    yield(rate)

    def bestContamination(self,varLimit=math.inf):
        for gene in self.genes():
            for rate in bestContaminationGenerator(gene,varLimit):
                if rate is not None:
                    yield(rate)

    def bestSpawnRates(self,varLimit=math.inf,rateLB=0,rateUB=math.inf): 
    #walks through database gene by gene, for each gene yields the (gene, pos, spawn rate) ensuring rate is not variable and is from the best run
        for gene in self.genes():
            for rate in bestSpawnRatesGenerator(gene,varLimit,rateLB,rateUB):
                if rate is not None:
                    yield(rate)

    def bestElongationRates(self,varLimit=math.inf,rateLB=0,rateUB=math.inf): #window index max is number of time points - 1
        for gene in self.genes():
            for rate in bestElongationRateGenerator(gene,varLimit,rateLB,rateUB):
                if rate is not None:
                    yield(rate)
    
    def bestTauBasedElongationRates(self,positionA,positionB,varLimit=math.inf,rateLB=0,rateUB=math.inf):
        for gene in self.genes():
            for rate in bestTauBasedElongationRateGenerator(gene,positionA,positionB,varLimit,rateLB,rateUB):
                if rate is not None:
                    yield(rate)
   
    def bestSplicingRates(self,varLimit=math.inf,rateLB=0,rateUB=math.inf):
        for gene in self.genes():
            for rate in bestSplicingRateGenerator(gene,varLimit,rateLB,rateUB):
                if rate is not None:
                    yield(rate)
                    
    def bestCleavageRates(self,varLimit=math.inf,rateLB=0,rateUB=math.inf):
        for gene in self.genes():
            for rate in bestCleavageGenerator(gene,varLimit,rateLB,rateUB):
                if rate is not None:
                    yield(rate)

def rateMatchGenerator(generator1,generator2): 
    """
    Function to match individual rates across SKaTERseq experiments. The data in each SKaTERseq database must be sorted by gene name, and the generators must be for the same rate type.  

    Parameters
    ----------
    generator1 and generator2: "best rate" generators from separate SKaTERseq experiments

    Returns
    ----------
    Tuple describing the matched rate: (gene,position,rate from gen1, rate from gen2)
    """
    try:
        x = next(generator1)
        y = next(generator2)
    except StopIteration:
        return(None)
    while True:
        try:
            if (x[0],x[1]) == (y[0],y[1]):
                result = (x[0],x[1],x[2],y[2]) #tuple of (gene,position,rate from gen1, rate from gen2)
                try:
                    x = next(generator1)
                    y = next(generator2)
                    yield(result)
                except StopIteration:
                    yield(result)
                    break
            # else:
            #     print(x,y)
            while (x[0],x[1]) < (y[0],y[1]):
                    x = next(generator1)
            while (x[0],x[1]) > (y[0],y[1]):
                    y = next(generator2)
        except StopIteration:
            break   

def rateMatchMatchGenerator(geneMatchGenerator1,geneMatchGenerator2):
    """
    Function to match MATCHED rates. 
    Typical use: match rates across replicates and across experiments so that average rates across replicates can be compared across experiments.
    This will walk through both generators and match GeneIDs only. 

    Parameters
    ----------
    generator1 and generator2: rate matched generators ((gene,position,rate from rep1,rate from rep2),(gene,position,rate from rep1,rate from rep2))

    Returns
    ----------
    Tuple describing the matched rates: ((gene,position,rate from gen1, rate from gen2),(gene,position,rate from gen1, rate from gen2))
    """
    x = next(geneMatchGenerator1)
    y = next(geneMatchGenerator2)
    while True:
        try:
            if (x[0],x[1]) == (y[0],y[1]):
                result = (x,y) #tuple of (gene,position,rate from gen1, rate from gen2)
                x = next(geneMatchGenerator1)
                y = next(geneMatchGenerator2)
                yield(result)
            while (x[0],x[1]) < (y[0],y[1]):
                    x = next(geneMatchGenerator1)
            while (x[0],x[1]) > (y[0],y[1]):
                    y = next(geneMatchGenerator2)
        except StopIteration:
            break   

def geneMatchGenerator(generator1,generator2): #update to accept a list of generators
    """
    Function to match GeneLevelData across SKaTERseq experiments. The data in each SKaTERseq database must be sorted by gene name

    Parameters
    ----------
    generator1 and generator2: "genes" generators from separate SKaTERseq experiments

    Returns
    ----------
    Tuple containing all of the information for the matched genes: (GeneLevelData from generator1, GeneLevelData from generator2)
    (rep1GeneLevelData(recordList), rep2GeneLevelData(recordList))
    """ 
    x = next(generator1)
    y = next(generator2) #use for loop to call next on each generator provided in list (for x in generators, next x)
    # items = [next x for x in generators]
    # if itertools.all 
    # find minimum, next on that one then recheck
    while True:
        try:
            if x.geneName == y.geneName:
                result = (x,y) # becomes tuple comprehension 
                x = next(generator1)
                y = next(generator2)
                yield(result)
            while x.geneName < y.geneName:
                x = next(generator1)
            while x.geneName > y.geneName:
                y = next(generator2)
        except StopIteration:
            break
        
def geneMatchMatchGenerator(geneMatchGenerator1,geneMatchGenerator2):
    """
    Function to match MATCHED genes. 
    Typical use: match genes across replicates and across experiments so that average rates across replicates can be compared across experiments.
    This will walk through both generators and match GeneIDs. 

    Parameters
    ----------
    generator1 and generator2: gene matched generators (rep1GeneLevelData(recordList), rep2GeneLevelData(recordList))

    Returns
    ----------
    Tuple describing the matched genes: ((rep1GeneLevelData(recordList), rep2GeneLevelData(recordList)),(rep1GeneLevelData(recordList), rep2GeneLevelData(recordList)))
    """
    x = next(geneMatchGenerator1)
    y = next(geneMatchGenerator2)
    while True:
        try:
            if x[0].geneName == y[0].geneName:
                result = (x,y)
                x = next(geneMatchGenerator1)
                y = next(geneMatchGenerator2)
                yield(result)
            while x[0].geneName < y[0].geneName:
                x = next(geneMatchGenerator1)
            while x[0].geneName > y[0].geneName:
                y = next(geneMatchGenerator2)
        except StopIteration:
            break

def intronAnnotationDict(intronAnnotationFile):
    #takes intron annotation file
    #returns dictionary: (gene, (intronStart, intronStop)) = (eventType, eventIndex)
    intronDict = dict()
    data = open(intronAnnotationFile)
    data.readline()
    for line in data:
        temp = line.strip().split("\t")
        intronDict[(temp[5],(int(temp[1]),int(temp[2])))] = (temp[3], temp[4])
    data.close()

def expressionDict(database,expressionData):
    geneSet = set()
    for i in dbReader(database).genes():
        if i.best().strand == "+":
            geneSet.add((i.best().geneID,i.best().tssPos[0],i.best().tesPos[0]))
        else:
            geneSet.add((i.best().geneID,i.best().tesPos[0],i.best().tssPos[0]))
    exp = dict()
    data = open(expressionData)
    for line in data:
        temp = line.strip().split("\t")
        gene = temp[3]  
        start = int(temp[1])
        stop = int(temp[2])          
        if (gene,start,stop) in geneSet:
            exp[gene] = float(temp[-1])
    data.close()
    return(exp)

def expressionCutoff(geneExpressionDict,fractionDenominator=1):
    geneExpression = []
    for k,v in geneExpressionDict.items():
        geneExpression.append((v,k))
    geneExpression.sort(reverse=True)
    #gene expression is sorted largest to smallest
    #length of list / denominator - take all ABOVE this number 
    return(geneExpression[int(len(geneExpression)/fractionDenominator)-1][0])
    #return gene expression 

def expressionCheck(gene,geneExpressionDict,cutoff,direction="GT"):
    try:
        if direction=="GT" and geneExpressionDict[gene] >= cutoff:
            return(True)
        elif direction=="LT" and geneExpressionDict[gene] < cutoff:
            return(True)
    except KeyError:
        return(False)

def relativePosition(tss,position,strand):
    if strand == "+":
        return(position - tss)
    else:
        return(tss - position)