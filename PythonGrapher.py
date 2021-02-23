import sys
import csv
import numpy as np
import PySide2
import PySide2.QtWidgets as Qt
import PySide2.QtGui as Gui
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.gridspec import GridSpec
from matplotlib.figure import Figure


#Same Curve should be set to true if the input data lies on the same curve as the
#extrapolation data. This is only needed in Nyquist mode, it is true always in
#Voigt mode regardless the value below.
SameCurve = True

class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent, mode, status, dpi=100):
        if status:
            self.initializeNormally(parent, mode, dpi=100)
        else:
            self.disableCaprinis(parent, mode, dpi=100)
        super(MplCanvas, self).__init__(self.fig)

    #initialize normally with caprinis.
    def initializeNormally(self, parent, mode, dpi=100):
        self.dataLabels = []
        self.lw = 1
        self.fig = Figure(dpi=dpi, constrained_layout=True)
        gs = GridSpec(2, 2, width_ratios=[3,2], figure=self.fig) #2x2 grid, more space alotted for left column.
        self.topCapriniAxes = self.fig.add_subplot(gs[0,1]) #top caprini is top left
        self.bottomCapriniAxes = self.fig.add_subplot(gs[1,1]) #bottom caprini is bottom left
        self.axisDict = {"TC":self.topCapriniAxes,
                        "BC":self.bottomCapriniAxes}
        if(mode == "Voigt"):
            self.topVoigtAxes = self.fig.add_subplot(gs[0,0])
            self.bottomVoigtAxes = self.fig.add_subplot(gs[1,0])
            self.axisDict["TV"] = self.topVoigtAxes
            self.axisDict["BV"] = self.bottomVoigtAxes
            self.setPlotInfo(self.topVoigtAxes, "Real Part", "Freq/Hz", "Z Real/Ohm")
            self.setPlotInfo(self.bottomVoigtAxes, "Imaginary Part", "Freq/Hz", "Z Imag/Ohm)")
            self.setPlotInfo(self.topCapriniAxes,"Caprini Function", "Freq/Hz", "")
            self.setPlotInfo(self.bottomCapriniAxes,"Zeros of the Caprini Function", "Freq/Hz", "")
        elif(mode == "Nyquist"):
            self.nyquistAxes = self.fig.add_subplot(gs[:,0]) #nyquist takes up whole left column
            self.axisDict["N"] = self.nyquistAxes
            self.setPlotInfo(self.nyquistAxes, "Nyquist Plot", "Re", "Im")
            self.setPlotInfo(self.topCapriniAxes,"Caprini Function", "t/S", "")
            self.setPlotInfo(self.bottomCapriniAxes,"Zeros of the Caprini Function", "t/S", "")

    #initialize without caprinis because they shouldn't be shown.
    def disableCaprinis(self, parent, mode, dpi=100):
        self.dataLabels = []
        self.lw = 1
        self.axisDict = {}
        self.fig = Figure(dpi=dpi, constrained_layout=True)
        if(mode == "Voigt"):
            gs = GridSpec(2, 1, figure=self.fig)
            self.topVoigtAxes = self.fig.add_subplot(gs[0,0])
            self.bottomVoigtAxes = self.fig.add_subplot(gs[1,0])
            self.axisDict["TV"] = self.topVoigtAxes
            self.axisDict["BV"] = self.bottomVoigtAxes
            self.setPlotInfo(self.topVoigtAxes, "Real Part", "Freq/Hz", "Z Real/Ohm")
            self.setPlotInfo(self.bottomVoigtAxes, "Imaginary Part", "Freq/Hz", "Z Imag/Ohm)")
        elif(mode == "Nyquist"):
            self.nyquistAxes = self.fig.add_subplot(111) #nyquist takes up whole left column
            self.axisDict["N"] = self.nyquistAxes
            self.setPlotInfo(self.nyquistAxes, "Nyquist Plot", "Re", "Im")


    def setGraphScale(self, axes, xMin, xMax, yMin, yMax):
        self.axisDict[axes].set_xlim(xMin, xMax)
        self.axisDict[axes].set_ylim(yMin, yMax)

    def drawSemiLogX(self, axes, x, y, dataLabel, curveColor):
        if not (dataLabel in self.dataLabels):
            self.dataLabels.append(dataLabel)
            self.axisDict[axes].semilogx(x,y, label=dataLabel, color=curveColor, linewidth=self.lw)
        else:
            self.axisDict[axes].semilogx(x,y,color=curveColor,linewidth=self.lw)
        if(axes == "TV" or axes == "N"):
            self.axisDict[axes].legend()

    def drawSemiLogXDots(self, axes, x, y, dataLabel, curveColor):
        if not (dataLabel in self.dataLabels):
            self.dataLabels.append(dataLabel)
            self.axisDict[axes].semilogx(x,y, '.',label=dataLabel, color=curveColor, linewidth=self.lw)
        else:
            self.axisDict[axes].semilogx(x,y,'.',color=curveColor,linewidth=self.lw)
        if(axes == "TV" or axes == "N"):
            self.axisDict[axes].legend()

    def drawLinLin(self, axes, x, y, dataLabel, curveColor):
        if not (dataLabel in self.dataLabels):
            self.dataLabels.append(dataLabel)
            self.axisDict[axes].plot(x,y, label=dataLabel, color=curveColor,linewidth=self.lw)
        else:
            self.axisDict[axes].plot(x,y, color = curveColor,linewidth=self.lw)
        if(axes == "TV" or axes == "N"):
            self.axisDict[axes].legend()

    def drawLinLinDots(self, axes, x, y, dataLabel, curveColor):
        if not (dataLabel in self.dataLabels):
            self.dataLabels.append(dataLabel)
            self.axisDict[axes].plot(x,y, '.',label=dataLabel, color=curveColor,linewidth=self.lw)
        else:
            self.axisDict[axes].plot(x,y, '.',color = curveColor,linewidth=self.lw)
        if(axes == "TV" or axes == "N"):
            self.axisDict[axes].legend()

    def setPlotInfo(self, axes, title, xTitle, yTitle):
        axes.set_title(title)
        axes.set_xlabel(xTitle)
        axes.set_ylabel(yTitle)


class StieltjesGraphWindow(Qt.QDialog):
    def __init__(self, sameCurveArg):
        super().__init__()#set up Qt
        self.topLevelLayout = Qt.QHBoxLayout()
        self.graphData = [[],[]] #set by populateData()
        self.mode, zetaFZeta = self.detectMode() #is set by populateData()
        self.sameCurve = sameCurveArg #set by user. whether the data lies on the extrapolation curve.
        self.status, fZeta = self.detectStatus() #whether to show caprini functions

        #setup graphs!
        self.graph = Qt.QFrame()
        graphLayout = Qt.QVBoxLayout()
        self.graphCanvas = MplCanvas(self,self.mode, self.status, dpi=100)
        toolbar = NavigationToolbar(self.graphCanvas, self)
        graphLayout.addWidget(toolbar)
        graphLayout.addWidget(self.graphCanvas)
        self.graph.setLayout(graphLayout)
        self.graph.setFrameShape(Qt.QFrame.Panel)
        self.graph.setFrameShadow(Qt.QFrame.Sunken)
        self.graph.setLineWidth(2)

        graphLayout = Qt.QVBoxLayout()  #update topLevelLayout to include graphs!
        graphLayout.addWidget(self.graph)
        self.topLevelLayout.addItem(graphLayout)
        self.populateData(zetaFZeta, fZeta)      #take in txt file data and add it to graphs.
        self.setLayout(self.topLevelLayout)
        self.setWindowState(PySide2.QtCore.Qt.WindowState.WindowMaximized)


    def populateData(self, zetaFZeta, fZeta): #figure out graph contents
        #populate data for correct graph
        if self.mode == "Nyquist":
            self.populateNyquist(zetaFZeta, fZeta)
        elif self.mode == "Voigt":
            self.populateVoigt(zetaFZeta, fZeta)
        if(self.status):
            self.populateCaprinis()


    def populateNyquist(self, zetaFZeta, fZeta): #update nyquist plot
        #Spectral rep and Extrapolation curves
        Z, W = [], []
        for row in zetaFZeta:
            Z.append(complex(row[0], row[1]))   #needed for spectral rep
            W.append(complex(row[2], row[3]))   #extrapolation
        fZ = self.computeSpectralRep(Z, fZeta)  #finish computing spectral rep
        #Uncertainty curves
        dataSizes = self.loadFileIntoNestedArrays("data_sizes.txt")
        NMC = int(dataSizes[2][0])
        nZeta = int(dataSizes[1][0])
        count = 0
        if NMC > 0:
            MonteCarlo = self.loadFileIntoNestedArrays("WMC.txt")
            WMC = []
            for k in range(nZeta): #iteration is not row->col to suit WMC.txt format.
                for j in range(NMC):
                    if(k==0):
                        WMC.append([])
                    WMC[j].append(complex(MonteCarlo[count][0], MonteCarlo[count][1]))
                    count=count+1
        #actually plot the results
        for j in range(NMC):
            self.graphCanvas.drawLinLin("N",np.real(WMC[j]),np.imag(WMC[j]),"Uncertainty", "silver")
        self.graphCanvas.drawLinLin("N",np.real(W), np.imag(W), "Extrapolation", "orange") #graph extrapolation
        self.graphCanvas.drawLinLin("N",np.real(fZ), np.imag(fZ), "Spectral Rep", "cyan") #graph spectral rep
        #plot original points if extrap is on same curve is on original data.
        if self.sameCurve:
            expData = self.loadFileIntoNestedArrays("exp_data.txt")
            if self.status:
                wfix = np.transpose(self.loadFileIntoNestedArrays("wfix.txt"))
                self.graphCanvas.drawLinLinDots("N", wfix[0], wfix[1], "Alternative Data", "red")
            wOrig = []
            for row in expData:
                wOrig.append(complex(row[2], row[3]))
            self.graphCanvas.drawLinLinDots("N", np.real(wOrig), np.imag(wOrig), "W Original", "black")


    def computeSpectralRep(self, Z, fZeta):
        sk = [] #sk is the 2nd column of fzeta.
        tk = [] #tk is the 1st column of fzeta, w/out the 1st element
        for i in range(len(fZeta)):
            sk.append(fZeta[i][1])
            if i != 0:
                tk.append(fZeta[i][0])
        M=[]  #M = 1/(tk - Z), tk is a column, Z is a row.
        for i in range(len(tk)):
            M.append([])
            for j in range(len(Z)):
                M[i].append( 1/(tk[i] - Z[j]) )
        Mt = np.transpose(M)
        fZ = np.matmul(Mt,sk[1:]) #fZ=[] #fZ=sk[0] + M.'* sk[1:]
        fZ = [sk[0] + i for i in fZ]
        return fZ


    def populateVoigt(self, zetaFZeta, fZeta):
        # fZ exp data [freq,re(Z),im(Z)]; wfix=alternative data
        # f=all frequencies, F=Z(f) is the model
        # Z(n_zeta,2) is the extrapolation
        # EIS(f)=F_spectral(f), WZ(n_zeta,Nr)=monte-carlo
        dataSizes=self.loadFileIntoNestedArrays('data_sizes.txt');
        fZ1Split = self.loadFileIntoNestedArrays("Voigt_data.txt")
        Nr = int(dataSizes[2][0])
        #unpack data
        freq, fZ1, f1, sk, tk, zetaFZetaCombined  = [],[],[],[],[],[]
        for row in fZ1Split:
            freq.append(row[0])
            fZ1.append(complex(row[1], row[2]))
        for row in zetaFZeta:
            f1.append(row[0])
        nZeta = len(f1)
        for i in range(len(fZeta)):
            if i > 0:
                tk.append(fZeta[i][0])
            sk.append(fZeta[i][1])
        for row in zetaFZeta:
            zetaFZetaCombined.append(complex(row[1], row[2]))
        #compute M and EIS - to show voigt circuit on graphs.
        M=[]  #M=1./(tk+2i*pi*f1');
        for i in range(len(tk)):
            M.append([])
            for j in range(len(f1)):
                M[i].append( 1/(complex(tk[i], f1[j]*2*np.pi)) )
        Mt = np.transpose(M)
        EIS = sk[0] + np.matmul(Mt, sk[1:])
        #uncertainty lines retrieved and computed
        if Nr > 0:
            WZ = []
            WMC = self.loadFileIntoNestedArrays("WMC.txt");
            count=0;
            for k in range(nZeta):
                for j in range(Nr):
                    if(k==0):
                        WZ.append([])
                    WZ[j].append( complex(WMC[count][0], WMC[count][1]) )
                    count+=1
        #uncertainty
        for i in range(Nr):
            self.graphCanvas.drawSemiLogX("TV", f1, np.real(WZ[i]), "Uncertainty", "silver")
            self.graphCanvas.drawSemiLogX("BV", f1, np.imag(WZ[i]), "Uncertainty", "silver")
        if self.status: #alternative data is only printed if status is true.
            wfix = self.loadFileIntoNestedArrays("wfix.txt")
            wfix = np.transpose(wfix)
            self.graphCanvas.drawSemiLogXDots("TV", freq, wfix[0], "Alternative Data","red")
            self.graphCanvas.drawSemiLogXDots("BV", freq, wfix[1], "Alternative Data","red")
        #upper graph - real part
        self.graphCanvas.drawSemiLogXDots("TV", freq, np.real(fZ1), "Data","black")
        self.graphCanvas.drawSemiLogX("TV", f1, np.real(zetaFZetaCombined), "Extrapolation", "black")
        self.graphCanvas.drawSemiLogX("TV", f1, np.real(EIS), "Voigt Circuit", "cyan")
        #lower graph - imaginary part
        self.graphCanvas.drawSemiLogXDots("BV", freq, np.imag(fZ1), "Data","black")
        self.graphCanvas.drawSemiLogX("BV", f1, np.imag(zetaFZetaCombined), "Extrapolation", "black")
        self.graphCanvas.drawSemiLogX("BV", f1, np.imag(EIS), "Voigt Circuit", "cyan")

    def populateCaprinis(self):
        capriniData = self.loadFileIntoNestedArrays("tC.txt") #The Caprini function [t(k),C(t(k))]
        tC = [] #data to be graphed - x
        C = []     #y axis data
        for row in capriniData:
            tC.append(row[0])
            C.append(row[1])
        self.graphCanvas.drawSemiLogX("TC",[tC[0],tC[len(tC)-1]],[0,0],"","red") #draw baseline
        self.graphCanvas.drawSemiLogX("TC", tC, C, "", "blue" ) #draw caprini
        self.graphCanvas.drawSemiLogX("BC",[tC[0],tC[len(tC)-1]],[0,0],"","red") #draw baseline
        self.graphCanvas.drawSemiLogX("BC",tC, C, "", "blue")   #draw caprini mins
        #scale the bottom graph to show zeros
        CMin, CMax = self.getgraphScale(C)
        self.graphCanvas.setGraphScale("BC",min(tC), max(tC), CMin, CMax)

    def getgraphScale(self, C): #find local minima, and figure out suitable graph scale.
        diffs = np.diff(C)
        CLocalMins = []
        lastNegative = False
        for i in range(len(diffs)):
            thisPositive = diffs[i] > 0
            if(lastNegative and thisPositive):
                CLocalMins.append(C[i])
            lastNegative = diffs[i] < 0
        biggestMin = max(CLocalMins)
        smallestMin = min(CLocalMins)
        CMax=max([-10*smallestMin,1.1*biggestMin,1.e-9*max(C)])
        CMin = -1.e-10*max(C)
        return CMin, CMax


    def detectMode(self):
        zetaFZeta=self.loadFileIntoNestedArrays('W_extr.txt')
        if len(zetaFZeta[0]) == 4:
            return "Nyquist", zetaFZeta
        elif len(zetaFZeta[0]) == 3:
            return "Voigt", zetaFZeta
        else:
            return "UnknownData", zetaFZeta

    def detectStatus(self):
        fZeta = self.loadFileIntoNestedArrays("spectral_measure.txt")
        if(fZeta[0][0] == 0): #determine if caprinis should be shown.
            return False, fZeta
        else:
            return True, fZeta

    def loadFileIntoNestedArrays(self, filename): #load a data file into a 2d array
        dataList = []
        with open(filename, 'r') as csvfile:
            csvreader = csv.reader(csvfile, delimiter=" ")
            #fields = next(csvreader)
            i = 0 #counter variable needed iteration over len(csvreader) not possible
            for row in csvreader: #fill in data
                dataList.append([]) #add next column of dataList
                for j in range(len(row)):
                    if(row[j] != ""):
                        dataList[i].append(float(row[j]))
                i+=1
        return dataList


if __name__ == '__main__':
    # Create the Qt Application
    app = Qt.QApplication(sys.argv)
    # Create and show the form
    form = StieltjesGraphWindow(SameCurve)
    form.show()
    # Run the main Qt loop
    sys.exit(app.exec_())
