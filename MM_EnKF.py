import matplotlib.pyplot as plt
from numpy import *
from random import choice
from copy import deepcopy

class Input(object):
    def __init__(self):
        self.EnNumber=100
        self.length=4.0
        self.laneNumber=3

        self.CTMNoiseMean=0.0
        self.CTMNoiseStd=1.0
        self.MeaNoiseDensity=3.0
        self.MeaNoiseVelocity=3.0

    def setLeftGhost(self,counter):
        self.leftGhost=5000+random.normal(0,150)
        return self.leftGhost
    
    def setLeftGhost_true(self,counter):
        self.leftGhost=5000
        return self.leftGhost
        
class KalmanEstimationModel(Input):
    def __init__(self):
        Input.__init__(self)
# define three parameters to determine the fundamental diagram
        self.Qmax=2000.0
        self.Vmax=65.0
        self.wf=28.0
# discretization 
        self.timeStep=20.0/3600
        self.spaceStep=self.Vmax*self.timeStep
        self.cellNumber=floor(self.length/self.spaceStep)
        self.spaceStep=self.length/self.cellNumber
        
        self.initializeDensity()
        self.initializeSpeed()

    def initializeDensity(self):
# initialize density for all cells. 
        self.cellDensity=self.laneNumber*self.Qmax/self.Vmax*ones((self.EnNumber,self.cellNumber))+random.normal(0,0.05*self.Qmax/self.Vmax,(self.EnNumber,self.cellNumber))
        self.updateCellDensity=zeros((self.EnNumber,self.cellNumber))

    def initializeSpeed(self):
# initialize speed for all cells
        self.cellSpeed=self.Vmax*ones((self.EnNumber,self.cellNumber))+random.normal(0,0.05*self.Vmax,(self.EnNumber,self.cellNumber))

    def setRightGhost(self):
# set the right boundary conditions
        self.rightGhost=inf
        return self.rightGhost

    def sending(self,density,mode):
# the sending function in CTM
        CriticalDensity=mode*self.Qmax/self.Vmax
        JamDensity=CriticalDensity+mode*self.Qmax/self.wf
        if density<0:
            density=.001
        if density<CriticalDensity:
            qSend=density*self.Vmax
        else:
            qSend=self.Qmax*mode
        if qSend<0:
            print 'qSend is smaller than 0, please check'
        return qSend

    def receiving(self,density,mode):
# the receiving function in CTM
        CriticalDensity=mode*self.Qmax/self.Vmax
        JamDensity=CriticalDensity+mode*self.Qmax/self.wf
        if density>JamDensity:
            density=JamDensity-0.001
        if density<CriticalDensity:
            qReceive=self.Qmax*mode
        else:
            qReceive=(JamDensity-density)*self.wf
        if qReceive<0:
            print 'qReceive is smaller than 0, please check'
        return qReceive

    def updateDensity(self, counter,mode):
# CTM forward prediction 
        for j in range(self.EnNumber):
            for i in range(0,int(self.cellNumber)):
                    if i==0:
                            self.updateCellDensity[j,i]=self.cellDensity[j,i]+(self.timeStep/self.spaceStep)*\
                                                   (min(self.setLeftGhost(counter),self.receiving(self.cellDensity[j,i],mode[i]))-\
                                                    min(self.sending(self.cellDensity[j,i],mode[i]),self.receiving(self.cellDensity[j,i+1],mode[i+1])))
                    elif i>0 and i<int(self.cellNumber-1):
                            self.updateCellDensity[j,i]=self.cellDensity[j,i]+(self.timeStep/self.spaceStep)*\
                                                   (min(self.sending(self.cellDensity[j,i-1],mode[i-1]),self.receiving(self.cellDensity[j,i],mode[i]))-\
                                                    min(self.sending(self.cellDensity[j,i],mode[i]),self.receiving(self.cellDensity[j,i+1],mode[i+1])))
                    else:
                            self.updateCellDensity[j,i]=self.cellDensity[j,i]+(self.timeStep/self.spaceStep)*\
                                                   (min(self.sending(self.cellDensity[j,i-1],mode[i-1]),self.receiving(self.cellDensity[j,i],mode[i]))-\
                                                    min(self.sending(self.cellDensity[j,i],mode[i]),self.setRightGhost()))
        self.cellDensity=self.updateCellDensity.copy()+random.normal(self.CTMNoiseMean, self.CTMNoiseStd,(self.EnNumber,self.cellNumber))

    def updateSpeed(self,mode):
# Velocity update
        for j in range(self.EnNumber):
            for i in range(0,int(self.cellNumber)):
                CriticalDensity=mode[i]*self.Qmax/self.Vmax
                JamDensity=CriticalDensity+mode[i]*self.Qmax/self.wf
                if self.cellDensity[j,i]<CriticalDensity:
                    self.cellSpeed[j,i]=self.Vmax
                else:
                    self.cellSpeed[j,i]=-self.wf*(1-JamDensity/self.cellDensity[j,i])


class TrueModel(Input):
    def __init__(self):
        Input.__init__(self)
# define three parameters to determine the fundamental diagram
        self.Qmax=1900.0
        self.Vmax=60.0
        self.wf=30.0
# discretization 
        self.timeStep=20.0/3600
        self.spaceStep=self.Vmax*self.timeStep
        self.cellNumber=floor(self.length/self.spaceStep)
        self.spaceStep=self.length/self.cellNumber
        self.mode=self.laneNumber*ones(self.cellNumber)
        
        self.initializeDensity()
        self.initializeSpeed()

    def initializeDensity(self):
# initialize density for all cells. 
        self.cellDensity=self.laneNumber*self.Qmax/self.Vmax*ones(self.cellNumber)
        self.updateCellDensity=zeros(self.cellNumber)

    def initializeSpeed(self):
# initialize speed for all cells
        self.cellSpeed=self.Vmax*ones(self.cellNumber)

    def setRightGhost(self):
# set the right boundary conditions
        self.rightGhost=inf
        return self.rightGhost

    def sending(self,density,mode):
# the sending function in CTM
        CriticalDensity=mode*self.Qmax/self.Vmax
        JamDensity=CriticalDensity+mode*self.Qmax/self.wf
        if density<0:
            density=.001
        if density<CriticalDensity:
            qSend=density*self.Vmax
        else:
            qSend=self.Qmax*mode
        if qSend<0:
            print 'qSend is smaller than 0, please check'
        return qSend

    def receiving(self,density,mode):
# the receiving function in CTM
        CriticalDensity=mode*self.Qmax/self.Vmax
        JamDensity=CriticalDensity+mode*self.Qmax/self.wf
        if density>JamDensity:
            density=JamDensity-0.001
        if density<CriticalDensity:
            qReceive=self.Qmax*mode
        else:
            qReceive=(JamDensity-density)*self.wf
        if qReceive<0:
            print 'qReceive is smaller than 0, please check'
        return qReceive

    def updateDensity(self, counter):
# CTM forward prediction
        if counter>60:
            self.mode[4]=1
        if counter>120:
            self.mode[4]=3
        for i in range(0,int(self.cellNumber)):
                if i==0:
                        self.updateCellDensity[i]=self.cellDensity[i]+(self.timeStep/self.spaceStep)*\
                                               (min(self.setLeftGhost_true(counter),self.receiving(self.cellDensity[i],self.mode[i]))-\
                                                min(self.sending(self.cellDensity[i],self.mode[i]),self.receiving(self.cellDensity[i+1],self.mode[i+1])))
                elif i>0 and i<int(self.cellNumber-1):
                        self.updateCellDensity[i]=self.cellDensity[i]+(self.timeStep/self.spaceStep)*\
                                               (min(self.sending(self.cellDensity[i-1],self.mode[i-1]),self.receiving(self.cellDensity[i],self.mode[i]))-\
                                                min(self.sending(self.cellDensity[i],self.mode[i]),self.receiving(self.cellDensity[i+1],self.mode[i+1])))
                else:
                        self.updateCellDensity[i]=self.cellDensity[i]+(self.timeStep/self.spaceStep)*\
                                               (min(self.sending(self.cellDensity[i-1],self.mode[i-1]),self.receiving(self.cellDensity[i],self.mode[i]))-\
                                                min(self.sending(self.cellDensity[i],self.mode[i]),self.setRightGhost()))
        self.cellDensity=self.updateCellDensity.copy()

    def updateSpeed(self):
# Velocity update
        for j in range(self.EnNumber):
            for i in range(0,int(self.cellNumber)):
                CriticalDensity=self.mode[i]*self.Qmax/self.Vmax
                JamDensity=CriticalDensity+self.mode[i]*self.Qmax/self.wf
                if self.cellDensity[i]<CriticalDensity:
                    self.cellSpeed[i]=self.Vmax
                else:
                    self.cellSpeed[i]=-self.wf*(1-JamDensity/self.cellDensity[i])



class GPSvehicle(TrueModel):
## This class simulates GPS vehicle trajectory
    def __init__(self):
        TrueModel.__init__(self)
        self.PR=40
        self.initializeGPSposition()
        self.initializeGPScellPosition()
        

## This function initialize the position of GPS vehicles 
    def initializeGPSposition(self):
        if self.PR==20:
                self.GPSposition=range(int(self.cellNumber))
        if self.PR==40:
                self.GPSposition=range(0,int(self.cellNumber),2)
        if self.PR==60:
                self.GPSposition=range(0,int(self.cellNumber),3)
        if self.PR==80:
                self.GPSposition=range(0,int(self.cellNumber),4)
        if self.PR==100:
                self.GPSposition=range(0,int(self.cellNumber),5)
        if self.PR==120:
                self.GPSposition=range(0,int(self.cellNumber),6)
        if self.PR==140:
                self.GPSposition=range(0,int(self.cellNumber),7)
        if self.PR==160:
                self.GPSposition=range(0,int(self.cellNumber),8)
        if self.PR==180:
                self.GPSposition=range(0,int(self.cellNumber),9)


## This function initilize the position of GPS vehicle in their cells
    def initializeGPScellPosition(self):
        self.GPScellPosition=zeros((len(self.GPSposition)))

## This function updates the position of GPS vehicle 
    def GPSpositionUpdate(self,counter):
            if counter != 0:
                    if counter%(self.PR/20)==0:
                            self.GPSposition.insert(0,0)
                            self.GPScellPosition=insert(self.GPScellPosition,0,0)
            self.GPSpositionCopy = deepcopy(self.GPSposition)
            self.GPScellPositionCopy = deepcopy(self.GPScellPosition)
            for i in range(len(self.GPSposition)-1,-1,-1):
                    if (self.cellSpeed[int(self.GPSposition[i])]*self.timeStep)>(self.spaceStep-self.GPScellPosition[i]):
                            if int(self.GPSposition[i]) != self.cellNumber-1:
                                    self.GPSpositionCopy[i]=self.GPSpositionCopy[i]+1
                                    self.GPScellPosition[i]=(self.timeStep-(self.spaceStep-self.GPScellPosition[i])/self.cellSpeed[int(self.GPSposition[i])])*self.cellSpeed[int(self.GPSposition[i]+1)]
                            else:
                                    self.GPSpositionCopy.remove(self.GPSposition[i])
                                    self.Xcell=delete(self.GPScellPosition,i)
                    else:
                            self.GPScellPosition[i]=self.GPScellPosition[i]+self.cellSpeed[int(self.GPSposition[i])]*self.timeStep
            self.GPSposition=deepcopy(self.GPSpositionCopy)

        

## Here starts the main code:

## Define the set of models
def SystemModel(cellNumber, laneNumber):
    ModelNumber=(cellNumber-4)*(laneNumber-1)+1
    Model=3*ones((ModelNumber,cellNumber))
    i=1
    for j in range(2,int(cellNumber-2)):
        for k in range(1,int(laneNumber)):
            Model[i,j]=k
            i=i+1
    return Model

def CalculateCov(EnNumber,cellNumber,cellDensity):
# compute the covariance matrix
    CovMatrix=zeros((cellNumber, cellNumber))
    for i in range(0,EnNumber):
        CovMatrix=CovMatrix+matrix(cellDensity[i,:]-average(cellDensity,axis=0)).T*matrix(cellDensity[i,:]-average(cellDensity,axis=0))
    CovMatrix=CovMatrix/(EnNumber-1)
    return CovMatrix


def OperatorC(cellNumber,GPSposition):
# linear operator to match the prediction to the measurements
    H=zeros(2*cellNumber)
    H[1]=1
    H[cellNumber-2]=1
    for i in GPS.GPSposition:
        H[cellNumber+i]=1
    C=zeros((int(sum(H)),int(2*cellNumber)))
    k=0
    for j in range(len(H)):
        if H[j] != 0:
            C[k,j]=1
            k=k+1
    return C

def computeKalmanGain(cellNumber,EnNumber,GPSposition,cellDensity,cellSpeed,MeaNoiseDensity,MeaNoiseVelocity):
# compute the kalman gain
        C=OperatorC(cellNumber,GPSposition)
        StatePredict=hstack((cellDensity, cellSpeed))
        StatePredictMean=average(StatePredict,axis=0)
        CovMatrixHPH=zeros((sum(C), sum(C)))
        CovMatrixPH=zeros((cellNumber,sum(C)))
        CovMatrixR=zeros((sum(C), sum(C)))
        for i in range(0,EnNumber):
            CovMatrixHPH=CovMatrixHPH+matrix(matrix(C)*matrix(StatePredict[i,:]).T-matrix(C)*matrix(StatePredictMean).T).T*matrix(matrix(C)*matrix(StatePredict[i,:]).T-matrix(C)*matrix(StatePredictMean).T)
            CovMatrixPH=CovMatrixPH+matrix(cellDensity[i,:]-average(cellDensity,axis=0)).T*matrix(matrix(C)*matrix(StatePredict[i,:]).T-matrix(C)*matrix(StatePredictMean).T).T
        CovMatrixHPH=CovMatrixHPH/(EnNumber-1)
        CovMatrixPH=CovMatrixPH/(EnNumber-1)
        for j in range(0,int(sum(C))):
            if j<2:
                CovMatrixR[j,j]=MeaNoiseDensity*MeaNoiseDensity
            else:
                CovMatrixR[j,j]=MeaNoiseVelocity*MeaNoiseVelocity
        kalmanGain=matrix(CovMatrixPH)*matrix(CovMatrixHPH+CovMatrixR).I
        return kalmanGain


def StateEstimation(cellDensity,cellSpeed,cellNumber,EnNumber,GPSposition,TMcellDensity,TMcellSpeed,kalmanGain,MeaNoiseDensity,MeaNoiseVelocity):
        StatePredict=hstack((cellDensity, cellSpeed))
        C=OperatorC(cellNumber,GPSposition)
        TrueMea=hstack((TMcellDensity, TMcellSpeed))
        meaNoise=zeros(sum(C))
        for i in range(0,int(sum(C))):
            if i<2:
                meaNoise[i]=MeaNoiseDensity*MeaNoiseDensity
            else:
                meaNoise[i]=MeaNoiseVelocity*MeaNoiseVelocity

        for i in range(0,EnNumber):
            EnKF.cellDensity[i,:]=EnKF.cellDensity[i,:]+array((kalmanGain*(matrix(C)*matrix(TrueMea).T-matrix(C)*matrix(StatePredict[i,:]).T+matrix(meaNoise).T)).T)[0]

def StateEstimationS(cellDensityCopy,cellDensity,cellSpeed,cellNumber,EnNumber,GPSposition,TMcellDensity,TMcellSpeed,kalmanGain,MeaNoiseDensity,MeaNoiseVelocity):
        StatePredict=hstack((cellDensity, cellSpeed))
        C=OperatorC(cellNumber,GPSposition)
        TrueMea=hstack((TMcellDensity, TMcellSpeed))
        meaNoise=zeros(sum(C))
        for i in range(0,int(sum(C))):
            if i<2:
                meaNoise[i]=MeaNoiseDensity*MeaNoiseDensity
            else:
                meaNoise[i]=MeaNoiseVelocity*MeaNoiseVelocity

        for i in range(0,EnNumber):
            cellDensityCopy[i,:]=cellDensityCopy[i,:]+array((kalmanGain*(matrix(C)*matrix(TrueMea).T-matrix(C)*matrix(StatePredict[i,:]).T+matrix(meaNoise).T)).T)[0]
        return cellDensityCopy

def ModelLikelihood(modeIndex):
    PosteriorModelDensity[modeIndex,:,:]=EnKF.cellDensity
    PosteriorModelSpeed[modeIndex,:,:]=EnKF.cellSpeed
    
    C=OperatorC(EnKF.cellNumber,GPS.GPSposition)
    TrueMea=hstack((TM.cellDensity, TM.cellSpeed))
    EstimationStateMean=hstack((average(EnKF.cellDensity,axis=0),average(EnKF.cellSpeed,axis=0)))
    
    Diff=(matrix(C)*matrix(TrueMea).T-matrix(C)*matrix(EstimationStateMean).T)

    for i in range(0,int(sum(C))):
        if i<2:
            ModelProb[modeIndex]=ModelProb[modeIndex]+log(1/(EnKF.MeaNoiseDensity*sqrt(2*pi))*exp(-(Diff[i])*(Diff[i])/(2*EnKF.MeaNoiseDensity*EnKF.MeaNoiseDensity)))
        else:
            ModelProb[modeIndex]=ModelProb[modeIndex]+log(1/(EnKF.MeaNoiseVelocity*sqrt(2*pi))*exp(-(Diff[i])*(Diff[i])/(2*EnKF.MeaNoiseVelocity*EnKF.MeaNoiseVelocity)))
    ModelProb[modeIndex]=exp(ModelProb[modeIndex])
    return (PosteriorModelDensity,PosteriorModelSpeed,ModelProb)


def ModelLikelihood(modeIndex):
    PosteriorModelDensity[modeIndex,:,:]=EnKF.cellDensity
    PosteriorModelSpeed[modeIndex,:,:]=EnKF.cellSpeed
    
    C=OperatorC(EnKF.cellNumber,GPS.GPSposition)
    TrueMea=hstack((TM.cellDensity, TM.cellSpeed))
    EstimationStateMean=hstack((average(EnKF.cellDensity,axis=0),average(EnKF.cellSpeed,axis=0)))
    
    Diff=(matrix(C)*matrix(TrueMea).T-matrix(C)*matrix(EstimationStateMean).T)

    for i in range(0,int(sum(C))):
        if i<2:
            ModelProb[modeIndex]=ModelProb[modeIndex]+log(1/(EnKF.MeaNoiseDensity*sqrt(2*pi))*exp(-(Diff[i])*(Diff[i])/(2*EnKF.MeaNoiseDensity*EnKF.MeaNoiseDensity)))
        else:
            ModelProb[modeIndex]=ModelProb[modeIndex]+log(1/(EnKF.MeaNoiseVelocity*sqrt(2*pi))*exp(-(Diff[i])*(Diff[i])/(2*EnKF.MeaNoiseVelocity*EnKF.MeaNoiseVelocity)))
    ModelProb[modeIndex]=exp(ModelProb[modeIndex])

    Sum=sum(ModelProb)
    for j in range(int(sum(C))):
        ModelProb[j]=ModelProb[j]/Sum
    return (PosteriorModelDensity,PosteriorModelSpeed,ModelProb)



                        
EnKF=KalmanEstimationModel()
TM=TrueModel()
GPS=GPSvehicle()


TimeStep=180
lag=1
pTran=0.99
Model=SystemModel(EnKF.cellNumber, EnKF.laneNumber)
ModelNumber=(EnKF.cellNumber-4)*(EnKF.laneNumber-1)+1


# define matrix to store the mean density and speed value
MeanDensity=zeros((ModelNumber,TimeStep,EnKF.cellNumber))
MeanSpeed=zeros((ModelNumber,TimeStep,EnKF.cellNumber))
TrueDensity=zeros((ModelNumber,TimeStep,EnKF.cellNumber))
TrueSpeed=zeros((ModelNumber,TimeStep,EnKF.cellNumber))

PosteriorModelDensity=zeros((ModelNumber,EnKF.EnNumber,EnKF.cellNumber))
PosteriorModelSpeed=zeros((ModelNumber,EnKF.EnNumber,EnKF.cellNumber))
ModelProb=zeros(ModelNumber)
EstimatedMode=zeros((TimeStep,EnKF.cellNumber))
EstimatedState=zeros((TimeStep,EnKF.cellNumber))
TrueMode=zeros((TimeStep,EnKF.cellNumber))
TrueState=zeros((TimeStep,EnKF.cellNumber))

ModelProbP=zeros((TimeStep,ModelNumber))

# Initialzation
MeanDensity[:,0,:]=average(EnKF.cellDensity,axis=0)
MeanSpeed[:,0,:]=average(EnKF.cellSpeed,axis=0)
TrueDensity[:,0,:]=TM.cellDensity
TrueSpeed[:,0,:]=TM.cellSpeed
EstimatedMode[0]=Model[0]
EstimatedState[0]=average(EnKF.cellDensity,axis=0)
TrueMode[0]=TM.mode
TrueState[0]=TM.cellDensity

# Initialize model likilihood
PModel=zeros(ModelNumber)
for i in range(int(ModelNumber)):
    if i==0:
        PModel[i]=1.0/ModelNumber
    else:
        PModel[i]=1.0/ModelNumber

ModelProbP[0]=PModel
for counter in range(1, TimeStep):
    print 'This is iteration',counter

    cellDensityCopy=EnKF.cellDensity.copy()
    cellSpeedCopy=EnKF.cellSpeed.copy()
    TcellDensityCopy=TM.cellDensity.copy()
    TcellSpeedCopy=TM.cellSpeed.copy()
    GPSpositionCopy=deepcopy(GPS.GPSposition)
    

    for modeIndex in range(0,int(ModelNumber)):
        mode=Model[modeIndex,:]
        
        EnKF.cellDensity=cellDensityCopy.copy()
        EnKF.cellSpeed=cellSpeedCopy.copy()
        TM.cellDensity=TcellDensityCopy.copy()
        TM.cellSpeed=TcellSpeedCopy.copy()
        GPS.GPSposition=deepcopy(GPSpositionCopy)
        
# forward prediction
        EnKF.updateDensity(counter,mode)
        EnKF.updateSpeed(mode)
        TM.updateDensity(counter)
        TM.updateSpeed()
        GPS.GPSpositionUpdate(counter)
    
# compute mean of ensembles
        MeanDensity[modeIndex,counter,:]=average(EnKF.cellDensity,axis=0)
        MeanSpeed[modeIndex,counter,:]=average(EnKF.cellSpeed,axis=0)

# compute Kalman gain
        kalmanGain=computeKalmanGain(EnKF.cellNumber,EnKF.EnNumber,GPS.GPSposition,EnKF.cellDensity,EnKF.cellSpeed,EnKF.MeaNoiseDensity, EnKF.MeaNoiseVelocity)

# state estimation
        StateEstimation(EnKF.cellDensity,EnKF.cellSpeed,EnKF.cellNumber,EnKF.EnNumber,GPS.GPSposition,TM.cellDensity,TM.cellSpeed,kalmanGain,EnKF.MeaNoiseDensity,EnKF.MeaNoiseVelocity)

# compute the likelihood of each model
        PosteriorModelDensity,PosteriorModelSpeed, ModelProb=ModelLikelihood(modeIndex)


        
    ModelProbSum=sum(ModelProb)
    ModelProb=ModelProb/ModelProbSum
    print ModelProb
    ModelProbP[counter]=ModelProb
    OptimalIndex=ModelProb.argmax()
    EnKF.cellDensity=PosteriorModelDensity[OptimalIndex]
    EnKF.cellSpeed=PosteriorModelSpeed[OptimalIndex]
    EstimatedMode[counter,:]=Model[OptimalIndex]
    EstimatedState[counter,:]=average(EnKF.cellDensity,axis=0)
#    print Model[OptimalIndex]

    TrueMode[counter,:]=TM.mode
    TrueState[counter,:]=TM.cellDensity
    
#### Plot section starts here
##save('EstimatedMode'+str(GPS.PR),EstimatedMode)
##save('EstimatedState'+str(GPS.PR),EstimatedState)
##save('TrueMode',TrueMode)
##save('TrueState',TrueState)
##
##plt.rc('xtick',labelsize=30)
##plt.rc('ytick',labelsize=30)
##imgplot=plt.imshow(EstimatedState,aspect='auto',origin='lower',interpolation='none')
##plt.ylabel('Time Step',fontsize=30)
##plt.xlabel('Cell Number',fontsize=30)
##plt.colorbar()
##plt.savefig('EstimatedState'+str(GPS.PR)+'.pdf',bbox_inches='tight')
##plt.show()
##plt.clf()
##
##
##
##plt.rc('xtick',labelsize=30)
##plt.rc('ytick',labelsize=30)
##imgplot1=plt.imshow(EstimatedMode,aspect='auto',origin='lower',interpolation='none')
##plt.ylabel('Time Step',fontsize=30)
##plt.xlabel('Cell Number',fontsize=30)
##plt.clim(1.0, 3.0)
##plt.colorbar()
##plt.savefig('EstimatedMode'+str(GPS.PR)+'.pdf',bbox_inches='tight')
##plt.show()
##plt.clf()

##plt.rc('xtick',labelsize=30)
##plt.rc('ytick',labelsize=30)
##imgplot=plt.imshow(TrueState,aspect='auto',origin='lower',interpolation='none')
##plt.ylabel('Time Step',fontsize=30)
##plt.xlabel('Cell Number',fontsize=30)
##plt.colorbar()
##plt.savefig('TrueState.pdf',bbox_inches='tight')
##plt.show()
##plt.clf()
##
##plt.rc('xtick',labelsize=30)
##plt.rc('ytick',labelsize=30)
##imgplot=plt.imshow(TrueMode,aspect='auto',origin='lower',interpolation='none')
##plt.ylabel('Time Step',fontsize=30)
##plt.xlabel('Cell Number',fontsize=30)
##plt.clim(1.0, 3.0)
##plt.colorbar()
##plt.savefig('TrueMode.pdf',bbox_inches='tight')
##plt.show()
##plt.clf()


plt.rc('xtick',labelsize=20)
plt.rc('ytick',labelsize=20)
plt.hold(True)
for i in range(int(ModelNumber)):
    if i==0:
        imgplot1=plt.plot(range(int(TimeStep)), ModelProbP[:,i],color='b',label='no incident')
    elif i==5:
        imgplot2=plt.plot(range(int(TimeStep)), ModelProbP[:,i],color='r',label='two lanes blocked at cell four')
    elif i==2:
        imgplot3=plt.plot(range(int(TimeStep)), ModelProbP[:,i],color='g',label='all other incidents')
    else:
        imgplot4=plt.plot(range(int(TimeStep)), ModelProbP[:,i],color='g')
plt.ylim([-0.05,1.05])
plt.legend(loc='center right')
plt.ylabel('Model Probability',fontsize=20)
plt.xlabel('Time Step',fontsize=20)
plt.savefig('1900nonIMMModelProb'+str(GPS.PR)+'.pdf',bbox_inches='tight')
plt.hold(False)
plt.show()
plt.clf()
