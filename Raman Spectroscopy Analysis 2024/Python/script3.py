#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 20:22:56 2023

Mat->CSV
@author: stephenevans



"""
import numpy as np
import sys 
import scipy.io as sio
from scipy import interpolate
import pandas as pd
#import pyqtgraph as pg
#from PyQt5.QtCore import Qt 
#from pyqtgraph.Qt import QtCore, QtGui,QtWidgets
from PyQt5 import uic
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import ( QFileDialog,QMainWindow)
from sklearn.metrics import auc
import matplotlib
matplotlib.use('Qt5Agg')

from matplotlib.figure import Figure 
from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT as NavigationToolbar)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas


# set-up graphs - 2 figures

class MplCanvas(FigureCanvas):
    def __init__(self,parent=None,width=5,height=4,dpi=150):
               
        fig = Figure(figsize=(width, height), dpi=dpi)         
        self.axes = fig.add_subplot(211)            
        self.axes2 = fig.add_subplot(212)         
        
                            
        fig.tight_layout(pad=5)
        #print ('in canvas')
        super().__init__(fig) 
        #super(MplCanvas(),self).__init__(fig)
        
class MainWindow(QMainWindow):      
            
    def __init__(self):     
        super().__init__()      
        QMainWindow.__init__(self)  
        self.ui=uic.loadUi('Mat_to_CSV_GUI.ui',self)   
        self.show()
        
        
        global fileCount,MessString,df2,df3,dframe,classDict, cols, cols2,color_by_name
        color_by_name=['tab:blue','tab:orange','tab:green','tab:red','tab:purple',
                  'tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
        classDict={} # used to store data and labels
               
        fileCount=1 # number of classes / different tyoes of cell / material
        MessString="Report:"
        self.spinBox.setSingleStep(1) # set in gui set-up next time compile
        self.setWindowTitle("PCA Analysis")
        self.pushButton_3.clicked.connect(self.getfiles) 
        self.pushButton_5.clicked.connect(self.end)
        self.pushButton_6.clicked.connect(self.clearAll)
        self.pushButton_9.clicked.connect(self.save)
    
        self.pushButton_8.setCheckable(True)
        self.pushButton_8.toggle()
        
        spc=self.spinBox.value()
        self.pushButton_8.clicked.connect(self.spectrum_to_display)
        self.pushButton_8.toggle()

        self.spinBox.valueChanged.connect(self.spectrum_to_display)
        #self.spinBox_3.valueChanged.connect(self.num_outliers)
        
        # Message box
        MessString="Messages\n"          
        self.textEdit.insertPlainText(MessString)
      

        # plot set up
        self.canvas=MplCanvas(self, width=8, height=9, dpi=100)        
        self.ui.verticalLayout_2.addWidget(self.canvas)      
        self.navi_toolbar = NavigationToolbar(self.canvas, self) #createa navigation toolbar for our plot canvas
        self.verticalLayout_2.addWidget(self.navi_toolbar)
        
        # it is functional but certain elements not yet implemented - eg dark mode, toggling grids, saving outputs
        Message="\n Step1. Open *.mat file\n Step2. Use spinbox to select outlier (red curve) on lower graph\n Step3. Hit the Zap button\
  button to remove outlier\n Step4. Repeat steps 2 & 3 \n Step5. Press OK to save, spline interpolated, spectra *.csv\n\n"
        self.textEdit.insertPlainText(Message)
        
        
    def save(self):     
         #print("old file name" ,fnam)
         temp=fnam.replace('.mat', '.csv')
         fname=temp
         print("new file name:",fname)
         print("df3=",len(df3.columns))
         #remove mean and std before saving
         if "mean" in df3.columns:
             del df3['mean']
             del df3['std']
         print("df3=",len(df3.columns),ncols3)
        
         col_name3=[]
         # for Raman spectra - then repeat for background
         #rename columns to account for those removed
         for spec in range(ncols3):
            if spec==0:
                col_name3=["Wavenumber"]
            else:
                colName="Data"+str(spec)            
                col_name3.append(colName)  
         df3.columns=col_name3
         df3.to_csv(fname,index=False)
         print("saved data")    
     
    def spectrum_to_display(self):  
        self.spinBox.setMaximum(ncols)         
        spc=self.spinBox.value()
        self.zap(spc)    
        
        
    def zap(self,spc):
        global df3, col_name3,ncols3
        
        #self.spinBox.setMaximum(ncols)         
        #spc=self.spinBox.value()
        
        
        
        
        if "mean" in df3.columns:
            del df3['mean']
            del df3['std']
        
        #
        x=df2["Wavenumber"]
        #df2['mean'] = df2.iloc[:, 1:ncols].mean(axis=1)
        #df2['std'] = df2.iloc[:, 1:ncols].std(axis=1)
        
        y=df2['mean']
        dy=df2["std"]
        #self.canvas.axes.cla()
        
        self.dataplot(x,y,dy,fileCount)
        
      
        num=spc-1
        
        
        outlier=sortedArr[[num]]
        outlier=int(np.take(outlier,0))
        outlierColumn="Data"+str(outlier)
        
        print("sb value=",spc," df2 column =",outlierColumn)
      
        #MessString="spectrum ("+str(spc)+")="+outlierColumn+"/"+str(ncols)+"\n"      
        #self.textEdit.insertPlainText(MessString)
        
        
        y2=df2[outlierColumn]
        self.dataplot2(x,y2,fileCount)
        
        
        ncols3=len(df3.columns)
        print("df3 length", ncols3,"\n")
        if self.pushButton_8.isChecked():
           
           self.pushButton_8.toggle()
           
           #print('zapping column',outlier)
          
           df3=df3.drop(columns=outlierColumn)
           ncols3=len(df3.columns)            
           print('zapping column',outlier," df3 length  = ",ncols3," columns")
           MessString="Zapping "+outlierColumn +"\n"     
           self.textEdit.insertPlainText(MessString)
           
           
          
        x=df3["Wavenumber"]
        #df3['mean'] = df3.iloc[:, 1:ncols3].mean(axis=1)
        #df3['std'] = df3.iloc[:, 1:ncols3].std(axis=1)
        
        #y=df3['mean']
        #dy=df3["std"]
        y= df3.iloc[:, 1:ncols3].mean(axis=1)
        dy=df3.iloc[:, 1:ncols3].std(axis=1)
        self.dataplot3(x,y,dy)
        
        
        
          
           
    
    
    def getfiles(self):
        global col_name,ncols,fnam, corrFactor, dataArray,dataSplineArray,sortedArr,df2,df3                   
        dlg = QFileDialog()  
        dlg.setNameFilter('(*.mat)')
        dlg.setFileMode(QFileDialog.AnyFile)
        
        df2 = pd.DataFrame() 
        if dlg.exec_():
             
             filenames = dlg.selectedFiles()
             fnam=filenames[0]
             print("open",fnam)
             
             sep="/"
             fnamLen=len(fnam)
             lastSep=fnam.rfind(sep)
             #pathName=fnam[:lastSep+1]  # identify path
             fileName=fnam[lastSep+1:fnamLen]  
             MessString=fileName+"\n"          
             self.textEdit.insertPlainText(MessString)
             
            
             
             self.pushButton_3.setEnabled(False) 
             mat_contents = sio.loadmat(fnam)
             print(mat_contents.keys())
        
             for key in mat_contents.keys():

                 if key=="__header__":
                     header=mat_contents['__header__']
                     print("____________header_____________")     
                     #print(header)    
                     print("     ")               
                 if key=="__version__":
                    version=mat_contents['__version__']
                    print("____________version_____________")   
                    #print(version)
                    print("     ")    
                 if key=="__globals__":
                    glob=mat_contents['__globals__']
                    print("____________globals____________")   
                    #print(glob)    
                    print("     ")                    
                 if key=='CalibrationPeakTemp':
                    cal_peak=mat_contents['CalibrationPeakTemp']
                    cal_peak_val=cal_peak[0,0]
                    print("____________calibration peak_____________")   
                    print(cal_peak_val)
                    print("     ")    
                 if key=='RamanBackgroundSpectra':
                     RamanBack=mat_contents[key]
                     print("____________Background____________")   
                     numBackSpectra=RamanBack.size
                     print("numberRamanBack",numBackSpectra)
                     print("     ")    
                 if key=='RamanSpectra':   
                     RamanSpectra=mat_contents[key]
                     numSpectra=RamanSpectra.size
                     
                     print("____________Spectra_____________")   
                     print("number of spectra=",numSpectra)
                     print("     ")    
            

             corrFactor=round((520.5-cal_peak_val),1)
             print("correction factor", corrFactor) 
             col_name=["Wavenumber"]
             # for Raman spectra - then repeat for background
             for spec in range(numSpectra):
                data=RamanSpectra[0,spec]
                lnData=len(data)
                colName="Data"+str(spec+1)
                
                col_name.append(colName)    
                data_x=data[:,0]+corrFactor
                data_y=data[:,1]
                
                #translate in y  so min value=0
                data_y_min=np.amin(data_y)
                data_y=data_y-data_y_min
                
                if data_x[0]>data_x[-1]:
                   #print("inverting data",spec) to make low to high wavenumber
                   data_x=data_x[::-1]
                   data_y=data_y[::-1]
                
               #take x data to 1 dp
                test_x=np.round(data_x,decimals=1)
                
                
                if spec==0:        
                #write x,y    
                    dataArray=np.column_stack((data_x,data_y))
                    xnew,ynew=self.spliInt(spec)
                    dataSplineArray=np.column_stack((xnew,ynew))
                else:
                #write y                    
                    dataArray=np.column_stack((dataArray,data_y))
                    xnew,ynew=self.spliInt(spec)
                    dataSplineArray=np.column_stack((dataSplineArray,ynew))

             #  create data frame at end from the arrays - this stops fragmentation
             df= pd.DataFrame(dataArray) 
             df2= pd.DataFrame(dataSplineArray) 
             df.columns=col_name
             df2.columns=col_name
             df2=df2.round(1)        
             #print(df)
             
             
             df3 = pd.DataFrame()
             df3=df2.copy(deep=True) 
             flag="showAll"
             ncols=numSpectra
             
             x=df2["Wavenumber"]
             df2['mean'] = df2.iloc[:, 1:ncols].mean(axis=1)
             df2['std'] = df2.iloc[:, 1:ncols].std(axis=1)
             
             y=df2['mean']
             dy=df2["std"]
             self.dataplot(x,y,dy,fileCount)
             print(df2)
             
             #determine outliers
             col_name=df2.columns
             #print(col_name)
             # calc area under curve using skilearn auc function
             meanArea=round(auc(x,y),2)
             MessString="Area under mean = "+str(meanArea)+"\n"          
             self.textEdit.insertPlainText(MessString)
             
             #find outlier spectra
             colArea=[]
             iCount=[]
             I=0
             for col in col_name:
                 if col=="Wavenumber":
                     x=df2.Wavenumber
                 elif col=="mean":
                     pass
                 elif col=="std":
                     pass
                 else:
                     I=I+1
                     y=df2[col]
                     newArea=abs(round(auc(x,y)-meanArea,2))
                     colArea.append(newArea)
                     iCount.append(I)
             areaArray=np.column_stack((iCount,colArea))         
             #print (areaArray)
             print(" ")
             sortedArr = areaArray[areaArray[:,1].argsort()]   
             
             print("______ Sorted_____")
             sortedArr=sortedArr[::-1]
             print (sortedArr)  
             
             #show 3 largest to zap
             
             out_first=sortedArr[[0]]
             out_first=int(np.take(out_first,0))
             
             print("first=",out_first)
             first="Data"+str(out_first)
             y2=df2[first]
             self.dataplot2(x,y2,fileCount)
             MessString="Outlier = "+str(out_first)+"\n"          
             self.textEdit.insertPlainText(MessString)
             


         
            
    def spliInt(self,spect):
         spect=spect+1
         #Spline Interpolation 
         x=dataArray[:,0]
         y=dataArray[:,spect]                 
         xStart=2*round(x[0]/2);xEnd=2*round(x[-1]/2)#1947
        # print("xstart,xstop",xStart,xEnd)          
         s=interpolate.InterpolatedUnivariateSpline(x, y)
         xnew=np.arange(xStart,xEnd,0.5) #+ corrFactor
         ynew=s(xnew)         
         xnew=xnew[::-1]
         ynew=ynew[::-1]
         return xnew,ynew       
           
#23456       
    def dataplot(self,x,y,dy,fileCount):
        #Average Data 
        
        if fileCount==1:
            self.canvas.axes.cla()
            self.canvas.axes.set_xlabel("Wavenumber (cm-1)",fontsize=12)
            self.canvas.axes.set_ylabel("Intensity (arb. units)",fontsize=12)
            self.canvas.axes.set_title("Raman Spectra",fontsize=14) 
        #offset_ratio=1.5   
        #y=y+(fileCount*offset_ratio)    # offset spectra
        pcol=color_by_name[fileCount-1]        
        self.canvas.axes.plot( x,y,color=pcol,linewidth=1.5,linestyle='-',alpha=1.0)  
        self.canvas.axes.fill_between(x,y+dy,y-dy,facecolor=pcol,linewidth=1.5,alpha=0.4)  
        
        cb2=self.checkBox_2.isChecked()
        if cb2==True:
            self.canvas.axes.grid(color = 'gray', linestyle = '--', linewidth = 0.5)
        #self.canvas.axes.legend()    
        self.canvas.draw () 
        
    def dataplot2(self,x,y,fileCount):
        #Average Data 
        
        if fileCount==1:
            self.canvas.axes2.cla()
            self.canvas.axes2.set_xlabel("Wavenumber (cm-1)",fontsize=12)
            self.canvas.axes2.set_ylabel("Intensity (arb. units)",fontsize=12)
            self.canvas.axes2.set_title("Raman Spectra",fontsize=14) 
        
          
        pcol=color_by_name[3]        
        self.canvas.axes.plot( x,y,color=pcol,linewidth=1.5,linestyle='-',alpha=1.0)  
       
        
        cb2=self.checkBox_2.isChecked()
        if cb2==True:
            self.canvas.axes2.grid(color = 'gray', linestyle = '--', linewidth = 0.5)
        #self.canvas.axes.legend()    
        self.canvas.draw ()   
        
    def dataplot3(self,x,y,dy):
        #Average Data 
                
        self.canvas.axes2.cla()
        self.canvas.axes2.set_xlabel("Wavenumber (cm-1)",fontsize=12)
        self.canvas.axes2.set_ylabel("Intensity (arb. units)",fontsize=12)
        self.canvas.axes2.set_title("Raman Spectra after Zapping",fontsize=14) 
        
        pcol=color_by_name[fileCount-1]        
        self.canvas.axes2.plot( x,y,color=pcol,linewidth=1.5,linestyle='-',alpha=1.0)  
        self.canvas.axes2.fill_between(x,y+dy,y-dy,facecolor=pcol,linewidth=1.5,alpha=0.4)  
        
        cb2=self.checkBox_2.isChecked()
        if cb2==True:
            self.canvas.axes2.grid(color = 'gray', linestyle = '--', linewidth = 0.5)
        #self.canvas.axes.legend()    
        self.canvas.draw ()    
        
        
        
        
        
        
#23456          
    def clearAll(self):
        print("Quit")
        app.quit()    
          
    def end(self):
        print("Quit")
        app.quit()
        
        
if __name__=="__main__":
    app=QtWidgets.QApplication(sys.argv)
    w=MainWindow()
    w.show()
    app.exec_()
