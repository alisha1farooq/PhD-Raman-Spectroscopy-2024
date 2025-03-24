#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 22/4/2023  - Raman Processing New 
based on old Raman processing - but modified to use mat plot lib rather than pyqt graphs.

to use Modpoly, IModpoly and Zhang baseline need to Install the BaselineRemoval 
library as pip install BaselineRemoval. 

used to  read a raw.csv and create an analysed.csv  - to be ready for PCA/LDA

@author: stephenevans

"""
import numpy as np

import sys 
import pandas as pd

from PyQt5 import uic
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import ( QFileDialog,QMainWindow)

import matplotlib
matplotlib.use('Qt5Agg')
from scipy.signal import savgol_filter
from matplotlib.figure import Figure 
from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT as NavigationToolbar)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

from BaselineRemoval import BaselineRemoval  
from scipy import sparse
from scipy.sparse.linalg import spsolve



class MplCanvas(FigureCanvas):
    def __init__(self,parent=None,width=5,height=4,dpi=300):
               
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
        self.ui=uic.loadUi('Raman_Processing_GUI.ui',self)   
        self.show()
        
        
        global fileCount,MessString,df2,dframe,classDict, cols, cols2,color_by_name
        color_by_name=['tab:blue','tab:orange','tab:green','tab:red','tab:purple',
                  'tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
        classDict={} # used to store data and labels
               
        fileCount=0 # number of classes / different tyoes of cell / material
        MessString="Report:"
       
        self.setWindowTitle("Raman Processing")
        self.pushButton_3.clicked.connect(lambda:self.getfiles()) 
        
        #limits for truncate
        
        self.pushButton_8.clicked.connect(lambda:self.truncate())        
        self.pushButton_7.clicked.connect(lambda:self.preview_baseline())
        self.pushButton_5.clicked.connect(lambda:self.apply_to_all())
        self.pushButton_6.clicked.connect(lambda:self.preview_smooth())
        self.pushButton_9.clicked.connect(lambda:self.smooth_all())        
        self.spinBox_3.valueChanged.connect(lambda:self.preview_smooth()) 
        self.spinBox_8.valueChanged.connect(lambda:self.preview_smooth())         
        self.pushButton.clicked.connect(lambda:self.normalise())
        self.pushButton_2.clicked.connect(lambda:self.mean_sd())        
        self.pushButton_4.clicked.connect(lambda:self.save())
      
    
        
        # Message box
        MessString="Files Loaded (and colour)\n"          
        self.textEdit.setPlainText(MessString)
        
        # plot set-up
        self.canvas=MplCanvas(self, width=8, height=9, dpi=100)        
        self.ui.verticalLayout_2.addWidget(self.canvas)      
        self.navi_toolbar = NavigationToolbar(self.canvas, self) #createa navigation toolbar for our plot canvas
        self.verticalLayout_2.addWidget(self.navi_toolbar)
        
        # it is functional but certain elements not yet implemented - eg dark mode, toggling grids, saving outputs
       
    
    def save(self):    
    #fname=fnam+"analysed"
        temp=fnam
        print("old file name",temp)
        temp2=temp.replace('raw', 'analysed')
        print("new file name",temp2)
        
        df4.to_csv(temp2,index=False)
        text=str(temp2)   
        self.textEdit.setPlainText(str(text))       
    


    
    
    def mean_sd(self):
        df4['mean'] = df4.iloc[:, 1:ncols].mean(axis=1)
        df4['SD1'] = df4.iloc[:, 1:ncols].std(axis=1)
        x=df4.Wavenumber
        y=df4['mean']
        dy=df4['SD1']
        flag=3
        self.dataplot3(x,y,dy,flag)
        
   
    def normalise (self):  
        # zoom in to noise and to sharp peak to optimize smooth without reducing peak height significantly
        print("in normalise")     
        xval=self.doubleSpinBox_2.value() # get value to normalise
        xpos=(df4[df4["Wavenumber"]==xval].index.values)
        xpos=np.take(xpos,0)
        
        text="Normalised against peak at "+ str(xval)+" (cm-1) \n"
        print(text)
        self.textEdit.setPlainText(text)
        
        
        for coln in col_name:
            if coln=='Wavenumber':
                pass
            else:
                y=df4[coln]
                ynorm=y[xpos]
                #print("x,ynorm",xpos,ynorm)
                
                ynew=y/ynorm
                df4[coln]=ynew
        self.dataplot(df4)
       
        
       
        
    def preview_smooth(self):  
        # zoom in to noise and to sharp peak to optimize smooth without reducing peak height significantly
        print("in smooth preview")  
        polynom=self.spinBox_3.value()     # polynomial
        wind=self.spinBox_8.value()         # window
        print("poly= ",polynom,"window= ",wind, "(w/p)=", str(wind/polynom))
        
        x=df3["Wavenumber"]
        y=df3["Data1"]
        flag=0
        self.dataplot2(x,y,flag) 
        
        Iny_SGS = savgol_filter(y, wind, polynom) # window size 51, polynomial order 3
        flag=1
        self.dataplot2(x,Iny_SGS,flag)           
        
    def smooth_all(self):      
        global df4
        print("in smooth all")  
        polynom=self.spinBox_3.value()     # polynomial
        wind=self.spinBox_8.value()  # window
        text="Smooth: poly= "+str(polynom)+" window= "+str(wind)+ " (w/p)= "+ str(wind/polynom)+"\n"
        print(text)
       
        self.textEdit.setPlainText(text)
        
        df4 = pd.DataFrame()# put smoothed data here
        
        
        for col in col_name:
           if col=="Wavenumber":
               x=df3.Wavenumber
               df4[col]=x
           else:
               y=df3[col]  
               Iny_SGS = savgol_filter(y, wind, polynom) # typical window size 15, polynomial order 2  note higher w/p ratio gives stronger smoothing
               df4[col]=Iny_SGS
        
        #show effect for smoothing on a single spectrum
        x=df4["Wavenumber"]
        y=df4["Data1"]
        flag=0
        self.dataplot2(x,y,flag) 
        flag=1
        x=df3["Wavenumber"]
        y=df3["Data1"]
        self.dataplot2(x,y,flag) 
        
        
        
        
        
        
        
    def preview_baseline(self):
        #first pick method from combo box
        choice=self.comboBox.currentIndex()
        choice_txt=self.comboBox.currentText()
        txt="baseline choice= "+ choice_txt+"\n"
        print(txt)            
        self.textEdit.insertPlainText(txt)
        
        
        # plot average and work out best base line then apply to all
        
        x=df2["Wavenumber"]       
        y = df2.iloc[:, 1:ncols].mean(axis=1)
        flag=0
        self.dataplot2(x,y,flag)
   
        if choice==0:
            #ALS fit
            #Asymmetric Least Squares Smoothing" by P. Eilers and H. Boelens in 2005. 
            
            txt="Asymmetric Least Square \n"
            print(txt)            
            self.textEdit.insertPlainText(txt)
           
            self.textEdit.insertPlainText("uses spin box values N, lambda, and p \n")
            
           
            lam=self.spinBox_6.value()
            p=self.doubleSpinBox.value()
            niter=self.spinBox_5.value()
            lam=(lam*1000)+1E3
            p=p/1000
            print("n. iterations=", niter," lambda=",lam, " p= ",p)
            y_base=self.baseline_als_optimized(y, lam, p, niter)
            y_res=y-y_base
            flag=1
            self.dataplot2(x,y_base,flag)
            flag=2
            self.dataplot2(x,y_res,flag)
        
        if choice==1:
            #IMOD poly
            #IMOD fit
            p=int(self.doubleSpinBox.value())
            if p<2:
                p=2
                #dont let p be less than 2
            niter=10*self.spinBox_5.value()   
            
            
            grad=0.001
            txt="IMOD Poly (degree, iterations, gradient) = "+str(p)+" "+str(niter)+" " +str(grad)+" \n"
            print(txt)            
            self.textEdit.insertPlainText(txt)
            baseObj=BaselineRemoval(y)
            y_res=baseObj.IModPoly(degree=p,repitition=niter,gradient=grad)
            y_base=y-y_res
            flag=1
            self.dataplot2(x,y_base,flag)
            flag=2
            self.dataplot2(x,y_res,flag)
        
        
        if choice==2:
             #Zhang
            txt="Zhang \n"
            print(txt)            
            self.textEdit.insertPlainText(txt)
           
            self.textEdit.insertPlainText("uses spin box values N, lambda \n") 
            
            niter=self.spinBox_5.value()
            lam=self.spinBox_6.value()
            
            baseObj=BaselineRemoval(y)
            y_res=baseObj.ZhangFit(lambda_=lam, porder=1, repitition=niter)
            y_base=y-y_res
            flag=1
            self.dataplot2(x,y_base,flag)
            flag=2
            self.dataplot2(x,y_res,flag)
            
            
            
  
        
        
#2345        
    def apply_to_all(self):  
        global df3
        temp_array=df2.to_numpy()
        x=df2["Wavenumber"]  
        
        choice=self.comboBox.currentIndex() 
        if choice==0: #ALS fit
            lam=self.spinBox_6.value()
            p=self.doubleSpinBox.value()
            niter=self.spinBox_5.value()
            lam=(lam*1000)+1E3
            p=p/1000
                        
            #x=df2["Wavenumber"]  
            temp_array[:,0]=x            
            for i in range(1,ncols):                    
                    y=temp_array[:,i]
                    yb=self.baseline_als_optimized(y, lam, p, niter)
                    yr=y-yb
                    temp_array[:,i]=yr
            #create new dataframe of baseline corrected spectra        
            df3= pd.DataFrame(temp_array) 
            df3.columns=col_name
            outText="Parameters: N;lambda; poly "+str(niter)+";"+ str(lam)+"; "+str(p)+" \n"
            
            self.textEdit.setPlainText(outText)
            self.dataplot(df3)
            
        if choice==1:
            #IMOD poly
            #IMOD fit
            p=int(self.doubleSpinBox.value())
            if p<2:
                p=2
                #dont let p be less than 2
            niter=10*self.spinBox_5.value()   
            
            
            grad=0.001
            txt="IMOD Poly (degree, iterations, gradient) = "+str(p)+" "+str(niter)+" " +str(grad)+" \n"
            print(txt)            
            self.textEdit.insertPlainText(txt)
            temp_array[:,0]=x
            for i in range(1,ncols):                    
                    y=temp_array[:,i]
                    baseObj=BaselineRemoval(y)
                    yr=baseObj.IModPoly(degree=p,repitition=niter,gradient=grad)
                    temp_array[:,i]=yr
                 #create new dataframe of baseline corrected spectra        
            df3= pd.DataFrame(temp_array) 
            df3.columns=col_name
            self.dataplot(df3)
        
        if choice==2: #Zhang
            niter=self.spinBox_5.value()
            lam=self.spinBox_6.value()
            temp_array[:,0]=x
            for i in range(1,ncols):                    
                    y=temp_array[:,i]
                    baseObj=BaselineRemoval(y)
                    yr=baseObj.ZhangFit(lambda_=lam, porder=1,repitition=niter)
                    temp_array[:,i]=yr
                 #create new dataframe of baseline corrected spectra        
            df3= pd.DataFrame(temp_array) 
            df3.columns=col_name
            outText="Parameters: N;lambda; poly "+str(niter)+";"+ str(lam)+"; "+str(1)+" \n"
            
            outText="Baseline: Zhang  \n"
            self.textEdit.setPlainText(outText)
            self.dataplot(df3)
        
        
        
        
        
        
        
        
        
        
        
    def baseline_als_optimized(self, y, lam, p, niter):
          L = len(y)
          D = sparse.diags([1,-2,1],[0,-1,-2], shape=(L,L-2))
          D = lam * D.dot(D.transpose()) # Precompute this term since it does not depend on `w`
          w = np.ones(L)
          W = sparse.spdiags(w, 0, L, L)
          for i in range(niter):
              W.setdiag(w) # Do not create a new matrix, just update diagonal values
              Z = W + D
              z = spsolve(Z, w*y)
              w = p * (y > z) + (1-p) * (y < z)
          return z           
        
        
       
    def truncate(self):     
         global df2,outText
         xmin=self.spinBox.value()
         xmax=self.spinBox_2.value()
         print("truncate between limits",xmin,xmax)          
         MessString="Truncated between "+str(xmin)+" and "+str(xmax)+"\n"                   
         self.textEdit.setPlainText(str(MessString))
        
         imin=(df[df["Wavenumber"]==xmin].index.values)
         imax=(df[df["Wavenumber"]==xmax].index.values)
         imin=np.take(imin,0)
         imax=np.take(imax,0)
         x=df['Wavenumber']
         newData=df.loc[imax:imin]   
         df2 = pd.DataFrame(newData) 
         df2.reset_index(drop=True, inplace=True)
         
         self.dataplot(df2)
        
        
             
        
    def getfiles(self):        
       global data,col_name,df,ncols,fnam,flag,outText
        
       dlg = QFileDialog()   
       dlg.setNameFilter('(*.csv)')
       dlg.setFileMode(QFileDialog.AnyFile)
       fileNameList=[]
    
     
       if dlg.exec_():
            filenames = dlg.selectedFiles()
            fnam=filenames[0]
            
            # get file & path name  (might only work for Mac i.e need to change separator?) i.e. if windows use "\"
            sep="/"
            fnamLen=len(fnam)
            lastSep=fnam.rfind(sep)
            #pathName=fnam[:lastSep+1]
            fileName=fnam[lastSep+1:fnamLen]  
            
            
            outText="File: "+str(fileName)+"\n"
            self.textEdit.setPlainText(outText)
            
            data= pd.read_csv (filenames[0])        
            ncols=len(data.columns)     
            print("number of columns=",ncols)
            col_name=[]
            df = pd.DataFrame(data) #create pandas dataframe called data
        
         # name the columns as wavenumber or data : first create a list of names then add to data frame
       for i in range(ncols):
            if i==0:
                col_name.append('Wavenumber')      
            else:
                col_name.append("Data"+str(i))
    
       df.columns=col_name  # add column names into df
       print(df)
      
       self.dataplot(df)
        
       x=df["Wavenumber"]
       xmax=int(max(x))
       xmin=int(min(x))
       print(xmin,xmax)
       
       self.spinBox.setValue(xmin)
       self.spinBox_2.setValue(xmax)
             
                
    
       
    
    def dataplot2(self,x,y,flag):
      
        pcol=color_by_name[flag]    
        if flag==0:
            print("in flag ==0")
            self.canvas.axes.cla()
            self.canvas.axes.set_xlabel("Wavenumber (cm-1)",fontsize=12)
            self.canvas.axes.set_ylabel("Intensity (arb. units)",fontsize=12)
            self.canvas.axes.set_title("Raman Spectra",fontsize=14) 
            pcol=color_by_name[flag]    
            self.canvas.axes.plot( x,y,color=pcol,linewidth=1.5,linestyle='-',alpha=1)  
                       
        elif flag==1:
            print("in flag ==1")
            pcol=color_by_name[flag]    
            self.canvas.axes.plot( x,y,color=pcol,linewidth=1,linestyle='--',alpha=1)  
        elif flag==2:
            pcol=color_by_name[flag-1]    
            self.canvas.axes.plot( x,y,color=pcol,linewidth=1.5,linestyle='-',alpha=1)  
            
            
        cb2=self.checkBox_2.isChecked()
        if cb2==True:
            self.canvas.axes.grid(color = 'gray', linestyle = '--', linewidth = 0.5)
        
        self.canvas.draw ()
        
       
    
    def dataplot(self,df):
        #Average Data 
              
        self.canvas.axes.cla()
        self.canvas.axes.set_xlabel("Wavenumber (cm-1)",fontsize=12)
        self.canvas.axes.set_ylabel("Intensity (arb. units)",fontsize=12)
        self.canvas.axes.set_title("Raman Spectra",fontsize=14) 
        
        pcol=color_by_name[0]    
        for col in col_name:
            if col=="Wavenumber":
                x=df.Wavenumber
            else:
                y=df[col]                   
                self.canvas.axes.plot( x,y,color=pcol,linewidth=1.5,linestyle='-',alpha=.5)  
          
        
        cb2=self.checkBox_2.isChecked()
        if cb2==True:
            self.canvas.axes.grid(color = 'gray', linestyle = '--', linewidth = 0.5)
        #self.canvas.axes.legend()    
        self.canvas.draw ()
        
        
    def dataplot3(self,x,y,dy,flag):
        #Average Data 
        pcol=color_by_name[flag]    
       
        self.canvas.axes2.cla()
        self.canvas.axes2.set_xlabel("Wavenumber (cm-1)",fontsize=12)
        self.canvas.axes2.set_ylabel("Intensity (arb. units)",fontsize=12)
        self.canvas.axes2.set_title("Raman Spectra ",fontsize=14)       
        self.canvas.axes2.plot( x,y,color=pcol,linewidth=1.5,linestyle='-',alpha=1.0)  
        self.canvas.axes2.fill_between(x,y+dy,y-dy,facecolor=pcol,linewidth=1.5,alpha=0.4)  
        
        
        cb2=self.checkBox_2.isChecked()
        if cb2==True:
            self.canvas.axes2.grid(color = 'gray', linestyle = '--', linewidth = 0.5)
        #self.canvas.axes.legend()    
        self.canvas.draw ()    
    
    
        
    
        
if __name__=="__main__":
    app=QtWidgets.QApplication(sys.argv)
    w=MainWindow()
    w.show()
    app.exec_()

