#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 15/4/2023  - PCA Analysis of Raman Spectra

Files input are analysed.csv files 

1) read in  multiple CSV files - one at a time (though could be changed to read all at once quite easily) these should be files that are "analyzed"
2) display data - the 'mean' data column for each data class displayed
3) truncate data so that all spectra cover the same wavenumber regime and cover the regions of interest
4) perform SNV 
5) PCA anaylsis using sk learn

6) display outputs
7) save in way suitable for publication - can produce tiff or other format output directly clicking on navigation bar
at the bottom of the plots - this also allows change colours. lines, etc 

Use PyQT & pyqtgrpah as framework for program & displays - could be altered for Matplot lib  which would give better histograms

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

from matplotlib.figure import Figure 
from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT as NavigationToolbar)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas


from sklearn.decomposition import PCA 
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.preprocessing import StandardScaler
#from sklearn import svm
#from sklearn.cluster import KMeans
from sklearn.model_selection import train_test_split, cross_val_score
#from collections import Counter





class MplCanvas(FigureCanvas):
    def __init__(self,parent=None,width=5,height=4,dpi=300):
               
        fig = Figure(figsize=(width, height), dpi=dpi)         
        self.axes = fig.add_subplot(221)    
        
        self.axes2 = fig.add_subplot(222)         
        self.axes3 = fig.add_subplot(223,projection='3d')    # for 3D plot 
        #self.axes3 = fig.add_subplot(223)     # for 2D plot
        self.axes4 = fig.add_subplot(224)  
                            
        fig.tight_layout(pad=5)
        #print ('in canvas')
        super().__init__(fig) 
        #super(MplCanvas(),self).__init__(fig)

class MainWindow(QMainWindow):      
            
    def __init__(self):     
        super().__init__()      
        QMainWindow.__init__(self)  
        self.ui=uic.loadUi('Raman_LDA_GUI.ui',self)   
        self.show()
        
        
        global fileCount,MessString,df2,dframe,classDict, cols, cols2,color_by_name
        color_by_name=['tab:blue','tab:orange','tab:green','tab:red','tab:purple',
                  'tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
        classDict={} # used to store data and labels
               
        fileCount=0 # number of classes / different tyoes of cell / material
        MessString="Report:"
       
        self.setWindowTitle("LDA Analysis")
        self.pushButton_3.clicked.connect(lambda:self.getfiles()) 
        self.pushButton_5.clicked.connect(lambda:self.LDA_fit())
        self.pushButton_7.clicked.connect(lambda:self.PCA_LDA_fit())
        self.pushButton_9.clicked.connect(lambda:self.save())
        self.pushButton_8.clicked.connect(lambda:self.truncate())
    
        #limits for truncate
        self.spinBox.valueChanged.connect(lambda:self.valueChange)
        self.spinBox_2.valueChanged.connect(self.valueChange)
        # number of ocmponents to use
        self.spinBox_8.valueChanged.connect(self.num_comp_PCA)
        self.spinBox_3.valueChanged.connect(self.num_comp_LDA)
        
        # changed axis so PC1 is on x-axis
        #self.spinBox_5.setValue(2)
        #self.spinBox_6.setValue(1)
        #self.spinBox_5.valueChanged.connect(self.pca_fit)
        #self.spinBox_6.valueChanged.connect(self.pca_fit)
        
        # Message box
        MessString="Files Loaded (and colour)\n"          
        self.textEdit.insertPlainText(MessString)
        
        # plot set-up
        self.canvas=MplCanvas(self, width=8, height=9, dpi=100)        
        self.ui.verticalLayout_2.addWidget(self.canvas)      
        self.navi_toolbar = NavigationToolbar(self.canvas, self) #createa navigation toolbar for our plot canvas
        self.verticalLayout_2.addWidget(self.navi_toolbar)
        
        # it is functional but certain elements not yet implemented - eg dark mode, toggling grids, saving outputs
        
        
     
        
     
     # need to add button to implement clear all   
    def clearAll(self):
        self.canvas.axes.cla()
        self.canvas.axes2.cla()
        self.canvas.axes3.cla()
        self.canvas.axes4.cla()
        self.canvas.draw ()
        self.pushButton_8.setEnabled(True)           
        
        # need to clear all relevant data - at present best to simply re-run program
       
    def truncate(self):     
         global XvaluesTrunc         
         xlowerLimit=self.spinBox.value()            
         xupperLimit=self.spinBox_2.value()   
         print("truncate between limits",xlowerLimit,xupperLimit) 
         
         MessString="Truncated between "+str(xlowerLimit)+" and "+str(xupperLimit)+"\n"          
         self.textEdit.insertPlainText(MessString)
         
         
         
         #data prep truncate for all data have same stop and start
         #print("Dictionary=\n",classDict.keys())
         df_temp=pd.DataFrame()
         
         for ky in classDict:
             df_temp = classDict[ky]
             alldataArray=[]
             
             col_name=df_temp.columns
             wnbr="Wavenumber"
             
             #find indices of upper and lower limits                
             imin=(df_temp[df_temp[wnbr]==xlowerLimit].index.values)
             imax=(df_temp[df_temp[wnbr]==xupperLimit].index.values)
         
             imin=np.take(imin,0)
             imax=np.take(imax,0)    
             print("")
             #print("data set:",ky,imin,imax)
              
             # create series / array of truncated x values 
             XvaluesTrunc=xFullSpectrum.truncate(imax,imin-1,copy=False);
             XvaluesTrunc.reset_index(drop=True, inplace=True)            
             newData=df_temp.iloc[imax:imin,:]           
             newData=newData.drop(['Wavenumber'], axis=1)
             # could use iloc[imax:imin,1:] instead of dropping Wavenumber
              
             newData.reset_index(drop=True, inplace=True)
             #print("new Data.head",newData.head)
         
             #replace the old dataframe with new truncated dataframe
             classDict[ky]=newData
             
         #plot new truncated dataframe - plots mean of the SNV normalised data for each class
         counter=0  
         for ky in classDict:
             counter=counter+1
             df_temp_new = classDict[ky]
             ncols=len(df_temp_new)
             x= XvaluesTrunc                        
             y = df_temp_new.iloc[:, 1:ncols].mean(axis=1)# calc mean for plotting  
             dy=df_temp_new.iloc[:, 1:ncols].std(axis=1)
             self.dataplot(x,y,dy,counter)    
         #self.pushButton_8.setEnabled(False)   # truncate only allowed once 
        # could add reset button and allow re-truncate?
        
    def valueChange(self):  
        # values for truncate
         global xlowerLimit,xupperLimit
         s=self.spinBox.value()       
         sb1_val=s    
         xlowerLimit=sb1_val
         sb2_val=self.spinBox_2.value()           
         xupperLimit=sb2_val
         #print("truncate region limits",xlowerLimit,xupperLimit)
         
    def num_comp_PCA(self):   
             global numPCA_components            
             numPCA_components=self.spinBox_8.value()      
    def num_comp_LDA(self):   
             global numLDA_components            
             numLDA_components=self.spinBox_3.value()               
      
#2345        
    def LDA_fit(self):
        print("In LDA")
        global XvaluesTrunc,alldataArray,lab,unique                
        numLDA_components=self.spinBox_3.value()        
        
        alldataArray=[]
        count=0;counter=0
        
        for ky in classDict:
            count=count+1
            df_temp_new = classDict[ky]
            col_name=df_temp_new.columns
           
            ncols=len(col_name)
            
            lbl=count*np.ones((ncols,1),dtype=int) #add labels to array to allow classification        
            temp_array=df_temp_new.to_numpy() # convert from df to np array            
            temp_array=temp_array.T   #transpose array so spectra are in rows             
            #add labels (nb (shape= (num spectra, no spectral point + 1 point for lable)))
            temp_array=np.hstack((lbl, temp_array))            
            #print("Label / class",ky," shape",temp_array.shape)            
            # add all data in one large array for PCA & LDA
            if count==1:
                alldataArray=temp_array
            else:
                alldataArray=np.concatenate((alldataArray, temp_array))
                
                lab = alldataArray[:,0].astype('uint8') # get labels
         # find unique labels (classes of data)
        unique=list(set(lab))        
        y=lab 
        
        
        # perform LDA with top n_components        
        X=alldataArray[:,1:]
        sk_lda=LinearDiscriminantAnalysis(n_components= numLDA_components)   
        X_r2=sk_lda.fit_transform(X,y)        
         
        scores = cross_val_score(sk_lda, X, y, cv=4)
        print("Accuracy LDA only but using cross validation: %0.4f (+/- %0.4f)" % (scores.mean(), scores.std() * 2))
        print("")
        
        Message="LDA Accuracy = "+str(round(scores.mean(),3))+"(+/- "+str(round(scores.std() * 2,3))+")\n"
        self.textEdit.insertPlainText(Message)
        
        lda_varExp=sk_lda.explained_variance_ratio_
        print("LDA variance explained=",lda_varExp,"\n")
       
        
        # plot variance explained (from LDA)        
        xpos=range(1,numLDA_components+1,1)
        ypos=lda_varExp
        
        counter=1# set colour 
        self.dataplot2(xpos,lda_varExp,counter) # plot
        ypos=np.cumsum(lda_varExp)  
        counter=2# set colour 
        self.dataplot2(xpos,ypos,counter) # plot
        
        xi=[];yi=[];zi=[]
        y_lbl="LD"+str(2)
        x_lbl="LD"+str(1)
        z_lbl="LD"+str(3)
        
        
        lda_z=(self.spinBox_7.value())-1  
        lda_y=(self.spinBox_6.value())-1  
        lda_x=(self.spinBox_5.value() )-1 
        print("lda x,y",lda_x,lda_y)
        
        
        
        
        
        for i,u in enumerate(unique):        
            xi = [X_r2[j,lda_x] for j in range(len(X_r2[:,0])) if lab[j] == u]
            yi = [X_r2[j,lda_y] for j in range(len(X_r2[:,1])) if lab[j] == u]
            zi = [X_r2[j,lda_z] for j in range(len(X_r2[:,2])) if lab[j] == u]   
            
            fig4_plot=self.comboBox.currentIndex()  
            if fig4_plot==0:
                pass
            elif fig4_plot==1:                
                freq1,x1 = np.histogram(xi)   
                x11=x1[:-1]                            
                self.histplot(x11,freq1,i)
            elif fig4_plot==2:
                freq1,x1 = np.histogram(yi)   
                x11=x1[:-1]  
                self.histplot(x11,freq1,i) 
            else:   
                print("whatever")
           
            x_lbl="PCA_LDA"+str(lda_x+1)
            y_lbl="PCA_LDA"+str(lda_y+1)
            z_lbl="PCA_LDA"+str(lda_z+1)
                
            self.scattplot(xi,yi,i,x_lbl,y_lbl,u)
            self.scattplot2(xi,yi,zi,i,x_lbl,y_lbl,z_lbl,u)
            
            xi=[];yi=[];zi=[]
        
        
        
        
     
#2345        
    def PCA_LDA_fit(self):        
        global XvaluesTrunc,alldataArray        
        numPCA_components=self.spinBox_8.value()  
        numLDA_components=self.spinBox_3.value()
        print("In PCA_LDA")
        print("number of PCA components=",numPCA_components)
        print("number of LDA components=",numLDA_components)
        
        alldataArray=[]
        count=0;counter=0
        
        for ky in classDict:
            count=count+1
            df_temp_new = classDict[ky]
            col_name=df_temp_new.columns           
            ncols=len(col_name)            
            lbl=count*np.ones((ncols,1),dtype=int) #add labels to array to allow classification
            temp_array=df_temp_new.to_numpy() # convert from df to np array            
            temp_array=temp_array.T   #transpose array so spectra are in rows                        
            temp_array=np.hstack((lbl, temp_array))            
            if count==1:
                alldataArray=temp_array
            else:
                alldataArray=np.concatenate((alldataArray, temp_array))                
                lab = alldataArray[:,0].astype('uint8') # get labels
        unique=list(set(lab))  # find unique labels (classes of data)      
        y=lab 
      
        #combined PCA / LDA
        X=alldataArray[:,1:]
        sk_lda=LinearDiscriminantAnalysis(n_components= numLDA_components)   
        #X_r2=sk_lda.fit_transform(X,y)        
        #print("shape of X_r2=", X_r2.shape)
        
        #LDA with cross validation on plain X data
        #scores = cross_val_score(sk_lda, X, y, cv=4)
        #print("Accuracy LDA only but using cross validation: %0.4f (+/- %0.4f)" % (scores.mean(), scores.std() * 2))
        #print("")
        
        
        #combined PCA / LDA
        sk_pca = PCA(n_components=numPCA_components) # or fixed value eg 15 components
        Xpc = sk_pca.fit_transform(X)
        scores = cross_val_score(sk_lda, Xpc, y, cv=4)
        
        print("Accuracy PCA-LDA (cv=4): %0.4f (+/- %0.4f)" % (scores.mean(), scores.std() * 2))
        
        Message="PCA_LDA Accuracy (cv=4) = "+str(round(scores.mean(),3))+"(+/- "+str(round(scores.std() * 2,3))+")\n"
        self.textEdit.insertPlainText(Message)
        
        lda_z=(self.spinBox_7.value())-1  
        lda_y=(self.spinBox_6.value())-1  
        lda_x=(self.spinBox_5.value() )-1 
        print("lda x,y",lda_x,lda_y)
        
        X=sk_lda.fit_transform(Xpc,y)    
        
        lda_pca_varExp=sk_lda.explained_variance_ratio_
        print("LDA-PCA varience explained=",lda_pca_varExp)
        print("")
        
        xpos=range(1,numLDA_components+1,1)
        ypos=lda_pca_varExp 
        counter=1# set colour 
        self.dataplot2(xpos,ypos,counter)
        ypos=np.cumsum(lda_pca_varExp)  
        counter=2# set colour 
        self.dataplot2(xpos,ypos,counter)
        
        
        
        xi=[];yi=[];zi=[]
        for i,u in enumerate(unique):
           
            xi = [X[j,lda_x] for j in range(len(X[:,0])) if lab[j] == u]
            yi = [X[j,lda_y] for j in range(len(X[:,1])) if lab[j] == u]
            zi = [X[j,lda_z] for j in range(len(X[:,2])) if lab[j] == u] 
            
            
            fig4_plot=self.comboBox.currentIndex()  
            if fig4_plot==0:
                pass
            elif fig4_plot==1:                
                freq1,x1 = np.histogram(xi)  
                x11=x1[:-1]                              
                self.histplot(x11,freq1,i)
            elif fig4_plot==2:
                freq1,x1 = np.histogram(yi) 
                x11=x1[:-1] 
                self.histplot(x11,freq1,i) 
            else:   
                print("whatever")
           
            
            x_lbl="PCA_LDA"+str(lda_x+1)
            y_lbl="PCA_LDA"+str(lda_y+1)
            z_lbl="PCA_LDA"+str(lda_z+1)
            
            # plot scatter plots 2D and 3D
            self.scattplot(xi,yi,i,x_lbl,y_lbl,u)
            self.scattplot2(xi,yi,zi,i,x_lbl,y_lbl,z_lbl,u)
            
                                 
            xi=[];yi=[];zi=[]
       
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
            
            
               
        
    def getfiles(self):        
        global x,y,dataArray,fileCount,MessString,df2,classDict,dfSNV,xFullSpectrum
        
        fileCount=fileCount+1  # number of different classes / sample types
        dfname="df"+str(fileCount) #create sequential names for dataframes as read in  eg df1,df2,...
        dfnameSNV="df"+str(fileCount) +"SNV"   # dataframe names for Standard Normal Variate eg, df1SNV,...
        print("file number",fileCount," dataframe name",dfname)  
        
        # open file (show .csv files)
        dlg = QFileDialog()   
        dlg.setNameFilter('(*.csv)')
        dlg.setFileMode(QFileDialog.AnyFile)
        fileNameList=[]
        
      
        minX_value=self.spinBox.value()   # low wavenumber limit
        maxX_value=self.spinBox_2.value() # high wavenumber limit
         
        
        if dlg.exec_():
           
            filenames = dlg.selectedFiles()
            fnam=filenames[0] # could put all data analysed files in a folder and open in one go            
            data= pd.read_csv (filenames[0])                   
            df = pd.DataFrame(data) #create pandas dataframe called df (reading analysed.csv files - as produced by Raman processing program)
            
    
            # get file & path name  (might only work for Mac i.e need to change separator?) i.e. if windows use "\"
            sep="/"
            fnamLen=len(fnam)
            lastSep=fnam.rfind(sep)
            #pathName=fnam[:lastSep+1]  # identify path
            fileName=fnam[lastSep+1:fnamLen]   # identify file name
            fileNameList.append((fileName))
           
            #print(df.columns)
            # set spin box values to values that are appropriate for the data
            x=df["Wavenumber"]
            if x.max()<maxX_value:
                maxX_value = int(x.max())
                self.spinBox_2.setValue(maxX_value) 
            if x.min()>minX_value:
                minX_value=int(x.min())
                self.spinBox.setValue(minX_value) 

            # plot x,y summary data just to show representation of loaded data (nb could have also plotted smoothed data)
            xFullSpectrum=df["Wavenumber"] # in case you need to know original spectral length
            y=df["mean"]  # get the mean value column from the data file (could also get the smoothed mean)
            y_sd=df["SD1"]
            
            #y=df["mean_smoothed"]  # get the mean value column from the data file (could also get the smoothed mean)
            #y_sd=df["SD_smoothed"]
            self.dataplot(x,y,y_sd,fileCount) # plot average spectrum from data analysed
            
            #put file names in message window - to show which files have been loaded
            pcol=color_by_name[fileCount-1]    #colours for plotting 
            MessString=fileName+" ("+ pcol+")"+"\n"          
            self.textEdit.insertPlainText(MessString)
            
            # drop unused columns
            df.drop('mean', axis=1, inplace=True)
            df.drop('SD1', axis=1, inplace=True)
            #df.drop('mean_smoothed', axis=1, inplace=True)
            #df.drop('SD_smoothed', axis=1, inplace=True)            
            #print(df.columns)
            
            
            
            # create SNV 
            col_name=df.columns
            col_nameSNV=["Wavenumber"]
            
            dataSNV=[]  #new dataframe
            #calc standard normal variance
            
            # create numpy array of data and column name header then create new SNV dataframe
            # this is a longwinded way and could be done using StandardScalar
            # feat_SNV = StandardScaler().fit_transform(feat)
          
            
            for col in col_name:
                if col=="Wavenumber":
                    x=df[col]
                elif col=="mean":
                    pass
                elif col=="SD1":
                    pass
                elif col=="mean_smoothed":
                    pass
                elif col=="SD_smoothed":
                    pass
                else:
                    newcolName=col+"_SNV"
                    
                    col_nameSNV.append(newcolName)
                    
                    #SNV is (data-average)/stdev  (column wise)
                    dataSNV=(df[col]-df[col].mean())/df[col].std()
                    if col=="Data1":
                        dataArray_SNV=np.column_stack((x,dataSNV))
                    else:
                        dataArray_SNV=np.column_stack((dataArray_SNV,dataSNV))
                    
            #create new data frame - without average and with SNV normalised data        
            dfSNV = pd.DataFrame(dataArray_SNV)
            dfSNV.columns=col_nameSNV
           
               
            #calc mean of SNV data as a check
            #ncols=len(dfSNV.columns)
            #meanSNV =dfSNV.iloc[:, 1:ncols].mean(axis=1) 
            #stdSNV= dfSNV.iloc[:, 1:ncols].std(axis=1)
            
            
           
            #add name data frame to dictionary of dataframes,  name=dfnameSNV           
            classDict[dfnameSNV]=dfSNV.copy() 
            #print("Dictionary=\n",classDict)
            #set spin box to be highest number allowed
           
            self.spinBox_3.setMinimum(1)
            if fileCount>1:
                self.spinBox_3.setValue(fileCount-1)
                self.spinBox_3.setMaximum(fileCount-1)
            else:
                self.spinBox_3.setValue(1)
                self.spinBox_3.setMaximum(1)
                
                
                
                
    def histplot(self,x,freq,colcount):
                  pcol=color_by_name[colcount]   
                  
                  if colcount==0:                     
                      self.canvas.axes4.cla()
                      self.canvas.axes4.set_xlabel("LD Score",fontsize=12)
                      self.canvas.axes4.set_ylabel("Frequency",fontsize=12)
                      #self.canvas.axes4.set_title("LD Score",fontsize=14)                     
                      self.canvas.axes4.bar(x,freq,width=0.5,bottom=None, align='edge',color=pcol,alpha=0.5, label=str(colcount+1))                    
                  else:  
                      self.canvas.axes4.set_xlabel("LD Score",fontsize=12)
                      self.canvas.axes4.bar(x,freq,width=0.5,bottom=None, align='edge',color=pcol,alpha=0.5, label=str(colcount+1))
                  
                  self.canvas.axes4.legend()       
                  cb2=self.checkBox_2.isChecked()
                  if cb2==True:
                      self.canvas.axes4.grid(color = 'gray', linestyle = '--', linewidth = 0.5)
                  self.canvas.axes4.legend()
                  self.canvas.draw ()
            
   
    def scattplot2(self,xi,yi,zi,colcount,x_lbl,y_lbl,z_lbl,u):
          # PCA Scores
         
         pcol=color_by_name[colcount]   
         if colcount==0:
             self.canvas.axes3.cla()
             self.canvas.axes3.set_xlabel(x_lbl,fontsize=12)
             self.canvas.axes3.set_ylabel(y_lbl,fontsize=12)
             self.canvas.axes3.set_zlabel(z_lbl,fontsize=12)
             #self.canvas.axes3.set_title("LDA Scores",fontsize=14) 
             
         self.canvas.axes3.scatter( xi,yi,zi,color=pcol, s=100,edgecolors='k',  alpha=0.4, label=str(u)) 
         self.canvas.axes3.legend()
         cb2=self.checkBox_2.isChecked()
         if cb2==True:
             self.canvas.axes3.grid(color = 'gray', linestyle = '--', linewidth = 0.5)
         self.canvas.draw ()  
    
              
            
    def scattplot(self,xi,yi,colcount,x_lbl,y_lbl,u):
         # LDA Scores
               
        pcol=color_by_name[colcount]   
        if colcount==0:
            self.canvas.axes2.cla()
            self.canvas.axes2.set_xlabel(x_lbl,fontsize=12)
            self.canvas.axes2.set_ylabel(y_lbl,fontsize=12)
            #self.canvas.axes2.set_title("LDA Scores",fontsize=14) 
            
        self.canvas.axes2.scatter( xi,yi,color=pcol, s=100,edgecolors='k',  alpha=0.4, label=str(u)) 
        self.canvas.axes2.legend()
        cb2=self.checkBox_2.isChecked()
        if cb2==True:
            self.canvas.axes2.grid(color = 'gray', linestyle = '--', linewidth = 0.5)
        self.canvas.draw ()
       
    
    def dataplot2(self,x,y,colcount):
        # plot varience and variance explained       
        self.canvas.axes4.set_xlabel("LD Component",fontsize=12)
        self.canvas.axes4.set_ylabel("Variance",fontsize=12)
        self.canvas.axes4.set_title("Variance Explained/ Cummulative Variance",fontsize=14) 
        self.canvas.axes4.legend()
        
        if colcount==1:
            self.canvas.axes4.cla()
            self.canvas.axes4.plot( x,y,marker='o',markerfacecolor="rosybrown", markersize=8,color='k',  linestyle='-',alpha=1,label='Var. Expl.') 
        else: 
            self.canvas.axes4.plot( x,y,marker='o',markerfacecolor="cornflowerblue", markersize=8,color='k',  linestyle='-',alpha=1,label='Cum. Var.') 
        
               
        cb2=self.checkBox_2.isChecked()
        if cb2==True:
            self.canvas.axes4.grid(color = 'gray', linestyle = '--', linewidth = 0.5)
        self.canvas.axes4.legend()
        self.canvas.draw ()
        
       
    
    def dataplot(self,x,y,dy,fileCount):
        #Average Data 
        
        if fileCount==1:
            self.canvas.axes.cla()
            self.canvas.axes.set_xlabel("Wavenumber (cm-1)",fontsize=12)
            self.canvas.axes.set_ylabel("Intensity (arb. units)",fontsize=12)
            self.canvas.axes.set_title("Raman Spectra",fontsize=14) 
        offset_ratio=1.5   
        y=y+(fileCount*offset_ratio)    # offset spectra
        pcol=color_by_name[fileCount-1]        
        self.canvas.axes.plot( x,y,color=pcol,linewidth=1.5,linestyle='-',alpha=1.0)  
        self.canvas.axes.fill_between(x,y+dy,y-dy,facecolor=pcol,linewidth=1.5,alpha=0.4)  
        
        cb2=self.checkBox_2.isChecked()
        if cb2==True:
            self.canvas.axes.grid(color = 'gray', linestyle = '--', linewidth = 0.5)
        #self.canvas.axes.legend()    
        self.canvas.draw ()
        
      
        
    def bgndChange(self):
        global styles, ax_col,cols
       
        
           
        
        
if __name__=="__main__":
    app=QtWidgets.QApplication(sys.argv)
    w=MainWindow()
    w.show()
    app.exec_()

