#!/usr/bin/env python
# encoding: utf-8
"""
processPatches.py
Loop through all patch files and evaluate significance of each one

Created by Mark Evans on 2012-01-11.
Copyright (c) 2012 __Monogram Biosciences__. All rights reserved.


07.16.2012 Revised to include R code from MinePatch-doAll-allSigSites.R 06.21.12
            Also modifed to output .csv
06.07.2012 Revised to include R code from MinePatch-doAll-SimpleOut.R_v4 06.06.12
02.28.2012 Revised to include R code from MinePatch-doAll-SimpleOut.R_v2 02.27.12
"""

import sys,os,glob,getopt,subprocess, random
from multiprocessing import Process, Manager, Queue
import multiprocessing
from operator import itemgetter

#####################################################################################
# glbl                                                                              #
# Creates a class containing global variables that can be accessed across processes # 
#####################################################################################
class glbl:
   jobs =[]
   data_dirs = {}

def processData(nid,results):
    path =''
    while 1:
        try: path = nid.get(False)                                                             # Get next query id from Queue generator, 
        except: 
            print " CPU not in use"
            break
        print path
        rcode = """
        ##################################################################
        # Neutralization - Patch Analysis
        # Multiple Linear Regression, Logistic Regression, SVM, and Fisher's Exact Test
        # Mojgan Haddad - 6-21-2012
        # listing all significant residues by MLR in patches with # of models
        ##################################################################


        ############################################
        ## Function for SVM

        doSVM<-function(myData, patch, maxIC, patchFile) {

         library(e1071) ## lib for SVM

         ## split into training and test set (even and odd rows)
         myData$train = as.numeric(rownames(myData)) %%%% 2  
         trainset = patch #patch[myData$train == 1,] 
         testset = patch #patch[myData$train == 0,]

         ##classification: sensitive (-1) or resistant (1)
         totrainClass <- ifelse(myData[,1] < maxIC, -1, 1) 
         totestClass <- ifelse(myData[,1] < maxIC, -1, 1)  
         #totrainClass <- ifelse(myData[myData$train == 1,1] < maxIC, -1, 1)  
         #totestClass <- ifelse(myData[myData$train == 0,1] < maxIC, -1, 1)  

         ## initialize performance characteristics
         SVM.Spec = 0
         SVM.Sens = 0

         # try to train SVM only if 2 classes (R/S) are available
         if (length(unique(totrainClass)) > 1) {

          ## SVM on predicting sensitive vs resistant
          mymodel <- svm(totrainClass  ~ ., data=trainset, 
                      type="C-classification", ## this is doing a classification
                      kernel="linear", ## options tried: "polynomial",radial basis, linear, sigmoid
                      cross=20, ## 20-fold cross validation was the best option
                      scale=FALSE)
                      ##cost=0.35) --> keep default

          ## Perform Prediction in SVMmodel is not empty
          if (length(summary(mymodel)) > 3) {

           mypred <- predict(mymodel, testset, decision.values=T)
           colnames(attr(mypred,"decision.values"))

           ## total accuracy from cross-val
           xValAccuracy = round(summary(mymodel)$tot.accuracy, 2)

           ##"1/-1" we are interested in the decision of classification of 1/-1
           tt <- table(Prediction=mypred, Sens=totestClass) 

           ## Check if the table has at least 2 rows 
           if (nrow(tt)>1) {
            totalS = tt[1,1]+tt[2,1]
            totalR = tt[2,2]+tt[1,2]

            ## make sure no DIV by zero; also
            ## if all samples are Sens (or R) 
            # then-> no SVM prediction needed
            if (totalR > 0 & totalS > 0) {

              SVM.Spec = round(tt[1,1]*100/ totalS, 2)
              SVM.Sens = round(tt[2,2]*100/ totalR, 2)
            }
           }
          } ## if SVM model not empty
         } ## end of if 2 classes

         ## If sensitivity or specificity > 70%% 
         ## -> and specificity/sensitivity > 30%%
         ## or total accuracy from cross-val > 70%%
         isSig = c(0,"-")
         if ((SVM.Sens > 70 & SVM.Spec > 30)  | (SVM.Sens > 30 & SVM.Spec > 70) | (xValAccuracy > 70)) {

             isSig = c(xValAccuracy,
                       max(SVM.Sens , SVM.Spec) )
          }

          ## return vector of accuracy from cross-val, and max of sens & spec
          return(isSig)

        } ## END of SVM 

        ############################################
        ## Function for MLR
        ## linear regression on IC50 versus patch composition

        doMLR<-function(myData, patch) {

          mylm = lm(myData[,1] ~ ., data=patch)  
          myModel <- summary(mylm)$coefficients

          ## get the p-values/significance for the patch parameters
          pvali = grep("Pr",colnames(myModel))
          pvalues = myModel[,pvali]

          ## get the coefficient
          coeffi = grep("Estimate",colnames(myModel))

          patchi <- names(pvalues)!="(Intercept)" & names(pvalues)!="ic50_val" & pvalues<0.05
          sigPVals <- round(pvalues[patchi],3)
          MLRCoeff <- round(myModel[patchi, coeffi],2)
          MLRResidues <- rownames(myModel)[patchi]

          isSig <- ifelse( length(sigPVals) > 0, 1, 0)

          if (isSig == 1) {  

            ## get the most significant coefficient & residue 
            # return as a vector
            sigVals = c(length(sigPVals),
                        MLRResidues[sigPVals==min(sigPVals)],
                        MLRCoeff[sigPVals==min(sigPVals)],
                        sigPVals[sigPVals==min(sigPVals)] )

          } else {
            sigVals = c(0,"-","-","-")
          }

          return(sigVals)

        } ## END of MLR


        ############################################
        ## Function for Logistic Regression
        ## logistic regression on IC50 R/S versus patch composition

        doLogR<-function(patch, myCall) {

          myLogr = glm(myCall ~ ., data=patch, family=binomial(link="logit") ) 
          myModel <- summary(myLogr)$coefficients

          ## get the p-values/significance for the patch parameters
          pvali = grep("Pr",colnames(myModel))
          pvalues = myModel[,pvali]

          ## get the coefficient
          coeffi = grep("Estimate",colnames(myModel))

          patchi <- names(pvalues)!="(Intercept)" & names(pvalues)!="ic50_val" & pvalues<0.05
          sigPVals <- round(pvalues[patchi], 3)
          LOG.REGCoeff <- round(myModel[patchi, coeffi], 2)
          LOG.REGResidues <- rownames(myModel)[patchi]

          ## after correction for multiple testing (~100 patches=> 0.05/100)
          ## Show trending results too: p-value<0.1
          isSig <- ifelse( length(sigPVals) > 0, 1, 0)

          if (isSig == 1) {  

            ## get the most significant coefficient & residue 
            # return as a vector
            sigVals = c( LOG.REGCoeff[sigPVals==min(sigPVals)],
                      LOG.REGResidues[sigPVals==min(sigPVals)])

          } else {
            sigVals = c(0,"-")
          }

          return(sigVals)

        } ## END of LogReg

        ############################################
        ## Functions for FET

        findUni<-function(posINpatch, myCall){
          mut_resistant = 0;
          mut_sensitive = 0;
          noMut_resistant = 0;
          noMut_sensitive = 0;

          ## generate matrix with:
          ### a) resistant-mut; b) resistant-noMut; c) sensitive-mut; d) sensitive-noMut

          mutated<-ifelse(posINpatch == 1, "Mut", "No")

          neutresp<-ifelse(myCall == 0, "NR", "R")

          tbl<-table(mut=mutated, resp=neutresp)

          ## only examine the table if # mutated > 0 and #NR > 0
          if (dim(tbl)[1] > 1) {
        	  mut_resistant<- tbl[1,2]
        	  mut_sensitive<- tbl[1,1]
        	  noMut_resistant<- tbl[2,2]
        	  noMut_sensitive<- tbl[2,1]
          }

          ## default to 1 if 0 to be able to calculate odds ratio (avoid DIV by 0)
          if (mut_resistant == 0) { mut_resistant = 1; }
          if (mut_sensitive == 0) { mut_sensitive = 1; }
          if (noMut_resistant == 0) { noMut_resistant = 1; }
          if (noMut_sensitive == 0) { noMut_sensitive = 1; }

         fisher_data<-matrix(c(mut_resistant,mut_sensitive,noMut_resistant,noMut_sensitive),nr=2,byrow=T,dimnames = list(c("Mut_Present", "Mut_Abscent"),c("Resistant","Sensitive")));
         #print(fisher_data);


         ftest <- c(fisher.test(fisher_data));
         ftest$mut_x = mut_resistant
         ftest$mut_r = mut_sensitive
         ftest$noMut_x = noMut_resistant
         ftest$noMut_r = noMut_sensitive
         #print(ftest);

         return(ftest);
        }


        doFET<-function(patch, myCall){
          names <- names(patch);
          pvalues<-c();
          use_names<-c();
          all_data<-c();
          ma<-c();
          mb<-c();
          aa<-c();
          ab<-c();

          for(i in 1:length(patch)) {
              # myCall is resistant/sensitive (-1/1)
              test_res <- c(findUni(patch[i], myCall));

              pvalues<-c(pvalues,test_res$p.value);
              #print(pvalues);

              use_names<-c(use_names,names[i]);
              ma <- c(ma,test_res$mut_x)
              mb <- c(mb,test_res$mut_r)
              aa <- c(aa,test_res$noMut_x)
              ab <- c(ab,test_res$noMut_r)
          }

          o<-c(order(pvalues,decreasing = TRUE));

          for(i in 1:length(o)) {
           if(mb[o[i]]*aa[o[i]] == 0){
        	OddsRatio = "NA";
           } else {
              OddsRatio = round( ((ma[o[i]]*ab[o[i]])/(mb[o[i]]*aa[o[i]])), 2)
           }

           isSig = c(0, "-")
           if (pvalues[o[i]] < 0.05) { 

              Residue = use_names[o[i]]

              isSig = c( OddsRatio[pvalues[o[i]] == min(pvalues[o[i]]) ], Residue)  
           }

          }

          return(isSig)

        } ############ End of FET #############

        ############################################
        ## Function for MLR
        ## linear regression on IC50 versus patch composition

        doMLRwOutput<-function(myData, patch, nrModels, patchFile) {

          mylm = lm(myData[,1] ~ ., data=patch)  
          myModel <- summary(mylm)$coefficients

          ## get the p-values/significance for the patch parameters
          pvali = grep("Pr",colnames(myModel))
          pvalues = myModel[,pvali]

          ## get the coefficient
          coeffi = grep("Estimate",colnames(myModel))

          ## after correction for multiple testing (~100 patches=> 0.05/100)
          #isSig <- ifelse( max(ifelse( sigPVals<0.005, 1, 0), na.rm=FALSE) == 1, 1, 0)

          ## Show significant results too: p-value<0.05
          patchi <- names(pvalues)!="(Intercept)" & names(pvalues)!="ic50_val" & pvalues<0.05
          sigPVals <- round(pvalues[patchi],3)
          MLRCoeff <- round(myModel[patchi, coeffi],2)
          MLRResidues <- rownames(myModel)[patchi]

          isSig <- ifelse( length(sigPVals) > 0, 1, 0)

          if (isSig == 1) {  

            patchSites = names(patch)
            ### make one string from patch sites -> then pass to print
            patchStr = paste(patchSites, collapse="-")

            nrRows = length(sigPVals)
            for (rw in 1:nrRows) {

              myTbl = c(patchFile, patchStr, nrModels, MLRResidues[rw], MLRCoeff[rw], sigPVals[rw], nrRows)
              toPrint = paste(myTbl, collapse=", ")
              print(toPrint)
            }

          }

        } ## END of MLRwithOutput
 
    # main
    #########################
    setwd("%s")
    filelist <- dir()
    
    ## Print column title
    print(paste (c("PatchFile", "Patch.30CHARs", "Nr-of-Models", "MLR-SigSite", "MLR", "MLR.p", "Nr-of-SigSites"), collapse=", ") )

    for(filei in filelist) # loop thru all files
    {
      mypatch <- read.delim(filei,header=T)

      ici = grep("ic50_val",colnames(mypatch))
      starti = ici+1
      lasti = length(mypatch)

      myData = mypatch[,ici:lasti]

      #handle if there is only one site in the patch
      if (starti == lasti) { 
      	patch = myData
      } else {
      	patch = mypatch[,starti:lasti]
      }

      ##classification: sensitive (0) or resistant (1)
      maxIC = 25
      myCall <- ifelse(myData[,1] < maxIC, 0, 1) 

      ## Perform linear regression on patch versus IC50
      myMLR = doMLR(myData, patch)

      ## Perform logistic regression on patch versus neut response
      LogR = doLogR(patch, myCall) 

      ## Perform SVM versus neut response
      mySVM = doSVM(myData, patch, maxIC, filei)

      ## Perform FET on mutations in the patch versus neut response
      myFET = doFET(patch, myCall) 

      nrModels = 0
      if (myMLR[1]!=0) { nrModels = nrModels + 1 }
      if (LogR[1]!=0) { nrModels = nrModels + 1 }
      if (mySVM[1]!=0) { nrModels = nrModels + 1 }
      if (myFET[1]!=0) { nrModels = nrModels + 1 }


      ## if any model significant => print patch
      if (nrModels > 0) {  

        ### redo MLR-> print all significant sites from MLR
        doMLRwOutput(myData, patch, nrModels, filei)

      }
    }
    ###-- end of main
 
""" % (path)
        tmp = str(random.random())[2:20]  # generate random temp filename
        tmp2 = str(random.random())[2:20]
        f = open(tmp,'w')
        f.write(rcode)
        f.close()
        f2 = open(tmp2,'w')
        try:
            process = subprocess.Popen(["R","--vanilla","--subordinate","-f",tmp], stdout=f2)
            process.wait()
        except:
            print "R subprocess failed"
            os.remove(tmp)
            f2.close()
            os.remove(tmp2)
      #      multiprocessing.current_process.terminate()
    
        os.remove(tmp)
        f2.close()
        if os.path.getsize(tmp2) != 0:
          f = open(tmp2,'r')
          results[path] = f.read()
          f.close()
          os.remove(tmp2)
        else:
          results[path] = 'None'
          os.remove(tmp2)
    return




def processFiles(patch_dir):
    root = os.getcwd()
    glbl.data_dirs = {}
    if root != patch_dir: working_path = root+"/"+patch_dir
    else: working_path = root

    for path, dirs, files in os.walk(working_path):
        if len(dirs) == 0: glbl.data_dirs[path] = ''
    

    # Multiprocessing Section
    #########################################
    Qids = glbl.data_dirs.keys()
    manager = Manager()                                      # creates shared memory manager object
    results = manager.dict()                                 # Add dictionary to manager, so it can be accessed across processes
    nextid = Queue()                                         # Create Queue object to serve as shared id generator across processes
    for qid in Qids: nextid.put(qid)                         # Load the ids to be tested into the Queue
    for x in range(0,multiprocessing.cpu_count()):           # Create one process per logical CPU
        p = Process(target=processData, args=(nextid,results)) # Assign process to processCBR function, passing in the Queue and shared dictionary
        glbl.jobs.append(p)                                   # Add the process to a list of running processes
        p.start()                                             # Start process running
    for j in glbl.jobs:
        j.join()                                              # For each process, join them back to main, blocking on each one until finished
    
    # write out results
    c = 1
    sets = results.keys()
    sets.sort()
    for x in sets:
        if results[x] != 'None':
            FINAL = open('result'+str(c)+'.csv','w')
            n = "\n************************************************************************************************\n"
            FINAL.write(n+"* "+x+'    *\n'+n+results[x]+"\n")
            FINAL.close()     
            c += 1
            

def main(argv):
    patch_dir=''
    #Read command line arguments
    try:
        opts, args = getopt.getopt(argv,"d:",["patch_dir"])
    except getopt.GetoptError:
        print "There was an error, no args\n"
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-d": 
            patch_dir = arg
            processFiles(patch_dir)
        else:
            print "There was an error, not right args\n"
            sys.exit()
       
   

if __name__ == '__main__':
    main(sys.argv[1:])