#!/usr/bin/env python
# encoding: utf-8
"""
processPatches.py
Loop through all patch files and evaluate significance of each one

Created by Mark Evans on 2012-01-11.
Copyright (c) 2012 __Monogram Biosciences__. All rights reserved.

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
############################################
## Function for SVM

doSVM<-function(myData, patch, maxIC, patchSites) 
{
  library(e1071) ## lib for SVM
  
  ## split into training and test set (even and odd rows)
  myData$train = as.numeric(rownames(myData)) %%%% 2  
  trainset = patch[myData$train == 1,] 
  testset = patch[myData$train == 0,]
  
  totrainClass <- ifelse(myData[myData$train == 1,1] < maxIC, 1, -1)  
  totestClass <- ifelse(myData[myData$train == 0,1] < maxIC, 1, -1)  
  
  
  ## SVM on predicting sensitive vs resistant
  mymodel <- svm(totrainClass  ~ ., data=trainset, 
              type="C-classification", ## this is doing a classification
              kernel="linear",
              cross=10, ## for 10-fold cross validation
              scale=T,
              cost=0.35)
  
  ## Prediction and Best model confusion matrix
  mypred <- predict(mymodel, testset, decision.values=T)
  colnames(attr(mypred,"decision.values"))
  
  ##"1/-1" we are interested in the decision of classification of 1/-1
  tt <- table(Prediction=mypred, Sens=totestClass) 
  
  SVM.Spec = 0
  SVM.Sens = 0
  if (nrow(tt)>1) 
  {
    SVM.Spec = round(tt[2,2]*100/(tt[2,1]+tt[2,2]), 2)
    SVM.Sens = round(tt[1,1]*100/(tt[1,1]+tt[1,2]), 2)
  }
  
  ## If sensitivity or specificity > 70%% 
  ## not 100%% -> all 0 muts would do that
  isSig = 0
  if (SVM.Sens > 70 | SVM.Spec > 70) 
  {   
     ### make one string from patch sites
     SVM.Site = paste(patchSites, collapse="-")
     
     SS = tt[2,2] 
     SR = tt[2,1] 
     RS = tt[1,2] 
     RR = tt[1,1] 
           
     myTbl = data.frame(SVM.Site, SVM.Sens, SVM.Spec, SS, SR, RS, RR)
     print(myTbl)
     
     isSig = 1
  }
  return(isSig)
} ## END of SVM 

############################################
## Function for MLR
## linear regression on IC50 versus patch composition

doMLR<-function(myData, patch) 
{
  mylm = lm(myData[,1] ~ ., data=patch)  
  myModel <- summary(mylm)$coefficients
  
  ## get the p-values/significance for the patch parameters
  pvali = grep("Pr",colnames(myModel))
  pvalues = myModel[,pvali]
  
  ## get the coefficient
  coeffi = grep("Estimate",colnames(myModel))
  
  patchi <- names(pvalues)!="(Intercept)" & names(pvalues)!="ic50_val" & pvalues<0.1
  sigPVals <- round(pvalues[patchi],3)
  MLRCoeff <- round(myModel[patchi, coeffi],2)
  MLRResidues <- rownames(myModel)[patchi]
  
  ## after correction for multiple testing (~100 patches=> 0.05/100)
  ## Show trending results too: p-value<0.1
  isSig <- ifelse( max(ifelse( sigPVals<0.1, 1, 0)) == 1, 1, 0)
  
  if (isSig == 1) 
  {  
    myTbl = data.frame(MLRResidues, MLRCoeff, sigPVals)
    print(myTbl)
  }
  return(isSig)
} ## END of MLR


############################################
## Function for Logistic Regression
## logistic regression on IC50 R/S versus patch composition

doLogR<-function(patch, myCall) 
{
  myLogr = glm(myCall ~ ., data=patch, family=binomial(link="logit") ) 
  myModel <- summary(myLogr)$coefficients
  
  ## get the p-values/significance for the patch parameters
  pvali = grep("Pr",colnames(myModel))
  pvalues = myModel[,pvali]
  
  ## get the coefficient
  coeffi = grep("Estimate",colnames(myModel))
  
  patchi <- names(pvalues)!="(Intercept)" & names(pvalues)!="ic50_val" & pvalues<0.1
  sigPVals <- round(pvalues[patchi], 3)
  LOG.REGCoeff <- round(myModel[patchi, coeffi], 2)
  LOG.REGResidues <- rownames(myModel)[patchi]
  
  ## after correction for multiple testing (~100 patches=> 0.05/100)
  ## Show trending results too: p-value<0.1
  isSig <- ifelse( max(ifelse( sigPVals<0.1, 1, 0)) == 1, 1, 0)
  
  if (isSig == 1) 
  {    
    myTbl = data.frame(LOG.REGResidues, LOG.REGCoeff, sigPVals)
    print(myTbl)
    
    #print(">>>LOG-REG----> Significant residues - Coefficients - pvalues")
    #print(paste(sigResidues, " - Coeff= ", LogRegCoeff, " - p=", sigPVals))
    #print(myModel)
  }  
  return(isSig)
} ## END of LogReg

############################################
## Functions for FET

findUni<-function(posINpatch, myCall)
{
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
  if (dim(tbl)[1] > 1) 
  {
    mut_resistant<- tbl[1,1]
    mut_sensitive<- tbl[1,2]
    noMut_resistant<- tbl[2,1]
    noMut_sensitive<- tbl[2,2]
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


doFET<-function(patch, myCall)
{
  names <- names(patch);
  pvalues<-c();
  use_names<-c();
  all_data<-c();
  ma<-c();
  mb<-c();
  aa<-c();
  ab<-c();
    
  for(i in 1:length(patch)) 
  {
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

  for(i in 1:length(o)) 
  {
    if(mb[o[i]]*aa[o[i]] == 0)
    {
      OddsRatio = "NA";
    } 
    else 
    {
      OddsRatio = round( ((ma[o[i]]*ab[o[i]])/(mb[o[i]]*aa[o[i]])), 2)
   }
   
   isSig = 0
   if (pvalues[o[i]] < 0.05) 
   {    
      isSig = 1
      
      Residue = use_names[o[i]]
      pVal = round(pvalues[o[i]], 3)
      MutR= ma[o[i]] 
      MutS= mb[o[i]]
      AbsR= aa[o[i]]
      AbsS= ab[o[i]]
      
      myTbl = data.frame(Residue, OddsRatio, pVal, MutR, MutS, AbsR, AbsS)
      print(myTbl)
      
      ### Results:
      #print(paste(">>>FET: Residue= ", use_names[o[i]], 
      #         " - pvalue= ", round(pvalues[o[i]], 3),
      #         " - MutR= ", ma[o[i]], " MutS= ", mb[o[i]], "NoneR= ", aa[o[i]], "NoneS= ", ab[o[i]], 
      #         " - ODDsRATIO= ", round(oddsratio,2) ) ) 
   }
  }
  return(isSig)
} ############ End of FET #############

##########################
# main

setwd("%s")
filelist <- dir()

for(filei in filelist) # loop thru all files
{
  mypatch <- read.delim(filei,header=T)
  
  ici = grep("ic50_val",colnames(mypatch))
  starti = ici+1
  lasti = length(mypatch)
  
  myData = mypatch[,ici:lasti]
  
  #handle if there is only one site in the patch
  if (starti == lasti) 
  { 
    patch = myData
  } 
  else 
  {
    patch = mypatch[,starti:lasti]
  }
  
  patchSites = names(patch)
  
  ##classification: sensitive (1) or resistant (0)
  maxIC = 25
  myCall <- ifelse(myData[,1] < maxIC, 1, 0) 
  
  ## Perform linear regression on patch versus IC50
  mlr = doMLR(myData, patch)
  
  ## Perform logistic regression on patch versus neut response
  logr = doLogR(patch, myCall) 
  
  ## Perform SVM versus neut response
  mySVM = doSVM(myData, patch, maxIC, patchSites)
  
  ## Perform FET on mutations in the patch versus neut response
  myFET = doFET(patch, myCall) 
  
  ## if any model significant => print patch
  if (mlr==1 || logr==1 || mySVM==1 || myFET==1) 
  {
  ## Print column title
    print(c("Site", "Algo", "Score", "SupScore", "Add1", "Add2", "Add3", "Add4") )
    print(patchSites)
    print(paste("----------------------------------Above results for Patch in FILE: ", filei,"-----------------------------") )
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
            process = subprocess.Popen(["R","--vanilla","--slave","-f",tmp], stdout=f2)
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
            FINAL = open('result'+str(c)+'.txt','w')
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