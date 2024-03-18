import os
import sys
import subprocess


def Make_CondorScr(sampleName, outputName, NumFileList, NumJob, isData, eraName, roccorUse, idisoUse, trigUse, massCut) :
    path = os.getcwd()
    os.system("mkdir -p %s/%s"%(outputName,sampleName))
    path=path.replace("/Batch", "/")
    os.getcwd()
    os.system ('pwd')
    os.system ('mkdir -p Log')
    os.system ('ls')
    print path
    numjobs = NumJob
    numfiles = NumFileList
    SampleFile = sampleName
    SubFile = SampleFile + "_"
    Lists = MakeSeparateList(numfiles,numjobs)

    for i in range(0,numjobs):
        idx_ = i +1
        runing = outputName + "/"+SampleFile + "_%s.sh" % (idx_)
        f = open( runing, 'w')
        f.write('#!/bin/tcsh \n')
        f.write('setenv SCRAM_ARCH slc6_amd64_gcc530 \n')
        f.write('source /cvmfs/cms.cern.ch/cmsset_default.csh \n')
        f.write('setenv LD_PRELOAD "/usr/lib64/libpdcap.so" \n')
        f.write('cd '+ path + " \n")
        f.write('mkdir -p ./output/%s/%s \n'%(outputName,sampleName))
        f.write('cmsenv \n')
        SampList = MakeSampleIdxList(SampleFile,Lists[i])
        print "SampList %s"%(SampList)
        f.write('set inputlists = (%s) \n'%(SampList) )
        f.write('foreach i ( $inputlists )\n')
        f.write('   mkdir -p output \n')
        f.write('   ./ssb_analysis %s/${i}.list %s/%s/${i}.root 0 %s %s %s %s %s %s \n'%(sampleName, outputName, sampleName, isData, eraName, roccorUse, idisoUse, trigUse, massCut))
        f.write('end \n')
        f.close()
        runchMod = "chmod 755 " + runing
        os.system (runchMod)

        ### making submit File ###
        f1 = open( "./%s/%s/condor_%s_%s.submit"%(outputName,sampleName,SampleFile,idx_) , 'w')
        f1.write('# Unix submit description file\n')
        f1.write('Universe = vanilla\n')
        f1.write('x509userproxy = $ENV(X509_USER_PROXY)\n')
        f1.write('Executable = ./%s/%s_%s.sh\n'%(outputName,SampleFile,idx_))
        f1.write('request_memory = 400 MB \n')
        f1.write('should_transfer_files   = Yes\n')
        f1.write('Output      = ./%s/%s/%s_%s.output\n'%(outputName,sampleName,SampleFile,idx_))
        f1.write('error       = ./%s/%s/errors_%s_%s.txt\n'%(outputName,sampleName,SampleFile,idx_))
        f1.write('log         = ./%s/%s/test_%s_%s.log\n'%(outputName,sampleName,SampleFile,idx_))
        f1.write('Queue 1\n')
        f1.close()
        subchMod = "condor_submit " + "./%s/%s/condor_%s_%s.submit"%(outputName,sampleName,SampleFile,idx_)
        os.system (subchMod)
    pass

def MakeSampleIdxList(Sample,Lists):
    SampleList =""
    for index_ in Lists:
        print "index_ %s in MakeSampleIdxList "%(index_)
        SampleList += '"%s_%s"'%(Sample,index_)
        SampleList += " "
        pass
        #SamleList
    print SampleList
    return SampleList
    pass

def MakeSeparateList(NumFiles, NumJob):
    quo = NumFiles/NumJob
    seplists=[]
    for inx_ in range(NumJob) :
        emptyarray =[]
        seplists.append(emptyarray)
    print "size of seplists : %s "%len(seplists)
    for inumfile in  range (1,NumFiles+1):
        #print "%s"%(inumfile%NumJob)
        seplists[inumfile%NumJob].append(inumfile)
    for idx_ in seplists:
        print "content of %s "%(idx_)
    return seplists
    pass


#Make_CondorScr("SingleNeutrino_v1",40,40)
if __name__ == '__main__':
    Study = "240213_200GeV_Roccor" 

    #              sampleName,         outputName, NumFileList, NumJob, isData, eraName, roccorUse, idisoUse, trigUse, massCut

    Make_CondorScr("Data_SingleMuon_Run2016B"   , Study, 1054 , 106 , "true", "Run2016B",   "true", "false", "false", 200)
    Make_CondorScr("Data_SingleMuon_Run2016C"   , Study, 348  , 35  , "true", "Run2016C",   "true", "false", "false", 200)
    Make_CondorScr("Data_SingleMuon_Run2016D"   , Study, 584  , 59  , "true", "Run2016D",   "true", "false", "false", 200)
    Make_CondorScr("Data_SingleMuon_Run2016E"   , Study, 496  , 50  , "true", "Run2016E",   "true", "false", "false", 200)
    Make_CondorScr("Data_SingleMuon_Run2016F"   , Study, 362  , 37  , "true", "Run2016F",   "true", "false", "false", 200)
    Make_CondorScr("Data_SingleMuon_Run2016G"   , Study, 854  , 86  , "true", "Run2016G",   "true", "false", "false", 200)
    Make_CondorScr("Data_SingleMuon_Run2016HV2" , Study, 925  , 93  , "true", "Run2016HV2", "true", "false", "false", 200)
    Make_CondorScr("Data_SingleMuon_Run2016HV3" , Study, 25   , 3   , "true", "Run2016HV3", "true", "false", "false", 200)

    Make_CondorScr("DYJetsToLL_M_10To50"        , Study, 178  , 18 , "false", "DY_10To50",  "true", "false", "false", 200)
    Make_CondorScr("DYJetsToLL_M_50"            , Study, 732  , 74 , "false", "DY_50",      "true", "false", "false", 200)
    Make_CondorScr("ST_tW_antitop"              , Study, 70   , 7  , "false", "tW_antitop", "true", "false", "false", 200)
    Make_CondorScr("ST_tW_top"                  , Study, 71   , 8  , "false", "tW_top",     "true", "false", "false", 200)
    Make_CondorScr("TTJets_Signal"              , Study, 621  , 63 , "false", "TTJets",     "true", "false", "false", 200)
    Make_CondorScr("WJetsToLNu"                 , Study, 120  , 12 , "false", "WJets",      "true", "false", "false", 200)
    Make_CondorScr("WW"                         , Study, 81   , 9  , "false", "WW",         "true", "false", "false", 200)
    Make_CondorScr("WZ"                         , Study, 41   , 5  , "false", "WZ",         "true", "false", "false", 200)
    Make_CondorScr("ZZ"                         , Study, 21   , 3  , "false", "ZZ",         "true", "false", "false", 200)

