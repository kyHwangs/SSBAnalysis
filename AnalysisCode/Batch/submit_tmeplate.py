import os
from os import path
import sys
import subprocess


def Make_CondorScr(sampleName, outputName, NumFileList, isData, eraName, processName, roccorUse, idisoUse, trigUse, massCut) :

    os.system("mkdir -p %s/log/%s"%(outputName, sampleName))
    os.system("mkdir -p ../output/%s/%s"%(outputName, sampleName))

    JobListFileName = outputName + "/joblist.txt"
    JobListFile = open(JobListFileName, 'a')

    for i in range(1, NumFileList + 1):
        line = "%s/%s_%s.list %s/%s/%s_%s.root 0 %s %s %s %s %s %s %s %s \n"%(sampleName, sampleName, i, outputName, sampleName, sampleName, i, isData, eraName, processName, roccorUse, idisoUse, trigUse, massCut, i)
        JobListFile.write(line)

    JobListFile.close()

def MakeCondorSubmit(outputName):
    CondorFile = open("%s/condor_submit.sub"%(outputName), 'w')
    CondorFile.write('universe              = vanilla \n')
    CondorFile.write('executable            = condor_wrapper.sh \n')
    CondorFile.write('getenv                = True \n')
    CondorFile.write('x509userproxy         = $ENV(X509_USER_PROXY) \n')
    CondorFile.write('arguments             = ./ssb_analysis $(inputlist) $(outputroot) $(genloop) $(isdata) $(era) $(procname) $(isroccor) $(isidiso) $(istrigg) $(masscut) \n')
    CondorFile.write('\n')    

    CondorFile.write('request_memory        = 200 MB  \n')
    CondorFile.write('should_transfer_files = YES \n')
    CondorFile.write('\n')    

    CondorFile.write('JobBatchName          = %s \n'%(outputName))
    CondorFile.write('\n')    

    CondorFile.write('output                = log/$(era)_$(procname)/$(era)_$(procname)_$(processnum).out \n')
    CondorFile.write('error                 = log/$(era)_$(procname)/$(era)_$(procname)_$(processnum).err \n')
    CondorFile.write('log                   = log/$(era)_$(procname)/$(era)_$(procname)_$(processnum).log \n')
    CondorFile.write('queue inputlist,outputroot,genloop,isdata,era,procname,isroccor,isidiso,istrigg,masscut,processnum from joblist.txt \n')
    CondorFile.close()

def MakeCondorWrapper(outputName):
    CondorFile = open("%s/condor_wrapper.sh"%(outputName), 'w')
    CondorFile.write('#!/bin/sh \n')
    CondorFile.write('\n')    

    CondorFile.write('export SCRAM_ARCH=slc6_amd64_gcc530 \n')
    CondorFile.write('source /cvmfs/cms.cern.ch/cmsset_default.sh \n')
    CondorFile.write('export LD_PRELOAD="/usr/lib64/libpdcap.so" \n')
    CondorFile.write('\n')    

    CondorFile.write('cd /u/user/kyhwang/WorkingDir/CMS/HighMassDY/ssb_analysis/240605_FullUL/CMSSW_8_0_26_patch1/src/SSBAnalysis/AnalysisCode/ \n')
    CondorFile.write('\n')    

    CondorFile.write('cmsenv \n')
    CondorFile.write('\n')    

    CondorFile.write('echo "$@" \n')
    CondorFile.write('eval "$@" \n')
    CondorFile.close()

if __name__ == '__main__':
    Study = "240724_10GeV_noCorrection" 

    os.system("mkdir -p %s"%(Study))

    if (not path.exists("%s/condor_submit.sub"%(Study))):
        MakeCondorSubmit(Study)

    if (not path.exists("%s/condor_wrapper.sh"%(Study))):
        MakeCondorWrapper(Study)

    doRoccoR = "false" 
    doIDISO = "false" 
    doTRIG = "false" 

    Make_CondorScr("UL2016APV_NNLO_10to50",       Study,  13, "false", "UL2016APV", "NNLO_10to50",        doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016APV_NNLO_10to50_v2",    Study,  72, "false", "UL2016APV", "NNLO_10to50_v2",     doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016APV_NNLO_100to200",     Study,   7, "false", "UL2016APV", "NNLO_100to200",      doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016APV_NNLO_200to400",     Study,   5, "false", "UL2016APV", "NNLO_200to400",      doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016APV_NNLO_400to500",     Study,   2, "false", "UL2016APV", "NNLO_400to500",      doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016APV_NNLO_500to700",     Study,   3, "false", "UL2016APV", "NNLO_500to700",      doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016APV_NNLO_700to800",     Study,   2, "false", "UL2016APV", "NNLO_700to800",      doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016APV_NNLO_800to1000",    Study,   3, "false", "UL2016APV", "NNLO_800to1000",     doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016APV_NNLO_1000to1500",   Study,   2, "false", "UL2016APV", "NNLO_1000to1500",    doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016APV_NNLO_1500to2000",   Study,   2, "false", "UL2016APV", "NNLO_1500to2000",    doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016APV_NNLO_2000toInf",    Study,   2, "false", "UL2016APV", "NNLO_2000toInf",     doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016APV_NNLO_inc",          Study, 443, "false", "UL2016APV", "NNLO_inc",           doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016APV_NNLO_tautau",       Study, 228, "false", "UL2016APV", "NNLO_tautau",        doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016APV_ST_s",              Study,  23, "false", "UL2016APV", "ST_s",               doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016APV_ST_t_AntiTop",      Study, 126, "false", "UL2016APV", "ST_t_AntiTop",       doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016APV_ST_t_Top",          Study, 225, "false", "UL2016APV", "ST_t_Top",           doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016APV_ST_tW_AntiTop",     Study,  10, "false", "UL2016APV", "ST_tW_AntiTop",      doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016APV_ST_tW_Top",         Study,  10, "false", "UL2016APV", "ST_tW_Top",          doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016APV_TTTo2L2Nu",         Study, 151, "false", "UL2016APV", "TTTo2L2Nu",          doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016APV_WJetsToLNu",        Study, 119, "false", "UL2016APV", "WJetsToLNu",         doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016APV_WW",                Study,  13, "false", "UL2016APV", "WW",                 doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016APV_WZ",                Study,  32, "false", "UL2016APV", "WZ",                 doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016APV_ZZ",                Study,   6, "false", "UL2016APV", "ZZ",                 doRoccoR, doIDISO, doTRIG, 10)

    Make_CondorScr("UL2016APV_Run2016B_v2_HIPM",  Study,  53, "true", "UL2016APV", "Run2016B_v2_HIPM",    doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016APV_Run2016C_HIPM",     Study,  18, "true", "UL2016APV", "Run2016C_HIPM",       doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016APV_Run2016D_HIPM",     Study,  30, "true", "UL2016APV", "Run2016D_HIPM",       doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016APV_Run2016E_HIPM",     Study,  25, "true", "UL2016APV", "Run2016E_HIPM",       doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016APV_Run2016F_HIPM",     Study,  16, "true", "UL2016APV", "Run2016F_HIPM",       doRoccoR, doIDISO, doTRIG, 10)

    Make_CondorScr("UL2016_NNLO_10to50",          Study,  15, "false", "UL2016", "NNLO_10to50",           doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016_NNLO_10to50_v2",       Study,  85, "false", "UL2016", "NNLO_10to50_v2",        doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016_NNLO_100to200",        Study,   8, "false", "UL2016", "NNLO_100to200",         doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016_NNLO_200to400",        Study,   6, "false", "UL2016", "NNLO_200to400",         doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016_NNLO_400to500",        Study,   2, "false", "UL2016", "NNLO_400to500",         doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016_NNLO_500to700",        Study,   2, "false", "UL2016", "NNLO_500to700",         doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016_NNLO_700to800",        Study,   3, "false", "UL2016", "NNLO_700to800",         doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016_NNLO_800to1000",       Study,   4, "false", "UL2016", "NNLO_800to1000",        doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016_NNLO_1000to1500",      Study,   2, "false", "UL2016", "NNLO_1000to1500",       doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016_NNLO_1500to2000",      Study,   3, "false", "UL2016", "NNLO_1500to2000",       doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016_NNLO_2000toInf",       Study,   3, "false", "UL2016", "NNLO_2000toInf",        doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016_NNLO_inc",             Study, 356, "false", "UL2016", "NNLO_inc",              doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016_NNLO_tautau",          Study, 194, "false", "UL2016", "NNLO_tautau",           doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016_ST_s",                 Study,  23, "false", "UL2016", "ST_s",                  doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016_ST_t_AntiTop",         Study, 124, "false", "UL2016", "ST_t_AntiTop",          doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016_ST_t_Top",             Study, 254, "false", "UL2016", "ST_t_Top",              doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016_ST_tW_AntiTop",        Study,  11, "false", "UL2016", "ST_tW_AntiTop",         doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016_ST_tW_Top",            Study,  11, "false", "UL2016", "ST_tW_Top",             doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016_TTTo2L2Nu",            Study, 175, "false", "UL2016", "TTTo2L2Nu",             doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016_WJetsToLNu",           Study, 116, "false", "UL2016", "WJetsToLNu",            doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016_WW",                   Study,  12, "false", "UL2016", "WW",                    doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016_WZ",                   Study,  31, "false", "UL2016", "WZ",                    doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016_ZZ",                   Study,   5, "false", "UL2016", "ZZ",                    doRoccoR, doIDISO, doTRIG, 10)

    Make_CondorScr("UL2016_Run2016F",             Study,   3, "true", "UL2016", "Run2016F",               doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016_Run2016G",             Study,  43, "true", "UL2016", "Run2016G",               doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2016_Run2016H",             Study,  48, "true", "UL2016", "Run2016H",               doRoccoR, doIDISO, doTRIG, 10)

    Make_CondorScr("UL2017_NNLO_10to50",          Study,  30, "false", "UL2017", "NNLO_10to50",           doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2017_NNLO_10to50_v2",       Study, 166, "false", "UL2017", "NNLO_10to50_v2",        doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2017_NNLO_100to200",        Study,  15, "false", "UL2017", "NNLO_100to200",         doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2017_NNLO_200to400",        Study,   5, "false", "UL2017", "NNLO_200to400",         doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2017_NNLO_400to500",        Study,   2, "false", "UL2017", "NNLO_400to500",         doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2017_NNLO_500to700",        Study,   2, "false", "UL2017", "NNLO_500to700",         doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2017_NNLO_700to800",        Study,   4, "false", "UL2017", "NNLO_700to800",         doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2017_NNLO_800to1000",       Study,   2, "false", "UL2017", "NNLO_800to1000",        doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2017_NNLO_1000to1500",      Study,   4, "false", "UL2017", "NNLO_1000to1500",       doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2017_NNLO_1500to2000",      Study,   3, "false", "UL2017", "NNLO_1500to2000",       doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2017_NNLO_2000toInf",       Study,   2, "false", "UL2017", "NNLO_2000toInf",        doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2017_NNLO_inc",             Study, 661, "false", "UL2017", "NNLO_inc",              doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2017_NNLO_tautau",          Study, 176, "false", "UL2017", "NNLO_tautau",           doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2017_ST_s",                 Study,  56, "false", "UL2017", "ST_s",                  doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2017_ST_t_AntiTop",         Study, 280, "false", "UL2017", "ST_t_AntiTop",          doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2017_ST_t_Top",             Study, 521, "false", "UL2017", "ST_t_Top",              doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2017_ST_tW_AntiTop",        Study,  23, "false", "UL2017", "ST_tW_AntiTop",         doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2017_ST_tW_Top",            Study,  23, "false", "UL2017", "ST_tW_Top",             doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2017_TTTo2L2Nu",            Study, 427, "false", "UL2017", "TTTo2L2Nu",             doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2017_WJetsToLNu",           Study, 110, "false", "UL2017", "WJetsToLNu",            doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2017_WW",                   Study,  29, "false", "UL2017", "WW",                    doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2017_WZ",                   Study,  32, "false", "UL2017", "WZ",                    doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2017_ZZ",                   Study,  11, "false", "UL2017", "ZZ",                    doRoccoR, doIDISO, doTRIG, 10)

    Make_CondorScr("UL2017_Run2017B",             Study,  25, "true", "UL2017", "Run2017B",               doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2017_Run2017C",             Study,  53, "true", "UL2017", "Run2017C",               doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2017_Run2017D",             Study,  28, "true", "UL2017", "Run2017D",               doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2017_Run2017E",             Study,  43, "true", "UL2017", "Run2017E",               doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2017_Run2017F",             Study,  59, "true", "UL2017", "Run2017F",               doRoccoR, doIDISO, doTRIG, 10)

    Make_CondorScr("UL2018_NNLO_10to50",          Study,  41, "false", "UL2018", "NNLO_10to50",           doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2018_NNLO_10to50_v2",       Study, 240, "false", "UL2018", "NNLO_10to50_v2",        doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2018_NNLO_100to200",        Study,  21, "false", "UL2018", "NNLO_100to200",         doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2018_NNLO_200to400",        Study,   5, "false", "UL2018", "NNLO_200to400",         doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2018_NNLO_400to500",        Study,   2, "false", "UL2018", "NNLO_400to500",         doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2018_NNLO_500to700",        Study,   2, "false", "UL2018", "NNLO_500to700",         doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2018_NNLO_700to800",        Study,   2, "false", "UL2018", "NNLO_700to800",         doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2018_NNLO_800to1000",       Study,   2, "false", "UL2018", "NNLO_800to1000",        doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2018_NNLO_1000to1500",      Study,   3, "false", "UL2018", "NNLO_1000to1500",       doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2018_NNLO_1500to2000",      Study,   2, "false", "UL2018", "NNLO_1500to2000",       doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2018_NNLO_2000toInf",       Study,   2, "false", "UL2018", "NNLO_2000toInf",        doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2018_NNLO_inc",             Study, 918, "false", "UL2018", "NNLO_inc",              doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2018_NNLO_tautau",          Study, 243, "false", "UL2018", "NNLO_tautau",           doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2018_ST_s",                 Study,  79, "false", "UL2018", "ST_s",                  doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2018_ST_t_AntiTop",         Study, 384, "false", "UL2018", "ST_t_AntiTop",          doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2018_ST_t_Top",             Study, 716, "false", "UL2018", "ST_t_Top",              doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2018_ST_tW_AntiTop",        Study,  32, "false", "UL2018", "ST_tW_AntiTop",         doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2018_ST_tW_Top",            Study,  33, "false", "UL2018", "ST_tW_Top",             doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2018_TTTo2L2Nu",            Study, 585, "false", "UL2018", "TTTo2L2Nu",             doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2018_WJetsToLNu",           Study, 120, "false", "UL2018", "WJetsToLNu",            doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2018_WW",                   Study,  41, "false", "UL2018", "WW",                    doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2018_WZ",                   Study,  33, "false", "UL2018", "WZ",                    doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2018_ZZ",                   Study,  15, "false", "UL2018", "ZZ",                    doRoccoR, doIDISO, doTRIG, 10)

    Make_CondorScr("UL2018_Run2018A",             Study,  57, "true", "UL2018", "Run2018A",              doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2018_Run2018B",             Study,  27, "true", "UL2018", "Run2018B",              doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2018_Run2018C",             Study,  27, "true", "UL2018", "Run2018C",              doRoccoR, doIDISO, doTRIG, 10)
    Make_CondorScr("UL2018_Run2018D",             Study, 127, "true", "UL2018", "Run2018D",              doRoccoR, doIDISO, doTRIG, 10)



