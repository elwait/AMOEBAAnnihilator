import os
import sys
import subprocess

def CallExternalAPI(annihilator,jobtolog,jobtojobpath,jobinfofilepath):
    annihilator.WriteToLog('Calling external API ')
    temp=open(jobinfofilepath,'w')
    for job,log in jobtolog.items():
        jobpath=jobtojobpath[job]
        temp.write('--job='+job+' '+'--outputlogpath='+log+' '+'--jobpath='+jobpath+'\n')
    temp.close()
    if annihilator.bashrcpath!=None:
        cmdstr='python'+' '+annihilator.externalapi+' '+'--bashrcpath='+annihilator.bashrcpath+' '+'--jobinfofilepath='+jobinfofilepath
    else:
        cmdstr='python'+' '+annihilator.externalapi+' '+'--jobinfofilepath='+jobinfofilepath

    print('cmdstr ',cmdstr,flush=True)
    call_subsystem(annihilator,cmdstr,wait=False,skiperrors=False)



def SubmitJobs(annihilator,jobtolog,jobtojobpath,jobinfofilepath):
    if annihilator.externalapi!=None:
        if len(jobtolog.keys())!=0:
            CallExternalAPI(annihilator,jobtolog,jobtojobpath,jobinfofilepath)
    else:
        CallJobsSeriallyLocalHost(annihilator,jobtolog)


def CallJobsSeriallyLocalHost(annihilator,jobtolog):
    for job in jobtolog.keys():
        call_subsystem(annihilator,job,wait=True)

def call_subsystem(annihilator,cmdstr,wait=False,skiperrors=False):
    annihilator.WriteToLog(" Calling: " + cmdstr+' '+'path'+' = '+os.getcwd())
    p = subprocess.Popen(cmdstr, shell=True,stdout=annihilator.logfh, stderr=annihilator.logfh)
    if wait==True:
        p.wait()
        if p.returncode != 0 and skiperrors==False:
            annihilator.WriteToLog("ERROR: " + cmdstr+' '+'path'+' = '+os.getcwd())
            raise ValueError("ERROR: " + cmdstr+' '+'path'+' = '+os.getcwd())


