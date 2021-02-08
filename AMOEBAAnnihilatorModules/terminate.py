import os

def CheckFileTermination(annihilator,f):
   term=False
   if os.path.exists(f):
       temp=open(f,'r')
       results=temp.readlines()
       temp.close()
       for line in results:
           if 'Performance:' in line or 'Normal Termination' in line or 'Potential Energy Values Written To :' in line or 'Free Energy Difference via BAR' in line:
               term=True
           if 'Dynamics Steps' in line and 'Average' in line:
               linesplit=line.split()
               stepnum=float(linesplit[8])
               steps=float(annihilator.proddynsteps)
               ratio=np.abs(steps-stepnum)/steps
               if ratio<=.01:       
                   term=True 

   return term

def CheckFilesTermination(annihilator,outputfilepathlist):
    finished=True
    for i in range(len(outputfilepathlist)):
        outputfilepath=outputfilepathlist[i]
        terminate=CheckFileTermination(annihilator,outputfilepath)
        if terminate==False:
            finished=False
            annihilator.nextfiletofinish=outputfilepath
            return finished
    return finished

