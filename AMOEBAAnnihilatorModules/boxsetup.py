import keyfilemodifications as keymods
import tables
import numpy as np
import itertools
import os
import submitjobs as submit
import numpy as np

def ComputeBoxSize(annihilator):
    longestdim=FindDimensionsOfMoleculeTinker(annihilator,annihilator.outputpath+annihilator.xyzfilename)
    if annihilator.fixedboxsize is None:
        annihilator.aaxis = round(longestdim + 2*annihilator.boxbufferlength +2*float(annihilator.vdwcutoff),1)
    else:
        annihilator.aaxis=annihilator.fixedboxsize
    annihilator.WriteToLog('Initial Box Length'+' '+str(annihilator.aaxis))
    annihilator.tabledict['Initial Box Length']=annihilator.aaxis
    tables.WriteTableUpdateToLog(annihilator)
    annihilator.volume=annihilator.aaxis**3
    keymods.AddKeyWord(annihilator,annihilator.outputpath+annihilator.configkeyfilename,'a-axis '+str(annihilator.aaxis)+'\n')

        
def ComputePhysiologicalIonNumber(annihilator):
    if annihilator.addphysioions==True:
        commasplit=annihilator.listofsaltcons.split(',')
        for ioncomplex in commasplit:
            ioncomplex=ioncomplex.lstrip().rstrip()
            equalsignsplit=ioncomplex.split('=')
            ioncomplexstring=equalsignsplit[0].lstrip().rstrip()
            conc=float(equalsignsplit[1].lstrip().rstrip()) # in mM
            complexnum=conc*6.022*10**-7*annihilator.volume
            Sum=0
            for el in annihilator.elementsymtotinktype.keys():
                if el in ioncomplexstring:
                    nextindex=ioncomplexstring.find(el)+len(el) # if K returns K index +1, if Cl returns Cl index +2, want to know if there is 2 or 3... after element, nextindex should always be defined because last elemt will be ]
                    if ioncomplexstring[nextindex].isdigit():
                        multfactor=int(ioncomplexstring[nextindex])
                    else:
                        multfactor=1
                    tinktype=annihilator.elementsymtotinktype[el]
                    if el not in annihilator.tabledict.keys():
                        annihilator.tabledict[el]=0
                    annihilator.tabledict[el]+=int(round(complexnum*multfactor))
                    Sum+=int(round(complexnum*multfactor))
                    if tinktype in annihilator.iontypetoionnumberphysio.keys():
                        annihilator.iontypetoionnumberphysio[tinktype]+=int(round(complexnum*multfactor))
                    else:
                        annihilator.iontypetoionnumberphysio[tinktype]=int(round(complexnum*multfactor))
                    

  
        annihilator.tabledict['Physio Counterions']=Sum
        tables.WriteTableUpdateToLog(annihilator)

def ComputeNeutralizingIonNumber(annihilator,systemcharge): # this can be just ligand charge or ligand and receptor charge depending on whether this is complexation or solvation    
    tinktypetoelementsym={v: k for k, v in annihilator.elementsymtotinktype.items()}
    if systemcharge>0:
        cltinktype=annihilator.elementsymtotinktype['Cl']
        annihilator.iontypetoionnumberneut[cltinktype]=np.abs(int(systemcharge))
        if 'Cl' not in annihilator.tabledict.keys():
            annihilator.tabledict['Cl']=0
        annihilator.tabledict['Cl']+=np.abs(int(systemcharge))


    elif systemcharge<0:
        ktinktype=annihilator.elementsymtotinktype['K']
        annihilator.iontypetoionnumberneut[ktinktype]=np.abs(int(systemcharge))
        if 'K' not in annihilator.tabledict.keys():
            annihilator.tabledict['K']=0
        annihilator.tabledict['K']+=np.abs(int(systemcharge))

    annihilator.tabledict['Neut Counterions']=np.abs(int(systemcharge))
    tables.WriteTableUpdateToLog(annihilator)


def ComputeWaterNumber(annihilator):
    annihilator.waternum=int(round(.0334*annihilator.volume)) # number of waters per cubic angstroms
    annihilator.WriteToLog('Water Number '+str(annihilator.waternum))

def RemoveBoxSizeTerms(annihilator,keyfile):
    keymods.RemoveKeyWord(annihilator,keyfile,'axis')


def FindDimensionsOfMoleculeTinker(annihilator,structurefilepath):
    veclist=[]
    temp=open(structurefilepath,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        linesplit=line.split()
        if len(linesplit)!=1 and '90.000000' not in line: # not line containing number of atoms
            vec=np.array([float(linesplit[2]),float(linesplit[3]),float(linesplit[4])])
            veclist.append(vec)

    pairs=list(itertools.combinations(veclist, 2))
    distlist=[]
    for pairidx in range(len(pairs)):
        pair=pairs[pairidx]
        progress=(pairidx*100)/len(pairs)
        dist=np.linalg.norm(np.array(pair[0])-np.array(pair[1]))
        distlist.append(dist)
    mindist=np.amax(np.array(distlist))
    return mindist


def TotalAtomNumber(annihilator,xyzfilename):
    atomnum=FindNumberTinkerXYZAtoms(annihilator,xyzfilename)
    annihilator.totalatomnumberxyzfilename=atomnum
    totalatomnum=atomnum+3*annihilator.waternum
    annihilator.tabledict['Total Atom Number']=totalatomnum
    tables.WriteTableUpdateToLog(annihilator)



def FindNumberTinkerXYZAtoms(annihilator,xyzfilename):
    temp=open(annihilator.outputpath+xyzfilename,'r')
    results=temp.readlines()
    temp.close()
    atomnum=int(results[0].replace('\n','').rstrip().lstrip())
    return atomnum

def GrabTypeNumber(annihilator,typedescrip):
    temp=open(annihilator.prmfilepath,'r')
    results=temp.readlines()
    temp.close()
    for line in results:
        if typedescrip in line:
            linesplit=line.split()
            typenum=linesplit[1]
    return typenum
   
def CreateWaterXYZ(annihilator):
    temp=open(annihilator.outputpath+'water.xyz','w')
    temp.write('3'+'\n')
    temp.write('1  O      8.019814    5.892935    0.449481   '+str(annihilator.waterOtypenum)+'     2     3'+'\n')
    temp.write('2  H      7.906681    4.942582    0.463215   '+str(annihilator.waterHtypenum)+'     1'+'\n')
    temp.write('3  H      8.149325    6.169963   -0.457590   '+str(annihilator.waterHtypenum)+'     1'+'\n')
    temp.close() 

      

def CreateSolventBox(annihilator):
    temp=open(annihilator.outputpath+'xyzedit.in','w')
    temp.write('water.xyz'+'\n')
    temp.write(annihilator.prmfilepath+'\n')
    temp.write('19'+'\n')
    temp.write(str(annihilator.waternum)+'\n')
    temp.write(str(annihilator.aaxis)+','+str(annihilator.aaxis)+','+str(annihilator.aaxis)+'\n')
    temp.write('Y'+'\n')
    temp.write(annihilator.prmfilepath+'\n')
    temp.close()
    cmdstr=annihilator.xyzeditpath+' '+'<'+' '+'xyzedit.in'
    submit.call_subsystem(annihilator,cmdstr,wait=True)    
 
def SoakMoleculeInSolventBox(annihilator,xyzfilename,keyfilename):
    cmdstr=annihilator.xyzeditpath+' '+xyzfilename+' '+'-k'+' '+keyfilename+' '+'20'+' '+'water.xyz_2'
    submit.call_subsystem(annihilator,cmdstr,wait=True)    

def AddIonToSolventBox(annihilator,solutexyzfilename,keyfilename,tinktype,ionnum,count):
    soluteatomnum=FindNumberTinkerXYZAtoms(annihilator,solutexyzfilename)
    num=2+count
    inputfile='xyzedit_'+str(num)+'.in'
    temp=open(annihilator.outputpath+inputfile,'w')
    temp.write('1'+' '+str(soluteatomnum)+'\n')
    string=''
    string+=str(tinktype)+' '+str(ionnum)+' '
    string+='\n'
    temp.write(string)
    temp.write('\n')
    temp.close()
    cmdstr=annihilator.xyzeditpath+' '+solutexyzfilename+'_'+str(num)+' '+'-k'+' '+keyfilename+' '+'21'+' '+' < '+inputfile
    submit.call_subsystem(annihilator,cmdstr,wait=True)    

def AddIonsToSolventBox(annihilator,solutexyzfilename,keyfilename):
    count=0
    for tinktype in annihilator.iontypetoionnumberneut.keys():
        ionnum=annihilator.iontypetoionnumberneut[tinktype]
        AddIonToSolventBox(annihilator,solutexyzfilename,keyfilename,tinktype,ionnum,count)
        count+=1
    for tinktype in annihilator.iontypetoionnumberphysio.keys():
        ionnum=annihilator.iontypetoionnumberphysio[tinktype]
        AddIonToSolventBox(annihilator,solutexyzfilename,keyfilename,tinktype,ionnum,count)
        count+=1
    os.rename(solutexyzfilename+'_'+str(count+2),annihilator.waterboxfilename)

def BoxSetupProtocol(annihilator):
    annihilator.WriteToLog('Computing volume ',prin=True)
    ComputeBoxSize(annihilator)
    if not os.path.isfile(annihilator.outputpath+annihilator.configkeyfilename):
        string='a-axis '+str(annihilator.aaxis)+'\n'
        keymods.AddKeyWord(annihilator,annihilator.configkeyfilename,string)
    ComputeWaterNumber(annihilator)
    ComputePhysiologicalIonNumber(annihilator)
    ComputeNeutralizingIonNumber(annihilator,annihilator.systemcharge)
    TotalAtomNumber(annihilator,annihilator.xyzfilename)
    if not os.path.isfile(annihilator.outputpath+'water.xyz'):
        CreateWaterXYZ(annihilator)
    if not os.path.isfile(annihilator.outputpath+annihilator.waterboxfilename):
        CreateSolventBox(annihilator)
        SoakMoleculeInSolventBox(annihilator,annihilator.xyzfilename,annihilator.configkeyfilename)
        AddIonsToSolventBox(annihilator,annihilator.xyzfilename,annihilator.configkeyfilename)
    annihilator.xyzfilesize=float(os.path.getsize(annihilator.waterboxfilename)) # in bytes
    annihilator.equilarcfilesize=annihilator.xyzfilesize*annihilator.equilframenum*10**-9 # in GB
    annihilator.singlepertubationfilesize=annihilator.xyzfilesize*int(annihilator.proddynframenum)*10**-9 # in GB
    annihilator.totalperturbationfilesize=len(annihilator.estatlambdascheme)*annihilator.singlepertubationfilesize # in GB
    annihilator.totalperturbationfilesize=len(annihilator.estatlambdascheme)*annihilator.singlepertubationfilesize # in GB
    annihilator.totalfilesize=annihilator.equilarcfilesize+annihilator.totalperturbationfilesize
    annihilator.tabledict['Prod MD Arc File Space']=annihilator.totalfilesize
    if annihilator.complexation==False:
        annihilator.ligandindices=[]
        for i in range(annihilator.totalatomnumberxyzfilename):
            index=i+1
            annihilator.ligandindices.append(index)
