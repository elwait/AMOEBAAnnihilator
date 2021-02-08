import csv
import os
import sys
import productiondynamics as prod
import terminate as term
import numpy as np
import plots

def GenerateAnnihilationProgressTable(annihilator):
    tableheader=['Dynamic Jobs','Modename','Writeout Freq (ps)','Dyntamic Time Step (fs)','Total Time (ns)','Total Dyn Extended Time','Total Dynamic Steps','Total Extended Dynamic Steps','Total ARC File Space Needed']
    tableheader2=[annihilator.simfoldname,annihilator.simname,annihilator.proddynwritefreq,annihilator.proddyntimestep,annihilator.proddyntime,annihilator.proddyntime,annihilator.proddynsteps,annihilator.proddynsteps,annihilator.totalfilesize]
    elelambda=['Ele-Lambda']
    vdwlambda=['Vdw-Lambda']
    restlambda=['Restraint-Lambda']
    outputname=['Outputname']
    percentdone=['Percent Complete']
    stepsleft=['Number of Steps Left']
    eta=['ETA Bounds (days)']
    percentdoneext=['Percent Complete Extended']
    spaceneeded=['Space Needed Left (GB)']

    for i in range(len(annihilator.vdwlambdascheme)):
        elelamb=annihilator.estatlambdascheme[i]
        vdwlamb=annihilator.vdwlambdascheme[i]
        if len(annihilator.restlambdascheme)==len(annihilator.estatlambdascheme):
            rest=annihilator.restlambdascheme[i]
        else:
            rest=0
        elelambda.append(elelamb)
        vdwlambda.append(vdwlamb)
        restlambda.append(rest)
        fold=annihilator.lambdafolderlist[i]
        outputfilepath=annihilator.proddynoutfilepath[i]
        if os.path.exists(fold):
            os.chdir(fold)
            outputname.append(outputfilepath)
            finish=term.CheckFileTermination(outputfilepath,annihilator.proddynsteps)
            if finish==True:
                 percent='100'
                 steps='0'
                 etastring='N/A'
                 spaceleft='0' 
            else:
                 stepstaken=prod.CheckLastNumberDynamicStepsCompleted(annihilator,outputfilepath) 
                 if stepstaken==None:
                     percent='N/A'
                     steps='N/A'
                     etastring='N/A'
                     spaceleft='0'
                 else:
                     percent=100*(stepstaken/int(annihilator.proddynsteps))
                     steps=float(int(annihilator.proddynsteps)-stepstaken)
                     timeleft=steps*float(annihilator.proddyntimestep)*10**-6
                     #lower,upper=annihilator.ComputeDynamicsTime(timeleft)
                     lower='N/A'
                     upper='N/A'
                     etastring=str(lower)+','+str(upper)
                     spaceleft=((100-percent)/100)*float(annihilator.totalfilesize)
            percentdone.append(percent)
            stepsleft.append(steps)
            eta.append(etastring)
            spaceneeded.append(spaceleft)
 

            os.chdir('..')
    tabledata=[tableheader,tableheader2,elelambda,vdwlambda,restlambda,outputname,percentdone,stepsleft,eta,spaceneeded]

    return tabledata

def ComputeDynamicsTime(annihilator,timeleft):
    upper=timeleft/annihilator.lowerperf
    lower=timeleft/annihilator.upperperf
    return lower,upper


def KeyLists(annihilator):
    annihilator.masterdict={}
    annihilator.masterdict['energy']={}
    annihilator.masterdict['freeenergy']={}
    annihilator.masterdict['summary']={}
    annihilator.masterdict['boxinfo']={}
    annihilator.masterdict['lambda']={}

    lambdakeylist=['Ele-Lambda','Vdw-Lambda','Restraint-Lambda']

    boxinfokeylist=['Total Atom Number','Average Box Size','Prod MD Ensemb','Prod MD Time','Prod MD Steps','Prod MD Arc File Space','Dynamic Writeout Frequency (ps)','Dynamic Time Step (fs)','Equil Time NPT','Equil Time NVT','Physio Counterions','Neut Counterions','Ligand Charge','Receptor Charge','Ligand Name','Receptor Name']
    energykeylist=[u'ΔGˢᵒˡᵛ',u'ΔGˢᵒˡᵛᵉʳʳ',u'ΔHˢᵒˡᵛ',u'ΔHˢᵒˡᵛᵉʳʳ',u'ΔSˢᵒˡᵛ',u'ΔSˢᵒˡᵛᵉʳʳ',u'ΔGᶜᵒᵐᵖᶜᵒʳʳ',u'ΔGᶜᵒᵐᵖᵘⁿᶜᵒʳʳ',u'ΔGᶜᵒᵐᵖᶜᵒʳʳᵉʳʳ',u'ΔHᶜᵒᵐᵖ',u'ΔHᶜᵒᵐᵖᵉʳʳ',u'ΔSᶜᵒᵐᵖ',u'ΔSᶜᵒᵐᵖᵉʳʳ',u'ΔGᵃⁿᵃᶜᵒᵐᵖᶜᵒʳʳ',u'ΔGᵇᶦⁿᵈᶜᵒʳʳ',u'ΔGᵇᶦⁿᵈᶜᵒʳʳᵉʳʳ',u'ΔHᵇᶦⁿᵈ',u'ΔHᵇᶦⁿᵈᵉʳʳ',u'ΔSᵇᶦⁿᵈ',u'ΔSᵇᶦⁿᵈᵉʳʳ']
    freeenergykeylist=[u'ΔGˢᵒˡᵛ',u'ΔGˢᵒˡᵛᵉʳʳ',u'ΔGᶜᵒᵐᵖᶜᵒʳʳ',u'ΔGᶜᵒᵐᵖᵘⁿᶜᵒʳʳ',u'ΔGᶜᵒᵐᵖᶜᵒʳʳᵉʳʳ',u'ΔGᵃⁿᵃᶜᵒᵐᵖᶜᵒʳʳ',u'ΔGᵇᶦⁿᵈᶜᵒʳʳ',u'ΔGᵇᶦⁿᵈᶜᵒʳʳᵉʳʳ']
    summarykeylist=[u'ΔGˢᵒˡᵛ',u'ΔGᶜᵒᵐᵖᶜᵒʳʳ',u'ΔGᵃⁿᵃᶜᵒᵐᵖᶜᵒʳʳ',u'ΔGᵇᶦⁿᵈᶜᵒʳʳ']
        
    if annihilator.averageenergies==True:
        energykeylist.extend([u'ΔGᵇᶦⁿᵈᶜᵒʳʳᵃᵛᵍ',u'ΔGᵇᶦⁿᵈᶜᵒʳʳᵃᵛᵍᵉʳʳ',u'ΔHᵇᶦⁿᵈᵃᵛᵍ',u'ΔSᵇᶦⁿᵈᵃᵛᵍ'])
        freeenergykeylist.extend([u'ΔGᵇᶦⁿᵈᶜᵒʳʳᵃᵛᵍ',u'ΔGᵇᶦⁿᵈᶜᵒʳʳᵃᵛᵍᵉʳʳ'])
        summarykeylist.extend([u'ΔGᵇᶦⁿᵈᶜᵒʳʳᵃᵛᵍ'])
    keylist=[]
    keylist.extend(lambdakeylist)
    keylist.extend(boxinfokeylist)
    keylist.extend(energykeylist)
    keylist.extend(freeenergykeylist)
    keylist.extend(summarykeylist)
    return lambdakeylist,boxinfokeylist,energykeylist,freeenergykeylist,summarykeylist,keylist


def CSVWriter(annihilator,tempname,progress=False):
    with open(tempname, mode='w') as energy_file:
        energy_writer = csv.writer(energy_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        grabbedfreeenergydict=annihilator.masterdict['freeenergy']
        grabbedenergydict=annihilator.masterdict['energy']
        grabbedboxinfodict=annihilator.masterdict['boxinfo']
        grabbeddictlist=[grabbedfreeenergydict,grabbedenergydict,grabbedboxinfodict]
        wrotecolumnheaderslist=[False,False,False]
        tablesummarylist=['Gibbs Free Energy Change Table','Enthalpy, Entropy, Gibbs Energy Change Table','Simulation Info Table']
        for dictidx in range(len(grabbeddictlist)):
            tabledict=grabbeddictlist[dictidx]
            for path in tabledict.keys():
                wrotecolumnheaders=wrotecolumnheaderslist[dictidx]
                table=tabledict[path]
                keylist=list(table.keys())
                valuelist=list(table.values())
                keyrowlist=['Name']
                keyrowlist.extend(keylist)
                emptyline=[None]*len(keyrowlist)
                summaryline=[None]*len(keyrowlist)
                summary=tablesummarylist[dictidx]
                summaryline[0]=summary
                head, tail = os.path.split(path)
                valuerowlist=[tail]
                valuerowlist.extend(valuelist)
                if wrotecolumnheaders==False:
                    energy_writer.writerow(emptyline)
                    energy_writer.writerow(summaryline)
                    energy_writer.writerow(keyrowlist)
                    wrotecolumnheaderslist[dictidx]=True
                energy_writer.writerow(valuerowlist) 
        if progress==True:
            progtable=GenerateAnnihilationProgressTable(annihilator)
            if progtable!=None:
                energy_writer.writerow(emptyline)
                for row in progtable:
                    energy_writer.writerow(row)
                energy_writer.writerow(emptyline)



def GenerateSimInfoTable(annihilator):
    
    curdir=os.getcwd()
    os.chdir(annihilator.outputpath)
    lambdakeylist,boxinfokeylist,energykeylist,freeenergykeylist,summarykeylist,keylist=KeyLists(annihilator)

    OrderTableData(annihilator,boxinfokeylist,lambdakeylist,energykeylist,freeenergykeylist,summarykeylist,annihilator.outputpath)  
    tempname='SimData.csv'
    CSVWriter(annihilator,tempname,True)
    os.chdir(curdir)


def OrderTableData(annihilator,boxinfokeylist,lambdakeylist,energykeylist,freeenergykeylist,summarykeylist,path):
    if path not in annihilator.masterdict['boxinfo'].keys():
        annihilator.masterdict['boxinfo'][path]={}
    for key in boxinfokeylist:
        if key in annihilator.tabledict.keys():
            annihilator.masterdict['boxinfo'][path][key]=annihilator.tabledict[key]
        else:
            annihilator.masterdict['boxinfo'][path][key]=None
    if path not in annihilator.masterdict['lambda'].keys():
        annihilator.masterdict['lambda'][path]={}

    for key in lambdakeylist:
        if key in annihilator.tabledict.keys():
            annihilator.masterdict['lambda'][path][key]=annihilator.tabledict[key]
        else:
            annihilator.masterdict['lambda'][path][key]=None

    if path not in annihilator.masterdict['energy'].keys():
        annihilator.masterdict['energy'][path]={}

    for key in energykeylist:
        if key in annihilator.tabledict.keys():
            annihilator.masterdict['energy'][path][key]=annihilator.tabledict[key]
        else:
            annihilator.masterdict['energy'][path][key]=None

    if path not in annihilator.masterdict['freeenergy'].keys():
        annihilator.masterdict['freeenergy'][path]={}

    for key in freeenergykeylist:
        if key in annihilator.tabledict.keys():
            annihilator.masterdict['freeenergy'][path][key]=annihilator.tabledict[key]
        else:
            annihilator.masterdict['freeenergy'][path][key]=None

    if path not in annihilator.masterdict['summary'].keys():
        annihilator.masterdict['summary'][path]={}

    for key in summarykeylist:
        if key in annihilator.tabledict.keys():
            annihilator.masterdict['summary'][path][key]=annihilator.tabledict[key]
        else:
            annihilator.masterdict['summary'][path][key]=None


def GenerateFreeEnergyMatrix(annihilator,ligandnames,receptornames,bindenergies):
    uniquereceptors=len(set(receptornames))
    uniqueligands=len(ligandnames)
    mat=np.array((uniquereceptors,uniqueligands),dtype=float)
    receptortoligands={}
    receptorligandtoenergy={}
    for i in range(len(receptornames)):
        receptorname=receptornames[i]
        ligandname=ligandnames[i] 
        free=bindenergies[i]
        if receptorname not in receptortoligands.keys():
            receptortoligands[receptorname]=[]
        receptortoligands[receptorname].append(ligandname)
        if receptorname not in receptorligandtoenergy.keys():
            receptorligandtoenergy[receptorname]={}
        receptorligandtoenergy[receptorname][ligandname]=free
    count=0
    yaxis=[]
    xaxis=[]
    for receptor,ligands in receptortoligands.items(): 
        yaxis.append(receptor)
        for i in range(len(ligands)):
            ligand=ligands[i]
            energy=receptorligandtoenergy[receptor][ligand]
            if len(ligands)>1:
                mat[count,i]=float(energy)
            else:
                mat[count]=float(energy)

            xaxis.append(ligand)
            
        count+=1
    mat=np.transpose(mat)
    return mat,xaxis,yaxis

def GrabSimDataFromPathList(annihilator):
    curdir=os.getcwd()
    lambdakeylist,boxinfokeylist,energykeylist,freeenergykeylist,summarykeylist,keylist=KeyLists(annihilator)
    ligandnames=[]
    receptornames=[]
    bindenergies=[]
    for pathls in annihilator.simpathlist:
        for path in pathls:
            os.chdir(path)
            with open('SimData.csv') as csv_file:
                rows=[]
                csv_reader = csv.reader(csv_file, delimiter=',')
                for row in csv_reader:
                    rows.append(row)
                for i in range(len(rows)):
                    row=rows[i]
                    for j in range(len(row)):
                        col=row[j]
                        if col in keylist:
                            nextrow=rows[i+1]
                            colvalue=nextrow[j]
                            annihilator.tabledict[col]=colvalue
       
                   
            OrderTableData(annihilator,boxinfokeylist,lambdakeylist,energykeylist,freeenergykeylist,summarykeylist,path)  
        comppath=pathls[0]
        solvpath=pathls[1]
        EnterMatchingData(annihilator,comppath,solvpath)      
        ligname=annihilator.masterdict['boxinfo'][comppath]['Ligand Name']
        receptorname=annihilator.masterdict['boxinfo'][comppath]['Receptor Name']
        freeenergy=annihilator.masterdict['energy'][comppath][u'ΔGᵇᶦⁿᵈᶜᵒʳʳ']
        ligandnames.append(ligname)
        receptornames.append(receptorname)
        bindenergies.append(freeenergy) 

    os.chdir(curdir)
    mat,xaxis,yaxis=GenerateFreeEnergyMatrix(annihilator,ligandnames,receptornames,bindenergies) 
    plots.PlotHeatmap(annihilator,mat,xaxis,yaxis)
    # now need to try and match existing solvation to complexation data


    if annihilator.averageenergies==True:
        delGbindmodelist=[]
        delGbindmodelisterr=[]
        delHbindmodelist=[]
        delSbindmodelist=[]
        for path in annihilator.masterdict['energy'].keys():
            if annihilator.masterdict['energy'][path][u'ΔGᵇᶦⁿᵈᶜᵒʳʳ']!=None:
                delGbindmodelist.append(annihilator.masterdict['energy'][path][u'ΔGᵇᶦⁿᵈᶜᵒʳʳ'])
            if annihilator.masterdict['energy'][path][u'ΔGᵇᶦⁿᵈᶜᵒʳʳᵉʳʳ']!=None:
                delGbindmodelisterr.append(annihilator.masterdict['energy'][path][u'ΔGᵇᶦⁿᵈᶜᵒʳʳᵉʳʳ'])
            if annihilator.masterdict['energy'][path][u'ΔHᵇᶦⁿᵈ']!=None:
                delHbindmodelist.append(annihilator.masterdict['energy'][path][u'ΔHᵇᶦⁿᵈ'])
            if annihilator.masterdict['energy'][path][u'ΔSᵇᶦⁿᵈ']!=None:
                delSbindmodelist.append(annihilator.masterdict['energy'][path][u'ΔSᵇᶦⁿᵈ'])

        AverageEnergyList(annihilator,delGbindmodelist,'G')
        AverageEnergyList(annihilator,delHbindmodelist,'H')
        AverageEnergyList(annihilator,delSbindmodelist,'S')
        annihilator.deltaGaverageerr=TotalBindingAffinityError(annihilator,delGbindmodelist,delGbindmodelisterr)
        for path in annihilator.masterdict['energy'].keys():
            annihilator.masterdict['freeenergy'][path][u'ΔGᵇᶦⁿᵈᶜᵒʳʳᵃᵛᵍ']=annihilator.deltaGaverage
            annihilator.masterdict['energy'][path][u'ΔGᵇᶦⁿᵈᶜᵒʳʳᵃᵛᵍ']=annihilator.deltaGaverage
            annihilator.masterdict['freeenergy'][path][u'ΔGᵇᶦⁿᵈᶜᵒʳʳᵃᵛᵍᵉʳʳ']=annihilator.deltaGaverageerr
            annihilator.masterdict['energy'][path][u'ΔGᵇᶦⁿᵈᶜᵒʳʳᵃᵛᵍᵉʳʳ']=annihilator.deltaGaverageerr
            annihilator.masterdict['summary'][path][u'ΔGᵇᶦⁿᵈᶜᵒʳʳᵃᵛᵍ']=annihilator.deltaGaverage
            annihilator.masterdict['energy'][path][u'ΔHᵇᶦⁿᵈᵃᵛᵍ']=annihilator.deltaHaverage
            annihilator.masterdict['energy'][path][u'ΔSᵇᶦⁿᵈᵃᵛᵍ']=annihilator.deltaSaverage


   
 
    tempname='GrabbedSimData.csv'
    CSVWriter(annihilator,tempname)

    return

def EnterMatchingData(annihilator,path,match):
    annihilator.masterdict['energy'][path][u'ΔSˢᵒˡᵛ']=annihilator.masterdict['energy'][match][u'ΔSˢᵒˡᵛ']
    annihilator.masterdict['energy'][path][u'ΔGˢᵒˡᵛ']=annihilator.masterdict['energy'][match][u'ΔGˢᵒˡᵛ']
    annihilator.masterdict['freeenergy'][path][u'ΔGˢᵒˡᵛ']=annihilator.masterdict['freeenergy'][match][u'ΔGˢᵒˡᵛ']
    annihilator.masterdict['summary'][path][u'ΔGˢᵒˡᵛ']=annihilator.masterdict['freeenergy'][match][u'ΔGˢᵒˡᵛ']
    annihilator.masterdict['energy'][path][u'ΔHˢᵒˡᵛ']=annihilator.masterdict['energy'][match][u'ΔHˢᵒˡᵛ']
    annihilator.masterdict['energy'][path][u'ΔSˢᵒˡᵛᵉʳʳ']=annihilator.masterdict['energy'][match][u'ΔSˢᵒˡᵛᵉʳʳ']
    annihilator.masterdict['energy'][path][u'ΔGˢᵒˡᵛᵉʳʳ']=annihilator.masterdict['energy'][match][u'ΔGˢᵒˡᵛᵉʳʳ']
    annihilator.masterdict['freeenergy'][path][u'ΔGˢᵒˡᵛᵉʳʳ']=annihilator.masterdict['freeenergy'][match][u'ΔGˢᵒˡᵛᵉʳʳ']
    annihilator.masterdict['energy'][path][u'ΔHˢᵒˡᵛᵉʳʳ']=annihilator.masterdict['energy'][match][u'ΔHˢᵒˡᵛᵉʳʳ']
    
    annihilator.masterdict['energy'][match][u'ΔSᶜᵒᵐᵖ']=annihilator.masterdict['energy'][path][u'ΔSᶜᵒᵐᵖ']
    annihilator.masterdict['energy'][match][u'ΔGᶜᵒᵐᵖᶜᵒʳʳ']=annihilator.masterdict['energy'][path][u'ΔGᶜᵒᵐᵖᶜᵒʳʳ']
    annihilator.masterdict['freeenergy'][match][u'ΔGᶜᵒᵐᵖᶜᵒʳʳ']=annihilator.masterdict['freeenergy'][path][u'ΔGᶜᵒᵐᵖᶜᵒʳʳ']
    annihilator.masterdict['summary'][match][u'ΔGᶜᵒᵐᵖᶜᵒʳʳ']=annihilator.masterdict['freeenergy'][path][u'ΔGᶜᵒᵐᵖᶜᵒʳʳ']
    annihilator.masterdict['energy'][match][u'ΔHᶜᵒᵐᵖ']=annihilator.masterdict['energy'][path][u'ΔHᶜᵒᵐᵖ']
    annihilator.masterdict['energy'][match][u'ΔSᶜᵒᵐᵖᵉʳʳ']=annihilator.masterdict['energy'][path][u'ΔSᶜᵒᵐᵖᵉʳʳ']
    annihilator.masterdict['energy'][match][u'ΔGᶜᵒᵐᵖᶜᵒʳʳᵉʳʳ']=annihilator.masterdict['energy'][path][u'ΔGᶜᵒᵐᵖᶜᵒʳʳᵉʳʳ']
    annihilator.masterdict['freeenergy'][match][u'ΔGᶜᵒᵐᵖᶜᵒʳʳᵉʳʳ']=annihilator.masterdict['freeenergy'][path][u'ΔGᶜᵒᵐᵖᶜᵒʳʳᵉʳʳ']
    annihilator.masterdict['energy'][match][u'ΔHᶜᵒᵐᵖᵉʳʳ']=annihilator.masterdict['energy'][path][u'ΔHᶜᵒᵐᵖᵉʳʳ']


    DelGBind=float(annihilator.masterdict['energy'][path][u'ΔGᶜᵒᵐᵖᶜᵒʳʳ'])-float(annihilator.masterdict['energy'][path][u'ΔGˢᵒˡᵛ'])
    DelSBind=float(annihilator.masterdict['energy'][path][u'ΔSᶜᵒᵐᵖ'])-float(annihilator.masterdict['energy'][path][u'ΔSˢᵒˡᵛ'])
    DelHBind=float(annihilator.masterdict['energy'][path][u'ΔHᶜᵒᵐᵖ'])-float(annihilator.masterdict['energy'][path][u'ΔHˢᵒˡᵛ'])
    solvGerr=float(annihilator.masterdict['energy'][path][u'ΔGˢᵒˡᵛᵉʳʳ'])
    compGerr=float(annihilator.masterdict['energy'][path][u'ΔGᶜᵒᵐᵖᶜᵒʳʳᵉʳʳ'])
    DelGBinderr=np.sqrt(solvGerr**2+compGerr**2)
    solvHerr=float(annihilator.masterdict['energy'][path][u'ΔHˢᵒˡᵛᵉʳʳ'])
    compHerr=float(annihilator.masterdict['energy'][path][u'ΔHᶜᵒᵐᵖᵉʳʳ'])
    DelHBinderr=np.sqrt(solvHerr**2+compHerr**2)
    solvSerr=float(annihilator.masterdict['energy'][path][u'ΔSˢᵒˡᵛᵉʳʳ'])
    compSerr=float(annihilator.masterdict['energy'][path][u'ΔSᶜᵒᵐᵖᵉʳʳ'])
    DelSBinderr=np.sqrt(solvSerr**2+compSerr**2)
    annihilator.masterdict['energy'][path][u'ΔGᵇᶦⁿᵈᶜᵒʳʳ']=str(DelGBind)
    annihilator.masterdict['energy'][path][u'ΔGᵇᶦⁿᵈᶜᵒʳʳᵉʳʳ']=str(DelGBinderr)
    annihilator.masterdict['energy'][path][u'ΔHᵇᶦⁿᵈ']=str(DelHBind)
    annihilator.masterdict['energy'][path][u'ΔHᵇᶦⁿᵈᵉʳʳ']=str(DelHBinderr)
    annihilator.masterdict['energy'][path][u'ΔSᵇᶦⁿᵈ']=str(DelSBind)
    annihilator.masterdict['energy'][path][u'ΔSᵇᶦⁿᵈᵉʳʳ']=str(DelSBinderr)
    annihilator.masterdict['freeenergy'][path][u'ΔGᵇᶦⁿᵈᶜᵒʳʳ']=str(DelGBind)
    annihilator.masterdict['summary'][path][u'ΔGᵇᶦⁿᵈᶜᵒʳʳ']=str(DelGBind)
    annihilator.masterdict['freeenergy'][path][u'ΔGᵇᶦⁿᵈᶜᵒʳʳᵉʳʳ']=str(DelGBinderr)


    annihilator.masterdict['energy'][match][u'ΔGᵇᶦⁿᵈᶜᵒʳʳ']=str(DelGBind)
    annihilator.masterdict['energy'][match][u'ΔGᵇᶦⁿᵈᶜᵒʳʳᵉʳʳ']=str(DelGBinderr)
    annihilator.masterdict['energy'][match][u'ΔHᵇᶦⁿᵈ']=str(DelHBind)
    annihilator.masterdict['energy'][match][u'ΔHᵇᶦⁿᵈᵉʳʳ']=str(DelHBinderr)
    annihilator.masterdict['energy'][match][u'ΔSᵇᶦⁿᵈ']=str(DelSBind)
    annihilator.masterdict['energy'][match][u'ΔSᵇᶦⁿᵈᵉʳʳ']=str(DelSBinderr)
    annihilator.masterdict['freeenergy'][match][u'ΔGᵇᶦⁿᵈᶜᵒʳʳ']=str(DelGBind)
    annihilator.masterdict['summary'][match][u'ΔGᵇᶦⁿᵈᶜᵒʳʳ']=str(DelGBind)
    annihilator.masterdict['freeenergy'][match][u'ΔGᵇᶦⁿᵈᶜᵒʳʳᵉʳʳ']=str(DelGBinderr)


def AverageEnergyList(annihilator,enlist,key):
    if 'G' in key:
        annihilator.deltaGaverage=TotalBindingAffinity(annihilator,enlist)
    elif 'H' in key:
        annihilator.deltaHaverage=BoltzmannAverage(annihilator,enlist)
    elif 'S' in key:
        annihilator.deltaSaverage=(annihilator.deltaHaverage-annihilator.deltaGaverage)/int(annihilator.roomtemp)

def TotalBindingAffinity(annihilator,enlist):
    kB=0.0019872041 # kcal/molK
    T=int(annihilator.roomtemp)
    sumboltzfactors=0
    for obs in enlist:
        sumboltzfactors+=np.exp(-obs/(kB*T))
    totalsum=-kB*T*np.log(sumboltzfactors)
    return totalsum

def TotalBindingAffinityError(annihilator,enlist,errlist):
    kB=0.0019872041 # kcal/molK
    T=int(annihilator.roomtemp)
    Sum=0
    for obsidx in range(len(enlist)):
        obs=enlist[obsidx]
        err=errlist[obsidx]
        Sum+=(err/obs)**2
    return np.sqrt(Sum)

                
def BoltzmannAverage(annihilator,enlist):
    kB=0.0019872041
    T=int(annihilator.roomtemp)
    sumboltzfactors=0
    for obs in enlist:      
        sumboltzfactors+=np.exp(-obs/(kB*T))
    boltzavg=0
    for obs in enlist:
        boltzfactor=np.exp(-obs/(kB*T))
        boltzavg+=obs*boltzfactor
    boltzavg=boltzavg/sumboltzfactors
    return boltzavg  


def WriteTableUpdateToLog(annihilator):
    for key in annihilator.tabledict.keys():
        value=annihilator.tabledict[key]
        if key not in annihilator.tabledictkeysused:
            annihilator.tabledictkeysused.append(key)
            annihilator.WriteToLog(key+' = '+str(value))

