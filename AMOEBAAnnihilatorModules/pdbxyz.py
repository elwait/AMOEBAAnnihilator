import os
import submitjobs as submit
import openbabel
import re

def GenerateProteinTinkerXYZFile(annihilator):
    if not os.path.isfile(annihilator.complexedxyzname):
        cmdstr=annihilator.pdbxyzpath+' '+annihilator.uncomplexedproteinpdbname+' '+annihilator.prmfilepath
        submit.call_subsystem(annihilator,cmdstr,wait=True)    
    atoms,coord,order,types,connections=readTXYZ(annihilator,annihilator.uncomplexedxyzname)
    uncomplexedatomnum=len(atoms) 
    indextocoordinates=GrabLigandCoordinates(annihilator,uncomplexedatomnum)
    annihilator.ligandindices=list(indextocoordinates.keys())
    if not os.path.isfile(annihilator.complexedxyzname):
        GenerateComplexedTinkerXYZFile(annihilator,annihilator.uncomplexedxyzname,indextocoordinates,uncomplexedatomnum) 
 
def readTXYZ(annihilator,TXYZ):
    temp=open(TXYZ,'r')
    lines = temp.readlines()[1:] #TINKER coordinate starts from second line
    atoms=[];coord=[]
    order=[];types=[];connections=[]
    for line in lines:
        data=line.split()
        order.append(data[0])
        types.append(data[5])
        connections.append(data[6:])
        atoms.append(data[1])
        coord.append([float(data[2]), float(data[3]), float(data[4])])
    return atoms,coord,order, types, connections


def GrabLigandCoordinates(annihilator,uncomplexedatomnum): # assumes appended to end of PDB file
    indextocoordinates={}
    pdbmol=openbabel.OBMol()
    obConversion = openbabel.OBConversion()
    obConversion.SetInFormat('pdb')
    obConversion.ReadFile(pdbmol,annihilator.complexedproteinpdbname)
    atomiter=openbabel.OBMolAtomIter(pdbmol)
    for atom in atomiter:
        atomidx=atom.GetIdx()
        if atomidx>uncomplexedatomnum:
            coords=[atom.GetX(),atom.GetY(),atom.GetZ()]
            indextocoordinates[atomidx]=coords
    if len(indextocoordinates.keys())==0:
        raise ValueError('Uncomplexed PDB and Complexed PDB missing atoms or ligand in Complexed PDB is not appended to the end of file')
    return indextocoordinates    

def GenerateComplexedTinkerXYZFile(annihilator,uncomplexedxyzname,indextocoordinates,uncomplexedatomnum):
    temp=open(uncomplexedxyzname,'r')
    results=temp.readlines()
    temp.close() 
    atoms,coord,order,types,connections=readTXYZ(annihilator,annihilator.ligandxyzfilename)
    temp=open(annihilator.complexedxyzname,'w')
    newatomnum=uncomplexedatomnum+len(indextocoordinates)
    for lineidx in range(len(results)):
        line=results[lineidx]
        linesplit=line.split()
        if lineidx==0:
            newline=str(newatomnum)+'\n' 
        else:
            newline=line
        temp.write(newline)
    for idx in range(len(atoms)):
        element=atoms[idx]
        oldindex=order[idx]
        newindex=int(oldindex)+uncomplexedatomnum
        coords=indextocoordinates[newindex]
        x=coords[0]
        y=coords[1]
        z=coords[2]
        typenum=types[idx]
        conns=connections[idx]
        conns=[int(i) for i in conns]
        conns=[i+uncomplexedatomnum for i in conns]
        newline='    '+str(newindex)+'  '+element+'     '+str(x)+'   '+str(y)+'   '+str(z)+'    '+str(typenum)+'     '
        for con in conns:
            newline+=str(con)+'     '
        newline+='\n'
        temp.write(newline)
    temp.close()
