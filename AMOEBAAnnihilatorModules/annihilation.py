import os
import sys
import time
import shutil
import numpy as np
import boxsetup as box
import keyfilemodifications as keymods
import minimization as mini     
import bar
import restraints as res
import plots
import tables
import productiondynamics as prod
import equilbriation as equil


def main(annihilator):
    shutil.copy(annihilator.keyfilename,annihilator.configkeyfilename)
    keymods.RemoveKeyWord(annihilator,annihilator.configkeyfilename,'parameters')
    keymods.RemoveKeyWord(annihilator,annihilator.configkeyfilename,'TARGET-DIPOLE')
    keymods.InsertKeyfileHeader(annihilator,annihilator.configkeyfilename)    
    box.BoxSetupProtocol(annihilator)
    mini.CheapMinimizationProtocol(annihilator)
    if annihilator.equiltimeNVT!=0 and annihilator.equiltimeNPT!=0:
        equil.EquilibriationProtocol(annihilator)
    if annihilator.proddyntime!=0:
        prod.ProductionDynamicsProtocol(annihilator) 
        bar.BARProtocol(annihilator)  
        tables.GenerateSimInfoTable(annihilator)
        plots.PlotEnergyData(annihilator)
