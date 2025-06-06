#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 11/6/22 3:15 PM

@author: alejandrobertolet
"""

import numpy as np
from copy import deepcopy, copy

from ChronoDNARepair.repair import tracking, output, processes
from ChronoDNARepair.induction.damage import DamageToDNA

import random
import os

DAMAGED = 1
REPAIRED = 2
MISREPAIRED = 3

# Outputs
DSB = 0
MISREPDSB = 1
SSB = 2
BD = 3
SURVIVAL = 4
FOCI = 5

class Simulator:
    def __init__(self, timeOptions=[], cellmodel='standard', diffusionmodel='free', dsbmodel='standard', ssbmodel='standard', bdmodel='standard', nucleusMaxRadius=None,
                 irradiationTime=0, doseratefunction=None, doseratefunctionargs=None, cellparams=None, diffusionparams=None, dsbparams=None, ssbparams=None, bdparams=None,
                 p_lethal_aberration=0.1, p_apoptosis=0.1, time_to_G1S_checkpoint=10 * 3600, time_to_G2M_checkpoint=24 * 3600):
        self.messages = []
        self.irradiationTime = irradiationTime
        self.doseratefunction = doseratefunction
        self.doseratefunctionargs = doseratefunctionargs
        self.runManager = RunManager(timeOptions=timeOptions, cellmodel=cellmodel, diffusionmodel=diffusionmodel, cellParameters=cellparams,
                                     dsbrepairmodel=dsbmodel, ssbrepairmodel=ssbmodel, bdrepairmodel=bdmodel, nucleusMaxRadius=nucleusMaxRadius,
                                     messages=self.messages, diffusionParameters=diffusionparams, dsbrepairParameters=dsbparams, ssbrepairParameters=ssbparams,
                                     bdrepairParameters=bdparams, p_lethal_aberration=p_lethal_aberration, p_apoptosis=p_apoptosis, time_to_G1S_checkpoint=time_to_G1S_checkpoint,
                                     time_to_G2M_checkpoint=time_to_G2M_checkpoint)

    def Run(self, nRuns, rereadDamageForNewRuns=True, basepath=None, maxDose=-1, version='2.0', plot=True, outputnorm=True, verbose=0, getVideo=False):
        self.nRuns = nRuns
        self.runManager.nRuns = nRuns
        self.runManager.maxDose = maxDose
        self.runManager.plotflag = False
        if plot:
            self.runManager.plotflag = True
        if not rereadDamageForNewRuns:
            self.runManager.Run(verbose=verbose, outputnorm=outputnorm, getVideo=getVideo)
        else:
            self.runManager.TotalRuns = self.nRuns
            self.runManager.plotflag = False
            for i in range(self.nRuns):
                self.runManager.nRuns = 1
                if i == self.nRuns - 1 and plot:
                    self.runManager.plotflag = True
                if i == 0:
                    self.runManager.Run(verbose=verbose, outputnorm=outputnorm, getVideo=getVideo)
                else:
                    self.ReadDamage(basepath, maxDose, version)
                    self.runManager.Run(verbose=verbose, outputnorm=outputnorm, getVideo=getVideo)
        self.avgRemainingDSBOverTime = self.runManager.runoutputDSB
        self.avgRemainFociOverTime = self.runManager.runoutputFoci
        self.avgMisrepairedDSBOverTime = self.runManager.runoutputMisrepairedDSB
        self.celloutput = self.runManager.cellcultureoutput

    def _get_numbered_dirs(self, path):
        """Get all numbered directories regardless of prefix"""
        numbered_dirs = []
        try:
            for d in os.listdir(path):
                full_path = os.path.join(path, d)
                if os.path.isdir(full_path):
                    # Extract number from directory name using last digits
                    number = ''.join(filter(str.isdigit, d))
                    if number:  # Only add if there are digits in the name
                        numbered_dirs.append((int(number), d))
        except OSError:
            print(f"Error accessing directory: {path}")
            return []
        
        # Sort by the numeric value
        numbered_dirs.sort(key=lambda x: x[0])
        return [d[1] for d in numbered_dirs]  # Return only directory names

    def ReadDamage(self, basepath, maxDose=2.0, version='2.0', recalculatePerEachTrack=False):
        damage = DamageToDNA(messages=self.messages)
        npaths = 1
        if type(basepath) == list:
            npaths = len(basepath)
            if type(maxDose) is not list or len(maxDose) != npaths:
                print('Doses have to match the number of radiations. The same dose will be used for all')
                if type(maxDose) is not list:
                    maxDose = [maxDose] * npaths
                else:
                    maxDose = [maxDose[0]] * npaths
        else:
            basepath = [basepath]
            maxDose = [maxDose]
        totalaccumulateddose = 0
        for ib, bpath in enumerate(basepath):
            # Section to get what directories actually contains both dose and SDD. Disregard others!
            listOfAvailableDirs = []
            for dir_name in self._get_numbered_dirs(bpath):
                newpath = os.path.join(bpath, dir_name) 
                files = os.listdir(newpath)
                if 'DNADamage_sdd.txt' in files and 'DNADamage.phsp' in files:
                    if os.path.getsize(os.path.join(newpath, 'DNADamage_sdd.txt')) > 0:  # only those with actual data
                        listOfAvailableDirs.append(dir_name)
        
            if not listOfAvailableDirs:
                raise ValueError(f"No valid damage files found in directory {bpath}")

            accumulatedose = 0  # Initialize before the loop
            neworder = random.sample(listOfAvailableDirs, len(listOfAvailableDirs))
            for i, e in enumerate(neworder):
                accumulatedose = damage.accumulateDose - totalaccumulateddose
                if accumulatedose > maxDose[ib]:
                    break
                time = self._getTimeForDose(accumulatedose)
                if 0 < self.irradiationTime < time:
                    time = 1e20
                path = os.path.join(bpath, e)
                damage.readSDDAndDose(os.path.join(path, ''), version=version, particleTime=time, lesionTime=time)
            totalaccumulateddose += accumulatedose
        
        damage.populateDamages(getVideo=False, stopAtDose=-1, stopAtTime=0, recalculatePerEachTrack=recalculatePerEachTrack)
        self.runManager.damage = damage

    def LoadDamageFromMGM(self, listOfDamageSites):
        self.runManager.mgmFlag = True
        self.runManager.mgmDamage = listOfDamageSites
        self.runManager.damage = None

    def _getTimeForDose(self, d):
        if d == 0:
            return 0
        if self.irradiationTime == 0:
            return 0
        if self.doseratefunction is None:
            return 0
        if self.doseratefunction == 'uniform':
            return d / self.doseratefunctionargs[0]
        if self.doseratefunction == 'linear':
            return (-self.doseratefunction[0] + np.sqrt(self.doseratefunction[0] ** 2 + 4 * self.doseratefunction[1] * d)) / (2 * self.doseratefunction[1])
        if self.doseratefunction == 'exponential':
            constant = np.log(2) / self.doseratefunctionargs[1]
            initialdoserate = self.doseratefunctionargs[0]
            return -1/constant * np.log(1-constant*d/initialdoserate)

    def GetSurvivalFraction(self):
        return self.runManager.GetSurvivalFraction()

class RunManager:
    def __init__(self, timeOptions = [], cellmodel='standard', diffusionmodel='free', dsbrepairmodel='standard', ssbrepairmodel='standard',
                 bdrepairmodel='standard', nucleusMaxRadius = None, outputs=[DSB, FOCI, MISREPDSB, SSB, BD, SURVIVAL], messages=[],
                 cellParameters=None, diffusionParameters=None, dsbrepairParameters=None, ssbrepairParameters=None, bdrepairParameters=None,
                 p_lethal_aberration=0.1, p_apoptosis=0.1, time_to_G1S_checkpoint=10 * 3600, time_to_G2M_checkpoint=24 * 3600):
        self.messages = messages
        self.maxDose = -1
        self._diffusionactivated = False
        self._dsbrepactivated = False
        self._ssbrepactivated = False
        self._bdrepactivated = False
        self._cellActivated = False
        self.trackid = 0
        self.mgmFlag = False
        self.mgmDamage = None
        self.p_lethal_aberration = p_lethal_aberration
        self.p_apoptosis = p_apoptosis
        self.time_to_G1S_checkpoint = time_to_G1S_checkpoint
        self.time_to_G2M_checkpoint = time_to_G2M_checkpoint
        if cellmodel is not None:
            if cellmodel.lower() != 'none':
                self.CellActivated = True
                self.cell = processes.CellEvolution(cellmodel, cellParameters)
        if diffusionmodel is not None:
            if diffusionmodel.lower() != 'none':
                self.DiffusionActivated = True
                self.diffusionModel = processes.Diffusion(diffusionmodel, diffusionParameters)
        if dsbrepairmodel is not None:
            if dsbrepairmodel.lower() != 'none':
                self.DSBRepairActivated = True
                self.dsbRepairModel = processes.DSBRepair(dsbrepairmodel, dsbrepairParameters)
        if ssbrepairmodel is not None:
            if ssbrepairmodel.lower() != 'none':
                self.SSBRepairActivated = True
                self.ssbRepairModel = processes.SSBRepair(ssbrepairmodel, ssbrepairParameters)
        if bdrepairmodel is not None:
            if bdrepairmodel.lower() != 'none':
                self.BDRepairActivated = True
                self.bdRepairModel = processes.BDRepair(bdrepairmodel, bdrepairParameters)
        if len(timeOptions) > 3:
            self.clock = Clock(timeOptions[0], timeOptions[1], timeOptions[2], timeOptions[3])
        else:
            self.clock = Clock(timeOptions[0], timeOptions[1], timeOptions[2])
        self.nucleusMaxRadius = nucleusMaxRadius
        self.outputs = outputs
        self.runoutputDSB = output.AverageTimeCurveOverRuns()
        self.runoutputFoci = output.AverageTimeCurveOverRuns()
        self.runoutputMisrepairedDSB = output.AverageTimeCurveOverRuns()
        self.cellcultureoutput = output.CellCulture(cellmodel)
        self.outputsurvival = []
        self.plotflag = True
        self.currentrun = 0

    def InitializeNewRun(self):
        self.betracks = []
        self.ssbdamages = []
        self.bdamages = []
        self.foci = []
        if self.CellActivated:
            self.cell.Initialize()
            self.celloutput = output.CellOutput()

    def InitializeNewTracks(self, dam):
        if not self.mgmFlag:
            for iCh in dam.DSBMap:
                for iBp in dam.DSBMap[iCh]:
                    for iCo in dam.DSBMap[iCh][iBp]:
                        if dam.DSBMap[iCh][iBp][iCo].type > 0:
                            pos = dam.DSBMap[iCh][iBp][iCo].position
                            time = dam.DSBMap[iCh][iBp][iCo].particletime
                            if time <= self.clock.CurrentTime and time > self.clock.CurrentTime - self.clock.CurrentTimeStep:
                                newBeStep = tracking.BeStep(pos, time, complexity=dam.DSBMap[iCh][iBp][iCo].complexity)
                                # Sample heterochromain vs euchromatin
                                if self.CellActivated:
                                    p = self.cell.GetCurrentProportionOfHeterochromatin()
                                    if np.random.rand() < p:
                                        newBeStep.IsHeterochromatin = True
                                newBeTrack = tracking.BeTrack(self.trackid, dam.DSBMap[iCh][iBp][iCo].dsbID)
                                newBeTrack.ChromosomeID = iCh
                                newBeTrack.BasePairID = iBp
                                newBeTrack.StrandID = iCo
                                newBeTrack.AddNewStep(newBeStep)
                                self.betracks.append(newBeTrack)
                                self.trackid += 1
            for iCh in dam.SSBMap:
                for iBp in dam.SSBMap[iCh]:
                    for iCo in dam.SSBMap[iCh][iBp]:
                        if dam.SSBMap[iCh][iBp][iCo].type > 0:
                            pos = dam.SSBMap[iCh][iBp][iCo].position
                            time = dam.SSBMap[iCh][iBp][iCo].particletime
                            if time <= self.clock.CurrentTime and time > self.clock.CurrentTime - self.clock.CurrentTimeStep:
                                newSSBDamage = tracking.DamageTrack(self.trackid)
                                newSSBDamage.Time = time
                                newSSBDamage.Position = pos
                                newSSBDamage.ChromosomeID = iCh
                                newSSBDamage.BasePairID = iBp
                                newSSBDamage.StrandID = iCo
                                newSSBDamage.Complexity = dam.SSBMap[iCh][iBp][iCo].complexity
                                self.ssbdamages.append(newSSBDamage)
                                self.trackid += 1
            for iCh in dam.BDMap:
                for iBp in dam.BDMap[iCh]:
                    for iCo in dam.BDMap[iCh][iBp]:
                        if dam.BDMap[iCh][iBp][iCo].type > 0:
                            pos = dam.BDMap[iCh][iBp][iCo].position
                            time = dam.BDMap[iCh][iBp][iCo].particletime
                            if time <= self.clock.CurrentTime and time > self.clock.CurrentTime - self.clock.CurrentTimeStep:
                                newBDDamage = tracking.DamageTrack(self.trackid)
                                newBDDamage.Time = time
                                newBDDamage.Position = pos
                                newBDDamage.ChromosomeID = iCh
                                newBDDamage.BasePairID = iBp
                                newBDDamage.StrandID = iCo
                                newBDDamage.Complexity = dam.BDMap[iCh][iBp][iCo].complexity
                                self.bdamages.append(newBDDamage)
                                self.trackid += 1
        else:
            if self.clock.CurrentTime == self.clock.InitialTime:
                for i in range(len(self.mgmDamage)):
                    # First break end. Move it randomly in a 3 nm radius
                    position = self.mgmDamage[i][0] + np.random.normal(0, 3, 3) * 1e-3
                    complexity = self.mgmDamage[i][1]
                    time = 0.0
                    newBeStep = tracking.BeStep(position, time, complexity=complexity)
                    # Sample heterochromain vs euchromatin
                    if self.CellActivated:
                        p = self.cell.GetCurrentProportionOfHeterochromatin()
                        if np.random.rand() < p:
                            newBeStep.IsHeterochromatin = True
                    newBeTrack = tracking.BeTrack(self.trackid, i)
                    newBeTrack.ChromosomeID = np.random.randint(0, 46)
                    newBeTrack.BasePairID = i
                    newBeTrack.StrandID = 1
                    newBeTrack.AddNewStep(newBeStep)
                    self.betracks.append(newBeTrack)
                    self.trackid += 1
                    # Second break end. Move it randomly in a 3 nm radius
                    position = self.mgmDamage[i][0] + np.random.normal(0, 3, 3) * 1e-3
                    complexity = self.mgmDamage[i][1]
                    time = 0.0
                    newBeStep = tracking.BeStep(position, time, complexity=complexity)
                    newBeTrack = tracking.BeTrack(self.trackid, i)
                    newBeTrack.ChromosomeID = np.random.randint(0, 46)
                    newBeTrack.BasePairID = i
                    newBeTrack.StrandID = 0
                    newBeTrack.AddNewStep(newBeStep)
                    self.betracks.append(newBeTrack)
                    self.trackid += 1

    def Run(self, verbose=0, outputnorm=True, getVideo=False):
        if not self.mgmFlag:
            self.originaldamage = deepcopy(self.damage)
        for i in range(self.nRuns):
            self.InitializeNewRun()
            self.InitializeNewTracks(self.damage)
            self.CountCurrentDSB()
            if getVideo:
                self._storeImages()
            self.repairedList = []
            self.misrepairedlist = []
            self.chromosomeAberrations = []
            self.currentrun += 1
            self.messages.append('Cell ' + str(self.currentrun) + ' of ' + str(self.TotalRuns) + '...')
            if verbose > 0:
                print(self.messages[-1])
            if DSB in self.outputs:
                self.outDSB = output.TimeCurveForSingleRun('Remaining DSB')
                self.misrepDSB = output.TimeCurveForSingleRun('Misrepaired DSB')
                self.DSBEvolution()
            if FOCI in self.outputs:
                self.outFoci = output.TimeCurveForSingleRun('Foci')
                self.FociEvolution()
            while self.clock.CurrentTime != self.clock.FinalTime:
                if verbose > 1:
                    if not self.mgmFlag:
                        print("Time " + str(round(self.clock.CurrentTime/3600,2)) + " h - Dose: " + str(round(self.damage.cumulativeDose, 2)) + " Gy. Number of DSB: " + str(self.damage.numDSB))
                    else:
                        print("Time " + str(round(self.clock.CurrentTime/3600,2)) + " h. Number of DSB: " + str(self.numDSB))
                self.clock.AdvanceTimeStep()
                self.DoOneStep()
                if not self.mgmFlag:
                    self.InitializeNewTracks(self.damage)
                if getVideo:
                    self._storeImages()
            if DSB in self.outputs:
                self.runoutputDSB.AddTimeCurveForSingleRun(self.outDSB)
                self.runoutputMisrepairedDSB.AddTimeCurveForSingleRun(self.misrepDSB)
                #self.outDSB.Plot()
                #self.outDSB.WriteCSV()
            if SURVIVAL in self.outputs:
                if not self.CellActivated:
                    self.outputsurvival.append(self.DetermineCellFateForRun())
                else:
                    self.outputsurvival.append(self.cell.Surviving)
                if self.outputsurvival[-1] == 0:
                    print('Cell died')
                else:
                    print('Cell survived')
            if FOCI in self.outputs:
                self.runoutputFoci.AddTimeCurveForSingleRun(self.outFoci)
            self.messages.append('Repaired: ' + str(len(self.repairedList)) + ' - Misrepaired: ' + str(len(self.misrepairedlist)))
            if verbose > 0:
                print(self.messages[-1])
            self.clock.Reset()
            self.resetDamage()
        if DSB in self.outputs:
            self.runoutputDSB.DoStatistics(outputnorm)
            self.runoutputMisrepairedDSB.DoStatistics(False)
            if self.plotflag:
                self.runoutputDSB.Plot()
                self.runoutputDSB.WriteCSV()
                self.runoutputMisrepairedDSB.Plot()
                self.runoutputMisrepairedDSB.WriteCSV()
        if FOCI in self.outputs:
            self.runoutputFoci.DoStatistics(outputnorm)
            if self.plotflag:
                self.runoutputFoci.Plot()
                self.runoutputFoci.WriteCSV()
        if SURVIVAL in self.outputs and self.CellActivated:
            self.AddCellOutput()

    def DetermineCellFateForRun(self):
        # Get number of misrepaired DSB
        numMisrepDSB = int(self.misrepDSB.yvalues[-1])
        # Check if cell is dead
        for i in range(numMisrepDSB):
            if np.random.rand() < self.p_lethal_aberration:
                return 0
        # Get number of unresolved DSB at the checkpoints
        numDSBcheckpoint1 = int(self.outDSB.GetValueForTimePoint(self.time_to_G1S_checkpoint))
        for i in range(numDSBcheckpoint1):
            if np.random.rand() < self.p_apoptosis:
                return 0
        numDSBcheckpoint2 = int(self.outDSB.GetValueForTimePoint(self.time_to_G2M_checkpoint))
        for i in range(numDSBcheckpoint2):
            if np.random.rand() < self.p_apoptosis:
                return 0
        # The cell survived!
        return 1

    def GetSurvivalFraction(self):
        return np.sum(self.outputsurvival)/len(self.outputsurvival)

    def resetDamage(self):
        if not self.mgmFlag:
            self.damage = deepcopy(self.originaldamage)

    def DoOneStep(self):
        if self.DiffusionActivated:
            self.DoDiffusion()
        if self.DSBRepairActivated:
            self.DoDSBRepair()
        if self.SSBRepairActivated:
            self.DoSSBRepair()
        if self.BDRepairActivated:
            self.DoBDRepair()
        if self.CellActivated:
            self.DoCellCycle()
        self.UpdateDamageMaps()
        if DSB in self.outputs:
            self.DSBEvolution()
        if FOCI in self.outputs:
            self.FociEvolution()

    def DoDiffusion(self):
        for i, t in enumerate(self.betracks):
            newpos = self.diffusionModel.Diffuse(t, self.clock.CurrentTimeStep)
            while self._checkPosWithinNucleus(newpos) is False:
                newpos = self.diffusionModel.Diffuse(t, self.clock.CurrentTimeStep)
            newstep = tracking.BeStep(newpos, self.clock.CurrentTime, complexity=self.betracks[i].GetLastStep().Complexity)
            newstep.Status = self.betracks[i].GetLastStep().Status
            self.betracks[i].AddNewStep(newstep)

    def DoDSBRepair(self):
        if self.CellActivated:
            repair, self.foci, fociformed, misrepaired = self.dsbRepairModel.Repair(self.betracks, self.clock.PreviousTime, self.clock.CurrentTimeStep, self.foci, self.cell.Stage, self.cell.OxygenConcentration)
        else:
            repair, self.foci, fociformed, misrepaired = self.dsbRepairModel.Repair(self.betracks, self.clock.PreviousTime, self.clock.CurrentTimeStep, self.foci)
        if self.foci is None:
            try: self.outputs.remove(FOCI)
            except: pass
        for i in range(repair.shape[0]):
            for j in range(i+1, repair.shape[1]):
                if fociformed is not None and fociformed[i, j]:
                    self.betracks[i].IsFixed = True
                    self.betracks[j].IsFixed = True
                if repair[i, j] and misrepaired is None:
                    if self.betracks[i].DSBid == self.betracks[j].DSBid:
                        self.betracks[i].GetLastStep().Status = REPAIRED
                        self.betracks[j].GetLastStep().Status = REPAIRED
                        self.repairedList.append([self.betracks[i], self.betracks[j]])
                    else:
                        self.betracks[i].GetLastStep().Status = MISREPAIRED
                        self.betracks[j].GetLastStep().Status = MISREPAIRED
                        self.misrepairedlist.append([self.betracks[i], self.betracks[j]])
                elif repair[i, j] and misrepaired is not None and misrepaired[i, j]:
                    self.betracks[i].GetLastStep().Status = MISREPAIRED
                    self.betracks[j].GetLastStep().Status = MISREPAIRED
                    self.misrepairedlist.append([self.betracks[i], self.betracks[j]])
                    if self.betracks[i].DSBid != self.betracks[j].DSBid:
                        self.chromosomeAberrations.append([self.betracks[i], self.betracks[j]])
                elif repair[i, j] and misrepaired is not None and not misrepaired[i, j]:
                    self.betracks[i].GetLastStep().Status = REPAIRED
                    self.betracks[j].GetLastStep().Status = REPAIRED
                    self.repairedList.append([self.betracks[i], self.betracks[j]])

    def DoSSBRepair(self):
        for i in range(len(self.ssbdamages)):
            if self.ssbdamages[i].Status == DAMAGED:
                if self.ssbRepairModel.Repair(self.ssbdamages[i], self.clock.CurrentTimeStep):
                    self.ssbdamages[i].Status = REPAIRED

    def DoBDRepair(self):
        for i in range(len(self.bdamages)):
            if self.bdamages[i].Status == DAMAGED:
                if self.bdRepairModel.Repair(self.clock.CurrentTimeStep):
                    self.bdamages[i].Status = REPAIRED

    def DoCellCycle(self):
        if self.cell.Surviving:
            stage = self.cell.Stage
            fate = self.EvaluateCheckpoint()
            if fate == 'Stalled':
                inCheckpoint = self.cell.EvolveOverCycle(self.clock.CurrentTime, True)
            elif fate == 'Apoptosis' or fate == 'Necrosis' or fate == 'MitoticCatastrophe':
                inCheckpoint = self.cell.EvolveOverCycle(self.clock.CurrentTime, True)
            else:
                inCheckpoint = self.cell.EvolveOverCycle(self.clock.CurrentTime, False)
            # if in checkpoint, fate can be definite
            if inCheckpoint:
                if fate == 'Apoptosis':
                    self.cell.Surviving = False
                    self.cell.Apoptosis = True
                    self.cell.TimeOfDeath = self.clock.CurrentTime + self.cell.SampleTimeToApoptosis()
                    self.celloutput.CauseOfDeath = 'Apoptosis'
                elif fate == 'Necrosis':
                    self.cell.Surviving = False
                    self.cell.Necrosis = True
                    self.cell.TimeOfDeath = self.clock.CurrentTime + self.cell.SampleTimeToNecrosis()
                    self.celloutput.CauseOfDeath = 'Necrosis'
                elif fate == 'MitoticCatastrophe':
                    self.cell.Surviving = False
                    self.cell.MitoticCatastrophe = True
                    self.cell.TimeOfDeath = self.clock.CurrentTime + self.cell.SampleTimeToMitoticCatastrophe()
                    self.celloutput.CauseOfDeath = 'MitoticCatastrophe'
                elif fate == 'Stalled':
                    self.celloutput.Senescence = True
                    self.celloutput.TimeOfSenescence = self.clock.CurrentTime
                elif fate == 'Continue':
                    if self.cell.Stage == 1:
                        self.celloutput.MadeCheckpointG1S = True
                    elif self.cell.Stage == 3:
                        self.celloutput.MadeCheckpointG2M = True

            newstage = self.cell.Stage
            if stage != newstage:
                print('Cell stage: ' + str(stage) + ' -> ' + str(newstage) + ' at time ' + str(self.clock.CurrentTime/3600) + ' hours')
                # Update heterochromatin status
                for i in range(len(self.betracks)):
                    # Sample heterochromain vs euchromatin
                    if self.CellActivated:
                        p = self.cell.GetCurrentProportionOfHeterochromatin()
                        if np.random.rand() < p:
                            self.betracks[i].GetLastStep().IsHeterochromatin = True

    def EvaluateCheckpoint(self):
        # Check the remaining damage and assign probabilities to stall the cell cycle
        self.CountCurrentDSB()
        if self.cell.Stage == 1:
            return self.cell.CheckpointG1S(self.numDSB)
        if self.cell.Stage == 3:
            return self.cell.CheckpointG2M(self.numDSB, len(self.misrepairedlist), len(self.chromosomeAberrations))

    def CountCurrentDSB(self):
        self.numDSB = self.damage.numDSB
        # for be in self.betracks:
        #     if be.GetLastStep().Status == DAMAGED:
        #         self.numDSB += 0.5
        # self.numDSB = int(self.numDSB)

    def UpdateDamageMaps(self):
        if not self.mgmFlag:
            for be in self.betracks:
                if be.GetLastStep().Status == REPAIRED or be.GetLastStep().Status == MISREPAIRED:
                    iCh = be.ChromosomeID
                    iBp = be.BasePairID
                    iCo = be.StrandID + 1
                    self.damage.damageMap[iCh][iBp][iCo].type = 0
            for ssb in self.ssbdamages:
                if ssb.Status == REPAIRED:
                    iCh = ssb.ChromosomeID
                    iBp = ssb.BasePairID
                    iCo = ssb.StrandID + 1
                    self.damage.damageMap[iCh][iBp][iCo].type = 0
            for bd in self.bdamages:
                if bd.Status == REPAIRED:
                    iCh = bd.ChromosomeID
                    iBp = bd.BasePairID
                    iCo = bd.StrandID
                    if iCo == 2:
                        iCo = 4
                    self.damage.damageMap[iCh][iBp][iCo].type = 0
            self.damage.recomputeDamagesFromReadSites(stopAtTime=self.clock.CurrentTime, stopAtDose=self.maxDose)

    def DSBEvolution(self):
        self.CountCurrentDSB()
        self.outDSB.AddTimePoint(self.clock.CurrentTime, self.numDSB)
        self.misrepDSB.AddTimePoint(self.clock.CurrentTime, len(self.misrepairedlist))

    def FociEvolution(self):
        unresolvedFoci = 0
        for f in self.foci:
            if f.status != 'Resolved' and f.status != 'Notformed':
                unresolvedFoci += 1
        self.outFoci.AddTimePoint(self.clock.CurrentTime, unresolvedFoci)

    def AddCellOutput(self):
        self.celloutput.Surviving = self.cell.Surviving
        self.celloutput.DoseExposure = self.damage.cumulativeDose
        self.cellcultureoutput.AddCell(self.celloutput)

    def _storeImages(self):
        # Get the root temp directory location
        repo_root = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../..'))
        temp_dir = os.path.join(repo_root, 'temp')
        folder_path = os.path.join(temp_dir, 'images2D')
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        folder_path_3D = os.path.join(temp_dir, 'images3D')
        if not os.path.exists(folder_path_3D):
            os.makedirs(folder_path_3D)
        time = str(np.round(self.clock.CurrentTime / 3600, 2))
        name = os.path.join(folder_path, '2D_' + str(self.currentrun) + '_' + time + '.png')
        name_3D = os.path.join(folder_path_3D, '3D_' + str(self.currentrun) + '_' + time + '.png')
        self.damage.produce2DImages(saveFile=name, onlyz=True, title= time + ' hours - Dose: ' + str(np.round(self.damage.cumulativeDose, 2)) + ' Gy - DSB: ' + str(self.numDSB))
        self.damage.produce3DImage(show=False, saveFile=name_3D, title= time + ' hours - Dose: ' + str(np.round(self.damage.cumulativeDose, 2)) + ' Gy - DSB: ' + str(self.numDSB))

    def _checkPosWithinNucleus(self, pos):
        if self.nucleusMaxRadius is None:
            return True
        else:
            pos = np.array(pos)
            if np.sqrt(np.sum(np.power(pos, 2))) > self.nucleusMaxRadius:
                return False

    @property
    def CellActivated(self):
        if self._cellActivated is False:
            return self._cellActivated
        else:
            return True
    @CellActivated.setter
    def CellActivated(self, b):
        self._cellActivated = b


    @property
    def DiffusionActivated(self):
        if self._diffusionactivated is False:
            return self._diffusionactivated
        else:
            return True
    @DiffusionActivated.setter
    def DiffusionActivated(self, b):
        self._diffusionactivated = b

    @property
    def DSBRepairActivated(self):
        if self._dsbrepactivated is False:
            return self._dsbrepactivated
        else:
            return True
    @DSBRepairActivated.setter
    def DSBRepairActivated(self, b):
        self._dsbrepactivated = b

    @property
    def SSBRepairActivated(self):
        if self._ssbrepactivated is False:
            return self._ssbrepactivated
        else:
            return True
    @SSBRepairActivated.setter
    def SSBRepairActivated(self, b):
        self._ssbrepactivated = b

    @property
    def BDRepairActivated(self):
        if self._bdrepactivated is False:
            return self._bdrepactivated
        else:
            return True
    @BDRepairActivated.setter
    def BDRepairActivated(self, b):
        self._bdrepactivated = b

    @property
    def TotalRuns(self):
        try:
            return self._totalruns
        except:
            return self.nRuns
    @TotalRuns.setter
    def TotalRuns(self, t):
        self._totalruns = t


class Clock:
    def __init__(self, initialTime, finalTime, nSteps, listOfTimePoints = None):
        self.CurrentIndex = 0
        if listOfTimePoints is None:
            self.timepoints = np.linspace(initialTime, finalTime, nSteps)
        else:
            self.timepoints = listOfTimePoints

    def AdvanceTimeStep(self):
        self.CurrentIndex = self._currentindex + 1

    def Reset(self):
        self.CurrentIndex = 0

    @property
    def CurrentIndex(self):
        return self._currentindex
    @CurrentIndex.setter
    def CurrentIndex(self, i):
        self._currentindex = i

    @property
    def CurrentTime(self):
        return self.timepoints[self.CurrentIndex]

    @property
    def CurrentTimeStep(self):
        if self.CurrentIndex < len(self.timepoints) - 1:
            return self.timepoints[self.CurrentIndex + 1] - self.timepoints[self.CurrentIndex]
        else:
            return self.timepoints[-1] - self.timepoints[-2]

    @property
    def PreviousTime(self):
        return self.timepoints[self.CurrentIndex - 1]

    @property
    def InitialTime(self):
        return self.timepoints[0]

    @property
    def FinalTime(self):
        return self.timepoints[-1]

