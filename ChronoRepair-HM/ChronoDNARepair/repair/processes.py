#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 11/6/22 3:00 PM

@author: alejandrobertolet
"""
import numpy as np
import scipy.stats as stats
from ChronoDNARepair.repair.models import DiffusionModel, DSBRepairModel, SSBRepairModel, BDRepairModel, CellModel
from ChronoDNARepair.repair.tracking import *

DAMAGED = 1
REPAIRED = 2
MISREPAIRED = 3

# Mechanisms
HR = 1
NHEJ = 2
MMEJ = 3

# Cell cycle stages
G0 = 0
G1 = 1
S = 2
G2 = 3
M = 4

class CellEvolution:
    def __init__(self, model='standard', pars=None):
        if model == 'standard':
            self.model = CellModel('standard', pars)
            self.Initialize()

    def Initialize(self):
        self.Stage = self.model.initialStage
        self.nextTransitionScheduledAt = self.SampleTimeForNextTransition()
        self.Stalled = False
        self.Apoptotic = False
        self.Necrotic = False
        self.MitoticCatastrophe = False
        self.Surviving = True
        self.OxygenConcentration = self.model.oxygen_concentration

    def GetCurrentProportionOfHeterochromatin(self):
        if self.Stage == 0:
            return self.model.propHeterochromatin_G1
        if self.Stage == 1:
            return self.model.propHeterochromatin_G1
        if self.Stage == 2:
            return self.model.propHeterochromatin_S
        if self.Stage == 3:
            return self.model.propHeterochromatin_G2
        if self.Stage == 4:
            return self.model.propHeterochromatin_M

    def SampleTimeForNextTransition(self):
        if self.Stage == 1:
            dist = TimeDistribution('exp', [self.model.mu_G1, self.model.unit_G1])
            return dist.SampleTime()
        if self.Stage == 2:
            dist = TimeDistribution('exp', [self.model.mu_S, self.model.unit_S])
            return dist.SampleTime()
        if self.Stage == 3:
            dist = TimeDistribution('exp', [self.model.mu_G2, self.model.unit_G2])
            return dist.SampleTime()
        if self.Stage == 4:
            dist = TimeDistribution('exp', [self.model.mu_M, self.model.unit_M])
            return dist.SampleTime()
        else:
            return np.inf

    def EvolveOverCycle(self, currenttime, stalled):
        self.Stalled = stalled
        if currenttime < self.nextTransitionScheduledAt:
            return 0
        elif self.Stalled:
            self.Stage = 0
            return 1
        else:
            if self.Stage == 4:
                self.Stage = 1
            else:
                self.Stage = self.Stage + 1
            self.nextTransitionScheduledAt = max(self.nextTransitionScheduledAt, currenttime) + self.SampleTimeForNextTransition()

    def CheckpointG1S(self, remainingdsb):
        if remainingdsb > self.model.damageThresholdForNecrosisAtG1S:
            return 'Necrosis'
        elif remainingdsb > self.model.damageThresholdForApoptosisAtG1S:
            return 'Apoptosis'
        elif remainingdsb > self.model.damageThresholdForArrestAtG1S:
            return 'Stalled'
        else:
            return 'Continue'

    def CheckpointG2M(self, remainingdsb, misrepaired, chromosomeAberrations):
        if chromosomeAberrations > 0:
            return 'Apoptosis'
        if remainingdsb > self.model.damageThresholdForNecrosisAtG2M:
            return 'Necrosis'
        elif remainingdsb > self.model.damageThresholdForApoptosisAtG2M:
            return 'Apoptosis'
        elif remainingdsb > self.model.damageThresholdForArrestAtG2M:
            return 'Stalled'
        elif misrepaired > self.model.misrepairThresholdForMitoticCatastrophe:
            return 'MitoticCatastrophe'
        else:
            return 'Continue'

    def SampleTimeToNecrosis(self):
        dist = TimeDistribution('exp', [self.model.mu_necrosis, self.model.unit_necrosis])
        return dist.SampleTime()

    def SampleTimeToApoptosis(self):
        dist = TimeDistribution('exp', [self.model.mu_apoptosis, self.model.unit_apoptosis])
        return dist.SampleTime()

    def SampleTimeToMitoticCatastrophe(self):
        dist = TimeDistribution('exp', [self.model.mu_mitoticCatastrophe, self.model.unit_mitoticCatastrophe])
        return dist.SampleTime()

    @property
    def Stage(self):
        return self._stage
    @Stage.setter
    def Stage(self, s):
        if s == 'G0' or s == 0:
            self._stage = 0
        if s == 'G1' or s == 1:
            self._stage = 1
        if s == 'S' or s == 2:
            self._stage = 2
        if s == 'G2' or s == 3:
            self._stage = 3
        if s == 'M' or s == 4:
            self._stage = 4

    @property
    def TimeOfDeath(self):
        return self._timeOfDeath
    @TimeOfDeath.setter
    def TimeOfDeath(self, t):
        self._timeOfDeath = t

    @property
    def DoseExposure(self):
        return self._doseExposure
    @DoseExposure.setter
    def DoseExposure(self, d):
        self._doseExposure = d

class Diffusion:
    def __init__(self, model='free', pars=None):
        if model == 'free':
            self.ActivateFreeDiffusion(pars)
        if model == 'subdiffusion':
            self.ActivateSubdiffusion(pars)

    def ActivateFreeDiffusion(self, pars):
        self.model = DiffusionModel('free', pars)

    def ActivateSubdiffusion(self, pars):
        self.model = DiffusionModel('subdiffusion', pars)

    def Diffuse(self, track, timestep):
        pos = np.array(track.GetLastStep().Position)
        if self.model.Model == 'free':
            Ns = np.random.normal(0, 1, 3)
            return pos + Ns * np.sqrt(self.model.diffusionCoefficient * timestep)
        elif self.model.Model == 'subdiffusion':
            if track.IsFixed:
                return pos
            else:
                pos = np.array(track.GetLastStep().Position)
                Ns = np.random.normal(0, 1, 3)
                return pos + Ns * np.sqrt(self.model.diffusionCoefficient * timestep)
        else:
            return pos

class TrackPairProcess:
    def __init__(self):
        pass

    @property
    def InteractionRadius(self):
        return self._intradius
    @InteractionRadius.setter
    def InteractionRadius(self, r):
        self._intradius = r

    def _getDistance(self, betrack1, betrack2):
        pos1 = betrack1.GetLastStep().Position
        pos2 = betrack2.GetLastStep().Position
        return np.sqrt(np.power(pos1[0] - pos2[0], 2) + np.power(pos1[1] - pos2[1], 2) + np.power(pos1[2] - pos2[2], 2))

class DSBRepair(TrackPairProcess):
    def __init__(self, model='standard', pars=None):
        if model == 'standard':
            self.ActivateDSBRepairStandard(pars)
        elif model == 'foci_nocycle':
            self.ActivateDSBRepairFoci(pars)
        elif model == 'foci_cycle':
            self.ActivateDSBRepairFociCycle(pars)

    def ActivateDSBRepairStandard(self, pars=None):
        self.model = DSBRepairModel('standard', pars)
        self.InteractionRadius = 1e8

    def ActivateDSBRepairFoci(self, pars=None):
        self.model = DSBRepairModel('foci_nocycle', pars)
        self.InteractionRadius = 1e8

    def ActivateDSBRepairFociCycle(self, pars=None):
        self.model = DSBRepairModel('foci_cycle', pars)
        self.InteractionRadius = 1e8

    def Repair(self, tracklist, currenttime, timestep, foci=None, currentcyclestage=None, oxygenperc=1.0):
        if self.model.Model == 'foci_cycle':
            return self.RepairThroughCycle(tracklist, currenttime, timestep, foci, currentcyclestage, oxygenperc)
        if self.model.Model == 'foci_nocycle':
            return self.RepairThroughFoci(tracklist, currenttime, timestep, foci)
        if self.model.Model == 'standard':
            return self.RepairStandard(tracklist, timestep)

    def RepairThroughCycle(self, tracklist, currenttime, timestep, focilist, currentcyclestage, oxygenperc=1.0):
        ntracks = len(tracklist)
        fociformmatrix = np.zeros([ntracks, ntracks])
        repairmechanisms = np.zeros([ntracks, ntracks])
        for i in range(ntracks):
            for j in range(i + 1, ntracks):
                if not tracklist[i].IsFixed and not tracklist[j].IsFixed:
                    # Pick the repair mechanism
                    mechanism = self.GetRepairMechanism(tracklist[i], tracklist[j], currentcyclestage)
                    if mechanism != 0:
                        fociformmatrix[i, j] = self.GetFociFormationProbability(tracklist[i], tracklist[j])
                        repairmechanisms[i, j] = mechanism
                else:
                    fociformmatrix[i, j] = 0.0
        rs = np.random.random([ntracks, ntracks])
        rs[fociformmatrix == 0.0] = 1.0
        fociformed = rs <= fociformmatrix
        # Avoid multiple uses of the same and: pick from maximum to minimum probability
        for i in range(ntracks):
            maxFori = 0.0
            repj = 0
            for j in range(i + 1, ntracks):
                if fociformed[i, j] and fociformmatrix[i, j] > maxFori:
                    maxFori = fociformmatrix[i, j]
                    repj = j
            if repj > 0:
                for k in range(ntracks):
                    fociformed[repj, k] = False
                for j in range(i + 1, ntracks):
                    if j != repj:
                        fociformed[i, j] = False
        # Sample the times of new foci formation
        for i in range(ntracks):
            for j in range(i + 1, ntracks):
                if fociformed[i, j]:
                    # Sample formation time
                    if repairmechanisms[i, j] == NHEJ:
                        dist = TimeDistribution('exp', [self.model.NHEJ.mu_form, self.model.NHEJ.unit_form])
                    elif repairmechanisms[i, j] == HR:
                        dist = TimeDistribution('exp', [self.model.HR.mu_form, self.model.HR.unit_form])
                    elif repairmechanisms[i, j] == MMEJ:
                        dist = TimeDistribution('exp', [self.model.MMEJ.mu_form, self.model.MMEJ.unit_form])
                    formationTime = currenttime + dist.SampleTime()
                    location = (tracklist[i].GetLastStep().Position + tracklist[j].GetLastStep().Position) / 2
                    complexity = int(np.round((tracklist[i].GetLastStep().Complexity + tracklist[j].GetLastStep().Complexity) / 2))
                    resolutionTime = -1
                    cellContext = currentcyclestage
                    tracks = [tracklist[i], tracklist[j]]
                    newFoci = RepairFoci(location, complexity, repairmechanisms[i, j], formationTime, resolutionTime, cellContext, tracks)
                    if formationTime <= currenttime + timestep:
                        newFoci.status = 'Ongoing'
                    focilist.append(newFoci)
        # For those foci formed, sample whether they are repaired
        repaired = np.zeros([ntracks, ntracks], dtype=bool)
        misrepaired = np.zeros([ntracks, ntracks], dtype=bool)
        for i, f in enumerate(focilist):
            if f.status == 'Notformed' and f.formationTime <= currenttime + timestep:
                f.status = 'Ongoing'
            if f.status == 'Ongoing':
                if f.repairPathway == NHEJ:
                    dist = TimeDistribution('exp', [self.model.NHEJ.mu_repair[min(f.complexity - 2, 18)], self.model.NHEJ.unit_repair])
                elif f.repairPathway == HR:
                    dist = TimeDistribution('exp', [self.model.HR.mu_repair[min(f.complexity - 2, 18)], self.model.HR.unit_repair])
                elif f.repairPathway == MMEJ:
                    dist = TimeDistribution('exp', [self.model.MMEJ.mu_repair[min(f.complexity - 2, 18)], self.model.MMEJ.unit_repair])
                prob = dist.SampleProbabilityForAGivenTime(timestep)
                if np.random.random() < prob:
                    indexi = tracklist.index(f.tracks[0])
                    indexj = tracklist.index(f.tracks[1])
                    repaired[indexi, indexj] = True
                    # Sample misrepair probability
                    if f.repairPathway == NHEJ:
                        pMisrepair = self.model.NHEJ.pMisrepair * f.complexity/2
                    elif f.repairPathway == HR:
                        pMisrepair = self.model.HR.pMisrepair * f.complexity/2
                    elif f.repairPathway == MMEJ:
                        pMisrepair = self.model.MMEJ.pMisrepair * f.complexity/2
                    if np.random.random() < pMisrepair:
                        misrepaired[indexi, indexj] = True
                    focilist[i].repairTime = currenttime
                    focilist[i].Repair()
                    # Sample a resolution time
                    if f.repairPathway == NHEJ:
                        dist = TimeDistribution('exp', [self.model.NHEJ.mu_resolve, self.model.NHEJ.unit_resolve])
                    elif f.repairPathway == HR:
                        dist = TimeDistribution('exp', [self.model.HR.mu_resolve, self.model.HR.unit_resolve])
                    elif f.repairPathway == MMEJ:
                        dist = TimeDistribution('exp', [self.model.MMEJ.mu_resolve, self.model.MMEJ.unit_resolve])
                    resolutionTime = currenttime + timestep + dist.SampleTime()
                    focilist[i].resolutionTime = resolutionTime
        # Resolve foci
        for i, f in enumerate(focilist):
            if f.status == 'Repaired' and f.resolutionTime <= currenttime + timestep:
                focilist[i].Resolve()
        return repaired, focilist, fociformed, misrepaired

    def GetRepairMechanism(self, track1, track2, stage):
        # Get heterochromatin status
        h = 0
        if track1.IsHeterochromatin: h += 1
        if track2.IsHeterochromatin: h += 1
        if h == 0: pAccessibility = 1
        if stage == G0 or stage == G1:
            # Try NHEJ
            pProteinAvailability = self.model.NHEJ.proteinAvailability
            if h == 2: pAccessibility = self.model.NHEJ.prob_NHEJ_heterochromatin
            if h == 1: pAccessibility = (1+self.model.NHEJ.prob_NHEJ_heterochromatin)/2
            p = pProteinAvailability * pAccessibility
            if np.random.random() < p:
                return NHEJ
            # Try MMEJ
            pProteinAvailability = self.model.MMEJ.proteinAvailability
            if h == 2: pAccessibility = self.model.MMEJ.prob_MMEJ_heterochromatin
            if h == 1: pAccessibility = (1+self.model.MMEJ.prob_MMEJ_heterochromatin)/2
            p = pProteinAvailability * pAccessibility
            if np.random.random() < p:
                return MMEJ
        if stage == S or stage == G2:
            if stage == S: pAvailableCopy = np.random.random()
            if stage == G2: pAvailableCopy = 1
            # Try HR
            pProteinAvailability = self.model.HR.proteinAvailability
            if h == 2: pAccessibility = self.model.HR.prob_HR_heterochromatin
            if h == 1: pAccessibility = (1+self.model.HR.prob_HR_heterochromatin)/2
            p = pProteinAvailability * pAccessibility * pAvailableCopy
            if np.random.random() < p:
                return HR
            # Try NHEJ
            pProteinAvailability = self.model.NHEJ.proteinAvailability
            if h == 2: pAccessibility = self.model.NHEJ.prob_NHEJ_heterochromatin
            if h == 1: pAccessibility = (1+self.model.NHEJ.prob_NHEJ_heterochromatin)/2
            p = pProteinAvailability * pAccessibility
            if np.random.random() < p:
                return NHEJ
            # Try MMEJ
            pProteinAvailability = self.model.MMEJ.proteinAvailability
            if h == 2: pAccessibility = self.model.MMEJ.prob_MMEJ_heterochromatin
            if h == 1: pAccessibility = (1+self.model.MMEJ.prob_MMEJ_heterochromatin)/2
            p = pProteinAvailability * pAccessibility
        if stage == 'M':
            return 0
        return 0

    def RepairThroughFoci(self, tracklist, currenttime, timestep, focilist):
        ntracks = len(tracklist)
        fociformmatrix = np.zeros([ntracks, ntracks])
        for i in range(ntracks):
            for j in range(i + 1, ntracks):
                if not tracklist[i].IsFixed and not tracklist[j].IsFixed:
                    fociformmatrix[i, j] = self.GetFociFormationProbability(tracklist[i], tracklist[j], self.model.proteinAvailability)
                else:
                    fociformmatrix[i, j] = 0.0
        rs = np.random.random([ntracks, ntracks])
        rs[fociformmatrix == 0.0] = 1.0
        fociformed = rs <= fociformmatrix
        # Avoid multiple uses of the same and: pick from maximum to minimum probability
        for i in range(ntracks):
            maxFori = 0.0
            repj = 0
            for j in range(i + 1, ntracks):
                if fociformed[i, j] and fociformmatrix[i, j] > maxFori:
                    maxFori = fociformmatrix[i, j]
                    repj = j
            if repj > 0:
                for k in range(ntracks):
                    fociformed[repj, k] = False
                for j in range(i + 1, ntracks):
                    if j != repj:
                        fociformed[i, j] = False
        # Sample the times of new foci formation
        for i in range(ntracks):
            for j in range(i + 1, ntracks):
                if fociformed[i, j]:
                    # Sample formation time
                    dist = TimeDistribution('exp', [self.model.mu_form, self.model.unit_form])
                    formationTime = currenttime + dist.SampleTime()
                    location = (tracklist[i].GetLastStep().Position + tracklist[j].GetLastStep().Position) / 2
                    complexity = int(np.round((tracklist[i].GetLastStep().Complexity + tracklist[j].GetLastStep().Complexity) / 2))
                    repairPathway = 'NHEJ'
                    resolutionTime = -1
                    cellContext = None
                    tracks = [tracklist[i], tracklist[j]]
                    newFoci = RepairFoci(location, complexity, repairPathway, formationTime, resolutionTime, cellContext, tracks)
                    focilist.append(newFoci)
        # For those foci formed, sample whether they are repaired
        repaired = np.zeros([ntracks, ntracks], dtype=bool)
        for i, f in enumerate(focilist):
            if f.status == 'Ongoing' and f.formationTime <= currenttime + timestep:
                dist = TimeDistribution('exp', [self.model.mu_repair[min(f.complexity - 2, 18)], self.model.unit_repair])
                prob = dist.SampleProbabilityForAGivenTime(timestep)
                if np.random.random() < prob:
                    indexi = tracklist.index(f.tracks[0])
                    indexj = tracklist.index(f.tracks[1])
                    repaired[indexi, indexj] = True
                    focilist[i].repairTime = currenttime
                    focilist[i].Repair()
                    # Sample a resolution time
                    dist = TimeDistribution('exp', [self.model.mu_resolve, self.model.unit_resolve])
                    resolutionTime = currenttime + timestep + dist.SampleTime()
                    focilist[i].resolutionTime = resolutionTime
        # Resolve foci
        for i, f in enumerate(focilist):
            if f.status == 'Repaired' and f.resolutionTime <= currenttime + timestep:
                focilist[i].Resolve()
        return repaired, focilist, fociformed, None

    def GetFociFormationProbability(self, track1, track2, proteinAvailability=1):
        distance = self._getDistance(track1, track2)
        if distance > self.InteractionRadius:
            return 0.0
        else:
            fd = np.exp(-distance ** 2 / (2 * self.model.sigmaDistance ** 2))
            return fd * proteinAvailability

    def RepairStandard(self, tracklist, timestep):
        ntracks = len(tracklist)
        probmatrix = np.zeros([ntracks, ntracks])
        for i in range(ntracks):
            for j in range(i+1, ntracks):
                if tracklist[i].GetLastStep().Status == DAMAGED and tracklist[j].GetLastStep().Status == DAMAGED:
                    probmatrix[i, j] = self.GetPairwiseProbabilityForStandardModel(tracklist[i], tracklist[j], timestep)
                else:
                    probmatrix[i, j] = 0.0
        rs = np.random.random([ntracks, ntracks])
        rs[probmatrix == 0.0] = 1.0
        repaired = rs <= probmatrix
        # Avoid multiple uses of the same and: pick from maximum to minimum probability
        for i in range(ntracks):
            maxFori = 0.0
            repj = 0
            for j in range(i+1, ntracks):
                if repaired[i, j] and probmatrix[i, j] > maxFori:
                    maxFori = probmatrix[i, j]
                    repj = j
            if repj > 0:
                for k in range(ntracks):
                    repaired[repj, k] = False
                for j in range(i+1, ntracks):
                    if j != repj:
                        repaired[i, j] = False
        return repaired, None, None, None

    def GetPairwiseProbabilityForStandardModel(self, betrack1, betrack2, timestep):
        distance = self._getDistance(betrack1, betrack2)
        if distance > self.InteractionRadius:
            return 0.0
        else:
            fd = np.exp(-distance**2 / (2*self.model.sigmaDistance**2))
            if self.model.competentInNHEJ:
                if betrack1.GetLastStep().Complexity <= 4 and betrack2.GetLastStep().Complexity <= 4:
                    return (1.0 - np.exp(-self.model.repairRateNCNC * timestep)) * fd
                else:
                    return (1.0 - np.exp(-self.model.repairRateComplex * timestep)) * fd
            else:
                return (1.0 - np.exp(-self.model.repairMMEJ * timestep)) * fd

    @property
    def CompetentInNHEJ(self):
        return self._nhej
    @CompetentInNHEJ.setter
    def CompetentInNHEJ(self, v):
        self._nhej = v

class SSBRepair:
    def __init__(self, model='standard', pars=None):
        if model == 'standard':
            self.ActivateSSBRepairStandard(pars)

    def ActivateSSBRepairStandard(self, pars):
        self.model = SSBRepairModel('standard', pars)

    def Repair(self, damtrack, timestep):
        if self.model.Model == 'standard':
            if damtrack.Complexity >= 2:
                prob = 1 - np.exp(-self.model.repairRateComplex * timestep)
            else:
                prob = 1 - np.exp(-self.model.repairRateNoComplex * timestep)
            r = np.random.random()
            return r <= prob

class BDRepair:
    def __init__(self, model='standard', pars=None):
        if model == 'standard':
            self.ActivateBDRepairStandard(pars)

    def ActivateBDRepairStandard(self, pars):
        self.model = BDRepairModel('standard', pars)

    def Repair(self, timestep):
        if self.model.Model == 'standard':
            prob = 1 - np.exp(-self.model.repairRate * timestep)
            r = np.random.random()
            return r <= prob

class TimeDistribution:
    def __init__(self, typeofdist='lognorm', pars=[0, 1, 'min']):
        if typeofdist == 'lognorm':
            mu = pars[0]
            sigma = pars[1]
            if len(pars) > 2:
                self.unit = pars[2]
            else:
                self.unit = 's'
            self.dist = stats.lognorm(sigma, scale=mu)
        if typeofdist == 'exp':
            mu = pars[0]
            if len(pars) > 1:
                self.unit = pars[1]
            else:
                self.unit = 's'
            self.dist = stats.expon(scale=mu)

    def SampleTime(self):
        if self.unit == 'min':
            return self.dist.rvs() * 60
        elif self.unit == 's':
            return self.dist.rvs()
        elif self.unit == 'h':
            return self.dist.rvs() * 3600

    def SampleProbabilityForAGivenTime(self, t):
        if self.unit == 'min':
            return self.dist.cdf(t / 60)
        elif self.unit == 's':
            return self.dist.cdf(t)
        elif self.unit == 'h':
            return self.dist.cdf(t / 3600)