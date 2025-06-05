#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 1/22/23 2:47 PM

@author: alejandrobertolet
"""

import numpy as np

class CellModel:
    def __init__(self, name, pars=None):
        self.Model = name
        self.SetParameters(pars)

    def SetParameters(self, pars=None):
        if self.Model == 'standard':
            # Cell cycle parameters
            if pars is None or 'cycling' not in pars: self.cycling = True
            else: self.cycling = pars['cycling']
            if pars is None or 'initialStage' not in pars: self.initialStage = 'G1'
            else: self.initialStage = pars['initialStage']
            if pars is None or 'mu_G1' not in pars: self.mu_G1 = 10.0
            else: self.mu_G1 = pars['mu_G1']
            if pars is None or 'sigma_G1' not in pars: self.sigma_G1 = 0.1 # Not necessary for exp distributions
            else: self.sigma_G1 = pars['sigma_G1']
            if pars is None or 'unit_G1' not in pars: self.unit_G1 = 'h'
            else: self.unit_G1 = pars['unit_G1']
            if pars is None or 'mu_S' not in pars: self.mu_S = 7.0
            else: self.mu_S = pars['mu_S']
            if pars is None or 'sigma_S' not in pars: self.sigma_S = 0.1 # Not necessary for exp distributions
            else: self.sigma_S = pars['sigma_S']
            if pars is None or 'unit_S' not in pars: self.unit_S = 'h'
            else: self.unit_S = pars['unit_S']
            if pars is None or 'mu_G2' not in pars: self.mu_G2 = 5.0
            else: self.mu_G2 = pars['mu_G2']
            if pars is None or 'sigma_G2' not in pars: self.sigma_G2 = 0.1 # Not necessary for exp distributions
            else: self.sigma_G2 = pars['sigma_G2']
            if pars is None or 'unit_G2' not in pars: self.unit_G2 = 'h'
            else: self.unit_G2 = pars['unit_G2']
            if pars is None or 'mu_M' not in pars: self.mu_M = 1.5
            else: self.mu_M = pars['mu_M']
            if pars is None or 'sigma_M' not in pars: self.sigma_M = 0.1 # Not necessary for exp distributions
            else: self.sigma_M = pars['sigma_M']
            if pars is None or 'unit_M' not in pars: self.unit_M = 'h'
            else: self.unit_M = pars['unit_M']
            if pars is None or 'propHeterochromatin_G1' not in pars: self.propHeterochromatin_G1 = 0.1
            else: self.propHeterochromatin_G1 = pars['propHeterochromatin_G1']
            if pars is None or 'propHeterochromatin_S' not in pars: self.propHeterochromatin_S = 0.25
            else: self.propHeterochromatin_S = pars['propHeterochromatin_S']
            if pars is None or 'propHeterochromatin_G2' not in pars: self.propHeterochromatin_G2 = 0.5
            else: self.propHeterochromatin_G2 = pars['propHeterochromatin_G2']
            if pars is None or 'propHeterochromatin_M' not in pars: self.propHeterochromatin_M = 0.9
            else: self.propHeterochromatin_M = pars['propHeterochromatin_M']
            # Damage at the checkpoints
            if pars is None or 'damageThresholdForNecrosisAtG1S' not in pars: self.damageThresholdForNecrosisAtG1S = 20
            else: self.damageThresholdForNecrosisAtG1S = pars['damageThresholdForNecrosisAtG1S']
            if pars is None or 'damageThresholdForApoptosisAtG1S' not in pars: self.damageThresholdForApoptosisAtG1S = 14
            else: self.damageThresholdForApoptosisAtG1S = pars['damageThresholdForApoptosisAtG1S']
            if pars is None or 'damageThresholdForArrestAtG1S' not in pars: self.damageThresholdForArrestAtG1S = 8
            else: self.damageThresholdForArrestAtG1S = pars['damageThresholdForArrestAtG1S']
            if pars is None or 'damageThresholdForNecrosisAtG2M' not in pars: self.damageThresholdForNecrosisAtG2M = 16
            else: self.damageThresholdForNecrosisAtG2M = pars['damageThresholdForNecrosisAtG2M']
            if pars is None or 'damageThresholdForApoptosisAtG2M' not in pars: self.damageThresholdForApoptosisAtG2M = 8
            else: self.damageThresholdForApoptosisAtG2M = pars['damageThresholdForApoptosisAtG2M']
            if pars is None or 'damageThresholdForArrestAtG2M' not in pars: self.damageThresholdForArrestAtG2M = 5
            else: self.damageThresholdForArrestAtG2M = pars['damageThresholdForArrestAtG2M']
            if pars is None or 'misrepairThresholdForMitoticCatastrophe' not in pars: self.misrepairThresholdForMitoticCatastrophe = 5
            else: self.misrepairThresholdForMitoticCatastrophe = pars['misrepairThresholdForMitoticCatastrophe']
            # Times for death times
            if pars is None or 'mu_necrosis' not in pars: self.mu_necrosis = 2.0
            else: self.mu_necrosis = pars['mu_necrosis']
            if pars is None or 'sigma_necrosis' not in pars: self.sigma_necrosis = 0.75 # Not necessary for exp distributions
            else: self.sigma_necrosis = pars['sigma_necrosis']
            if pars is None or 'unit_necrosis' not in pars: self.unit_necrosis = 'h'
            else: self.unit_necrosis = pars['unit_necrosis']
            if pars is None or 'mu_apoptosis' not in pars: self.mu_apoptosis = 3.0
            else: self.mu_apoptosis = pars['mu_apoptosis']
            if pars is None or 'sigma_apoptosis' not in pars: self.sigma_apoptosis = 0.5 # Not necessary for exp distributions
            else: self.sigma_apoptosis = pars['sigma_apoptosis']
            if pars is None or 'unit_apoptosis' not in pars: self.unit_apoptosis = 'h'
            else: self.unit_apoptosis = pars['unit_apoptosis']
            if pars is None or 'mu_mitoticcat' not in pars: self.mu_mitoticCatastrophe = 10.0
            else: self.mu_mitoticCatastrophe = pars['mu_mitoticcat']
            if pars is None or 'sigma_mitoticcat' not in pars: self.sigma_mitoticCatastrophe = 0.75 # Not necessary for exp distributions
            else: self.sigma_mitoticCatastrophe = pars['sigma_mitoticcat']
            if pars is None or 'unit_mitoticcat' not in pars: self.unit_mitoticCatastrophe = 'h'
            else: self.unit_mitoticCatastrophe = pars['unit_mitoticcat']
            if pars is None or 'oxygen_concentration' not in pars: self.oxygen_concentration = 1.0
            else: self.oxygen_concentration = pars['oxygen_concentration']

class DiffusionModel:
    def __init__(self, name, pars=None):
        self.Model = name
        self.SetParameters(pars)

    def SetParameters(self, pars=None):
        # pars is a dictionary with the parameters
        if self.Model == 'free' or self.Model == 'subdiffusion':
            if pars is None or 'D' not in pars.keys(): self.diffusionCoefficient = 5.0e-8
            else: self.diffusionCoefficient = pars['D']
            if pars is None or 'Dunits' not in pars.keys(): self.diffusionCoefficientUnits = 'um^2/s'
            else: self.diffusionCoefficientUnits = pars['Dunits']


class DSBRepairModel:
    def __init__(self, name, pars=None):
        self.Model = name
        self.SetParameters(pars)

    def SetParameters(self, pars=None):
        # pars is a dictionary with the parameters
        if self.Model == 'standard':
            if pars is None or 'NHEJ' not in pars.keys(): self.competentInNHEJ = True
            else: self.competentInNHEJ = pars['NHEJ']
            if pars is None or 'rNCNC' not in pars.keys(): self.repairRateNCNC = 2.0e-4
            else: self.repairRateNCNC = pars['rNCNC']
            if pars is None or 'rNCNCunits' not in pars.keys(): self.repairRateNCNCUnits = 'rep/s'
            else: self.repairRateNCNCUnits = pars['rNCNCunits']
            if pars is None or 'rComplex' not in pars.keys(): self.repairRateComplex = 1.0e-5
            else: self.repairRateComplex = pars['rComplex']
            if pars is None or 'rComplexunits' not in pars.keys(): self.repairRateComplexUnits = 'rep/s'
            else: self.repairRateComplexUnits = pars['rComplexunits']
            if pars is None or 'rMMEJ' not in pars.keys(): self.repairMMEJ = 2.361e-7
            else: self.repairMMEJ = pars['rMMEJ']
            if pars is None or 'rMMEJunits' not in pars.keys(): self.repairMMEJUnits = 'rep/s'
            else: self.repairMMEJUnits = pars['rMMEJunits']
            if pars is None or 'sigma' not in pars.keys(): self.sigmaDistance = 0.25
            else: self.sigmaDistance = pars['sigma']
            if pars is None or 'sigmaunits' not in pars.keys(): self.sigmaDistanceUnits = 'um'
            else: self.sigmaDistanceUnits = pars['sigmaunits']
            if pars is None or 'competentInNHEJ' not in pars.keys(): self.competentInNHEJ = True
            else: self.competentInNHEJ = pars['NHEJ']
        if self.Model == 'foci_nocycle':
            if pars is None or 'mu_form' not in pars.keys(): self.mu_form = 1.0
            else: self.mu_form = pars['mu_form']
            if pars is None or 'sigma_form' not in pars.keys(): self.sigma_form = 0.5 # Not necessary for exp distributions
            else: self.sigma_form = pars['sigma_form']
            if pars is None or 'unit_form' not in pars.keys(): self.unit_form = 'min'
            else: self.unit_form = pars['unit_form']
            if pars is None or 'mu_repair' not in pars.keys(): self.mu_repair = np.linspace(15, 600, 19)
            else: self.mu_repair = pars['mu_repair']
            if pars is None or 'sigma_repair' not in pars.keys(): self.sigma_repair = 0.5 # Not necessary for exp distributions
            else: self.sigma_repair = pars['sigma_repair']
            if pars is None or 'unit_repair' not in pars.keys(): self.unit_repair = 'min'
            else: self.unit_repair = pars['unit_repair']
            if pars is None or 'mu_resolve' not in pars.keys(): self.mu_resolve = 90.0
            else: self.mu_resolve = pars['mu_resolve']
            if pars is None or 'sigma_resolve' not in pars.keys(): self.sigma_resolve = 0.5 # Not necessary for exp distributions
            else: self.sigma_resolve = pars['sigma_resolve']
            if pars is None or 'unit_resolve' not in pars.keys(): self.unit_resolve = 'min'
            else: self.unit_resolve = pars['unit_resolve']
            if pars is None or 'proteinAvailability' not in pars.keys(): self.proteinAvailability = 1.0
            else: self.proteinAvailability = pars['proteinAvailability']
            if pars is None or 'sigma' not in pars.keys(): self.sigmaDistance = 0.
            else: self.sigmaDistance = pars['sigma']
            if pars is None or 'sigmaunits' not in pars.keys(): self.sigmaDistanceUnits = 'um'
            else: self.sigmaDistanceUnits = pars['sigmaunits']
        if self.Model == 'foci_cycle':
            if pars is None or 'sigma' not in pars.keys(): self.sigmaDistance = 0.25
            else: self.sigmaDistance = pars['sigma']
            if pars is None or 'sigmaunits' not in pars.keys(): self.sigmaDistanceUnits = 'um'
            else: self.sigmaDistanceUnits = pars['sigmaunits']
            self.NHEJ = NHEJ(pars)
            self.HR = HR(pars)
            self.MMEJ = MMEJ(pars)


class SSBRepairModel:
    def __init__(self, name, pars=None):
        self.Model = name
        self.SetParameters(pars)

    def SetParameters(self, pars=None):
        # pars is a dictionary with the parameters
        if self.Model == 'standard':
            if pars is None or 'rNC' not in pars.keys(): self.repairRateNoComplex = 1.774e-3
            else: self.repairRateNoComplex = pars['rNC']
            if pars is None or 'rC' not in pars.keys(): self.repairRateComplex = 2.247e-16
            else: self.repairRateComplex = pars['rC']
            if pars is None or 'rNCunits' not in pars.keys(): self.repairRateNoComplexUnits = 'rep/s'
            else: self.repairRateNoComplexUnits = pars['rNCunits']
            if pars is None or 'rCunits' not in pars.keys(): self.repairRateComplexUnits = 'rep/s'
            else: self.repairRateComplexUnits = pars['rCunits']

class BDRepairModel:
    def __init__(self, name, pars=None):
        self.Model = name
        self.SetParameters(pars)

    def SetParameters(self, pars=None):
        # pars is a dictionary with the parameters
        if self.Model == 'standard':
            if pars is None or 'r' not in pars.keys(): self.repairRate = 1.774e-3
            else: self.repairRate = pars['r']
            if pars is None or 'runits' not in pars.keys(): self.repairRateUnits = 'rep/s'
            else: self.repairRateUnits = pars['runits']

class RepairMechanism:
    def __init__(self, pars=None):
        if pars is None or 'prob_NHEJ_heterochromatin' not in pars.keys(): self.prob_NHEJ_heterochromatin = 0.6
        else: self.prob_NHEJ_heterochromatin = pars['prob_NHEJ_heterochromatin']
        if pars is None or 'prob_HR_heterochromatin' not in pars.keys(): self.prob_HR_heterochromatin = 0.3
        else: self.prob_HR_heterochromatin = pars['prob_HR_heterochromatin']
        if pars is None or 'prob_MMEJ_heterochromatin' not in pars.keys(): self.prob_MMEJ_heterochromatin = 0.8
        else: self.prob_MMEJ_heterochromatin = pars['prob_MMEJ_heterochromatin']

class NHEJ(RepairMechanism):
    def __init__(self, pars=None):
        super().__init__(pars)
        if pars is None or 'mu_repair_NHEJ' not in pars.keys(): self.mu_repair = np.linspace(10, 600, 19)  # Mu of the lognormal dist as a function of complexity
        else: self.mu_repair = np.linspace(pars['mu_repair_NHEJ_min'], pars['mu_repair_NHEJ_max'], 19)
        if pars is None or 'sigma_repair_NHEJ' not in pars.keys(): self.sigma_repair = 0.5  # Sigma of the lognormal dist as a function of complexity  # Not necessary for exp distributions
        else: self.sigma_repair = pars['sigma_repair_NHEJ']
        if pars is None or 'unit_repair_NHEJ' not in pars.keys(): self.unit_repair = 'min'
        else: self.unit_repair = pars['unit_repair_NHEJ']
        if pars is None or 'mu_form_NHEJ' not in pars.keys(): self.mu_form = 10.0
        else: self.mu_form = pars['mu_form_NHEJ']
        if pars is None or 'sigma_form_NHEJ' not in pars.keys(): self.sigma_form = 0.5 # Not necessary for exp distributions
        else: self.sigma_form = pars['sigma_form_NHEJ']
        if pars is None or 'unit_form_NHEJ' not in pars.keys(): self.unit_form = 'min'
        else: self.unit_form = pars['unit_form_NHEJ']
        if pars is None or 'mu_resolve_NHEJ' not in pars.keys(): self.mu_resolve = 90.0
        else: self.mu_resolve = pars['mu_resolve_NHEJ']
        if pars is None or 'sigma_resolve_NHEJ' not in pars.keys(): self.sigma_resolve = 0.5 # Not necessary for exp distributions
        else: self.sigma_resolve = pars['sigma_resolve_NHEJ']
        if pars is None or 'unit_resolve_NHEJ' not in pars.keys(): self.unit_resolve = 'min'
        else: self.unit_resolve = pars['unit_resolve_NHEJ']
        if pars is None or 'proteinAvailability_NHEJ' not in pars.keys(): self.proteinAvailability = 1.0
        else: self.proteinAvailability = pars['proteinAvailability_NHEJ']
        if pars is None or 'pMisrepair_NHEJ' not in pars.keys(): self.pMisrepair = 0.01
        else: self.pMisrepair = pars['pMisrepair_NHEJ']

class MMEJ(RepairMechanism):
    def __init__(self, pars=None):
        super().__init__(pars)
        if pars is None or 'mu_repair_MMEJ' not in pars.keys(): self.mu_repair = np.linspace(0.1, 90, 19)  # Mu of the lognormal dist as a function of complexity
        else: self.mu_repair = np.linspace(pars['mu_repair_MMEJ_min'], pars['mu_repair_MMEJ_max'], 19)
        if pars is None or 'sigma_repair_MMEJ' not in pars.keys(): self.sigma_repair = 0.5  # Sigma of the lognormal dist as a function of complexity  # Not necessary for exp distributions
        else: self.sigma_repair = pars['sigma_repair_MMEJ']
        if pars is None or 'unit_repair_MMEJ' not in pars.keys(): self.unit_repair = 'min'
        else: self.unit_repair = pars['unit_repair_MMEJ']
        if pars is None or 'mu_form_MMEJ' not in pars.keys(): self.mu_form = 1.0
        else: self.mu_form = pars['mu_form_MMEJ']
        if pars is None or 'sigma_form_MMEJ' not in pars.keys(): self.sigma_form = 0.5 # Not necessary for exp distributions
        else: self.sigma_form = pars['sigma_form_MMEJ']
        if pars is None or 'unit_form_MMEJ' not in pars.keys(): self.unit_form = 'min'
        else: self.unit_form = pars['unit_form_MMEJ']
        if pars is None or 'mu_resolve_MMEJ' not in pars.keys(): self.mu_resolve = 90.0
        else: self.mu_resolve = pars['mu_resolve_MMEJ']
        if pars is None or 'sigma_resolve_MMEJ' not in pars.keys(): self.sigma_resolve = 0.5 # Not necessary for exp distributions
        else: self.sigma_resolve = pars['sigma_resolve_MMEJ']
        if pars is None or 'unit_resolve_MMEJ' not in pars.keys(): self.unit_resolve = 'min'
        else: self.unit_resolve = pars['unit_resolve_MMEJ']
        if pars is None or 'proteinAvailability_MMEJ' not in pars.keys(): self.proteinAvailability = 1.0
        else: self.proteinAvailability = pars['proteinAvailability_MMEJ']
        if pars is None or 'pMisrepair_MMEJ' not in pars.keys(): self.pMisrepair = 0.1
        else: self.pMisrepair = pars['pMisrepair_MMEJ']

class HR(RepairMechanism):
    def __init__(self, pars=None):
        super().__init__(pars)
        if pars is None or 'mu_repair_HR' not in pars.keys(): self.mu_repair = np.linspace(50, 900, 19)  # Mu of the lognormal dist as a function of complexity
        else: self.mu_repair = np.linspace(pars['mu_repair_HR_min'], pars['mu_repair_HR_max'], 19)
        if pars is None or 'sigma_repair_HR' not in pars.keys(): self.sigma_repair = 0.5  # Sigma of the lognormal dist as a function of complexity  # Not necessary for exp distributions
        else: self.sigma_repair = pars['sigma_repair_HR']
        if pars is None or 'unit_repair_HR' not in pars.keys(): self.unit_repair = 'min'
        else: self.unit_repair = pars['unit_repair_HR']
        if pars is None or 'mu_form_HR' not in pars.keys(): self.mu_form = 20.0
        else: self.mu_form = pars['mu_form_HR']
        if pars is None or 'sigma_form_HR' not in pars.keys(): self.sigma_form = 0.5 # Not necessary for exp distributions
        else: self.sigma_form = pars['sigma_form_HR']
        if pars is None or 'unit_form_HR' not in pars.keys(): self.unit_form = 'min'
        else: self.unit_form = pars['unit_form_HR']
        if pars is None or 'mu_resolve_HR' not in pars.keys(): self.mu_resolve = 90.0
        else: self.mu_resolve = pars['mu_resolve_HR']
        if pars is None or 'sigma_resolve_HR' not in pars.keys(): self.sigma_resolve = 0.5 # Not necessary for exp distributions
        else: self.sigma_resolve = pars['sigma_resolve_HR']
        if pars is None or 'unit_resolve_HR' not in pars.keys(): self.unit_resolve = 'min'
        else: self.unit_resolve = pars['unit_resolve_HR']
        if pars is None or 'proteinAvailability_HR' not in pars.keys(): self.proteinAvailability = 1.0
        else: self.proteinAvailability = pars['proteinAvailability_HR']
        if pars is None or 'pMisrepair_HR' not in pars.keys(): self.pMisrepair = 0.001
        else: self.pMisrepair = pars['pMisrepair_HR']

