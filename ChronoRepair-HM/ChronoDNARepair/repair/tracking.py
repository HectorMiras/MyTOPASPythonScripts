#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 11/6/22 2:30 PM

@author: alejandrobertolet
"""

class BeStep:
    def __init__(self, initialPosition, initialTime = 0, complexity = 0):
        self.Position = initialPosition
        self.Time = initialTime
        self.Complexity = complexity
        self.Status = 'Damaged'
        self.RepairProteinStatus = 'NotPresent'
        self._isHeterochromatin = False

    @property
    def Time(self):
        return self._time
    @Time.setter
    def Time(self, t):
        self._time = t

    @property
    def Position(self):
        return self._position
    @Position.setter
    def Position(self, p):
        self._position = p

    @property
    def Complexity(self):
        return self._complexity
    @Complexity.setter
    def Complexity(self, c):
        self._complexity = c

    @property
    def Status(self):
        return self._status
    @Status.setter
    def Status(self, s):
        if s == 'Damaged' or s == 1:
            self._status = 1
        if s == 'Repaired' or s == 2:
            self._status = 2
        if s == 'Misrepaired' or s == 3:
            self._status = 3

    @property
    def IsHeterochromatin(self):
        return self._isHeterochromatin

    @IsHeterochromatin.setter
    def IsHeterochromatin(self, h):
        self._isHeterochromatin = h

class DamageTrack:
    def __init__(self, trackid = -1):
        self.TrackID = trackid
        self.Status = 'Damaged'

    def GetTrackID(self):
        return self.TrackID

    @property
    def TrackID(self):
        return self._trackid

    @TrackID.setter
    def TrackID(self, id):
        self._trackid = id

    @property
    def ChromosomeID(self):
        return self._chrId

    @ChromosomeID.setter
    def ChromosomeID(self, c):
        self._chrId = c

    @property
    def BasePairID(self):
        return self._bpId

    @BasePairID.setter
    def BasePairID(self, b):
        self._bpId = b

    @property
    def StrandID(self):
        return self._strandid

    @StrandID.setter
    def StrandID(self, s):
        self._strandid = s

    @property
    def Time(self):
        return self._time
    @Time.setter
    def Time(self, t):
        self._time = t

    @property
    def Position(self):
        return self._position
    @Position.setter
    def Position(self, p):
        self._position = p

    @property
    def Complexity(self):
        return self._complexity
    @Complexity.setter
    def Complexity(self, c):
        self._complexity = c

    @property
    def Status(self):
        return self._status
    @Status.setter
    def Status(self, s):
        if s == 'Damaged' or s == 1:
            self._status = 1
        if s == 'Repaired' or s == 2:
            self._status = 2
        if s == 'Misrepaired' or s == 3:
            self._status = 3

class BeTrack(DamageTrack):
    def __init__(self, trackid = -1, dsbid = 0):
        super().__init__(trackid)
        self.DSBid = dsbid
        self.Steps = []
        self.inFoci = False

    def AddNewStep(self, step):
        self.Steps.append(step)

    def GetLastStep(self):
        if len(self.Steps) > 0:
            return self.Steps[-1]
        else:
            print("No steps were found for this break end track.")
            return -1

    def GetStepAtIndex(self, i):
        if i >=0 and i < len(self.Steps):
            return self.Steps[i]
        else:
            print("Error at accessing at step ", str(i), " in this break end track.")
            return -1

    @property
    def IsFixed(self):
        if self.inFoci:
            return True
        else:
            return False
    @IsFixed.setter
    def IsFixed(self, f):
        self.inFoci = f

    @property
    def OriginTime(self):
        if len(self.Steps) > 0:
            return self.Steps[0].Time
        else:
            print ("Error at getting origin time, this break end has no steps.")
            return -1

    @property
    def IsRepaired(self):
        if len(self.Steps) > 0:
            for s in self.Steps:
                if s.Status == 2:
                    return True
            return False
        else:
            return -1

    @property
    def IsMisrepaired(self):
        if len(self.Steps) > 0:
            for s in self.Steps:
                if s.Status == 3:
                    return True
            return False
        else:
            return -1

    @property
    def RepairTime(self):
        if len(self.Steps) > 0 and self.IsRepaired:
            for s in self.Steps:
                if s.Status == 2:
                    return s.Time
        else:
            return -1

    @property
    def MisrepairTime(self):
        if len(self.Steps) > 0 and self.IsMisrepaired:
            for s in self.Steps:
                if s.Status == 3:
                    return s.Time
        else:
            return -1

    @property
    def IsHeterochromatin(self):
        if len(self.Steps) > 0:
            return self.Steps[-1].IsHeterochromatin
        else:
            return -1

class RepairFoci:
    def __init__(self, location, complexity, repairPathway, formationTime, resolutionTime, cellContext, tracks):
        self.location = location
        self.complexity = complexity
        self.repairPathway = repairPathway
        self.formationTime = formationTime
        self.resolutionTime = resolutionTime
        self.cellContext = cellContext
        self.tracks = tracks
        self.proteins = []
        self.status = 'Notformed'

    def Resolve(self):
        self.status = 'Resolved'

    def Stall(self):
        self.status = 'Stalled'

    def Repair(self):
        self.status = 'Repaired'

    def AddGammaH2AX(self):
        self.proteins.append('GammaH2AX')

    def Add53BP1(self):
        self.proteins.append('53BP1')

    @property
    def Time(self):
        return self._time

    @Time.setter
    def Time(self, t):
        self._time = t

    @property
    def RepairTime(self):
        return self._originTime

    @RepairTime.setter
    def RepairTime(self, t):
        self._originTime = t

