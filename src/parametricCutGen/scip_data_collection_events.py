from pyscipopt import Model, SCIP_EVENTTYPE
import logging

# when a primal dual gap is updated recorded when  a seperator is called.
# record gap, BinvArow, Binvbrow, Binvc
class CutGapDataRecording(Eventhdlr):
    def __init__(self, model, cut_name, number_of_cuts):
        Eventhdlr.__init__(model)
        self.number_of_cuts = number_of_cuts
        self.gap_data = {cut_name}

    def eventinit(self):
        self.model.catchEvent(SCIP_EVENTTYPE.ROWADDEDSEPA, self)
        self.model.catchEvent(SCIP_EVENTTYPE.GAPUPDATED, self)

    # def eventexit(self):
        # self.model.dropEvent(SCIP_EVENTTYPE.ROWDELETEDSEPA , self)

    def eventexec(self, event):
        name = self.model.Event.getName()
        if name not in self.gap_data.keys():
            self.gap_data[name] = [abs(self.model. self.model.getPrimalbound()- self.model.getDualbound())/abs(self.model.getPrimalbound())]
        else:
            self.gap_data[name].append(abs(self.model. self.model.getPrimalbound()- self.model.getDualbound())/abs(self.model.getPrimalbound()))
        if self.gap_data[self.tracked_cut_name] == self.number_of_cuts:
            model.interruptSolve()

    def write_data(self, path):
        return self.gap_data
