from pyscipopt import Model, SCIP_EVENTTYPE
import logging

# when a primal dual gap is updated recorded when  a seperator is called.
# record gap, BinvArow, Binvbrow, Binvc
class CutGapDataRecording(Eventhdlr):
    def __init__(self, model, cut_name, number_of_cuts):
        Eventhdlr.__init__(model)
        self.gap_data = {}

    def eventinit(self):
        self.model.catchEvent(SCIP_EVENTTYPE.ROWADDEDSEPA, self)
        self.model.catchEvent(SCIP_EVENTTYPE.GAPUPDATED, self)

    # def eventexit(self):
        # self.model.dropEvent(SCIP_EVENTTYPE.ROWDELETEDSEPA , self)

    def eventexec(self, event):
        name = self.model.Event.getName()
        if name not in self.gap_data.keys():
            #
            self.gap_data[name] = [abs(self.model. self.model.getPrimalbound()- self.model.getDualbound())/abs(self.model.getPrimalbound())]
        else:
            self.gap_data[name].append(abs(self.model. self.model.getPrimalbound()- self.model.getDualbound())/abs(self.model.getPrimalbound()))
        if cut_name in self. gap_data.keys():
            if len(self.gap_data[cut_name]) == number_of_cuts:
                model.interruptSolve() # Once we've collected enough data, exit.
