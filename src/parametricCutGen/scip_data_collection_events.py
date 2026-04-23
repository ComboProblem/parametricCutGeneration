from pyscipopt import Model, SCIP_EVENTTYPE, Eventhdlr
import logging

# when a primal dual gap is updated recorded when  a seperator is called.
# record gap, BinvArow, Binvbrow, Binvc
class CutGapDataRecording(Eventhdlr):
    def __init__(self, model, cut_name, number_of_cuts):
        """
        >>> from parametricCutGen.optimal_cut_generation import OptimalCut
        >>> from parametricCutGen.scip_data_collection_events import CutGapDataRecording
        >>> from pyscipopt import Model
        >>> import logging; logging.disable()
        >>> model = Model()
        >>> sepa = OptimalCut()
        >>> model.includeSepa(sepa, "optimal_cut", "gmic equiv", priority=1000, freq=1)
        >>> data_record = CutGapDataRecording(model, "optimal_cut", 10)
        >>> model.includeEventhdlr(data_record, "record_gap_data", "Records dual gap data when optimal_cut is called")
        >>> model.readProblem("/home/acadia/Downloads/flugpl.mps")
        >>> model.optimize()
        """
        Eventhdlr.__init__(model)
        self.number_of_cuts = number_of_cuts
        self._cut_name = cut_name
        self.gap_data = {cut_name:[]}

    def eventinit(self):
        self.model.catchEvent(SCIP_EVENTTYPE.ROWADDEDSEPA, self)

    def eventexit(self):
        self.model.dropEvent(SCIP_EVENTTYPE.ROWDELETEDSEPA , self)

    def eventexec(self, event):
        name = event.getName()
        try:
            self.gap_data[name].append(abs(self.model.getPrimalbound()- self.model.getDualbound())/abs(self.model.getPrimalbound()))
        except KeyError:
            self.gap_data[name] = [abs(self.model.getPrimalbound()- self.model.getDualbound())/abs(self.model.getPrimalbound())]
        if len(self.gap_data[self._cut_name]) == self.number_of_cuts:
            model.interruptSolve()

    def write_data(self, path):
        return self.gap_data
