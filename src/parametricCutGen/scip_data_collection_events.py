from pyscipopt import Model, SCIP_EVENTTYPE, Eventhdlr
import logging

# when a primal dual gap is updated recorded when  a seperator is called.
# record gap, BinvArow, Binvbrow, Binvc

data_collection_logger = logging.getLogger(__name__)
data_collection_logger.setLevel(logging.DEBUG)
# for now, a test class to count number of sepa added and interupt after some have been added
# how to haldel inf gap?
class CutGapDataRecording(Eventhdlr):
    def __init__(self, model, max_number_of_sepa):
        """
        TESTS::
        >>> from parametricCutGen.optimal_cut_generation import OptimalCut
        >>> from parametricCutGen.scip_data_collection_events import CutGapDataRecording
        >>> from pyscipopt import Model
        >>> import logging; logging.disable()
        >>> model = Model()
        >>> sepa = OptimalCut()
        >>> model.includeSepa(sepa, "optimal_cut", "gmic equiv", priority=1000, freq=1)
        >>> data_record = CutGapDataRecording(model, 10)
        >>> model.includeEventhdlr(data_record, "record_gap_data", "Records dual gap data when optimal_cut is called")
        >>> model.readProblem("/home/acadia/Downloads/markshare_4_0.mps")
        >>> model.optimize()
        >>> data_record.get_data()
        """
        Eventhdlr.__init__(model)
        self.max_number_of_sepa = max_number_of_sepa
        self.number_of_sepa_added = 0
        self.gap_data = []

    def eventinit(self):
        self.model.catchEvent(SCIP_EVENTTYPE.ROWADDEDSEPA, self)

    def eventexit(self):
        data_collection_logger.debug(f"Here")
        self.model.dropEvent(SCIP_EVENTTYPE.GAPUPDATED, self)

    def eventexec(self, event):
        self.gap_data.append(self.model.getGap())
        self.number_of_sepa_added += 1
#        if self.number_of_sepa_added == self.max_number_of_sepa:
#            self.model.interruptSolve()

    def get_data(self, path=''):
        return self.gap_data, self.number_of_sepa_added 
