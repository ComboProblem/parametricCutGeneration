from pyscipopt import Model, SCIP_EVENTTYPE, SCIP_RESULT, Eventhdlr
import logging

# Data collection moved to optimal cut generation.

data_collection_logger = logging.getLogger(__name__)
data_collection_logger.setLevel(logging.DEBUG)
# for now, a test class to count number of sepa added and interupt after some have been added
# how to haldel inf gap?
# type s of data; how long to solve; ect. 
#class completedOptimalCuts(Eventhdlr):
#    def __init__(self, model, max_number_of_optimal_cut, optimalCutSepa):
#        """
#        TESTS::
#        >>> from parametricCutGen.scip_data_collection_events import CutGapDataRecording
#        >>> from pyscipopt import Model
#        >>> import logging; logging.disable()
#        >>> model = Model()
#        >>> data_record = CutGapDataRecording(model, 100)
#        >>> model.includeEventhdlr(data_record, "record_gap_data", "Records dual gap data when optimal_cut is called")
#        """
#        Eventhdlr.__init__(model)
#        self.optimalCutSepa = optimalCutSepa
#
#    def eventinit(self):
#        # do a check at every event
#    def eventexec(self, event):
#        # stop if our 
#        if self.optimalCutSepa.get_nseperating_cuts() == self.max_number_of_optimal_cut:
#            self.model.interruptSolve()
#
