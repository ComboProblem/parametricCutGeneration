from pyscipopt import Model, SCIP_EVENTTYPE, SCIP_RESULT, Eventhdlr, SCIP_PARAMSETTING
import logging

# Data collection moved to optimal cut generation.

data_collection_logger = logging.getLogger(__name__)
data_collection_logger.setLevel(logging.DEBUG)
# for now, a test class to count number of sepa added and interupt after some have been added
# how to haldel inf gap?
# type s of data; how long to solve; ect. 

# get all nodes which 
class CollectGapNodes(Eventhdlr):
    def __init__(self, model, max_number_of_gap_updates,  cut_name="", write_node_data=False, write_file_name_style=None, write_file_dir=""):
        """
        TESTS::
        >>> from parametricCutGen.scip_data_collection_events import CutGapDataRecording
        >>> from pyscipopt import Model
        >>> import logging; logging.disable()
        >>> model = Model()
        >>> data_record = CollectGapNodes(model, 10)
        >>> model.includeEventhdlr(data_record, "record_gap_data", "Records dual gap data when optimal_cut is called")
        """
        pass

class checkCutsAdded(Eventhdlr):
    def __init__(self, model):
        Eventhdlr.__init__(model)
        self.count = 0 
        
    def eventinit(self):
        self.model.catchEvent(SCIP_EVENTTYPE.LPSOLVED, self)

    def eventexec(self, event):
        data_collection_logger.debug(f"event: {event.getName()}")       
        if self.model.getCurrentNode().getDepth() > 1:
            self.count += 1
            self.model.writeMIP(f"node_lp_{self.count}.lp")
        if self.count == 10:
            self.model.interruptSolve()
#


class disableCuts(Eventhdlr):
    def __init__(self, model):
        self.model = model
    def eventinit(self):
        self.model.catchEvent(SCIP_EVENTTYPE.NODEFOCUSED, self) # the event is called whenever a node is about to be solved

    def eventexec(self, event):
        if self.model.getCurrentNode().getDepth() > 1: # if we aren't in the root node
            self.model.setSeparating(SCIP_PARAMSETTING.OFF) # disable separators
        else:
            self.model.setSeparating(SCIP_PARAMSETTING.DEFAULT)
            self.model.readParams("./src/Experiments/paramFiles/scip_disable_other_cuts.set") # Disable scip avilaible 

class GapData(Eventhdlr):
    def __init__(self, model):
        Eventhdlr.__init__(model)
        self.lp_count = 0
        self.gap_update_count = 0
        self.ndual_bound_changes = 0
        self.measured_depths = [0, 1, 2, 3, 4, 5]
        
    def eventinit(self):
        # self.model.catchEvent(SCIP_EVENTTYPE.GAPUPDATED, self)
        self.model.catchEvent(SCIP_EVENTTYPE.LPSOLVED, self)
        # self.model.catchEvent(SCIP_EVENTTYPE.DUALBOUNDIMPROVED, self)
        # self.model.catchRowEvent(row, eventtype, eventhdlr)
    def eventexec(self, event):
        # LP solution
        data_collection_logger.debug(f"event: {event.getName()}")
        if event.getName() is "DUALBOUNDIMPROVED":
            data_collection_logger.debug(event.getOldBound())
            data_collection_logger.debug(event.getNewBound())
            self.ndual_bound_changes += 1
            # self.model.constructLP()
            # self.model.writeLP(f"update_lp_{self.gap_update_count}.lp") # node record            
        if event.getName() is "LPSOLVED":
            data_collection_logger.debug(f"number of sepa rounds: {self.model.getNCutsApplied()}")
            if self.model.getNCutsApplied() == 1:
                node = event.getNode()
#            model_node = self.model.getCurrentNode()
#            data_collection_logger.debug(f"event node:{node.getNumber()}, depth: {node.getDepth()}, added cons: {node.getNAddedConss()}")
#            for con in node.getAddedConss():
#                data_collection_logger.debug(f"con handlr: {con.getConshdlrName()}")
#                data_collection_logger.debug(con.isStickingAtNode())
#                data_collection_logger.debug(f"con dual sol: {self.model.getDualsolLinear(con)}")
#            if node.getDepth() in self.measured_depths:
                self.lp_count += 1
                self.model.constructLP()
                self.model.writeMIP(f"node_lp_{self.lp_count}.lp") # node record
                for con in node.getAddedConss():
                    data_collection_logger.debug(f"con handlr: {con.getConshdlrName()}")
                    data_collection_logger.debug(con.isStickingAtNode())
                    data_collection_logger.debug(f"con dual sol: {self.model.getDualsolLinear(con)}")
    # gap updated
#        if event.getName is "GAPUPDATED":
#            data_collection_logger.debug(event.getOldBound())
#            data_collection_logger.debug(event.getNewBound())
#            self.gap_update_count += 1
#            self.model.writeLP(f"update_lp_{self.gap_update_count}.lp") # node record
        if self.gap_update_count > 10 or self.lp_count > 10 or self.ndual_bound_changes>100 :
            data_collection_logger.debug(f"counts: {self.gap_update_count}, {self.lp_count}, {self.ndual_bound_changes}")
            self.model.interruptSolve()


