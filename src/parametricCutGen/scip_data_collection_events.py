from pyscipopt import Model, SCIP_EVENTTYPE, SCIP_RESULT, Eventhdlr, SCIP_PARAMSETTING
import os
import logging

data_collection_logger = logging.getLogger(__name__)
data_collection_logger.setLevel(logging.DEBUG)

def record_data(model: Model, snaphot_frequency=10000,  write_path=""):
    """
    Attaches an event handler to a given SCIP model that collects primal and dual solutions,
    along with the solving time when they were found.
    The data is saved in model.data["primal_log"] and model.data["dual_log"]. They consist of
    a list of tuples, each tuple containing the solving time and the corresponding solution.

    A usage example can be found in examples/finished/plot_primal_dual_evolution.py. The
    example takes the information provided by this recipe and uses it to plot the evolution
    of the dual and primal bounds over time. 
    """
    class GapEventhdlr(Eventhdlr):
        def eventinit(self): # we want to collect best primal solutions and best dual solutions
            self.model.catchEvent(SCIP_EVENTTYPE.BESTSOLFOUND, self)
            self.model.catchEvent(SCIP_EVENTTYPE.DUALBOUNDIMPROVED, self)

        def eventexec(self, event):
            # if a new best primal solution was found, we save when it was found and also its objective
            if event.getType() == SCIP_EVENTTYPE.BESTSOLFOUND:
                self.model.data["primal_log"].append([self.model.getSolvingTime(), self.model.getPrimalbound()])
            
            if event.getType() == SCIP_EVENTTYPE.DUALBOUNDIMPROVED:
                self.model.data["dual_log"].append([self.model.getSolvingTime(), self.model.getDualbound()])

    class SnapShothdlr(Eventhdlr):
        def __init__(self, snaphot_frequency, write_path):
            self.snaphot_frequency = snaphot_frequency
            self.write_path = write_path
            self.count = 0
        def eventinit(self):
            self.model.catchEvent(SCIP_EVENTTYPE.LPSOLVED, self)

        def eventexec(self, event):
            if self.count % self.snaphot_frequency == 0:
                self.model.writeMIP(os.path.join(self.write_path, f"node_lp_{self.count}.lp"))
                self.model.data["snaphot_count"].append(self.count)
                self.model.data["snaphot_n_cuts"].append(self.model.getNCutsApplied())
                self.model.data["snaphot_node_depth"].append(self.model.getCurrentNode().getDepth())
            self.count += 1
            

             

    if not hasattr(model, "data") or model.data==None:
        model.data = {}

    model.data.update({
            'primal_log': [],
            'dual_log': [],
            'snaphot_count': [],
            'snaphot_node_depth': [],
            'snaphot_n_cuts': []
            })

    gap_hdlr = GapEventhdlr()
    snapshot_hdlr = SnapShothdlr(snaphot_frequency=snaphot_frequency, write_path=write_path)
    model.includeEventhdlr(gap_hdlr, "gapEventHandler", "Event handler which collects primal and dual solution evolution")
    model.includeEventhdlr(snapshot_hdlr, "snapShothdlr", "Event handler which takes LP snapshots and verifies experimental paramteters")
    return model

