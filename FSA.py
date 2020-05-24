#!/usr/bin/env python
"""
NOTE: I'm not sure if it's faster line-by-line
or as a chunked file. This implementation is line-by-line

Given a BLAST output, extract the necessary lines from it.

The objective is:

key: query
value: [PDBid, chain, 
        (start_chain, stop_chain), 
        (start_alignseq, stop_alignseq), 
        seqID_proteome]


NOTE - ignore "\n" only commands
Query-> Sbject name = Maybe you want the length?

TODO:
Transitions that require rules:
1. Start-> Query; Grab the Query name
2. Query-> SbjectName; Parse the subject name
3. SbjectName-AlignQuery; Get the query start
4. AlignQuery-> AlignSbject; Get Query end, Get the sbject start
5. AlignSbject -> SbjectName; get the sbject end, get sbject name
6. AlignSbject -> Query; get sbject end, get query end


"""
import gzip
import pickle as pkl
import time

class fsm:
    """
    General Finite state machine (FSM) architecture.
    
    This will be initialized with per-line processing rules.
    """

    def __init__(self, states=[]):
        self._states=states # List of states possible
        self.oldState = None
        self.currentState = None #Initialize the current state
        self.running = False # Initialize if you've started the FSM

    def start(self,startState=None):
        """ Start the finite state machine """
        if not startState or not (startState in [x[0] for x in self._states]):
            raise ValueError(startState, " not a valid starting state.")

        print("FSM initialized with rules added.", flush=True)
        self.oldState = startState
        self.currentState = startState
        self.running = True 

    def stop(self):
        """ Terminate the FSM """
        self.currentState = None
        self.running = False
        print("FSM has completed", flush=True)

    def addRule(self, fromState, toState, condition, operation=None):
        """ 
        (Optional to add more rules)
        Transition rules between states.
        Provide the following:

        1. Name of the starting state "fromState"
        2. Name of the ending state "toState"
        3. Condition - some function (ideally anonymous) to evaluate
        4. Perform an action (operation)
        """
        if self.running:
            raise ValueError("The FSM has been started. Cannot add rules.")

        # add a transition to the state table
        self._states.append( (fromState, 
                              toState, 
                              condition, 
                              operation))

    def event(self, value):
        """ 
        Evaluate the transition
        """
        if not self.currentState:
            raise ValueError("StateMachine not Started - cannot process event")

        # Get all the valid next states, given the current      
        self.nextStates = list(filter(lambda x: x[0]==self.currentState and 
                                        (x[2]==True or (callable(x[2]) and x[2](value))), 
                                        self._states))

        if not self.nextStates: 
            raise ValueError("No Transition defined from state {0} with value '{1}'".format(self.currentState, value))
        
        elif len(self.nextStates) > 1:
            raise ValueError("Too many transitions {0} with value '{1}' ->  New states defined {2}".format(self.currentState, value, [x[0] for x in self.nextStates]))
        
        else:
            if len(self.nextStates[0]) == 4:
                currstate, nextstate, condition, operation = self.nextStates[0]
            else:
                currstate, nextstate, condition = self.nextStates[0]
                callback = None

            # Check if the current state has changed to a new one
            self.oldState = self.currentState
            self.currentState, changed = (nextstate, True) if self.currentState != nextstate else (nextstate, False)
            
            # Execute the callback if defined
            if callable(operation):
                ret = operation(value)
            else:
                ret = None

            return self.currentState, changed, ret

    def CurrentState(self):
        """ 
        Return the current state of FSM
        """
        return self.currentState


class ProcessBLAST:
    """
    Uses the FSM model to process the BLASTP results
    """
    def __init__(self, fname):
        """ Initializes the filename for the model """
        self.fname = fname

    def init_FSM(self, transitions):
        """Initialize the FSM model"""
        model = fsm(states=[])

        for rule in transitions:
            if len(rule) < 4:
                model.addRule(rule[0], rule[1], rule[2])
            else:
                model.addRule(rule[0], rule[1], rule[2], rule[3])

        self.model = model


    def run_model(self, reportlines=1000000):
        """ Begin the FSM calculation """

        if self.model.currentState is None:
            self.model.start("Start")
        else:
            raise ValueError("Reinitialize the model with self.init_FSM")

        # Model initializations
        self.queryname = None
        self.sname = None
        self.alignQ = []
        self.alignS = []

        proteins = {}
        countr = 0
        print("Beginning model", flush=True)
        with open(self.fname, 'r') as f:

            for line in f:
                parsed = self.model.event(line)
                
                countr += 1
                if (countr % reportlines == 0):
                    print("Processed ", countr, "lines", flush=True)
                    
                # Condition 1: Query Name
                if (parsed[0] == "Query") and parsed[1]:
                    self.queryname = parsed[2]
                    proteins.update({self.queryname: {}})

                elif (parsed[0] == "SubjectName") and parsed[1]:
                    self.alignQ = []
                    self.alignS = []
                    self.sname = parsed[2]

                elif (parsed[0] == "AlignQ") and parsed[1]:
                    self.alignQ += parsed[2]

                elif (parsed[0] == "AlignS") and parsed[1]:
                    self.alignS += parsed[2]

                elif (self.queryname is not None and
                      self.sname is not None and
                      len(self.alignQ) > 0 and
                      len(self.alignS) > 0):
                    aQ = (self.alignQ[0], self.alignQ[-1])
                    aS = (self.alignS[0], self.alignS[-1])
                    proteins[self.queryname].update({self.sname: [aQ, aS]})

                else:
                    continue
                    
        self.model.stop()
        return proteins

# ----------- #
# Define the BLAST processing functions for the rules

def fxn_align_startstop(line):
    """ Split an alignment string for start/end state """
    line = line.strip(">\n").split()
    return [int(line[1]), int(line[-1])]

fxn_sbjct_name = lambda x: x.strip("\n>").split()[0]
fxn_query_name = lambda x: x.strip("\nQuery= ")

#Transition Rules
rules = [("Start", "Start", lambda x: not (x.startswith("Query") or x.startswith("Sbjct")) ),
         ("Start", "Query", lambda x: x.startswith("Query= "), fxn_query_name), 
         ("Query", "Query", lambda x: not x.startswith(">") ), 
         ("Query", "SubjectName", lambda x: x.startswith(">"), fxn_sbjct_name),
         ("SubjectName", "SubjectName", lambda x: not x.startswith("Query ") ),
         ("SubjectName", "AlignQ", lambda x: x.startswith("Query "), fxn_align_startstop ),
         ("AlignQ", "AlignS", lambda x: x.startswith("Sbjct "), fxn_align_startstop ),
         ("AlignQ", "AlignQ", lambda x: not x.startswith("Sbjct ") ),              
         ("AlignS", "AlignQ", lambda x: x.startswith("Query "), fxn_align_startstop),
         ("AlignS", "SubjectName", lambda x: x.startswith(">"), fxn_sbjct_name),
         ("AlignS", "Query", lambda x: x.startswith("Query= "), fxn_query_name),
         ("AlignS", "AlignS", lambda x: not(x.startswith("Query") or
                                            x.startswith(">")) )
         ]

if __name__ == "__main__":

    filename = "blastp.out"

    parser = ProcessBLAST(filename)
    parser.init_FSM(transitions=rules)

    t = time.time()
    proteins = parser.run_model()
    print("Elapsed time", time.time()-t)

    with gzip.open('proteins.pkl.gzip', 'wb') as f:
        pkl.dump(proteins, f)
