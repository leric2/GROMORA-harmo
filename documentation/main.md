# Pist for the main script

Modules vs toolchain methods

Add some tests

## Suggestion
Take a sort of "object oriented" approach for launching the retrievals. 

We have an object (let's call it retrievalTool) that will contain all the relevant information for the retrieval we want to perform (date, name of the instruments, etc...). We manually define this object for the operation we want to perform (for instance, for running a retrieval between 2 given dates in Payerne, without quality checking anything) and then launch the retrieval through a "run" function.

As launching a retrieval is a step-by-step "linear" operation, we do not need to edit the script effectively doing the step-by-step operation but we only need to edit the retrievalTool object. 

On of the nice thing is that we can store function in this retrievalTool object, so that the retrievals can use different functions for different instruments in the case where it would be needed. This would be defined and documented in the retrievalTool structure and implemented in the main "run" function.

CAUTION: switch cases inside functions should be avoided because it would be painfull to edit them all when a new instruments is added (or a new spectrometers for instance).


