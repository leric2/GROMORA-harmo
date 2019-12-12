# Pist for the main script

Modules vs toolchain methods

Add some tests

## Suggestion
Take a sort of "object oriented" approach for launching the retrievals. 

We have an object (let's call it retrievalTool) that will contain all the relevant information for the retrieval we want to perform (date, name of the instruments, etc...). We manually define this object for the operation we want to perform (for instance, for running a retrieval on between 2 given dates in Payerne, without quality checking anything) and then launch the retrieval through a "run" function.

As launching a retrieval is a step-by-step operation, we do not need to edit the script effectively doing the step-by-step operation but we only need to edit the retrievalTool object. 

On of the nice thing is that we can store function in this retrievalTool object, so that the retrievals can use different functions for different instruments in the case where it would be needed (see when ??). When possible, only one function with switch case shall be prefered. 

CAUTION: what about adding an instrument ? would be nice not to edit all switch cases !!!!!!!

