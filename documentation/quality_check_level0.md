# Quality check level0 data

## Objective and role of this function
perform a quick check of the raw data to spot serious problem

### Called from
Main routine

### Calling
None

### Previously
Reading routine level0

### Next
Calibration

## Inputs
* 1 standardized matlab structure containing the **raw** data.
* 1 standardized matlab structure containing the **log** data.

## Outputs
* Error code depending on the spotted problem

## Structure
At this level, there are not many check to do but here some ideas:
* Number of cycle per days ?
* Size of the raw data corresponding to its log ?

## How it was done in the past
### GROMOS

Not done

### SOMORA

Not done

## Questions remaining
...
