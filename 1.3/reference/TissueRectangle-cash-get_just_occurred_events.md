# Getting the events that just occurred

This method returns the list of the just occurred events.

## Value

A list describing the events that occurred at the current simulated
time. Each element is a named list. The names describe the event type
and are among `\"rate update\"`, `\"mutant emerged\"`, and
`\"sampling\"`. The type of values, instead, depend on the name. When
the name is `\"rate update\"`, the value is a data frame containing the
rate updates. When the name is `\"mutant emerged\"`, the value is the
name of just emerged mutant. Finally, when the name is `\"sampling\"`,
the value is the name of the sample.
