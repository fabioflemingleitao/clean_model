# clean_model
This repository is about a "fatal error" that occurs in an OpenSeesPy model

# Hello everyone,
Sorry for resurrecting this topic.  (https://www.facebook.com/share/p/1DuNexzfdK/) https://www.facebook.com/groups/opensees
However, I haven’t been able to find a solution. I’ve been running analyses 
on frame structures subjected to seismic accelerations. Some analyses converge, 
and others do not (which is normal). However, in some specific situations 
(usually with high acceleration signals due to a high scale factor), a fatal error 
occurs and the analysis is automatically interrupted (crashing the analysis).
I’ve been trying every day to solve this error, but without success. 


By debugging this code, the "fatal error" occurs in an "ops.analyze" command 
around time = 42.421s of the earthquake. Perhaps by debugging inside OpenSees (in C++) 
it might be possible to identify what’s happening. But I don’t know how to do that.

# Today, I decided to create a simple model where I could replicate this "fatal error" 
so you can try to help me.

Note: I noticed that in this analysis, when I remove the reinforcement from my 
beam section, the fatal error does not occur.


Thank you for the help.
Fábio Leitão
