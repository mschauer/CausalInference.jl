# Developer comments

## Whether to use Meek's or Chickerings completion algorithm.

Meek's approach we currently use in PC and GES is: Take the PDAG and apply R1-R4 repeatedly until no rule applies. To obtain the CPDAG from a pattern the rules R1-R3 suffice, R4 is not needed.

Chickering's approach is to take the PDAG and find a DAG with same skeleton and v-structures. (This can be done e.g. with this algorithm by Dor and Tarsi: https://ftp.cs.ucla.edu/pub/stat_ser/r185-dor-tarsi.pdf)
From the DAG construct the CPDAG (which is possible in linear-time, we also have this implemented cpdag(g)).

Implementing Chickerings approach would be an option because the complicated part (from an implementation perspective) is DAG-to-CPDAG which we already have. A straightforward implementation of PDAG to DAG by Dor and Tarsi are just a few lines. For speed of PC and GES this will likely not matter too much overall, other parts are more costly in case of PC the pattern learned from data might not have a DAG with same skeleton and v-structures. Then, it is not clear how to use Chickerings approach. 

## On Meek's rule 4

R4 is only needed in case of background knowledge, i.e. some edge orientations which do not follow from v-structures. For soundness of the rule, it is necessary that v and l are adjacent. For completeness of the Meek rules it suffices to apply it only if there is an undirected edge.