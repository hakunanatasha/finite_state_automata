# finite_state_automata
Text parsing semi-structured data via an FSA

2020 May 22

The following is the completed FSA that I created for parsing BLAST outputs. BLAST is a sequence alignment tool, where provided a database (or the canonical db), it will tell you how similar a sequence is to a list of your own target sequences based on various alignment heuristics.

Blast will output several possible candidates based on particular scores. When the dB is huge or if the target is huge, you could potentially get many matches. Moreover, if your threshold of score is low (for example, in an exploratory sense), then you may find that you have to subsequently find better quality candidates one the calculation has been run.

From a computational perspective, BLAST outputs are semi-structured. While the details themselves may be different (i.e different sequence lengths), the format is particular.

The following code allows you to identify these text transitions and scrape what is required from them.

The program is written quite generally, with the parser expecting certain transition rules of when a new (query, subject) pair has been identified.

The following program operates using a 'blastp.out' file (a typical blast file). 