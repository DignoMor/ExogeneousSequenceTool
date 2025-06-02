# ExogeneousSequenceTool
Tools to work on Exogeneous Sequences

## Exogeneous Sequence

An Exogeneous Sequence is a data structure similar 
to a [Genomic Element](https://github.com/DignoMor/GenomicElementTool).
However, Exogeneous Sequences do not 
have a region file, thus it has only 2 components.

- sequence file 
- annotations

### Exogeneous Sequence file

Just like in Genomic Element, Exogeneous Sequence 
has sequence file in fasta format. The difference 
is that instead of being immutable genomic sequence, 
functionalities of ExogeneousSequenceTool can actually 
change the content of the fasta sequence, thus generate 
new Exogeneous Sequence instances.

### Annotation file

Annotation files are the same as Genomic Element.

