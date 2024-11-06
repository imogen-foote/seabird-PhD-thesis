#Chapter 3
Variant calling from population samples of Antipodean and Gibson's albatross (Diomedea antipodensis antipodensis & D. a. gibsoni) aligned to a reference genome (see albatross-refgenomes and albatross-popgenomics for assembly and alignment), and multiple sequence alignment.

Variant calling uses the Genome Analysis ToolKit (GATK) joint genotyping approach. Reads are filtered with BCFtools. Consensus FASTA sequences for all individuals in the populations are created using GATK's FastaAlternateReferenceMaker. All sequences are aligned using MAFFT multiple sequence alignment.
