var documenterSearchIndex = {"docs":
[{"location":"man/act_guide/#Activity-Guide-1","page":"Activity Guide","title":"Activity Guide","text":"","category":"section"},{"location":"man/spec_guide/#Specificity-Guide-1","page":"Specificity Guide","title":"Specificity Guide","text":"","category":"section"},{"location":"man/spec_guide/#","page":"Specificity Guide","title":"Specificity Guide","text":"How to use specificity module of siRNATools","category":"page"},{"location":"man/act_index/#Activity-Functions-1","page":"Activity Functions","title":"Activity Functions","text":"","category":"section"},{"location":"#siRNATools-1","page":"Home","title":"siRNATools","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"A Package with tools for siRNA","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Package is sorted into several submodules","category":"page"},{"location":"#Specificity-1","page":"Home","title":"Specificity","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Pages = [\n    \"man/spec_guide.md\",\n    \"man/spec_index.md\"\n]","category":"page"},{"location":"#Activity-1","page":"Home","title":"Activity","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Pages = [\n    \"man/act_guide.md\",\n    \"man/act_index.md\"\n]","category":"page"},{"location":"man/spec_index/#Specificity-Functions-1","page":"Specificity Functions","title":"Specificity Functions","text":"","category":"section"},{"location":"man/spec_index/#","page":"Specificity Functions","title":"Specificity Functions","text":"","category":"page"},{"location":"man/spec_index/#Reference-Sequence-1","page":"Specificity Functions","title":"Reference Sequence","text":"","category":"section"},{"location":"man/spec_index/#","page":"Specificity Functions","title":"Specificity Functions","text":"siRNATools.Specificity.ReferenceSequence\nsiRNATools.Specificity.get_refseq_pos\nsiRNATools.Specificity.encode_refseq\nsiRNATools.Specificity.decode_refseq\nsiRNATools.Specificity.decode_refseq_partial\nsiRNATools.Specificity.download_RefSeq\nsiRNATools.Specificity.process_RefSeq\nsiRNATools.Specificity.save_RefSeq","category":"page"},{"location":"man/spec_index/#siRNATools.Specificity.ReferenceSequence","page":"Specificity Functions","title":"siRNATools.Specificity.ReferenceSequence","text":"Structure for reference sequences.  Compresses RNA data into 2 bits of information from the 8 of a normal character string.  Can only use bases A, C, G, U.\n\nA => 00\nC => 01\nG => 10\nU => 11\n\n\n\n\n\n","category":"type"},{"location":"man/spec_index/#siRNATools.Specificity.get_refseq_pos","page":"Specificity Functions","title":"siRNATools.Specificity.get_refseq_pos","text":"get_refseq_pos(::ReferenceSequence, ::Int)\n\nFunction to retrieve bit value of base at specificied position in a reference sequence.\n\n\n\n\n\n","category":"function"},{"location":"man/spec_index/#siRNATools.Specificity.encode_refseq","page":"Specificity Functions","title":"siRNATools.Specificity.encode_refseq","text":"encode_refseq(::String) :: ReferenceSequence\n\nFunction that takes a string containing only A, C, G, and U and encodes it into a ReferenceSequence type\n\n\n\n\n\n","category":"function"},{"location":"man/spec_index/#siRNATools.Specificity.decode_refseq","page":"Specificity Functions","title":"siRNATools.Specificity.decode_refseq","text":"decode_refseq(::ReferenceSequence) :: String\n\nFunction takes a ReferenceSequence type and returns a String containing the RNA bases.\n\n\n\n\n\n","category":"function"},{"location":"man/spec_index/#siRNATools.Specificity.decode_refseq_partial","page":"Specificity Functions","title":"siRNATools.Specificity.decode_refseq_partial","text":"decode_refseq_partial(::ReferenceSequence, ::UnitRange) :: String\n\nFunction which takes a ReferenceSequence type and outputs specificied range of bases as a String.\n\n\n\n\n\n","category":"function"},{"location":"man/spec_index/#siRNATools.Specificity.download_RefSeq","page":"Specificity Functions","title":"siRNATools.Specificity.download_RefSeq","text":"download_RefSeq(::UnitRange{Int64}=1:8, ::String=PATH)\n\nDownloads mRNA reference sequence from ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/ to the PATH folder.  Defaults to downloading 8 files.\n\n\n\n\n\n","category":"function"},{"location":"man/spec_index/#siRNATools.Specificity.process_RefSeq","page":"Specificity Functions","title":"siRNATools.Specificity.process_RefSeq","text":"process_RefSeq(::UnitRange{Int64}=1:8, ::String=PATH)\n\nProcesses raw gzipped fasta files into DataFrame and saves it as a CSV in the PATH folder\n\n\n\n\n\n","category":"function"},{"location":"man/spec_index/#siRNATools.Specificity.save_RefSeq","page":"Specificity Functions","title":"siRNATools.Specificity.save_RefSeq","text":"save_RefSeq(::String=PATH)\n\nSaves relevant data structures from the processes mRNA reference sequence for use in searches.\n\nTranscriptGene => dictionary of Transcripts to Genes\nGeneTranscripts => dictionary of Genes to Transcripts\nallT => dictionary of Transcript name to base sequence as String\nallRefSeq => dictionary of Transcript name to base sequence as ReferenceSequence\n\n\n\n\n\n","category":"function"},{"location":"man/spec_index/#Genome-Search-1","page":"Specificity Functions","title":"Genome Search","text":"","category":"section"},{"location":"man/spec_index/#","page":"Specificity Functions","title":"Specificity Functions","text":"siRNATools.Specificity.Calculate_Specificity\nsiRNATools.Specificity.reverse_complement\nsiRNATools.Specificity.calculate_Peq\nsiRNATools.Specificity.motif_to_transcript_match\nsiRNATools.Specificity.find_match_sequences\nsiRNATools.Specificity.mismatch_positions\nsiRNATools.Specificity.find_genome_matches\nsiRNATools.Specificity.compress_genome_matches\nsiRNATools.Specificity.final_calc","category":"page"},{"location":"man/spec_index/#siRNATools.Specificity.Calculate_Specificity","page":"Specificity Functions","title":"siRNATools.Specificity.Calculate_Specificity","text":"Calculate_Specificity(patterns, excluded_gene=\"\", rg=2:18, verbose=true) :: DataFrame\n\nFunction takes as input the pattern to be searched against the genome, the excluded gene to be ignored, the range of pattern being used, and a boolean which controls whether progress bars will be shown.  Output is a DataFrame with a column for the pattern, the number of genes with minimum mismatch distance of 0-4 and the specificity score.\n\n\n\n\n\n","category":"function"},{"location":"man/spec_index/#siRNATools.Specificity.reverse_complement","page":"Specificity Functions","title":"siRNATools.Specificity.reverse_complement","text":"reverse_complement(::String) :: String\n\nTakes the reverse complement of an RNA string\n\n\n\n\n\n","category":"function"},{"location":"man/spec_index/#siRNATools.Specificity.calculate_Peq","page":"Specificity Functions","title":"siRNATools.Specificity.calculate_Peq","text":"calculate_Peq(::String) :: Array{UInt64, 1}\n\nTakes as input an RNA strand and translates it into an array of Unsigned Integers each one representing one letter and each bit referring to a position.   Used as a pre-processing step for Myers searching algorithm.  Maximum length of searched pattern is 64.\n\nExample\n\nAACGCU becomes [110000, 001010, 000100, 000001]\n\n\n\n\n\n","category":"function"},{"location":"man/spec_index/#siRNATools.Specificity.motif_to_transcript_match","page":"Specificity Functions","title":"siRNATools.Specificity.motif_to_transcript_match","text":"motif_to_transcript_match(Peq::Array{Unit64, 1}, m::Int64, refseq::ReferenceSequence, min\\_K::Int) :: Array{Tuple{UInt64, UInt64}}\n\nTakes as input Peq calculated in calculate_Peq, m = the length of the pattern to search for, the ReferenceSequence to search in, and mink =  the maximum distance to output.  Implements the string matching algorithm described in _A Fast Bit-Vector Algorithm for Approximate String Matching Based on Dynamic Programming by Gene Myers.  Output is a tuple of all Levenshtein distances less than min_k and their positions in refseq.   \n\n\n\n\n\n","category":"function"},{"location":"man/spec_index/#siRNATools.Specificity.find_match_sequences","page":"Specificity Functions","title":"siRNATools.Specificity.find_match_sequences","text":"find_match_sequences(motif::String, sequence::String, mismatches::Int) :: Array{String, 1}\n\nFunction takes as input a motif, sequence, and number of mismatches to search for.  Output is an array of all substrings of sequence which have a  Hamming distance of exactly mismatches to motif. \n\n\n\n\n\n","category":"function"},{"location":"man/spec_index/#siRNATools.Specificity.mismatch_positions","page":"Specificity Functions","title":"siRNATools.Specificity.mismatch_positions","text":"mismatch_positions(::String, ::String) :: Array{Int, 1}\n\nFunction returns an array of all positions in which strings used as input differ.  Strings must be the same length.\n\n\n\n\n\n","category":"function"},{"location":"man/spec_index/#siRNATools.Specificity.find_genome_matches","page":"Specificity Functions","title":"siRNATools.Specificity.find_genome_matches","text":"find_genome_matches(pattern::String, excluded_gene::String = \"\", verbose::Bool = true, minimum_matches = 5) :: Array{Tuple{String, Int64}}\n\nFunction takes as input a pattern to search the genome for, excludedgene to exclude the gene of interest from the search, verbose set to true displays a progress bar showing progress of search, and minimummatches which is the amount of mismatches searched for + 1.  Output is an array of tuples of transcript names and the lowest Hamming distance found to the pattern within that transcript.\n\n\n\n\n\n","category":"function"},{"location":"man/spec_index/#siRNATools.Specificity.compress_genome_matches","page":"Specificity Functions","title":"siRNATools.Specificity.compress_genome_matches","text":"compress_genome_matches(::Array{Tuple{String, Int64}}) :: Dict{String, Array{Int64, 1}}\n\nFunction takes as input the output of find_genome_matches and collapses it into a dictionary of Gene names and an array of all lowest Hamming mismatch  of the transcripts of that gene.\n\n\n\n\n\n","category":"function"},{"location":"man/spec_index/#siRNATools.Specificity.final_calc","page":"Specificity Functions","title":"siRNATools.Specificity.final_calc","text":"final_calc(pattern::String, raw_data::Array{Tuple{String, Int64}}, compressed_data::Dict{String, Array{Int64, 1}})\n\nFunction takes as input the pattern being searched for, the rawdata from [find\\genome_matches](@ref), and the compressed data from compress_genome_matches. Ouput is a Tuple containing the mismatch_counts dictionary which adds up the number of genes with a minimum Hamming distance of 0-4, and the specificity score.\n\n\n\n\n\n","category":"function"}]
}
