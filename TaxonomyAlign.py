import argparse
import logging
import time
import shutil
import os
import subprocess as sp
from Bio import SeqIO
from Bio.Seq import Seq

def get_args():
    """
    Get arguments from command line.
    """
    parser = argparse.ArgumentParser(
            description="""Taxonomy-Guided-Alignment: Use to perform multiple sequence alignment for a highly variable gene 
            and a large number of divergent sequences from many species. This program automatically divides the sequences 
            into smaller subsets based on user-supplied taxonomic groupings, and builds sub-alignments for the groups using MAFFT 
            (several algorithm options available). The defined groupings should be related to both taxonomy and phylogeny (e.g., genus, 
            family, superfamily, etc). The resulting sub-alignments are then merged together with MAFFT, with or without an optional 
            polishing step. The original use-case was for creating high-quality alignments from many 12S and 16S mtDNA sequences 
            across a large phylogenetic scale. However, it should prove useful for other genes with many species and divergent sequences 
            to align.""")
    
    parser.add_argument("-f", "--fasta",
                            required=True,
                            help="REQUIRED: The full path to the input fasta file of unaligned sequences.")
        
    parser.add_argument("-m", "--map",
                            required=True,
                            help="REQUIRED: The full path to a map file containing sequence names "
                            "and corresponding groupings/taxonomy.")
        
    parser.add_argument("-o", "--outdir",
                            required=True,
                            help="REQUIRED: The full path to an existing directory "
                            "to write output files.")
    
    parser.add_argument("-a", "--algorithm",
                            required=False,
                            default="auto",
                            choices=["auto", "FFT-NS-i", "E-INS-i", "L-INS-i", "G-INS-i"],
                            help="OPTIONAL: The MAFFT algorithm to use: 1) auto-select based on input sequences, "
                            "2) FFT-NS-i is a progressive and general algorithm, 3) E-INS-i is best for several "
                            "conserved core regions and variable regions (e.g., 12S and 16S mtDNA), "
                            "4) L-INS-i is best for a single conserved core region with variable flanking regions, and 5) "
                            "G-INS-i is best for full-length single core sequences. Defaults to auto.")
    
    parser.add_argument("-t", "--threads",
                            default=None,
                            help="OPTIONAL: Specifies number of threads to use.")
        
    parser.add_argument("--polish",
                            action='store_true',
                            help="OPTIONAL: Performs additional distance measure refinement on final merged alignment, "
                            "at cost of increasing run time (sometimes substantially).")
    
    parser.add_argument("--output_seqnames",
                            action='store_true',
                            help="OPTIONAL: Write sorted names of all sequences in fasta file to stdout and quit. "
                            "Must still provide -f, -m, and -o required flags above (though value only matters for -f).")

    parser.add_argument("--verbose",
                            action='store_true',
                            help="OPTIONAL: Show mafft progress during final merge step instead of redirecting to log file. "
                            "May be useful for monitoring the --polish option, which can be slow for larger datasets.")

    return parser.parse_args()

def setup_logging(outdir):
    """
    Set up logging to log file and on-screen.

    Args:
    outdir - full path to existing directory to write outputs
    """
    # set up logging to file
    logging.basicConfig(filename=os.path.join(outdir, "Taxonomy-Align_{}.log".format(time.strftime("%b-%d-%Y_%I.%M.%S"))),
                            format="%(levelname)s: %(asctime)s: %(message)s",
                            datefmt='%d-%b-%y %H:%M:%S',
                            level=logging.DEBUG)
    
    # set up logging to console 
    console = logging.StreamHandler() 
    console.setLevel(logging.INFO) 
    formatter = logging.Formatter("%(asctime)s: %(levelname)s: %(message)s",
                                      datefmt='%H:%M:%S')
    console.setFormatter(formatter) 
    logging.getLogger().addHandler(console) 
    
def make_output_paths(outdir):
    """
    Function to create output directories in directory specified by user.

    Args:
    outdir - full path to existing directory to write outputs

    Returns:
    dir_groups - full path to directory 1-Group-Fasta-Files
    dir_alns - full path to directory 2-Group-Alignments
    dir_merged - full path to directory 3-Merged-Alignments
    """
    os.chdir(outdir)
    curpath = os.getcwd()
    dir_groups = os.path.join(curpath, "1-Group-Fasta-Files")
    dir_alns = os.path.join(curpath, "2-Group-Alignments")
    dir_merged = os.path.join(curpath, "3-Merged-Alignments")

    for path in [dir_groups, dir_alns, dir_merged]:
        if not os.path.exists(path):
            os.mkdir(path)

    return dir_groups, dir_alns, dir_merged

def fasta_to_dict(fasta_file):
    """
    Function to convert any fasta file (f) into
    a dictionary structure with taxon as key and 
    sequence as value.

    Args:
    fasta_file - full path to the input fasta file of unaligned sequences

    Returns:
    fasta_dict - dictionary of seq names as keys and seqs as values
    """
    logging.info("Parsing fasta file.")
    
    fasta_dict = {}
    with open(fasta_file, 'r') as fh:
        lines = [l.strip() for l in fh if l.strip()]
        
    for line in lines:
        if line.startswith(">"):
            new_key = line.replace(">",'')
            fasta_dict[new_key] = ""
        else:
            fasta_dict[new_key] += line.upper()
            
    logging.info("Found {:,} sequences in fasta file.\n\n".format(len(fasta_dict)))
        
    return fasta_dict

def map_to_dict(map_file):
    """
    Processes map file to obtain a dictionary of group names as keys and 
    a list of taxa as the values. Returns this dictionary. 

    The map file should be tab-delimited and contain exactly two columns. Column one 
    is the name of a sequence (ideally in format like 'Hyperolius_guttulatus') and 
    column two is the group assignment for that sequence ('Hyperoliidae'). It is 
    expected that multiple sequences will be assigned to the same group, and that there 
    must be at least two groups in order to run this program. 

    Args:
    map_file - full path to the input map file of seq names and group names

    Returns:
    map_dict - dictionary of group names as keys and list of seq names as values
    """
    logging.info("Parsing map file.")
    
    map_dict = {}
    with open(map_file, 'r') as fh:
        for line in fh:
            if line.strip() and len(line.split()) == 2:
                if line.strip().split()[1] not in map_dict:
                    map_dict[line.strip().split()[1]] = [line.strip().split()[0]]
                else:
                    map_dict[line.strip().split()[1]].append(line.strip().split()[0])
                    
    logging.info("Found {:,} groups in map file.\n\n".format(len(map_dict)))
                    
    return map_dict
            
def check_map_fasta_compatibility(fasta_dict, map_dict):
    """
    Perform a quick check to see which sequences in the fasta file 
    can be matched to those provided in the mapping file. Reports all 
    unmatched sequences in the fasta file. 

    Args:
    fasta_dict - dictionary of seq names as keys and seqs as values (produced by fasta_to_dict())
    map_dict - dictionary of group names as keys and list of seq names as values (produced by map_to_dict())
    """
    fasta_taxa = set(list(fasta_dict.keys()))
    
    map_taxa = set()
    for k, v in map_dict.items():
        for i in v:
            map_taxa.add(i)
            
    unique_fasta_taxa = fasta_taxa - map_taxa
    
    if len(unique_fasta_taxa) == 0:
        logging.info("All sequences in fasta file have an assigned grouping.\n\n")
    else:
        logging.info("The following {:,} sequences are not assigned to a group:\n{}".format(len(unique_fasta_taxa), "\n\t".join(sorted(unique_fasta_taxa))))
        logging.info("They will be excluded from the analysis.\n\n")

def write_group_fastas(fasta_dict, map_dict, dir_groups):
    """
    Write group-specific fasta files to output directory.

    Args:
    fasta_dict - dictionary of seq names as keys and seqs as values (produced by fasta_to_dict())
    map_dict - dictionary of group names as keys and list of seq names as values (produced by map_to_dict())
    dir_groups - full path to full path to directory 1-Group-Fasta-Files
    """
    logging.info("Writing fasta files for groups:")
    
    for group in sorted(map_dict.keys()):
        # get intersection (common vals) of the taxon sets
        matched_taxa = set(map_dict[group]) & set(list(fasta_dict.keys()))
        
        # if one sequence present for group, write singleton fasta
        if len(matched_taxa) == 1:
            # remove any existing version of this group fasta file
            out_fasta = os.path.join(dir_groups, "{}.singleton.fasta".format(group))
            if os.path.exists(out_fasta):
                logging.warning("Removing prior existing file: {}".format(out_fasta))
                os.remove(out_fasta)
            with open(out_fasta, 'a') as fh:
                for taxon in sorted(matched_taxa):
                    fh.write(">{}\n{}\n".format(taxon, fasta_dict[taxon]))
            logging.info("\tGroup {} is represented by {} fasta sequences.".format(group, len(matched_taxa)))

        # if 2 or more sequences present for group, write multiples fasta            
        elif len(matched_taxa) > 1:
            # remove any existing version of this group fasta file
            out_fasta = os.path.join(dir_groups, "{}.multiples.fasta".format(group))
            if os.path.exists(out_fasta):
                logging.warning("Removing prior existing file: {}".format(out_fasta))
                os.remove(out_fasta)
            with open(out_fasta, 'a') as fh:
                for taxon in sorted(matched_taxa):
                    fh.write(">{}\n{}\n".format(taxon, fasta_dict[taxon]))
            logging.info("\tGroup {} is represented by {} fasta sequences.".format(group, len(matched_taxa)))
            
        # if no sequences present for group, skip fasta writing
        else:
            logging.info("\tGroup {} is NOT represented in fasta sequences, skipping.".format(group, len(matched_taxa)))
            
    logging.info("Finished writing group fasta files.\n\n")

def reformat_mafft_output_fasta(temp_fasta, out_fasta):
    """
    Checks the alignment output file of mafft to ensure worked, and if so 
    will cleanup the format and write to the final alignment file.

    Args:
    temp_fasta - the full path to the alignment file produced directly by mafft 
    out_fasta - the full path to the final alignment file to write after cleaning
    """
    if not os.stat(temp_fasta).st_size != 0:
        logging.info("ERROR: Alignment failed. Please check log file for details.\n\n")
        
    else:
        if os.path.exists(out_fasta):
            logging.warning("Removing prior existing file: {}".format(out_fasta))
            os.remove(out_fasta)

        fdict = SeqIO.index(temp_fasta, "fasta")
        with open(out_fasta, 'a') as fh:
            for record in fdict:
                newseq = fdict[record].seq.upper()
                fh.write( ">{}\n{}\n".format(fdict[record].description, newseq))           
        os.remove(temp_fasta)
        logging.info("Finished alignment!\n\n")

def align_individual_group_fastas(dir_groups, dir_alns, threads, algorithm):
    """
    Runs mafft with the specified algorithm for all group fasta files containing 
    at least two sequences. Will clean the outputs from mafft to write a nicely 
    formatted output alignment fasta file.

    Args:
    dir_groups - full path to directory 1-Group-Fasta-Files
    dir_alns - full path to directory 2-Group-Alignments
    threads - optional number of threads to use for mafft
    algorithm - user choice for alignment algorithm to run in mafft
    """
    # show information about singleton groups, if they exist
    singleton_fasta_list = [os.path.join(dir_groups, f) for f in os.listdir(dir_groups) if f.endswith(".singleton.fasta")]
    if singleton_fasta_list:
        logging.info("There are {:,} group fasta files containing a single sequence:".format(len(singleton_fasta_list)))
        for f in singleton_fasta_list:
            logging.info("\t{}".format(f.split('/')[-1]))
        logging.info("These will be individually aligned to the merged group alignment.\n\n".format(len(singleton_fasta_list)))

    # show information about multiples groups
    multiples_fasta_list = sorted([os.path.join(dir_groups, f) for f in os.listdir(dir_groups) if f.endswith(".multiples.fasta")])
    logging.info("Found {:,} group fasta files suitable for alignment.\n\n".format(len(multiples_fasta_list)))
    # if all groups are single sequence, kill the program
    if not multiples_fasta_list:
        message = "ERROR: No group fasta files contain >1 sequence!"
        logging.info("{}".format(message))
        raise ValueError(message)

    if threads:
        thread_str = "--thread {}".format(threads)
    else:
        thread_str = ""

    algo_dict = {"FFT-NS-i":"--maxiterate 1000",
                     "E-INS-i":"--genafpair --maxiterate 1000",
                     "L-INS-i":"--localpair --maxiterate 1000",
                     "G-INS-i":"--globalpair --maxiterate 1000",
                     "auto":"--auto"}
        
    # iterate over group fasta files to perform individual alignments
    for input_fasta in multiples_fasta_list:
        name = input_fasta.split('/')[-1].split(".")[0]
        logging.info("Beginning alignment of {} using {} algorithm.".format(name, algorithm))

        # set names for all outputs
        temp_fasta = os.path.join(dir_alns, "{}.mafft.temp.fasta".format(name))
        if os.path.exists(temp_fasta):
            logging.warning("Removing prior existing file: {}".format(temp_fasta))
            os.remove(temp_fasta)
        out_fasta = os.path.join(dir_alns, "{}.mafft-{}.fasta".format(name, algorithm))
        log = os.path.join(dir_alns, "{}.mafft-{}.log".format(name, algorithm))

        # generate mafft command and execute
        mafft_call = "mafft {0} {1} {2} > {3} 2> {4}".format(thread_str, algo_dict[algorithm], input_fasta, temp_fasta, log)
        logging.info("{}".format(mafft_call))
        sp.call(mafft_call, shell=True)

        # reformat output file if alignment successful - will also notify if alignment failed
        reformat_mafft_output_fasta(temp_fasta, out_fasta)
        
    logging.info("Finished all group-specific alignments.\n\n")

def prepare_subalignments_for_merge(dir_groups, dir_alns, dir_merged, algorithm):
    """
    This function prepares a concatenated sub-alignments fasta file and a 
    sequence assignment table. These are mafft-specific input files required 
    for the --merge algorithm. 

    The sub-alignments are simply concatenated in fasta format, but the numerical sequence order 
    is used to specify subalignment sequence sets in the table file. 

    Example table:
    1 2 3 4 5  # this is comment. subMSA1
    6 7 8      # you can write anything after #
    9 10       # subMSA3

    In the above, sequences 1-5 are aligned, 6-8 are aligned, 9-10 are aligned, and any other additional sequences in 
    the concatenated sub-alignments fasta file (e.g., 10, 11, 12, etc.) are considered unaligned. Those unaligned 
    sequences will get aligned to the merged profile alignment during this final step. 

    The output table produced here will have the group name commented after the sequence numbers:
    1 2 3 4 5 # Group = Hyperoliidae

    Args:
    dir_groups - full path to directory 1-Group-Fasta-Files
    dir_alns - full path to directory 2-Group-Alignments
    dir_merged - full path to directory 3-Merged-Alignments
    algorithm - user choice for alignment algorithm to run in mafft

    Returns:
    table - full path to the sequence assignment table
    cataln - full path to the concatenated sub-alignments fasta file
    """

    logging.info("Preparing sub-alignment file and sequence table.\n\n")
    
    # start sequence counter
    seq_number = int(1)

    # create table filename
    table = os.path.join(dir_merged, "Merge-Table.txt")
    if os.path.exists(table):
        logging.warning("Removing prior existing file: {}".format(table))
        os.remove(table)

    # create concatenated msa filename
    cataln = os.path.join(dir_merged, "Sub-Alignments.fasta")
    if os.path.exists(cataln):
        logging.warning("Removing prior existing file: {}".format(cataln))
        os.remove(cataln)

    # iterate over aligned fasta files and write required input files
    aligned_fasta_list = sorted([os.path.join(dir_alns, f) for f in os.listdir(dir_alns) if f.endswith(".mafft-{}.fasta".format(algorithm))])
    if not aligned_fasta_list:
        message = "ERROR: No aligned fasta files found!"
        logging.info("{}".format(message))
        raise ValueError(message)
    else:
        logging.info("Adding sub-alignments:")
        
    for aln_file in aligned_fasta_list:
        temp_seq_numbers = []
        with open(cataln, 'a') as fhout, open(aln_file, 'r') as fhin:
            for line in fhin:
                if line.startswith('>'):
                    temp_seq_numbers.append(str(seq_number))
                    seq_number += 1
                fhout.write(line)
        with open(table, 'a') as fh:
            fh.write("{} # Group = {}\n".format(" ".join(temp_seq_numbers), aln_file.split('/')[-1].split('.')[0]))
            
        logging.info("\tsub-alignment {}: wrote {:,} aligned sequences ({:,} - {:,})".format(aln_file.split('/')[-1],
                                                                                         len(temp_seq_numbers),
                                                                                         int(temp_seq_numbers[0]), int(temp_seq_numbers[-1])))
        
    logging.info("Finished adding sub-alignments.\n\n")
    
    # add singleton sequences to end of fasta file, no numbers added to table because these are unaligned
    singleton_fasta_list = [os.path.join(dir_groups, f) for f in os.listdir(dir_groups) if f.endswith(".singleton.fasta")]
    
    if singleton_fasta_list:
        logging.info("Writing singleton group sequences to sub-alignment file:")
        for fasta_file in singleton_fasta_list:
            with open(cataln, 'a') as fhout, open(fasta_file, 'r') as fhin:
                for line in fhin:
                    if line.startswith('>'):
                        logging.info("\tAdding single sequence from {}: {}".format(fasta_file.split('/')[-1].split('.')[0], line.strip().strip('>')))
                    fhout.write(line)
        logging.info("Finished adding singleton group sequences.\n\n")

    else:
        logging.info("No singleton group sequences to add to sub-alignment file.\n\n")
        
    logging.info("Finished preparing sub-alignment file and sequence table.\n\n")
    
    return table, cataln

def run_merge(table, cataln, dir_merged, threads, polish, verbose):
    """
    Performs the merge algorithm in mafft based on the sequence assignment table 
    and concatenated sub-alignments fasta file produced. 

    Options for running mafft merge algorithm:
    1) For basic merge run:
    mafft --merge table input > output

    2) A more rigorous distance measure can be used for small data:
    mafft --localpair --merge table input > output

    3) Iterative refinement is supported in version 7.040 and higher:
    mafft --localpair --maxiterate 10 --merge table input > output

    The absence of --polish invokes option 1 above, the presence of --polish 
    will invoke option 3. 

    Args:
    table - full path to the sequence assignment table
    cataln - full path to the concatenated sub-alignments fasta file
    dir_merged - full path to directory 3-Merged-Alignments
    algorithm - user choice for alignment algorithm to run in mafft
    threads - optional number of threads to use for mafft
    polish - a true/false toggle for the distance measure + interative refinement step after merge
    """
            
    logging.info("Starting the alignment merging step.")

    # set names for all outputs
    temp_fasta = os.path.join(dir_merged, "temp.fasta")
    out_fasta = os.path.join(dir_merged, "Merged-alignment.fasta")
    log = os.path.join(dir_merged, "Merged-alignment.log")

    if threads:
        thread_str = "--thread {} ".format(threads)
    else:
        thread_str = ""

    if polish is True:
        polish_str = "--localpair --maxiterate 10"
    else:
        polish_str = ""

    if verbose is True:
        err_str = ""
    else:
        err_str = "2> {}".format(log)
    
    # generate mafft command and execute
    mafft_call = "mafft {0} {1} --merge {2} {3} > {4} {5}".format(thread_str, polish_str, table, cataln, temp_fasta, err_str)
    logging.info("\t{}".format(mafft_call))
    sp.call(mafft_call, shell=True)

    # reformat output file if alignment successful - will also notify if alignment failed
    reformat_mafft_output_fasta(temp_fasta, out_fasta)
        
    logging.info("Process complete!")
    logging.info("Please inspect Merged-alignment.fasta file available in: {}.\n\n".format(dir_merged))
        
def main():
    args = get_args()

    if args.output_seqnames is True:
        # turn input fasta file into dictionary structure (keys = seq names, vals = seqs)
        fasta_dict = fasta_to_dict(args.fasta)
        for name in sorted(fasta_dict.keys()):
            print(name)

    else:
        # setup log file and on-screen settings 
        setup_logging(args.outdir)

        # create output directories and return their paths
        dir_groups, dir_alns, dir_merged = make_output_paths(args.outdir)
        
        # turn input fasta file into dictionary structure (keys = seq names, vals = seqs)
        fasta_dict = fasta_to_dict(args.fasta)

        # turn mapping file into dictionary structure (keys = group names, vals = lists of seq names)
        map_dict = map_to_dict(args.map)

        # check whether all sequences in fasta file can be assigned to groups
        check_map_fasta_compatibility(fasta_dict, map_dict)

        # write the group-specific fasta files
        write_group_fastas(fasta_dict, map_dict, dir_groups)

        # align each group-specific fasta file, e.g. prepare sub-alignments
        align_individual_group_fastas(dir_groups, dir_alns, args.threads, args.algorithm)
        
        # prepare two required input files for mafft --merge based on sub-alignments produced, return paths to those input files
        table, cataln = prepare_subalignments_for_merge(dir_groups, dir_alns, dir_merged, args.algorithm)
        
        # run final merge algorithm in mafft
        run_merge(table, cataln, dir_merged, args.threads, args.polish, args.verbose)
    
if __name__ == '__main__':
    main()
