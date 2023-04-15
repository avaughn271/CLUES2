import sys
import os
import argparse
import logging
import msprime
from mcmc import *
'''
command line interface for arginfer
'''
logger = logging.getLogger(__name__)
log_format = "%(asctime)s %(levelname)s %(message)s"

def setup_logging(args):
    log_level = "WARN"
    if args.verbose:
        log_level = "DEBUG"#"INFO"
    logging.basicConfig(level=log_level, format=log_format)

def arginfer_cli_parser():
    high_parser = argparse.ArgumentParser(prog="arginfer",
        description="This is the command line interface for arginfer, "
                    "a probabilistic method to infer the Ancestral Recombination Graph.")
    subparsers = high_parser.add_subparsers(dest="subcommand")
    subparsers.required = True
    parser = subparsers.add_parser(
        "infer",
        help=(
            "Takes the data or the ARG in tree sequence full_ARG format and "
            "returns MCMC sampled ARGs."
        ),
    )
    parser.add_argument('--tsfull', type=argparse.FileType('r', encoding='UTF-8'), default=None,
                                            help='an msprime .args file.'
                                                 ' If None, build an ARG from haplotype data')
    parser.add_argument('--input_path',type=str,
                        default=os.getcwd()+"/data", help='The path to input data, '
                                    'this is the path to haplotype, ancestral allele, and snp_pos ')
    parser.add_argument('--haplotype_name' , type = str,
                        default= None, help='the haplotype file name',#"haplotype_ready.txt"
                        required=False)
    parser.add_argument('--ancAllele_name' , type = str,
                        default= "ancestral_allele_ready.txt",
                        help='a txt file of ancestral allele for each snp',
                        required=False)
    parser.add_argument('--snpPos_name' , type = str,
                        default= "ancestral_allele_ready.txt",
                        help='a txt file of SNP chrom position',
                        required=False)
    parser.add_argument('--iteration','-I', type=int, default=20,
                        help= 'the number of mcmc iterations')
    parser.add_argument('--thin', type=int, default= 10, help=' thining steps')
    parser.add_argument('--burn', '-b', type=int, default= 0, help=' The burn-in')
    parser.add_argument('--sample_size', '-n', type=int, default= 5, help=' sample size')
    parser.add_argument('--seq_length','-L', type=float, default=1e4,help='sequence length')
    parser.add_argument('--Ne', type=int, default= 5000, help=' effective population size')
    parser.add_argument('--recombination_rate', '-r', type=float, default=1e-8,
                        help=' the recombination rate per site per generation ')
    parser.add_argument('--mutation_rate', '-mu', type=float, default=1e-8,
                        help='the mutation rate per site per generation')
    parser.add_argument('--outpath', '-O',type=str,
                        default=os.getcwd()+"/output", help='The output path')
    parser.add_argument( '-p','--plot', help="plot the output", action="store_true")
    parser.add_argument("--random-seed", "-s", type = int, default=1)
    parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
    parser.set_defaults(runner=run_mcmc)
    #if you need any other subparsers, they are added here
    return high_parser

def run_mcmc(args):
    input_data_path = args.input_path
    haplotype_data_name = args.haplotype_name
    ancAllele_data_name = args.ancAllele_name
    snpPos_data_name= args.snpPos_name
    iteration = args.iteration
    thin = args.thin
    burn = args.burn
    n = args.sample_size
    seq_length = args.seq_length
    mu = args.mutation_rate
    r= args.recombination_rate
    Ne= args.Ne
    outpath = args.outpath
    tsfull = None
    if args.tsfull !=None:#else real data
        try:
            tsfull = msprime.load(args.tsfull.name) #trees is a fh
        except AttributeError:
            tsfull = msprime.load(args.tsfull)
    mcmc = MCMC(tsfull, n, Ne, seq_length, mu, r,
                 input_data_path,
                 haplotype_data_name,
                 ancAllele_data_name,
                 snpPos_data_name, outpath, args.verbose)
    mcmc.run(iteration, thin, burn)
    if args.plot:
        p= Trace(outpath)
        p.arginfer_trace()
    if args.verbose:
        mcmc.print_state()

def arginfer_main(arg_list=None):
    parser = arginfer_cli_parser()
    args = parser.parse_args(arg_list)
    setup_logging(args)
    args.runner(args)
