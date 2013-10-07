#!/usr/bin/env python
'''
Generate a parallelizable (either via GNU parallel or Nick Carriero's SimpleQueue) set of commands for gen_structure_coms

Usage:
    gen_structure_coms.py [options] -f FILE -l INT -i INT

'''
import sys
from os import path
from string import Template
import random
import argparse

MAINP = Template('''
KEY PARAMETERS FOR THE PROGRAM structure.  YOU WILL NEED TO SET THESE
IN ORDER TO RUN THE PROGRAM.  VARIOUS OPTIONS CAN BE ADJUSTED IN THE
FILE extraparams.

"(int)" means that this takes an integer value.
"(B)"   means that this variable is Boolean 
        (ie insert 1 for True, and 0 for False)
"(str)" means that this is a string (but not enclosed in quotes!) 

Basic Program Parameters

#define BURNIN    $b   // (int) length of burnin period
#define NUMREPS   $n   // (int) number of MCMC reps after burnin

Data file format
#define INFILE    $d // (str) in file name
#define NUMINDS    $i    // (int) number of diploid individuals in data file
#define NUMLOCI    $l    // (int) number of loci in data file
#define PLOIDY       2    // (int) ploidy of data
#define MISSING     -9    // (int) value given to missing genotype data
#define ONEROWPERIND 0    // (B) store data for individuals in a single line
#define LABEL     1     // (B) Input file contains individual labels
#define POPDATA   1     // (B) Input file contains a population identifier
#define POPFLAG   0     // (B) Input file contains a flag which says 
                              whether to use popinfo when USEPOPINFO==1
#define LOCDATA   0     // (B) Input file contains a location identifier
#define PHENOTYPE 0     // (B) Input file contains phenotype information
#define EXTRACOLS 0     // (int) Number of additional columns of data 
                             before the genotype data start.
#define MARKERNAMES      1  // (B) data file contains row of marker names
#define RECESSIVEALLELES 0  // (B) data file contains dominant markers (eg AFLPs)
                            // and a row to indicate which alleles are recessive
#define MAPDISTANCES     0  // (B) data file contains row of map distances 
                            // between loci

Advanced data file options

#define PHASED           0 // (B) Data are in correct phase (relevant for linkage model only)
#define PHASEINFO        0 // (B) the data for each individual contains a line
                                  indicating phase (linkage model)
#define MARKOVPHASE      0 // (B) the phase info follows a Markov model.
#define NOTAMBIGUOUS  -999 // (int) for use in some analyses of polyploid data
''')

EXTRAP = Template('''
EXTRA PARAMS FOR THE PROGRAM structure.  THESE PARAMETERS CONTROL HOW THE
PROGRAM RUNS.  ATTRIBUTES OF THE DATAFILE AS WELL AS K AND RUNLENGTH ARE 
SPECIFIED IN mainparams.

"(int)" means that this takes an integer value.
"(d)"   means that this is a double (ie, a Real number such as 3.14).
"(B)"   means that this variable is Boolean 
        (ie insert 1 for True, and 0 for False).

PROGRAM OPTIONS

#define NOADMIX     0 // (B) Use no admixture model (0=admixture model, 1=no-admix)
#define LINKAGE     0 // (B) Use the linkage model model 
#define USEPOPINFO  0 // (B) Use prior population information to pre-assign individuals
                             to clusters
#define LOCPRIOR    0 //(B)  Use location information to improve weak data

#define FREQSCORR   1 // (B) allele frequencies are correlated among pops
#define ONEFST      0 // (B) assume same value of Fst for all subpopulations.

#define INFERALPHA  1 // (B) Infer ALPHA (the admixture parameter)
#define POPALPHAS   0 // (B) Individual alpha for each population
#define ALPHA     1.0 // (d) Dirichlet parameter for degree of admixture 
                             (this is the initial value if INFERALPHA==1).

#define INFERLAMBDA 0 // (B) Infer LAMBDA (the allele frequencies parameter)
#define POPSPECIFICLAMBDA 0 //(B) infer a separate lambda for each pop 
                    (only if INFERLAMBDA=1).
#define LAMBDA    $lambd // (d) Dirichlet parameter for allele frequencies 

PRIORS

#define FPRIORMEAN 0.01 // (d) Prior mean and SD of Fst for pops.
#define FPRIORSD   0.05  // (d) The prior is a Gamma distribution with these parameters

#define UNIFPRIORALPHA 1 // (B) use a uniform prior for alpha;
                                otherwise gamma prior
#define ALPHAMAX     10.0 // (d) max value of alpha if uniform prior
#define ALPHAPRIORA   1.0 // (only if UNIFPRIORALPHA==0): alpha has a gamma 
                            prior with mean A*B, and 
#define ALPHAPRIORB   2.0 // variance A*B^2.  


#define LOG10RMIN     -4.0   //(d) Log10 of minimum allowed value of r under linkage model
#define LOG10RMAX      1.0   //(d) Log10 of maximum allowed value of r
#define LOG10RPROPSD   0.1   //(d) standard deviation of log r in update
#define LOG10RSTART   -2.0   //(d) initial value of log10 r
                         
USING PRIOR POPULATION INFO (USEPOPINFO)

#define GENSBACK    2  //(int) For use when inferring whether an indiv-
                         idual is an immigrant, or has an immigrant an-
                         cestor in the past GENSBACK generations.  eg, if 
                         GENSBACK==2, it tests for immigrant ancestry 
                         back to grandparents. 
#define MIGRPRIOR 0.01 //(d) prior prob that an individual is a migrant 
                             (used only when USEPOPINFO==1).  This should 
                             be small, eg 0.01 or 0.1.
#define PFROMPOPFLAGONLY 0 // (B) only use individuals with POPFLAG=1 to update P.
                                  This is to enable use of a reference set of 
                                  individuals for clustering additional "test" 
                                  individuals.

LOCPRIOR MODEL FOR USING LOCATION INFORMATION

#define LOCISPOP      1    //(B) use POPDATA for location information 
#define LOCPRIORINIT  1.0  //(d) initial value for r, the location prior
#define MAXLOCPRIOR  20.0  //(d) max allowed value for r

OUTPUT OPTIONS

#define PRINTNET     1 // (B) Print the "net nucleotide distance" to screen during the run
#define PRINTLAMBDA  1 // (B) Print current value(s) of lambda to screen
#define PRINTQSUM    1 // (B) Print summary of current population membership to screen

#define SITEBYSITE   0  // (B) whether or not to print site by site results. 
                       (Linkage model only) This is a large file!
#define PRINTQHAT    0  // (B) Q-hat printed to a separate file.  Turn this 
                           on before using STRAT.
#define UPDATEFREQ   100  // (int) frequency of printing update on the screen.
                                 Set automatically if this is 0.
#define PRINTLIKES   0  // (B) print current likelihood to screen every rep
#define INTERMEDSAVE 0  // (int) number of saves to file during run

#define ECHODATA     1  // (B) Print some of data file to screen to check
                              that the data entry is correct.
(NEXT 3 ARE FOR COLLECTING DISTRIBUTION OF Q:)
#define ANCESTDIST   0  // (B) collect data about the distribution of an-
                              cestry coefficients (Q) for each individual
#define NUMBOXES   1000 // (int) the distribution of Q values is stored as 
                              a histogram with this number of boxes. 
#define ANCESTPINT 0.90 // (d) the size of the displayed probability  
                              interval on Q (values between 0.0--1.0)

MISCELLANEOUS

#define COMPUTEPROB 1     // (B) Estimate the probability of the Data under 
                             the model.  This is used when choosing the 
                             best number of subpopulations.
#define ADMBURNIN  500    // (int) [only relevant for linkage model]: 
                             Initial period of burnin with admixture model (see Readme)
#define ALPHAPROPSD 0.025 // (d) SD of proposal for updating alpha
#define STARTATPOPINFO 0  // Use given populations as the initial condition 
                             for population origins.  (Need POPDATA==1).  It 
                             is assumed that the PopData in the input file 
                             are between 1 and k where k<=MAXPOPS.
#define RANDOMIZE      0  // (B) use new random seed for each run 
#define METROFREQ    10   // (int) Frequency of using Metropolis step to update
                             Q under admixture model (ie use the metr. move every
                             i steps).  If this is set to 0, it is never used.
                             (Proposal for each q^(i) sampled from prior.  The 
                             goal is to improve mixing for small alpha.)
#define REPORTHITRATE 0 //   (B) report hit rate if using METROFREQ

''')

def main():
    parser = argparse.ArgumentParser(description="Generate a parallelizable (either via GNU parallel or Nick Carriero's SimpleQueue) set of commands for gen_structure_coms")
    parser.add_argument("-b", "--burnin",
                        help = "Burnin to discard [default: 500000]",
                        type = int, default = 500000)
    parser.add_argument("-n", "--reps",
                        help = "Number of MCMC reps [default: 1000000]",
                        type = int, default = 1000000)
    parser.add_argument("-r", "--runs",
                        help = "Number of replicate runs of each K [default: 20]",
                        type = int, default = 20)
    parser.add_argument("--kstart", help = "First value of K to try [default: 1]",
                        type = int, default = 1)
    parser.add_argument("--kend", help = "Last value of K to try [default: 20]",
                        type = int, default = 20)
    parser.add_argument("--lambd", help = "Value for lambda [default: 1.0]",
                        type = float, default = 1.0)
    parser.add_argument("--pref", help="The prefix for results files [default: data file name without .txt]", default='')
    parser.add_argument("--wd", help = "Working directory [default: ./]",
                        default = './')
    parser.add_argument("--outdir", help = "Output directory [default: ./]",
                        default = './')
    parser.add_argument("--bindir", help = "Directory where structure can be found",
                        default = '')
    parser.add_argument("-f", "--datafile", help = "Structure data file",
                        required = True)
    parser.add_argument("-l", "--loci", help = "Number of loci",
                        required = True, type = int)
    parser.add_argument("-i", "--indiv", help = "Number of individuals ",
                        required = True, type = int)
    args = parser.parse_args()

    wd = path.abspath(args.wd)
    outd = path.abspath(args.outdir)
    if args.bindir != '':
        bind = path.abspath(args.bindir)
    else:
        bind = ''
    data = path.abspath(args.datafile)
    if args.pref == '':
            pref = path.basename(data)
    else:
        pref = args.pref

    mainparams = open(path.join(wd,'mainparams'), 'w')
    mainparams.write(MAINP.substitute(b = args.burnin,
                                      n = args.reps,
                                      d = data,
                                      i = args.indiv,
                                      l = args.loci))
    mainparams.close()
    extrap = open(path.join(wd,'extraparams'), 'w')
    extrap.write(EXTRAP.substitute(lambd = args.lambd))
    extrap.close()

    kstart = args.kstart
    kend = args.kend+1
    runs = args.runs
    total = (kend - kstart) * runs

    seeds = random.sample(range(100000*total), total)
    snum = 0
    for r in range(1,runs+1):
        for k in range(kstart, kend):
            dname = data
            out = path.join(outd, '{0}_K{1}_run{2}'.format(
                            pref, k, r))
            print('cd {0}; {1}/structure -K {2} -o {3} -D {4} > {3}.log 2>&1'.format(
                  wd, bind, k, out, seeds[snum]))
            snum += 1

if __name__ == '__main__':
    main()
