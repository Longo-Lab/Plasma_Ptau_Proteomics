#!/usr/bin/env python
#!interpreter [optional-arg]
# -*- coding: utf-8 -*-

"""
Pipeline for generating a gene-set analysis using MAGMA a set of gene lists.
Runs with a range of annotation windows surrounding the gene, and can
incorporate gene-set covariate files as defined by magma
"""

# Futures
from __future__ import print_function
# […]

# Built-in/Generic Imports
import os
import sys
# […]

# Libs
import logging
import argparse
import datetime
import subprocess as sp
# […]

# Own modules
# […]

# global variables
__author__ = 'Robert R Butler III, William A. Johnson'
__copyright__ = 'Copyright 2023, Longo Lab'
__version__ = '0.0.12'
__maintainer__ = 'Robert R Butler III'
__email__ = 'rrbutler@stanford.edu'


def getargs():
    """
    argparse function

    Returns
    -------
    Parsed arguments for use in the rest of the script
    """
    parser = argparse.ArgumentParser(prog='magma_geneset.py',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     description='')
    # parser.add_argument('--log', action='store_true', help='create logfile')
    parser.add_argument('--threads', default=1, type=int,
                        help='number of threads. If 0, threads is the '
                             'number of GWAS')
    parser.add_argument('--long_log', action='store_true',
                        help='create detailed logfile')
    parser.add_argument('--grch37', action='store_true',
                        help='Use GRCh37 snps and genes. Requires:\n'
                             '  *_magma.GRCh37.sumstats for each GWAS\n'
                             '  ENSG_coord_vXXX.GRCh37.txt for gene-locs')
    parser.add_argument('-o', metavar='log prefix', dest='outprefix',
                        default='magma_geneset',
                        help='choose alternate prefix for log files')
    parser.add_argument('-c', '--covars', nargs='*',
                        help='list of covariates to condition on \n'
                             '  (in a batch_cov.geneset file)')
    parser.add_argument('--interaction', action='store_true',
                        help='run interaction-pairs on genesets')
    parser.add_argument('--gene_model', default='snp-wise=mean',
                        choices=['snp-wise=mean', 'snp-wise=top',
                                 'snp-wise=top,P', 'multi=snp-wise'],
                        help='model to run for gene level analysis')
    parser.add_argument('--geneset-only', action='store_true',
                        help='skip annotation and gene level analysis?\n'
                             '  (requires previous run in folder)')
    parser.add_argument("--version", action='version',
                        version='\n'.join([
                            'Magma_geneset v' + __version__, __doc__
                        ]))
    parser.add_argument('-a', "--annot_prefix",
                        help="Filename prefix for custom annot files.")
    parser.add_argument('-p', '--protein_coding', action='store_true',
                        help="Only assess protein coding genes.")
    parser.add_argument('--gwas-dir', metavar='gwas directory', dest='gwas_dir',
                        default='/labs/flongo/reference/gwas_sumstats',
                        help='directory containing gwas files')
    parser.add_argument('--gref', metavar='genome reference', dest='gref',
                        default='/labs/flongo/reference/bed/g1000_eur/g1000_eur',
                        help='filepath/fileprefix of reference plink files '
                             'default is the 1000G EUR genome reference')
    parser.add_argument('--ensembl-ver', metavar='Enseml version', dest='ensver',
                        default='110',
                        help='Ensembl database build version (for gene IDs)')
    # parser.add_argument("input", metavar=('file'), nargs='+',
    # help="input file(s) (returns outfile for each)")
    required_named = parser.add_argument_group('required arguments')
    required_named.add_argument('-u', '--up', required=True,
                                help='upstream annotation window')
    required_named.add_argument('-d', '--down', required=True,
                                help='downstream annotation window')
    required_named.add_argument('-b', '--batches', nargs='*',
                                default=['magma'],
                                help='list of batch prefixes on geneset files\n'
                                     'of the format batch_set.geneset or'
                                     ' batch_cov.geneset')
    required_named.add_argument('-g', '--gwas', nargs='*',
                                help='list of GWAS datasets to analyze files\n'
                                     'located in lab_flongo/reference/gwas_sumstats')
    return parser.parse_args()


def error_handling():
    """
    config for error reporting

    Returns
    -------
    str
        Human-readable error output
    """
    return ' {}. {}, line: {}'.format(sys.exc_info()[0], sys.exc_info()[1],
                                      sys.exc_info()[2].tb_lineno)


def log_args(long_log, outprefix):
    """
    Choosing log options from argparse

    Parameters
    ----------
    long_log : bool
        True writes all levels of log events to file
        False only writes INFO and higher level events to file

    outprefix : str
        prefix for output filename

    Returns
    -------
    None
        Writes logged errors to file defined by {outprefix}.log
    """
    if long_log:
        logging.basicConfig(level=logging.DEBUG,
                            filename="{}.log".format(outprefix))
    else:
        logging.basicConfig(level=logging.INFO,
                            filename="{}.log".format(outprefix))


def infiles(batch):
    """
    Reads in a list of single column geneset files by prefix and converts to dict

    Parameters
    ----------
    batch : str
        prefix for geneset filenames: {batch}_{geneset}.geneset

    Returns
    -------
    Dict
        keys are the {geneset} name, and values are list[str] of gene_ids
    """
    fdict = {}
    for r, d, files in os.walk('.'):
        for file in files:
            prefix = batch + '_'
            if file.startswith(prefix) and file.endswith('.geneset'):
                setname = os.path.splitext(file[len(prefix):])[0]
                setlist = open(os.path.join(r, file), 'r').read().splitlines()
                fdict.update({setname: setlist})
    return fdict


def print_sideways(name, a_dict):
    """
    Writes a batch of geneset files as a single set-annot file for MAGMA 
    with genesets horizontally per row

    Parameters
    ----------
    name : str
        batch filename prefix for geneset files

    a_dict : dict
        dictionary of genesets to be written to file w/ keys set names and
        values as lists[str] of genes

    Returns
    -------
    None
        writes file with filename {name}.geneset
    """
    logging.info(
        '{} gene sets: {} containing {} genes'.format(
            name, len(a_dict), [len(x) for x in a_dict.values()]
        )
    )
    with open('{}.geneset'.format(name), 'w') as out:
        out.writelines(
            ('\t'.join([k] + v) + '\n' for k, v in a_dict.items())
        )
    return


def generate_geneset(batch):
    """
    Generates MAGMA's set-annot file given a batchname prefix for genelist files
    Calls infiles() to read and print_sideways() to write

    Parameters
    ----------
    batch : str
        batchname provided in script call

    Returns
    -------
    None
        writes {batch}.geneset file for the given batchname
    """
    fdict = {}
    try:
        fdict = infiles(batch)
        # We don't know what the below two lines are trying to catch...
        #     commented out for now.
        # logging.fatal(error_handling())
        # print(f'Unable to open file, {error_handling()}')
        assert (all(len(x) >= 10 for x in fdict.values()))
        print_sideways(batch, fdict)
        return

    except AssertionError:
        logging.warning(error_handling())
        for k, v in fdict.items():
            if len(v) < 10:
                print(f'{k} in {batch} set too small, skipping '
                      f'{error_handling()}')
                del fdict[k]
        print_sideways(batch, fdict)
        return

    except OSError:
        logging.fatal(error_handling())
        print('Unable to open file, {}'.format(error_handling()))

    except BaseException:
        logging.fatal(error_handling())


def get_37(g37, protein_coding, ensver, gwas):
    """
    getting grch37 loci

    Parameters
    ----------
    g37 : bool
        True if geneset comes from g37 ENSG
        False if geneset comes from g38 ENSG

    protein_coding : bool
        True if geneset comes from Protein Coding Only set
        False if geneset does not come from Protein Coding Only set

    ensver : str
        String specifying version of the Ensembl database used for gene names

    gwas : list[str]
        List of string prefixes for gwas sets to test

    Returns
    -------
    list[str] of two:
        list[str] of formatted gwas (from -g, --gwas) with appropriate 
            filename suffix and extension.
        str filename of Ensembl gene names to use for analysis.
    """
    if g37:
        if protein_coding:
            return [
                       x + '_magma.GRCh37' for x in gwas
                   ], 'ENSG_coord_vGRCh37_PC.GRCh37.filtered.txt'
        else:
            return [
                       x + '_magma.GRCh37' for x in gwas
                   ], 'ENSG_coord_vGRCh37.GRCh37.filtered.txt'
    else:
        if protein_coding:
            return [
                       x + '_magma' for x in gwas
                   ], 'ENSG_coord_v' + ensver + '_PC.GRCh38.filtered.txt'
        else:
            return [
                       x + '_magma' for x in gwas
                   ], 'ENSG_coord_v' + ensver + '.GRCh38.filtered.txt'


def geneset_model(covars, interaction):
    """
    Generate non-standard models for covariates and interactions

    Parameters
    ----------
    covars : list[str]
        List of covariates for MAGMA to incorporate into the analysis

    interaction : bool
        True include the interaction-pairs setting in the model
        False no interactions

    Returns
    -------
    list[str]
        list of two strings to add to the MAGMA command to specify the model
    """
    if covars and not interaction:
        covar_str = ','.join(covars)
        return ['--model', '='.join(['condition-hide', covar_str])]
    elif covars and interaction:
        covar_str = ','.join(covars)
        return ['--model', '='.join(['condition-hide', covar_str]),
                'interaction-pairs']
    elif not covars and interaction:
        return ['--model', 'interaction-pairs']
    else:
        return None


def gene_annot(up, down, outfile, sumstats, genloc, thr_string, gwas_dir):
    """ MAGMA Step 1
    Build gene_annotation files if annot_prefix argument is not used

    Parameters
    ----------
    up : int
        Length of upstream window (0 for hMAGMA)

    down : int
        Length of downstream window (0 for hMAGMA)

    outfile : list[str]
        List of output file prefixes (one per gwas)

    sumstats : list[str]
        List of sumstats filename prefixes for each gwas

    genloc : str
        Filename of Ensembl gene names to use for analysis

    thr_string : str
        Number of GNU Parallel threads

    gwas_dir : str
        Filepath for the location of gwas sumstats, genloc file

    Returns
    -------
    None
        MAGMA writes annotation file for each gwas as: 
            ${gwas}_${up}_${down}.genes.annot
    """
    cmd_list = ['parallel', thr_string, '--link', 'magma',
                '--annotate', ''.join(['window=', up, ',', down]),
                '--snp-loc', '/'.join([gwas_dir, '{2}.sumstats']),
                '--gene-loc', '/'.join([gwas_dir, genloc]),
                '--out', '{1}',
                ':::', *outfile,
                ':::', *sumstats]
    logging.debug('Command to run: {}'.format(cmd_list))
    cmd = sp.Popen(cmd_list, stdout=sp.PIPE, stderr=sp.STDOUT)
    stdout, stderr = cmd.communicate()
    logging.info(stdout.decode())
    return


def gene_level(outfile, sumstats, model, thr_string, gref, gwas_dir):
    """ MAGMA Step 2
    Gene-level analysis on SNP summary statistics (p-values).

    Parameters
    ----------
    outfile : list[str]
        List of strings serving as output file prefixes (one per gwas)

    sumstats : list[str]
        List of sumstats filename prefixes for each gwas

    model : str
        Model type to be used with (h)MAGMA (see MAGMA manual)

    thr_string : str
        String specifying number of GNU Parallel threads

    gref : str
        String specifying the filepath and prefix of the plink genome reference

    gwas_dir : str
        String filepath for the location of gwas sumstats, genloc file

    Returns
    -------
    None
        MAGMA writes various outfiles for (h)MAGMA gene-level analysis
    """
    cmd_list = ['parallel', thr_string, '--link', 'magma',
                '--gene-model', model,
                '--bfile', gref,
                '--gene-annot', '{1}.genes.annot',
                # '--gene-settings adap-permp',
                '--pval', '/'.join([gwas_dir, '{2}.sumstats']),
                'ncol=N',
                '--out', '{1}',
                ':::', *outfile,
                ':::', *sumstats]
    # run
    logging.debug('Command to run: {}'.format(cmd_list))
    cmd = sp.Popen(cmd_list, stdout=sp.PIPE, stderr=sp.STDOUT)
    stdout, stderr = cmd.communicate()
    logging.info(stdout.decode())
    return


def gene_set(outfile, batches, covars, interaction, thr_string):
    """ MAGMA Step 3
    MAGMA Gene-set analysis

    Parameters
    ----------
    outfile : list[str]
        List output file prefixes (one per gwas)

    batches : list[str]
        List of batchnames provided in script call (args.batches)

    covars : list[str]
        List of covariates for MAGMA to incorporate into the analysis

    interaction : bool
        True include the interaction-pairs setting in the model
        False no interactions

    thr_string : str
        Number of GNU Parallel threads

    Returns
    -------
    None
        MAGMA writes various outfiles for (h)MAGMA gene-set analysis
    """
    cmd_list = ['parallel', thr_string, 'magma',
                '--gene-results', '{1}.genes.raw',
                '--set-annot', '{2}.geneset',
                '--out', '{1}_{2}',
                ':::', *outfile,
                ':::', *batches]
    # check if model
    model = geneset_model(covars, interaction)
    if model:
        cmd_list[3:3] = model
    # run
    logging.debug('Command to run: {}'.format(cmd_list))
    cmd = sp.Popen(cmd_list, stdout=sp.PIPE, stderr=sp.STDOUT)
    stdout, stderr = cmd.communicate()
    logging.info(stdout.decode())
    return


def main():
    # setup args and logging
    args = getargs()
    log_args(args.long_log, args.outprefix)
    logging.info('Run date: {}'
                 .format(datetime.date.today().strftime('%Y-%m-%d')))
    logging.info('CLI inputs are up: {} down: {} outprefix: {} batches: {} '
                 'covars: {} gwas: {} annots: {} ensver: {}'
                 .format(args.up, args.down, args.outprefix, args.batches,
                         args.covars, args.gwas, args.annot_prefix, args.ensver))
    # defining reference data
    gwas_dir, gref = args.gwas_dir, args.gref
    logging.info('References for magma are in gwas: {} genome reference: {}'
                 .format(gwas_dir, gref))
    # configuring for reference data
    sumstats, genloc = get_37(args.grch37, args.protein_coding, args.ensver, args.gwas)
    logging.debug('sumstats list: {}'.format(sumstats))
    logging.debug('gene-loc file {}'.format(genloc))
    # set threading
    if args.threads == 0:
        thr_string = '-j' + str(len(args.gwas))
    else:
        thr_string = '-j' + str(args.threads)
    logging.debug('thread string: {}'.format(thr_string))
    # define output filenames
    if args.annot_prefix is None:
        outfile = ['_'.join([x, args.up, args.down]) for x in args.gwas]
    else:
        outfile = ['_'.join([x, args.annot_prefix]) for x in args.gwas]
    logging.debug('outfile list: {}'.format(outfile))
    # make genesets
    batches = args.batches
    [generate_geneset(batch) for batch in batches]
    # run analyses
    if args.geneset_only:
        pass
    elif args.annot_prefix:
        gene_level(outfile, sumstats, args.gene_model, thr_string, gref, gwas_dir)
    else:
        gene_annot(args.up, args.down, outfile, sumstats, genloc, thr_string, gwas_dir)
        gene_level(outfile, sumstats, args.gene_model, thr_string, gref, gwas_dir)
    gene_set(outfile, batches, args.covars, args.interaction, thr_string)


if __name__ == "__main__":
    main()
