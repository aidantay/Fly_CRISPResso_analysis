#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports
import argparse
import os
import sys
import re
from pathlib import Path

# External imports
import typer
from typing import List
import numpy as np
import pandas as pd
from Bio.Seq import Seq

# Internal imports

#------------------- Constants ------------------------------#

app = typer.Typer()

#------------------- Public Classes & Functions -------------#

@app.command()
def freqToSam(
    i: str = typer.Argument(..., help='Allele frequency table (i.e., Alleles_frequency_table.txt)'),
    n: str = typer.Option('Reference_Sequence', help='Reference sequence name'),
    o: str = typer.Option('output.sam', help='Output file')):

    ## Read file
    freq_df = readAlleleFreq(i)

    ## Convert allele frequencies into SAM alignments
    freq_sam = toSAM(freq_df.copy(), n)

    ## Create DIR if they don't exist & write allele frequency SAM file
    p = Path(o)
    p.parent.mkdir(parents=True, exist_ok=True)
    freq_sam.to_csv(o, sep='\t', index=False, header=None)

@app.command()
def extractRef(
    i: str = typer.Argument(..., help='Allele frequency table (i.e., Alleles_frequency_table.txt)'),
    n: str = typer.Option('Reference_Sequence', help='Reference sequence name'),
    o: str = typer.Option('reference.fasta', help='Output file')):

    ## Read file
    freq_df = readAlleleFreq(i)

    ## Get the the reference sequence (nucleotide & amino acid)
    nt_ref_seq = getNtRef(freq_df)

    ## Create DIR if they don't exist & write reference FASTA file
    p = Path(o)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(o, 'w') as file:
        file.write('>{}\n'.format(n))
        file.write(nt_ref_seq)

@app.command()
def grnaToSam(
    i: str = typer.Argument(..., help='Allele frequency table (i.e., Alleles_frequency_table.txt)'),
    g: List[str] = typer.Argument(..., help='gRNAs'),
    n: str = typer.Option('Reference_Sequence', help='Reference sequence name'),
    o: str = typer.Option('gRNA.sam', help='Output file')):

    ## Read file
    freq_df = readAlleleFreq(i)

    ## Get the the reference sequence (nucleotide & amino acid)
    nt_ref_seq = getNtRef(freq_df)

    ## Convert gRNA sequences into SAM alignments
    grna_sam = convertGrnaToSam(g, nt_ref_seq, n)

    ## Create DIR if they don't exist & write gRNA SAM file
    p = Path(o)
    p.parent.mkdir(parents=True, exist_ok=True)
    grna_sam.to_csv(o, sep='\t', index=False, header=None)

@app.command()
def countStopCodons(
    i: str = typer.Argument(..., help='Allele frequency table (i.e., Alleles_frequency_table.txt)'),
    g: List[str] = typer.Argument(..., help='gRNAs'),
    o: str = typer.Option('stop_codon_counts.tsv', help='Output file')):

    ## Read file
    freq_df = readAlleleFreq(i)

    ## Get the the reference sequence (nucleotide & amino acid)
    nt_ref_seq = getNtRef(freq_df)

    ## Convert gRNA sequences into SAM alignments
    grna_sam = convertGrnaToSam(g, nt_ref_seq, 'Reference_Sequence')

    ## Find the position of stop codons in each alignment
    ## and count how many of them occur in each gRNA.
    stop_codon_count_df = getStopCodons(grna_sam, freq_df)

    ## Create DIR if they don't exist & Write stop codon counts
    p = Path(o)
    p.parent.mkdir(parents=True, exist_ok=True)
    stop_codon_count_df.to_csv(o, sep='\t', index=False)

@app.command()
def countmutation(
    i: str = typer.Argument(..., help='Allele frequency table (i.e., Alleles_frequency_table.txt)'),
    g: List[str] = typer.Argument(..., help='gRNAs'),
    o: str = typer.Option('mutation_counts.tsv', help='Output file')):

    ## Read file
    freq_df = readAlleleFreq(i)

    ## Get the the reference sequence (nucleotide & amino acid)
    nt_ref_seq = getNtRef(freq_df)

    ## Convert gRNA sequences into SAM alignments
    grna_sam = convertGrnaToSam(g, nt_ref_seq, 'Reference_Sequence')

    ## Count nucleotide coverage at each gRNA position
    mutation_count_df = getMutations(freq_df, grna_sam)

    ## Create DIR if they don't exist & Write stop codon counts
    p = Path(o)
    p.parent.mkdir(parents=True, exist_ok=True)
    mutation_count_df.to_csv(o, sep='\t', index=False)

@app.command()
def countframeshift(
    i: str = typer.Argument(..., help='Allele frequency table (i.e., Alleles_frequency_table.txt)'),
    q: str = typer.Argument(..., help='Query nucleotide sequence'),
    o: str = typer.Option('frameshift_counts.tsv', help='Output file')):

    ## Read file
    freq_df = readAlleleFreq(i)

    ## Find sequences containing frameshift mutations within the query sequence
    frameshift_df = getFrameshifts(freq_df, q)

    ## Create DIR if they don't exist & Write stop codon counts
    p = Path(o)
    p.parent.mkdir(parents=True, exist_ok=True)
    frameshift_df.to_csv(o, sep='\t', index=False)

def getFrameshifts(freq_df, q):
    def convertRefIdxToAlnIdx(seq, ref_start_idx, ref_end_idx):
        pos = {}
        ref_idx = 0
        for idx, c in enumerate(seq):
            pos[ref_idx] = idx
            if (seq[idx] == '-'):
                ref_idx -= 1

            ref_idx += 1

        aln_start_idx = pos[ref_start_idx]
        aln_end_idx   = pos[ref_end_idx]
        return pd.Series([aln_start_idx, aln_end_idx])

    ## Get the the reference sequence (nucleotide & amino acid)
    nt_ref_seq = getNtRef(freq_df)

    ## Get the start and end positions of the query in the reference
    nt_ref_start_idx = nt_ref_seq.find(q)
    nt_ref_end_idx   = nt_ref_start_idx + len(q) - 1

    ## Map the start and end positions of the reference to the alignment
    ## These positions can be different due to indels
    df = freq_df[['qName', 'Aligned_Sequence', 'Reference_Sequence']].copy()
    f = lambda x: convertRefIdxToAlnIdx(x['Reference_Sequence'], nt_ref_start_idx, nt_ref_end_idx)
    df[['aln_start_idx', 'aln_end_idx']] = df.apply(f, axis=1)

    ## Extract the query sequence
    df['ref_nt_seq'] = df.apply(lambda x: x['Reference_Sequence'][x['aln_start_idx']:x['aln_end_idx'] + 1], axis=1)
    df['aln_nt_seq'] = df.apply(lambda x: x['Aligned_Sequence'][x['aln_start_idx']:x['aln_end_idx'] + 1], axis=1)

    ## Count the number of indels in each query sequence
    df['indel_count'] = df['aln_nt_seq'].str.count('-') + df['ref_nt_seq'].str.count('-')

    ## Determine whether translating gives us amino acids in a different frame
    df['aln_aa_seq'] = df['aln_nt_seq'].str.replace('-', '').apply(getAa, frame=0)
    df['aln_nt_seqlen'] = df['aln_nt_seq'].str.replace('-', '').str.len()
    df['is_frameshift'] = (df['aln_nt_seqlen'] % 3 != len(q) % 3)

    ## Select columns
    df = df[['qName', 'aln_nt_seq', 'aln_aa_seq', 'aln_nt_seqlen', 'indel_count', 'is_frameshift']]
    return df

#------------------- Private Classes & Functions ------------#

def readAlleleFreq(f):
    freq_df = pd.read_csv(f, sep='\t')
    # freq_df = freq_df[freq_df['%Reads'].cumsum() <= 95] ## Keep reads within 95 percentile. Still quite a lot of reads...
    # freq_df = freq_df.head(100)
    # freq_df = freq_df[(freq_df['Aligned_Sequence'].str.contains('-')) & (freq_df['Reference_Sequence'].str.contains('-'))]
    # freq_df = freq_df.head(150).tail(20)
    freq_df = insertQNameCol(freq_df)
    return freq_df

def insertQNameCol(df):
    f_len = len(str(df.iloc[-1].name + 1))
    df['qName'] = ['Seq_{}'.format(str(n).zfill(f_len)) for n in df.index + 1]
    if ('%Reads' in df.columns):
        df['qName'] = ['{}_{}_{}'.format(a, b, c) for a,b,c in zip(df['qName'], df['#Reads'], df['%Reads'])]
    return df

def toSAM(df, rName):
    def getCigarStr(ref, qry):
        if (len(ref) != len(qry)):
            raise Exception('unequal length')

        cigar = []
        for i in range(len(ref)):
            r, q = ref[i], qry[i];
            if (r == '-' and q == '-'):
                raise Exception('both gaps')

            op = 'M' if r == q else 'I' if r == '-' else 'D' if q == '-' else 'X';
            if (len(cigar) > 0 and cigar[-1][1] is op): # add to the last operation
                cigar[-1][0] += 1

            else:
                cigar.append([1, op])              # a new operation

        return "".join(map(lambda x: str(x[0]) + x[1], cigar)); # turn to string

    df['flag']  = 0
    df['rName'] = rName
    df['pos']   = 1
    df['mapq']  = 60
    df['cigar'] = df.apply(lambda x: getCigarStr(x['Reference_Sequence'], x['Aligned_Sequence']), axis=1)
    df['rnext'] = '*'
    df['pnext'] = 0
    df['tlen']  = 0
    df['seq']   = df['Aligned_Sequence'].str.replace('-', '')
    df['qual']  = [''.join(['I'] * len(s)) for s in df['Aligned_Sequence'].str.replace('-', '')]

    ## Select columns
    cols = [
        'qName', 'flag', 'rName', 'pos', 'mapq', 'cigar',
        'rnext', 'pnext', 'tlen', 'seq', 'qual'
    ]
    df = df[cols].copy()
    return df

def getNtRef(freq_df):
    nt_ref_seq = freq_df['Reference_Sequence'].str.replace('-', '').unique()
    if (len(nt_ref_seq) != 1):
        raise Exception('More than 1 reference sequence')
    return nt_ref_seq[0]

def convertGrnaToSam(grnas, nt_ref_seq, rName):
    def grnaStrToSeries(nt_ref_seq, grna, i):
        pfxLen = nt_ref_seq.find(grna)
        pfx = ''.join(['-'] * pfxLen)
        sfxLen = len(nt_ref_seq) - (pfxLen + len(grna))
        sfx = ''.join(['-'] * sfxLen)
        a = pfx + grna + sfx

        d = {'Aligned_Sequence':a, 'Reference_Sequence': nt_ref_seq}
        s = pd.Series(d, name=i)
        return s

    grna_s = [grnaStrToSeries(nt_ref_seq, grna, i) for i, grna in enumerate(grnas)]
    grna_s = pd.DataFrame(grna_s)
    grna_s = insertQNameCol(grna_s)
    grna_sam = toSAM(grna_s, rName)
    grna_sam['pos']   = grna_sam['cigar'].str.extract('([0-9]+)[A-Z][0-9]+[A-Z][0-9]+[A-Z]')
    grna_sam['pos']   = grna_sam['pos'].astype(int) + 1
    grna_sam['cigar'] = grna_sam['cigar'].str.extract('[0-9]+[A-Z]([0-9]+[A-Z])[0-9]+[A-Z]')
    return grna_sam

def getStopCodons(grna_sam, freq_df):
    def insertQryStopIdxsCol(df):
        ## Idxs refer to the first base in the codon
        df['nt_qry_seq']       = df['Aligned_Sequence'].str.replace('-', '')
        df['aa_qry_seq']       = df['nt_qry_seq'].str.replace('-', '').apply(getAa, frame=1)
        df['nt_qry_stop_idxs'] = df['aa_qry_seq'].apply(getNtStopIdxs, frame=1).apply(tuple)
        df = df.drop(columns=['nt_qry_seq', 'aa_qry_seq'])
        return df

    def insertStopCodonsCols(df):
        df['nt_qry_seq']         = df['Aligned_Sequence'].str.replace('-', '')
        df['aa_qry_seq']         = df['nt_qry_seq'].str.replace('-', '').apply(getAa, frame=1)
        df['nt_qry_stop_codons'] = df.apply(lambda x: getNtStopCodon(x['aa_qry_seq'], x['nt_qry_seq'], frame=1), axis=1)
        df['n_qry_stop_codons']  = df['nt_qry_stop_codons'].apply(len)
        df = df.drop(columns=['nt_qry_seq', 'aa_qry_seq'])
        return df

    def explodeStopIdxsCol(df):
        df = df.explode('nt_qry_stop_idxs')
        df = df.rename(columns={'nt_qry_stop_idxs':'nt_qry_stop_idx'})

        ## If the sequence doesn't contain any stop codons, then we'll
        ## set the index to -1
        cond = (df['nt_qry_stop_idx'].isna())
        df.loc[cond, 'nt_qry_stop_idx'] = -1
        return df

    def convertQryStopIdxsToRefStopIdxs(df):
        ## Idxs for aligned qry and ref sequences are the same
        f = lambda x: toAlnStopIdx(x['Aligned_Sequence'], x['nt_qry_stop_idx'])
        df['nt_aln_qry_stop_idx'] = df.apply(f, axis=1)

        f = lambda x: toRefStopIdx(x['Reference_Sequence'], x['nt_aln_qry_stop_idx'])
        df['nt_ref_stop_idx'] = df.apply(f, axis=1)

        ## Adjust indices so that we start at 1 instead of 0
        df['nt_qry_stop_idx'] = df['nt_qry_stop_idx'] + 1
        df['nt_aln_qry_stop_idx'] = df['nt_aln_qry_stop_idx'] + 1
        df['nt_ref_stop_idx'] = df['nt_ref_stop_idx'] + 1
        return df

    def toAlnStopIdx(nt_aln_qry_seq, nt_qry_stop_idx):
        if (nt_qry_stop_idx == -1):
            return -1

        seq = np.array(list(nt_aln_qry_seq))
        nt_aln_qry_stop_idx = np.nan * np.ones(len(seq))
        nt_aln_qry_stop_idx[seq != '-'] = np.arange(np.sum(seq != '-'))
        nt_aln_qry_stop_idx = np.where(nt_aln_qry_stop_idx == nt_qry_stop_idx)[0][0]

        qry_prefix = nt_aln_qry_seq[:nt_aln_qry_stop_idx + 1]
        qry_n_gaps = qry_prefix.count('-')
        return nt_aln_qry_stop_idx

    def toRefStopIdx(nt_aln_ref_seq, nt_aln_ref_stop_idx):
        if (nt_aln_ref_stop_idx == -1):
            return -1

        ref_prefix = nt_aln_ref_seq[:nt_aln_ref_stop_idx + 1]
        ref_n_gaps = ref_prefix.count('-')
        nt_ref_stop_idx = nt_aln_ref_stop_idx - ref_n_gaps
        return nt_ref_stop_idx

    def getGrnaPositions(grna_sam):
        ## Find the start and end idxs of each gRNA
        grna_pos_df = grna_sam[['qName', 'seq', 'pos']]
        grna_pos_df = grna_pos_df.rename(columns={'pos':'pos_start'})
        grna_pos_df['pos_end'] = grna_pos_df['pos_start'] + grna_pos_df['seq'].str.len()
        grna_pos_df['pos'] = list(zip(grna_pos_df['pos_start'], grna_pos_df['pos_end']))
        grna_pos_df['qName'] = (
            grna_pos_df['seq'] 
            + '_' + grna_pos_df['pos_start'].astype(str)
            + '_' + grna_pos_df['pos_end'].astype(str)
        )

        ## Convert the df to a dict
        grna_pos = grna_pos_df[['qName', 'pos']]
        grna_pos = grna_pos.set_index('qName')['pos']
        grna_pos = grna_pos.to_dict()
        return grna_pos

    ## Find stop codons in each sequence
    df = freq_df[['qName', 'Aligned_Sequence', 'Reference_Sequence']].copy()
    df = insertQryStopIdxsCol(df)
    df = insertStopCodonsCols(df)
    df = explodeStopIdxsCol(df)

    ## Before we can match them, we need to convert the idxs
    ## from Qry -> AlnQry -> AlnRef -> Ref.
    ## Idxs in the Qry are shifted in the Ref due to insertions/deletions.
    df = convertQryStopIdxsToRefStopIdxs(df)

    ## Count how many stop codons occur in each gRNAs
    grna_pos = getGrnaPositions(grna_sam)
    for k in grna_pos.keys():
        df[k] = 0
        cond1 = (df['nt_ref_stop_idx'] >= grna_pos[k][0])
        cond2 = (df['nt_ref_stop_idx'] <= grna_pos[k][1])
        df.loc[cond1 & cond2, k] = 1

    ## Sum totals
    cols = ['qName', 'n_qry_stop_codons']
    agg  = {k:'sum' for k in grna_pos.keys()}
    agg['nt_ref_stop_idx'] = lambda x: tuple(x)
    agg['nt_qry_stop_idx'] = lambda x: tuple(x)
    df = df.groupby(cols).agg(agg)
    df = df.reset_index()

    ## Add/rename columns
    df['#Reads'] = df['qName'].str.split('_', expand=True)[2]
    df['%Reads'] = df['qName'].str.split('_', expand=True)[3]
    df = df.rename(columns={'n_qry_stop_codons':'#Stop_codons'})

    ## Select columns
    cols = ['qName', '#Reads', '%Reads', '#Stop_codons', 'nt_qry_stop_idx', 'nt_ref_stop_idx']
    cols = cols + list(grna_pos.keys())
    df = df[cols]
    return df

def getAa(nt_seq, frame=0):
    aa_seq = Seq(nt_seq[frame:]).translate()
    return str(aa_seq)

def getAaStopIdxs(aa_seq):
    aa_stop_idxs = [i for i, aa in enumerate(aa_seq) if aa == '*']
    return aa_stop_idxs

def getNtStopIdxs(aa_seq, frame=0):
    aa_stop_idxs = getAaStopIdxs(aa_seq)
    nt_stop_idxs = [frame + (i * 3) for i in aa_stop_idxs]
    return nt_stop_idxs

def getNtStopCodon(aa_seq, nt_seq, frame=0):
    nt_stop_idxs = getNtStopIdxs(aa_seq, frame=frame)
    nt_codons = [nt_seq[i:i+3] for i in nt_stop_idxs]
    return nt_codons

def getMutations(freq_df, grna_sam):
    def getGrnaIdxs(grna_sam):
        grna_pos = grna_sam[['pos', 'seq']].copy()
        grna_pos['min'] = grna_pos['pos']
        grna_pos['max'] = grna_pos['pos'] + grna_pos['seq'].str.len()
        grna_pos['range'] = grna_pos.apply(lambda x: list(range(x['min'], x['max'])), axis=1)
        grna_pos['seq'] = grna_pos['seq'].apply(list)
        grna_pos = grna_pos[['range', 'seq']].apply(lambda x: list(zip(x['range'], x['seq'])), axis=1)
        grna_pos = [v for l in grna_pos.values for v in l]
        return grna_pos

    ## Calculate the gRNA position
    grna_pos = getGrnaIdxs(grna_sam)

    ## Construct a multiple sequence alignment
    freq_sam = toSAM(freq_df.copy(), 'Reference_Sequence')
    freq_msa = toMSA(freq_sam)

    ## Count nucleotide coverage at each gRNA position
    df = freq_msa[[i[0] for i in grna_pos]].copy()

    ## This preserves the qName and displays the nucleotides in each read
    ## at each position, so that the actual number of reads can be calculated
    # df = df.reset_index()
    # df['#Reads'] = df['qName'].str.split('_', expand=True)[2]
    # df['%Reads'] = df['qName'].str.split('_', expand=True)[3]
    # cols = ['qName', '#Reads', '%Reads'] + [i[0] for i in grna_pos]
    # df = df[cols]
    # df = df.rename(columns={i[0]:'{}_{}'.format(str(i[0]), i[1]) for i in grna_pos})

    ## Rather than displaying the nucleotides in each read at each position,
    ## this cauclates the actual number of reads with a given nucleotide
    df = df.rename(columns={i[0]:'{}_{}'.format(str(i[0]), i[1]) for i in grna_pos})
    df = df.reset_index()
    df['#Reads'] = df['qName'].str.split('_', expand=True)[2]
    df['#Reads'] = df['#Reads'].astype(int)
    df = df.drop(columns=['qName'])

    mutation_count_df = [df.groupby(c)['#Reads'].sum() for c in df.columns if c != '#Reads']
    mutation_count_df = [c.rename(c.index.name) for c in mutation_count_df]
    mutation_count_df = pd.concat(mutation_count_df, axis=1).fillna(0).astype(int)
    mutation_count_df.index.name = 'Nucleotide'
    mutation_count_df = mutation_count_df.reset_index()
    df = mutation_count_df.copy()
    return df

def toMSA(sam_df):
    ## Remove insertions because they mess up the position of SNPs across
    ## reads in the MSA. However, we keep the deletions so that the sequence
    ## lengths are consistent across sequences
    sam_df['seq'] = sam_df.apply(lambda x: insertIndels(x['seq'], x['cigar']), axis=1)

    ## Construct a table. Each cell represnts a base at a particular
    ## position (column) in a particular sequence (row)
    df = sam_df['seq']
    df = pd.DataFrame(df.values.tolist())
    df.index = sam_df['qName']
    df.columns = [c + 1 for c in df.columns]
    return df

def insertIndels(seq, cigar):
    def splitCigarStr(cigar):
        cigar = re.split('([MDIX])', cigar)
        cigar = cigar[:-1]
        cLens = list(map(int, cigar[::2]))
        cOps  = cigar[1::2]
        return (cLens, cOps)

    (cLens, cOps) = splitCigarStr(cigar)

    seq = '@' + seq
    pos = 1
    seqParts = []
    for cLen, cOp in zip(cLens, cOps):
        ## Convert deletions into '-''
        if (cOp == 'D'):
            s = '-' * cLen
            pos = pos - cLen
            s = list(s)

        ## Convert insertions into '+''
        elif (cOp == 'I'):
            s = '+' + seq[pos:pos + cLen]
            s = [s]

        else:
            s = seq[pos:pos + cLen]
            s = list(s)

        seqParts = seqParts + s
        pos = pos + cLen

    idx = [i for i, x in enumerate(seqParts) if '+' in x]
    for i in reversed(idx):
        if (i + 1 >= len(seqParts)):
            pass

        else:
            m = seqParts.pop(i + 1)
            s = seqParts.pop(i)
            s = s + m
            seqParts.insert(i, s)

    return seqParts

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    app()

#------------------------------------------------------------------------------
