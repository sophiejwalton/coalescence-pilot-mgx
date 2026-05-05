import glob 
import pandas as pd
import argparse
from os import path,mkdir
from Bio import SeqIO
import glob
import io, os, stat, sys, resource, gzip, platform, bz2, Bio.SeqIO

#### FROM MIDAS 
def iopen(inpath, mode='r'):
	""" Open input file for reading regardless of compression [gzip, bzip] or python version """
	ext = inpath.split('.')[-1]
	# Python2
	if sys.version_info[0] == 2:
		if ext == 'gz': return gzip.open(inpath, mode)
		elif ext == 'bz2': return bz2.BZ2File(inpath, mode)
		else: return open(inpath, mode)
	# Python3
	elif sys.version_info[0] == 3:
		if ext == 'gz': return io.TextIOWrapper(gzip.open(inpath, mode))
		elif ext == 'bz2': return bz2.BZ2File(inpath, mode)
		else: return open(inpath, mode)

### ALTERED FORM MIDAS 
def build_genome_db(outdir, genomedir, species_list):
    """ Build FASTA and BT2 database of representative genomes """
    outfile = open('/'.join([outdir, 'genomes.fa']), 'w')
    db_stats = {'total_length':0, 'total_seqs':0, 'species':0}
    contig_names = []
    contig_lengths = []
    species_names = []
    for sp in species_list:
        db_stats['species'] += 1
        infile = iopen(sp)
        for r in Bio.SeqIO.parse(infile, 'fasta'):
            outfile.write('>%s\n%s\n' % (r.id, str(r.seq).upper()))
            db_stats['total_length'] += len(r.seq)
            db_stats['total_seqs'] += 1
            sp_name = sp.split('/')[-1][:-len('.gff.gz')]
            species_names.append(sp_name)
            contig_lengths.append(len(r.seq))
            contig_names.append(r.id)
        infile.close()
    outfile.close()
    print("  total genomes: %s" % db_stats['species'])
    print("  total contigs: %s" % db_stats['total_seqs'])
    print("  total base-pairs: %s" % db_stats['total_length'])
    return pd.DataFrame(data = {'species_id':species_names, 'contig_name': contig_names, 'contig_length':contig_lengths})

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Make genome fasta file')

    # add arguments
    parser.add_argument('--outdir', action='store',
                    help='where to store fasta file')
   # parser.add_argument('--genomedir', action='store',
    #                help='where genomes stored at')

    genome_dir = '/home/groups/bhgood/UHGG_download_SJW/'
    species_list = glob.glob(f'{genomedir}/*.gff.gz')
 #   genome_dir = 
    species_list = [i.split('/')]
    args = parser.parse_args()
    if not path.isdir(args.outdir):
        mkdir(args.outdir)

    db_info = build_genome_db(args.outdir, args.genomedir, species_list)
    db_info.to_csv(f'{args.outdir}/genomes_info.csv')
   # db_stats.to_csv(f'{args.outdir}/genomes_stats.csv')


