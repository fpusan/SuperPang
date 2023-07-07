from datetime import datetime
import resource


def read_fasta(fasta, ambigs = 'ignore', Ns = 'ignore', split_name = True):
    assert Ns in ('ignore', 'split')
    assert ambigs in ('ignore', 'as_Ns')
    seqDict = {}
    for seq in open(fasta).read().strip().lstrip('>').split('>'):
        name, seq = seq.split('\n',1)
        if split_name:
            name = name.split(' ')[0]
        seq = seq.upper().replace('\n','').replace('.','').replace('-','').replace('U','A')
        if ambigs == 'as_Ns': # just translate the
            seq = fix_Ns(seq)
        s = 0
        if name in seqDict:
            raise Exception(f'Sequence "{name}" is duplicated in your input file')
        if Ns == 'ignore':
            seqDict[name] = seq
        else:
            if 'N' not in seq:
                seqDict[name] = seq
            else:
                i = 0
                for seq in seq.split('N'):
                    if seq:
                        seqDict[f'{name}_Nsplit_{i}'] = seq
                        i += 1
    return seqDict


def read_fastq(fastq, ambigs = 'ignore'):
    assert ambigs in ('ignore', 'as_Ns')
    seqDict = {}
    with open(fastq) as infile:
        while True:
            name = infile.readline().strip().lstrip('@')
            seq = infile.readline().strip()
            if ambigs == 'as_Ns':
                seq = fix_Ns(seq)
            sep = infile.readline().strip()
            qual = infile.readline().strip()
            if not name:
                break
            if name in seqDict:
                raise Exception(f'Sequence "{name}" is duplicated in your input file')
            seqDict[name] = seq
    return seqDict


ambig2N = {ord(a): ord('N') for a in ('R', 'Y', 'K', 'M', 'S', 'W', 'B', 'D', 'H', 'V')} # ord bc str.translate wants ascii codes
allowedChars = {'A', 'C', 'T', 'G', 'N'}
def fix_Ns(seq):
    seq = seq.translate(ambig2N)
    if not set(seq).issubset(allowedChars):
        badChars = set(seq) - allowedChars
        raise ValueError(f'Illegal chars in seq "{name}": {badChars}')
    return seq


def write_fasta(seqDict, fasta):
    with open(fasta, 'w') as outfile:
        for name, seq in seqDict.items():
            outfile.write(f'>{name}\n{seq}\n')
        


def write_fastq(seqDict, fastq, qualDict = {}, default_qual = 'I'):
    with open(fastq, 'w') as outfile:
        for name, seq in seqDict.items():
            qual = qualDict[name] if qualDict else default_qual*len(seq)
            outfile.write(f'@{name}\n{seq}\n+\n{qual}\n')


def fasta2fastq(fasta, fastq):
    write_fastq(read_fasta(fasta), fastq)


def fastq2fasta(fastq, fasta):
    write_fasta(read_fastq(fastq), fasta)

base_for = 'ACGTN'
base_rev = 'TGCAN'
comp_tab = str.maketrans(base_for, base_rev)

def reverse_complement(seq):
    return seq.translate(comp_tab)[::-1]



def print_time(msg, end = None):
    dt = datetime.now().strftime('%d/%m/%Y %H:%M:%S')
    mem = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    print(f'{dt}\t{mem}MB\t{msg}', end = end)


