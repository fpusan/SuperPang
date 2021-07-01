from datetime import datetime
import resource


def read_fasta(fasta):
    seqDict = {}
    for seq in open(fasta).read().strip().lstrip('>').split('>'):
        name, seq = seq.split('\n',1)
        seq = seq.upper().replace('\n','').replace('.','').replace('-','').replace('N','')
        assert name not in seqDict
        seqDict[name] = seq
    return seqDict


def read_fastq(fastq):
    seqDict = {}
    with open(fastq) as infile:
        while True:
            name = infile.readline().strip().lstrip('@').split(' ')[0]
            seq = infile.readline().strip()
            sep = infile.readline().strip()
            qual = infile.readline().strip()
            if not name:
                break
            assert name not in seqDict
            seqDict[name] = seq
    return seqDict


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


def reverse_complement(seq):
    complementarity_matrix = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N', 
                              'W':'W', 'S':'S', 'R':'Y', 'Y':'R', 'M':'K', 
                              'K':'M', 'B':'V', 'V':'B', 'D':'H', 'H':'D',
                              '-':'-', '.':'.'}
    return ''.join([complementarity_matrix[b] for b in seq[::-1]])



def print_time(msg, end = None):
    dt = datetime.now().strftime('%d/%m/%Y %H:%M:%S')
    mem = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    print(f'{dt}\t{mem}MB\t{msg}', end = end)


