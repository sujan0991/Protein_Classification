from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import pandas as pd
import os


# Create two sequence files
seq1 = SeqRecord(Seq("evlefddgtpatmsqvakdvctflrwaae"),
                   id="d1qcrd2")
seq2 = SeqRecord(Seq("sdlelhppsypwshrgllssldhtsirrgfqvykqvcsschsmdyvayrhlvgvcytedeakalaeevevqdgpnedgemfmrpgklsdyfpkpypnpeaaraanngalppdlsyivrarhggedyvfslltgycepptgvslreglyfnpyfpgqaigmappiynevlefddgtpatmsqvakdvctflrwaae"),
                   id="d1ppjd1")
SeqIO.write(seq1, "protein1.fasta", "fasta")
SeqIO.write(seq2, "protein2.fasta", "fasta")



## output directoy
fwd_out = os.path.join('/home/bioinf/mdho200b', 'fwd-results.tab')
rev_out = os.path.join('/home/bioinf/mdho200b', 'rev-results.tab')


# Run BLAST
fwd_output = NcbiblastpCommandline(query="protein1.fasta", subject="protein2.fasta",out =fwd_out, outfmt="\'6 qseqid sseqid pident qcovs qlen slen length bitscore evalue\'")()[0]
rev_output = NcbiblastpCommandline(query="protein2.fasta", subject="protein1.fasta",out =rev_out, outfmt="\'6 qseqid sseqid pident qcovs qlen slen length bitscore evalue\'")()[0]


fwd_results = pd.read_csv(fwd_out, sep="\t", header=None)
fwd_results.columns = ['query','subject','pident','qcovs','qlen','slen','length','bitscore','evalue']
rev_results = pd.read_csv(rev_out, sep="\t", header=None)
rev_results.columns = ['query','subject','pident','qcovs','qlen','slen','length','bitscore','evalue']


print(fwd_results)
print(rev_results)

##command line
##blastp -query seq1.fasta -subject seq2.fasta -out blast.csv