import sys
import pyskani
import Bio.SeqIO

database = pyskani.Database()

records = Bio.SeqIO.parse(str(sys.argv[2]).strip(), "fasta")
ref = (bytes(record.seq) for record in records)
database.sketch("ref", *ref)

query = Bio.SeqIO.parse(str(sys.argv[1]), "fasta")
query_contigs = (bytes(record.seq) for record in query)
hits = database.query("query", *query_contigs)


with open(str(sys.argv[4]).strip(), "a") as outfile:
    for hit in hits:
        print("\t".join([str(sys.argv[1]).strip(), str(sys.argv[2]).strip(), str(hit.identity).strip(), str(hit.query_fraction).strip(), str(hit.reference_fraction).strip()]), file = outfile)
