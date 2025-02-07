import sys
import pyfastani
import Bio.SeqIO

sketch = pyfastani.Sketch()

ref = list(Bio.SeqIO.parse(str(sys.argv[2]).strip(), "fasta"))
sketch.add_draft("ref", (bytes(record.seq) for record in ref))

mapper = sketch.index()
query = Bio.SeqIO.parse(str(sys.argv[1]), "fasta")
hits = mapper.query_draft((bytes(record.seq) for record in query), int(sys.argv[3]))


with open(str(sys.argv[4]).strip(), "a") as outfile:
    for hit in hits:
        print("\t".join([str(sys.argv[1]).strip(), str(sys.argv[2]).strip(), str(hit.name).strip(), str(hit.identity).strip(), str(hit.matches).strip(), str(hit.fragments).strip()]), file = outfile)
