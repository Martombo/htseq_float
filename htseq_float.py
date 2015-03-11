import HTSeq
import collections

gtf_file = HTSeq.GFF_Reader( 'prova.gtf' )
exons = HTSeq.GenomicArrayOfSets( 'auto', stranded=True )

for feature in gtf_file:
    if feature.type == 'exon':
        exons[ feature.iv ] += feature.attr['gene_id']

counts = collections.Counter( )

almnt_file = HTSeq.SAM_Reader( 'prova.sam' )
for almnt in almnt_file:
   if not almnt.aligned:
      counts[ '_unmapped' ] += 1
      continue

gene_ids = set()
for cigop in almnt.cigar:
   if cigop.type != "M":
      continue
   for iv, val in exons[ cigop.ref_iv ].steps():
      gene_ids |= val
   if len(gene_ids) == 1:
      gene_id = list(gene_ids)[0]
      counts[ gene_id ] += 1
   elif len(gene_ids) == 0:
      counts[ '_no_feature' ] += 1
   else:
      counts[ '_ambiguous' ] += 1

for gene_id in counts:
   print gene_id, counts[ gene_id ]
