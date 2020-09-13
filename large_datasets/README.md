## Large genomics files

**Need to add DOI links once the data are published**

To run the scripts in this repository, download the following files and add
them to this folder.  These are all from a set of _Miscanthus sacchariflorus_
variant calls available at [Illinois Data Bank](https://doi.org).

`180208Msa_filtered.vcf.bgz` and `180208Msa_filtered.vcf.bgz.tbi`, containing
allelic read depths from SNPs called in _M. sacchariflorus_.

`180110alignedtags.sam`, containing
alignment locations of sequence tags from the above _M. sacchariflorus_ dataset.

`180813tagtaxadist.txt`, containing read depth for each tag in each sample.

Additional files from a separate repository on
[Illinois Data Bank](https://doi.org)
specific to the testing of the Hind/He statistic:

`Msi_split_1_align.csv` and `Msi_split_1_depths.csv`, containing read depths
and alignment locations output by `process_sam_multi.py` in the polyRAD
variant calling pipeline, used by `scripts/get_inbreeding.R`.
