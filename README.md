# umi

Preprocessing tools for unique molecular index (UMI) sequencing reads

See the [wiki](https://github.com/aryeelab/umi/wiki) for documentation




# UCSF Francis Lab UMI Updates


demultiplex.py has been updated by SP to use a 3 column tab-delimited barcodes file with the columns `sampleid, i7barcode, i5barcode`.

```BASH
> head ../sampleindexes.txt 
1	GTCAACCA	GCATGCAA
2	TTCTGGTT	AGGCCTTG
3	CCGCAACG	CTCGAACC
4	CTAGGACC	TATATTGA
5	GCATGGTT	AACCAACG
6	AGGCCAAG	GCATTGGT
```


SP has also added a hamming distance check when matching barcodes to samples and uses the `--max_hamming` parameter to more flexibly assign reads to samples.

