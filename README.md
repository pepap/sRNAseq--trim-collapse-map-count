# How do I trim the adapters, collapse the trimmed reads, map them and count the mature miRNAs expression

## <ins>Adapter trimming & read collapsing</ins>
The trimming was performed using the [cutAdapt](https://github.com/marcelm/cutadapt) software. The trimmed reads of the same sequence were collapsed using the [fastx_toolkit](https://github.com/agordon/fastx_toolkit) :
```
bash 01.cutAdapt.parallel.NEW.bash
```

## <ins>Mapping<ins>
Collapsed-fasta files were mapped to the DL genome using the [STAR](https://github.com/alexdobin/STAR/releases/tag/2.7.7a) software :
```
bash 02.map.bash
```

## <ins>Mature miRNA counting<ins>
The final counting is performed by a home-mmade R-script :
```
R --no-save < 03.counts.mirAnnot.mm10.R
```
