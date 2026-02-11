# How do I trim the adapters, collapse the trimmed reads, map them and count the mature miRNAs expression

## <ins>Adapter trimming & read collapsing</ins>
The trimming was performed using the cutAdapt software. The trimmed reads of the same sequence were collapsed using the fastx_toolkit :
```
bash 01.cutAdapt.parallel.NEW.bash
```

## <ins>Mapping<ins>
Collapsed-fasta files were mapped to the DL genome using the [STAR](https://github.com/alexdobin/STAR/releases/tag/2.7.7a) program :
```
bash 02.map.bash
```

## <ins>Mature miRNA counting<ins>
```
R --no-save < 03.counts.mirAnnot.mm10.R
```
