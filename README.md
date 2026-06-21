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
The final counting is performed by a home-made R-script :
```
R --no-save < 03.counts.mirAnnot.mm10.R
```


## *<ins>Reference<ins>*
 + **Cutadapt** : Martin M. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal. 2011;17(1):10. doi: 10.14806/ej.17.1.200
 + **FASTX-Toolkit** : https://github.com/agordon/fastx_toolkit
 + **STAR** : Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013 Jan 1;29(1):15-21. doi: 10.1093/bioinformatics/bts635. Epub 2012 Oct 25. PMID: 23104886; PMCID: PMC3530905.
 + **miRBase** : Kozomara A, Birgaoanu M, Griffiths-Jones S. miRBase: from microRNA sequences to function. Nucleic Acids Res. 2019 Jan 8;47(D1):D155-D162. doi: 10.1093/nar/gky1141. PMID: 30423142; PMCID: PMC6323917.

