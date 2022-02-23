This is a general outline of the procedure used to calculate the CMI between chromosome 8 gene expression values and 8q24.3 region copy number variants. The selection of the data is done with the bash commands in linux and with Libre Office, the calculations with the R programming language and in particular with the infotheo library. We must have three fundamental files, the biomart file, which gives characteristics of the genes of the human genome, called here biomart.tsv; the gene expression file for luminal cancer variant B, called here luminalB.tsv, and the file with copy number variants of the 8q24.3 region, called here cnvs-lumB.tsv, which actually contains genome-wide cnvs.

The biomart file, biomart.tsv, contains nine columns that respectively have the Ensemble identifier, the gene name, the GC percentage, the gene type, the band it belongs to, the chromosome, the gene description, the start base and the end base.

1) The first thing is to select the genes of chromosome 8 with the following instruction:
```
awk '$6 == 8 {print $0}' biomart.tsv
```
To indicate that the field separator is the tabulator, use:
```
awk 'BEGIN {FS="\t"}; $6 == 8 {print $0}' biomart.tsv
```
From the previous result we are left with those that in column 4 have the type 'protein_coding'. Together these two things can be done like this:
```
awk '$6 == 8 {print $0}' | grep protein_coding > chr8CodingGenes.tsv
```
The final result, the coding genes for chromosome 8, are stored in the file chr8CodingGenes.tsv

2) Then the chr8CodingGenes.tsv file is sorted in ascending order according to column 8. This is simple in Libre Office, you have to extend the sort to all the columns of the file.

3) We obtain the Ensemble identifiers of the coding genes of chromosome 8 in ascending order.
```
awk '{print $1}' chr8CodingGenes.tsv > geneIDsChr8.txt
```
4) We have, in addition to the biomart.tsv file, the expression file for luminalB cancer type, luminalB.tsv. From this file we take all those genes contained in geneIDsChr8.txt in the order in which they appear. We do it with this bash script:

```
#!/bin/bash

file1='geneIDsChr8.txt'
file2='luminalB.tsv'

while read r;
do
grep $r $file2;
done < $file1
```

We save the result in:

geneExpressionCr8LumB.tsv

5) By a procedure similar to the previous one, we find the values ​​of the cnvs of the 8q24.3 region:
```
grep 8q24.3 crh8CodingGenes.tsv > genes8q24-3.txt
```
and we order the last file with Libre Office.

6) We keep only the first column of this file to get the Ensembl identifiers.
```
awk '{print $1}' genes8q24-3.txt > geneIDs8q24-3.txt
```
7) With the file geneIDs8q24-3.txt and the file containing the CNVs, cnvs-lumB.tsv, you can find the file with the cnvs for the ordered genes of the region
8q24.3. We make use of the script used in step 4. We obtain a file that we will call geneCNVs8q24-3.tsv

8) The two important files to perform our calculations are the one obtained in step 4) geneExpressionCr8LumB.tsv and the one obtained in step 7) geneCNVs8q24-3.tsv
with them it is possible to use the program to calculate the conditional mutual information CMI.R. An important point is that both must have the same number of columns. In case we want, say, 200 columns, we use the following instruction:
```
awk '{for(i=1;i<=200;i++) printf $i" "; print " "}' file.tsv
```
where file.tsv is the file to trim.


9) From the CMI.R program you get a file called IMC-lumB.Rcompatible, IMC-lumB.rds is a layer of it.

10) With IMC-lumB.Rcompatible and the circosPlot.R program, the ciros plot are created.
