# CapiceQuickFilter
A simple tool to prioritize candidate pathogenic variants using mainly CAPICE 
scores (see: https://www.medrxiv.org/content/10.1101/19012229v1). Report 
variants of potential clinical interest for a specific case sample identifier 
in a VCF, while allowing for additional control samples to drop non-relevant 
case genotypes.

## Features
- Simple to operate with only a few command-line arguments.
- Allows any number of control samples to filter out non-relevant genotypes.
- Drops 99.9% tot 99.99% of variants even under relatively sensitive settings.
- Output is (very close to proper) VCF but suitable for human interpretation.
- It doesn't take ages to run (~20,000 variants/sec).

## Limitations
- For now, genome build GRCh37/hg19 only.
- Requires input VCF to be pre-annotated with CAPICE and GnomAD, which may 
not always be a trivial task. See details below.
- Popular VCF annotation tools don't always report missing annotation values 
(e.g. CAPICE, GnomAD) for multi-allelic variants correctly, we therefore cannot 
rely on  exact allele matching of these criteria. A more relaxed strategy is 
 applied, matching any allele annotation, leading to overdetection.
- Similarly, due to complexity of matching multi-allelic variant genotypes, 
any non-reference, non-missing genotype is considered as alternative genotype, 
leading to overdetection.
- Control samples provided do not need to be parents, therefore potential 
compound heterozygote detection is simplified to be defined as two or more 
heterozygous candidate variants within the same gene for the case, while 
controls are also allowed to have heterozygous genotypes for this variant. 
This leads to overdetection.
- No SVs or CNVs are considered even if present in the VCF file, leading to 
perhaps some heterozygous variants compounded with a potential SV/CNV to be 
missed. Of course, if this variant is de novo, it will be detected anyhow.
- Non-autosomal detection is simplified to avoid problems with hemizygous 
calls and no sex information is considered. Any interesting variant located 
on an allosome is reported, leading to overdetection.
 - Variants with missing annotations (e.g. CAPICE, GnomAD) cannot be excluded
  based on score or frequency and therefore pass these checks. This leads to 
  overdetection. Of course, genotype matching is done as usual.

## Compiling
Compile using Java 8+ with these dependencies:
```
com.github.samtools:htsjdk:2.22.0
org.molgenis:vcf-io:1.1.1
```

## Demonstration

Download the JAR, a demo file, and run.
```
wget https://github.com/joerivandervelde/capice-quick-filter/releases/download/v0.0.1/capice-quick-filter-0.0.1.jar
wget https://github.com/joerivandervelde/capice-quick-filter/releases/download/demo-data-v0.0.1/CapiceQuickFilter_Demo_1000G.vep.vcfanno.vcf.gz
java -jar capice-quick-filter-0.0.1.jar CapiceQuickFilter_Demo_1000G.vep.vcfanno.vcf.gz CapiceQuickFilter_Demo_HG00096.vcf 0.2 0.05 HG00096 HG00171,HG00403
```

## How to run on your own data

This section explains how to prepare any VCF in order to be analysed by CapiceQuickFilter.

#### Step 1: Annotate your VCF with Ensembl VEP

The input VCF file must be annotated by [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html). The GnomAD allele 
frequencies and gene symbols are used by CapiceQuickFilter. Shown here is how
 to run VEP version 94 with these recommended settings:
```
vep \
--offline \
--cache \
--dir_cache /my/installation/Ensembl/VEP/94 \
--fasta /my/installation/Ensembl/VEP/94/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
--i MyGenomes.vcf.gz --format vcf \
-o MyGenomes.vep.vcf.gz --vcf --compress_output bgzip --force_overwrite \
--species homo_sapiens \
--assembly GRCh37 \
--use_given_ref \
--merged \
--hgvs \
--af_gnomad \
--verbose
```

In order to ensure that the CSQ fields per allele/transcript contain exactly:
```
Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|REFSEQ_MATCH|SOURCE|HGVS_OFFSET|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|CLIN_SIG|SOMATIC|PHENO]
```
It is essential that that ``gnomAD_AF`` is at index 26 and 
``SYMBOL`` is at index 3 (0 based). The same annotations may also be produced by the VEP [web service](http://grch37.ensembl.org/Homo_sapiens/Tools/VEP), but this is untested.


#### Step 2: Annotate your VCF with CAPICE scores

The input VCF file must be annotated by CAPICE. This can be locally installed
 via [Easybuild](https://github.com/molgenis/easybuild-easyconfigs/blob/master/easybuild/easyconfigs/c/CAPICE/CAPICE-v1.0-foss-2018b.eb) or [source code](https://github.com/molgenis/capice).
 The easiest way (though perhaps most labor-intensive in the long run) is to 
 run the openly available CAPICE web service, found at:
```
https://molgenis.org/capice
```

The input it expects has the 5-column VCF layout:
```
1	216497582	.	C	A
5	157217708	.	G	GT
```

Beware that multi allelics get discarded, so check 'discardedInput' and run 
these again.

CAPICE web service output looks like this:


```
1	216497582	C	A	NON_SYNONYMOUS	0.9883755445480347
5	157217708	G	GT	INTRONIC	5.081834751763381e-05
```

This output is  easily converted to VCF using a regex replace:
```
(.+?)\t(.+?)\t(.+?)\t(.+?)\t(.+?)\t(.+)
```
to
```
\1\t\2\t.\t\3\t\4\t.\t.\tCAPICE=\6
```

Of course, prepend a suitable header, so the final VCF file looks like this:

```
##fileformat=VCFv4.3
##fileDate=20200120
##CapiceVersion="1.0"
##contig=<ID=1,length=249250621,assembly=b37>
##contig=<ID=2,length=243199373,assembly=b37>
##contig=<ID=3,length=198022430,assembly=b37>
##contig=<ID=4,length=191154276,assembly=b37>
##contig=<ID=5,length=180915260,assembly=b37>
##contig=<ID=6,length=171115067,assembly=b37>
##contig=<ID=7,length=159138663,assembly=b37>
##contig=<ID=8,length=146364022,assembly=b37>
##contig=<ID=9,length=141213431,assembly=b37>
##contig=<ID=10,length=135534747,assembly=b37>
##contig=<ID=11,length=135006516,assembly=b37>
##contig=<ID=12,length=133851895,assembly=b37>
##contig=<ID=13,length=115169878,assembly=b37>
##contig=<ID=14,length=107349540,assembly=b37>
##contig=<ID=15,length=102531392,assembly=b37>
##contig=<ID=16,length=90354753,assembly=b37>
##contig=<ID=17,length=81195210,assembly=b37>
##contig=<ID=18,length=78077248,assembly=b37>
##contig=<ID=19,length=59128983,assembly=b37>
##contig=<ID=20,length=63025520,assembly=b37>
##contig=<ID=21,length=48129895,assembly=b37>
##contig=<ID=22,length=51304566,assembly=b37>
##contig=<ID=X,length=155270560,assembly=b37>
##contig=<ID=Y,length=59373566,assembly=b37>
##contig=<ID=MT,length=16569,assembly=b37>
##INFO=<ID=CAPICE,Number=1,Type=Float,Description="CAPICE score">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	216497582	.	C	A	.	.	CAPICE=0.9883755445480347
5	157217708	.	G	GT	.	.	CAPICE=5.081834751763381e-05
```

And be sure to sort the chromosomes and positions correctly using:
```
cat MyUnsortedCAPICE_Scores.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1V -k2,2n"}' > MyCAPICE_Scores.vcf
```

And then compress this file:
```
bgzip MyCAPICE_Scores.vcf
```

Now it is ready to be used by [VCFAnno](https://github.com/brentp/vcfanno). A suitable configuration may look 
like this: 
```
[[annotation]]
file="MyCAPICE_Scores.vcf.gz"
fields = ["CAPICE"]
ops=["self"]
names=["CAPICE"]
```

Finally, run VCFAnno using:
```
vcfanno ../vcfanno/CAPICE_conf.toml MyGenomes.vep.vcf.gz 2> vcfanno.log | bgzip > MyGenomes.vep.vcfanno.vcf.gz
```

#### Step 3: Run CapiceQuickFilter

When running CapiceQuickFilter, please supply 5 or 6 arguments:
1. File location of your input .VCF.GZ file.
2. Output file location. May not exist yet.
3. CAPICE score threshold. Lower scoring variants are dropped. Suggesting 0.2 
for 90% sensitivity.
4. GnomAD allele frequency threshold. Higher frequency variants are dropped. 
Suggesting 0.05 to be safe.
5. Case sample ID (ie. the proband, or index).
6. [optional] Control sample ID(s), comma-separated if multiple.

So, in case of ``MyGenomes.vep.vcfanno.vcf.gz``, containing for example, 
perhaps a sample quartet of unaffected mother, unaffected father, affected 
child, and an unaffected sibling:

```
java -jar capice-quick-filter-0.0.1.jar MyGenomes.vep.vcfanno.vcf.gz MyGenomes_AffChild01.vcf 0.2 0.05 AffChild01 Father01,Mother01,Sib01
```

## To do
- Unit and integration testing
- Proper dependency management
- Proper cmdline option parsing
- Address issues mentioned in manual
- Address issues mentioned in code
