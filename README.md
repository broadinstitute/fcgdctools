# fcgdctools
A python project containing utilities to support intereoperability between FireCloud and the [NIH/NCI Genomics Data Commons](https://gdc.cancer.gov/).  The project currently contains the tool `genFcWsLoadFiles` for creating FireCloud workspace load files (TSV-formatted files) from a file manifest downloaded from the GDC.  The manifest may have been downloaded from either the GDC's main or legacy archive portal.  Additional tools will be added over time.  

## Requirements

fcgdctools requires: 

* python 3 
* python 3 `requests` package.  Install command: `pip3 install requests` 


## Installation

To get started, from a command line terminal clone the repo, build, and install:

```
	% git clone https://github.com/broadinstitute/fcgdctools.git
	% cd fcgdctools
	% python3 setup.py build
	% python3 setup.py install
```
Note that if you are installing to a protected location, you may need to preface the `python3 setup.py install` command with `sudo`.  

We plan on making fcgdctools installable via pip3 (from [PyPI](https://pypi.python.org/pypi)), and will update these instructions when the fcgdctools package is submitted to the PyPI repository.

## Description
Following installation you should be able to run the `genFcWsLoadFiles` command from the command line

```
	% cd
	% genFcWsLoadFiles -h
	usage: genFcWsLoadFiles [-h] [-r RESOLVE_UUIDS] [-l] manifest

	create FireCloud workspace load files from GDC manifest

	positional arguments:
  	manifest                manifest file from the GDC Data Portal

	optional arguments:
  	  -h, --help            show this help message and exit
  	  -r RESOLVE_UUIDS, --resolve_uuids RESOLVE_UUIDS
                            resolve uuids in the manifest to urls
  	  -l, --legacy          point to GDC Legacy Archive
  ```
By default, the tool assumes the manifest references harmonized data from the GDC's principal portal.  For each file listed in the manifest, the tool queries the GDC for file metadata (e.g., the cases and samples it is associated with, the file's access status (open or controlled), data category and data type, etc.). After assembling the files' metadata, the tool creates FireCloud Workspace Load Files for populating a FireCloud workspace with participant, sample and pair entities containing attributes whose contents reference the listed files.  For each entity type, an attribute is defined for each type of file associated with that entity type.  Attribute names are derived as follows:

```
    [<experimental strategy abbrev>__][<workflow type abbrev>__]<data type abbrev>__uuid_and_filename
```
Here are some examples of attribute names:

```
    biospecimen_data__uuid_and_name
    miRNAseq__BCGSC__isoform_expression_quantification__uuid_and_filename
    RNA_seq__STAR2Pass__aligned_reads__uuid_and_filename
    WXS__BWAMDupCoClean__aligned_reads__uuid_filename
    WXS__MuTect2__raw_simple_somatic_mutation__uuid_and_filename
    
```

Attribute values are file references, consisting of a concatenation a uuid and filename:

```
   <file uuid>/<filename>
```

Here are some examples of attribute values corresponding to the above attribute names:

```
    40a28c9d-317d-4589-8071-d7ae1ac6430d/nationwidechildrens.org_biospecimen.TCGA-ZH-A8Y3.xml
    3d1f143a-c4e4-4c38-a99b-357f370ebaa5/isoforms.quantification.txt
    19fecfdf-d758-44d7-acaa-d25afc4c45c6/1e8db45c-5078-4985-ac7e-09082f2b2297_gdc_realn_rehead.bam
    9ee608a6-975b-4c69-b558-f37b56edd657/7075a31b8c33ef962d8a336d0ad83abd_gdc_realn.bam
    d014e1db-b77d-42ec-9db7-96c7a4bfd23e/d014e1db-b77d-42ec-9db7-96c7a4bfd23e.vcf.gz
    
```
The tool also creates load files for defining sets of participants, samples and pairs.  An entity set is defined for each file attribute attached to that entity type; the set consists of those entities that have a non-empty value for that attribute.  The set may be used to target workflows that operate on that file type.  In particular, the set may be use to run a workflow that retrieves from the GDC the files referenced by the corresponding attribute.  

The sets' have the following identifier naming convention:

```
	(OA|CA)__<attribute basename>
	where
	<attribute basename> = [<experimental strategy abbrev>__][<workflow type abbrev>__]<data type abbrev>
``` 

The OA/CA prefix specify whether the files reference by the corresponding attribute are open access (OA) or controlled access (CA).

Here are some example set identifiers:

```
	OA__biospecimen_data
	OA__miRNAseq__BCGSC__isoform_expression_quantification
	CA__RNAseq__STAR2Pass__aligned_reads
	CA__WXS__BWAMDupCoClean__aligned_reads
	CA__WXS__MuTect2__raw_simple_somatic_mutation
	
```

The optional input `RESOLVE_UUIDS` is a TSV file containing mappings of file uuids to urls of the locations of the files on cloud storage.  If this optional input is provided, `genFcWsLoadFiles` will add to the load files it generates attributes with suffix `__url`, which contain the url mapped to the uuid.

If your manifest was downloaded from the GDC Legacy Archive, you must use the `-l` option.