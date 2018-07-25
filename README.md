# fcgdctools
A python project containing utilities to support intereoperability between FireCloud and the [NIH/NCI Genomics Data Commons](https://gdc.cancer.gov/).  The project currently contains the tool `genFcWsLoadFiles` for creating FireCloud workspace load files (TSV-formatted files) from a file manifest downloaded from the GDC.  The manifest may have been downloaded from either the GDC's main or legacy archive portal.  Additional tools will be added over time.  

## Requirements

fcgdctools requires: 

* python 3 
* python 3 `requests` package.  Install command: `pip install requests` 


## Installation

Install this package via pip:

```
	% pip install fcgdctools
```

Or download the source and install with setup.py:

```
	% git clone https://github.com/broadinstitute/fcgdctools.git
	% cd fcgdctools
	% python setup.py build
	% python setup.py install
```
Note that if you are installing to a protected location, you may need to preface the `python setup.py install` command with `sudo`.  

## Description
Following installation you should be able to run the `genFcWsLoadFiles` command from the command line

```
	% genFcWsLoadFiles -h
	usage: genFcWsLoadFiles [-h] [-r RESOLVE_UUIDS] [-c] manifest

	create FireCloud workspace load files from GDC manifest

	positional arguments:
	  manifest              manifest file from the GDC Data Portal

	optional arguments:
	  -h, --help            show this help message and exit
	  -r RESOLVE_UUIDS, --resolve_uuids RESOLVE_UUIDS
                        TSV file mapping GDC UUIDs to URLs
	  -c, --all_cases       create participant entities for all referenced cases
  ```
By default, the tool assumes the manifest references harmonized data from the GDC's principal portal.  For each file listed in the manifest, the tool queries the GDC for file metadata (e.g., the cases and samples it is associated with, the file's data category, data type, etc.). After assembling the files' metadata, the tool creates FireCloud Workspace Load Files for populating a FireCloud workspace with participant, sample and pair entities containing attributes whose contents reference the listed files.  For each entity type, an attribute is defined for each type of file associated with that entity type.  Attribute names are derived as follows:

```
    [<experimental strategy abbrev>__][<workflow type abbrev>__]<data type abbrev>__<data format abbrev>__uuid_and_filename
```
Here are some examples of attribute names:

```
    biospecimen_supplement__bcr_ssf_xml__uuid_and_filename
    WXS__BWAMDupCoClean__aligned_reads__bam__uuid_filename
    RNAseq__STAR2Pass__aligned_reads__bam__uuid_and_filename
    WXS__MuTect2__annotated_simple_somatic_mutation__vcf__uuid_and_filename
```

Attribute values are file references, consisting of a concatenation a uuid and filename:

```
   <file uuid>/<filename>
```

Here are some examples of attribute values corresponding to the above attribute names:

```
    6266607c-5caa-4bea-be2e-74271846c171/nationwidechildrens.org_ssf.TCGA-C8-A137.xml
    c61bb1ab-688f-4d58-8388-60ae77c28840/TCGA-BH-A0HA-11A-31D-A12Q-09_IlluminaGA-DNASeq_exome_gdc_realn.bam
    5282c243-8d44-485a-beeb-ffb62201a60b/98cfb9c2-7c1e-4bc3-bea7-33b1ecb3ec0d_gdc_realn_rehead.bam
    2bb1f8bb-8834-4095-b1e0-028212d26731/2bb1f8bb-8834-4095-b1e0-028212d26731.vep.vcf.gz
```

Slide images (Tissue Slides and Diagnostic Slides) are handled a bit differently than the genomic data files.  Frequently a single biospecimen sample has several slide images associated with it; for example, top and bottom tissue slides or multiple diagnostic slides. Slide images cannot be distiguished from one another via the GDC's file metadata and researchers may want to include a sample's multiple slide images in the workspace.  Encoded in the file names is image metadata (e.g., the TCGA Slide barcode); this can be used to distinguish between slide images, and we incorporate the slide barcode's slide ID into the attribute name.

The tool also creates load files for defining sets of participants, samples and pairs.  An entity set is defined for each file attribute attached to that entity type; the set consists of those entities that have a non-empty value for that attribute.  The set may be used to target workflows that operate on that file type.  In particular, the set may be use to run a workflow that retrieves from the GDC the files referenced by the corresponding attribute.  

The sets have the following identifier naming convention:

```
	<attribute basename>
	where
	<attribute basename> = [<experimental strategy abbrev>__][<workflow type abbrev>__]<data type abbrev>__<data format abbrev>
``` 

Here are some example set identifiers:

```
	biospecimen_data__bcr_ssf__xml
	WXS__BWAMDupCoClean__aligned_reads__bam
	RNAseq__STAR2Pass__aligned_reads__bam
	WXS__MuTect2__raw_simple_somatic_mutation__vcf
	
```
This tool DOES NOT support manifests downloaded from the GDC Legacy Archive.

The optional input `RESOLVE_UUIDS` is a TSV file containing mappings of file uuids to urls of the locations of the files on cloud storage.  If this optional input is provided, `genFcWsLoadFiles` will add to the load files it generates attributes with suffix `__url`, which contain the url mapped to the uuid.

Finally, the tool creates a .tsv file with general workflow attributes.
Right now, the two attributes that are created are:

`legacy_flag` - a flag that indicates if the manifest was downloaded from the legacy archive. The flag is needed for other scripts to know where to get more file information from, e.g. the size of a file. The flag is boolean and will always be set to "false".

`workspace-column-defaults` - the default order in which the attribute columns should be shown in the table.  

Please note that there are instances where multiple files map to the same attribute name.  In these situations, fcgdctools attempts to select the "best" file based on metadata stored in the aliquot submitter id (for TCGA, the aliquot barcode).  In cases where the aliquot submitter ids are identical fcgdctools makes an arbitrary selection and prints a warning to stdout.  Users should search stdout for these warnings and adjust their loadfiles if fcgdctools' choice is incorrect.