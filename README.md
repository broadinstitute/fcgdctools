# fcgdctools
A python project containing utilities to support interoperability between FireCloud and the [NIH/NCI Genomics Data Commons](https://gdc.cancer.gov/).  The project currently contains the tool `genFcWsLoadFiles` for creating FireCloud workspace load files (TSV-formatted files) from a file manifest downloaded from the GDC.  The manifest may have been downloaded from either the GDC's main or legacy archive portal.  Additional tools will be added over time.  

## Requirements

fcgdctools requires: 
* python 3 


## Installation

Install this package via PyPI:

```
	% pip install fcgdctools
```

Or download the source and install:

```
	% git clone https://github.com/broadinstitute/fcgdctools.git
	% cd fcgdctools
	% pip install .
```
Note that if you are installing to a protected location, you may need to preface the `pip install .` command with `sudo`.  

## Description

Following installation you should be able to run the `gdcLoadFiles` and `gdcWorkspace` commands from the command line

```
    % gdcLoadFiles -h
    usage: gdcLoadFiles [-h] [-c] [-a {main,legacy,awg}] [-t TOKEN] manifest
    
    create FireCloud workspace load files from GDC manifest
    
    positional arguments:
      manifest              manifest file from the GDC Data Portal
    
    optional arguments:
      -h, --help            show this help message and exit
      -c, --all_cases       create participant entities for all referenced cases
      -a {main,legacy,awg}, --api {main,legacy,awg}
                            Select API endpoint. awg requires a valid GDC download
                            token.
      -t TOKEN, --token TOKEN
                            GDC download token to access controlled APIs
```
```
    % gdcWorkspace -h
    usage: gdcWorkspace [-h] [-d AUTH_DOMAIN] [-a] [-t TOKEN]
                        project_name cohort_name billing_project ws_suffix
    
    create FireCloud workspace for a specific project + cohort + access type
    
    positional arguments:
      project_name          the name of the project. e.g: TCGA
      cohort_name           name_of_cancer_cohort. e.g: LUAD
      billing_project       name of billing project to create the workspace under.
                            e.g: broad-firecloud-tcga
      ws_suffix             descriptive suffix to add to the workspace auto-
                            generated name. e.g: ControlledAccess_hg38_V1-0_DATA
    
    optional arguments:
      -h, --help            show this help message and exit
      -d AUTH_DOMAIN, --auth_domain AUTH_DOMAIN
                            authorization domain. for dbGaP controlled access the
                            domain name is TCGA-dbGaP-Authorized.
      -a, --awg             use AWG api endpoint. Requires a valid GDC download
                            token.
      -t TOKEN, --token TOKEN
                            GDC download token to access controlled APIs. Required
                            for AWG API.
```

By default, the tool assumes the manifest references harmonized data from the GDC's principal portal.  For each file listed in the manifest, the tool queries the GDC for file metadata (e.g., the cases and samples it is associated with, the file's data category, data type, etc.). After assembling the files' metadata, the tool creates FireCloud Workspace Load Files for populating a FireCloud workspace with participant, sample and pair entities containing attributes whose contents reference the listed files.  For each entity type, an attribute is defined for each type of file associated with that entity type.  Attribute names are derived as follows:

```
    [FFPE__][<experimental strategy abbrev>__][<workflow type abbrev>__]<data type abbrev>__<data format abbrev>__{drs_url,gdc_url}
```
Here are some examples of attribute names:

```
    biospecimen_supplement__bcr_ssf_xml__drs_url
    WXS__BWAMDupCoClean__aligned_reads__bam__drs_url
    RNAseq__STAR2Pass__aligned_reads__bam__drs_url
    WXS__MuTect2__annotated_simple_somatic_mutation__vcf__drs_url
```

Attribute values are file references, consisting of a DRS or GDC API URI:

```
   {drs://dataguids.org,https://<GDC_API_root>/data}/<file uuid>
```

Here are some examples of attribute values corresponding to the above attribute names:

```
    drs://dataguids.org/6266607c-5caa-4bea-be2e-74271846c171
    drs://dataguids.org/c61bb1ab-688f-4d58-8388-60ae77c28840
    drs://dataguids.org/5282c243-8d44-485a-beeb-ffb62201a60b
    drs://dataguids.org/2bb1f8bb-8834-4095-b1e0-028212d26731
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

Finally, the tool creates a .tsv file with general workflow attributes.
Right now, the two attributes that are created are:

`legacy_flag` - a flag that indicates if the manifest was downloaded from the legacy archive. The flag is needed for other scripts to know where to get more file information from, e.g. the size of a file. The flag is boolean and will always be set to "false".

`workspace-column-defaults` - the default order in which the attribute columns should be shown in the table.  

Please note that there are instances where multiple files map to the same attribute name.  In these situations, fcgdctools attempts to select the "best" file based on metadata stored in the aliquot submitter id (for TCGA, the aliquot barcode).  In cases where the aliquot submitter ids are identical fcgdctools makes an arbitrary selection and prints a warning to stdout.  Users should search stdout for these warnings and adjust their loadfiles if fcgdctools' choice is incorrect.
