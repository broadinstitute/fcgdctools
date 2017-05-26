# Version: 6
# Date: March 26, 2017
# Comments: combine uuid and filename into single attributeincorporated uuid-to-url resolution
# ToDo: should probably use gdctools for querying of GDC.
#       have multiple side images per case - we are not accounting for that

import csv
import requests
import argparse
import pprint
import os.path
import sys
import time

from fcgdctools import gdc_uuidresolver 

GDC_API_ROOT = "https://gdc-api.nci.nih.gov"
GDC_LEGACY_API_ROOT = "https://gdc-api.nci.nih.gov/legacy"

#access category
class GDC_FileAccessType:
    OPEN = 'open'
    CONTROLLED = 'controlled'

#data categories
class GDC_DataCategory:
    CLINICAL = 'Clinical'
    BIOSPECIMEN = "Biospecimen"
    DNA_METHYLATION = "DNA Methylation"
    RAW_SEQUENCING_DATA = "Raw Sequencing Data"
    COPY_NUMBER_VARIATION = "Copy Number Variation"
    TRANSCRIPTOME_PROFILING = "Transcriptome Profiling"
    SNV = "Simple Nucleotide Variation"

    #legacy data categories
    RAW_MICROARRAY_DATA = "Raw microarray data"
    
#data types
class GDC_DataType:
    #data types associated with Clinical data category
    CLINICAL_SUPPLEMENT = "Clinical Supplement"
    #legacy data types associated with Clinical data category
    TISSUE_SLIDE_IMAGE = "Tissue slide image"
    DIAGNOSTIC_IMAGE = "Diagnostic image"
    PATHOLOGY_REPORT = "Pathology report"
    CLINICAL_DATA = "Clinical data"
    BIOSPECIMEN_DATA = "Biospecimen data"

    #associated with Biospecimen data category
    BIOSPECIMEN_SUPPLEMENT = "Biospecimen Supplement"

    #assocaited with DNA Methylation data category
    METHYLATION_BETA_VALUE = "Methylation Beta Value"

    #associated with RAW Sequencing Data data category
    ALIGNED_READS = "Aligned Reads"

    #associated with Copy Number Variation data category
    COPY_NUMBER_SEGMENT = 'Copy Number Segment'
    MASKED_COPY_NUMBER_SEGMENT = 'Masked Copy Number Segment'

    #associated with Transcriptome Profiling data category
    MIRNA_EXPRESSION_QUANTIFICATION = "miRNA Expression Quantification"
    ISOFORM_EXPRESSION_QUANTIFICATION = "Isoform Expression Quantification"
    GENE_EXPRESSION_QUANTIFICATION = "Gene Expression Quantification"

    #associated with Simple Nucleotide Variation data category
    RAW_SIMPLE_SOMATIC_MUTATION = "Raw Simple Somatic Mutation"
    ANNOTATED_SOMATIC_MUTATION = "Annotated Somatic Mutation"
    AGGREGATED_SOMATIC_MUTATION = "Aggregated Somatic Mutation"
    MASKED_SOMATIC_MUTATION = "Masked Somatic Mutation"


class DataSource:
    ABBREV_TRANSLATE_TABLE = ''.maketrans({'.' : '', '-' : '', ' ' : '', '_' :''})    
    def __init__(self, abbreviations):
        self.abbreviations = abbreviations

    def getAbbreviation(self, type):
        if type in self.abbreviations:
            return self.abbreviations[type]
        else:
            return type.translate(DataSource.ABBREV_TRANSLATE_TABLE)
        
EXP_STRATEGY_ABBREVIATIONS = {
        'WXS' : 'WXS',
        'RNA-Seq' : 'RNAseq',
        'miRNA-Seq' : 'miRNAseq',
        'Genotyping Array' : 'GeneArray',
        'Methylation Array' : 'MethArray'}

EXP_STRATEGY = DataSource(EXP_STRATEGY_ABBREVIATIONS)
    
WORKFLOW_ABBREVIATIONS = {
        'SomaticSniper': 'SomSnip',
        'MuTect2' : 'MuTect2',
        'VarScan2' : 'VarScan2',
        'MuSE' : 'MuSE',
        'SomaticSniper Annotation' : 'SomSnipAnnot',
        'MuTect2 Annotation' : 'MuTect2Annot',
        'VarScan2 Annotation' : 'VarScan2Annot',
        'MuSE Annotation' : 'MuSEAnnot',
        'MuSE Variant Aggregation and Masking' : 'MuSEAggrMask',
        'MuTect2 Variant Aggregation and Masking' : 'MuTect2AggrMask',
        'SomaticSniper Variant Aggregation and Masking' : 'SomSnipAggrMask',
        'VarScan2 Variant Aggregation and Masking' : 'VarScan2AggrMask',
        'BCGSC miRNA Profiling' : 'BCGSCmiRNA',
        'HTSeq - Counts' : 'HTSeqCounts',
        'HTSeq - FPKM' : 'HTSeqFPKM',
        'HTSeq - FPKM-UQ' : 'HTSEQFPKMUQ',
        'DNACopy' : 'DNACopy',
        'BWA with Mark Duplicates and Cocleaning' : 'BWAMDupCoClean',
        'STAR 2-Pass' : 'STAR2Pass',
        'BWA-aln' : 'BWAaln',
        'Liftover' : 'Lift'}

WORKFLOW = DataSource(WORKFLOW_ABBREVIATIONS)
        
class Platform(DataSource):

    def __init__(self, abbreviations):
        DataSource.__init__(self, abbreviations)
        

PLATFORM_ABBREVIATIONS = {
        'Affymetrix SNP 6.0' : 'AffySNP6',
        'Illumina' : 'Illum',
        'Illumina Human Methylation 450' : 'IllumHuMeth450',
        'Illumina Human Methylation 27' : 'IllumHuMeth27'}

PLATFORM = Platform(PLATFORM_ABBREVIATIONS)    
    
class GDC_FileAccess:

    def __init__(self):

        self.accessTypePrefix = {
            GDC_FileAccessType.OPEN : 'OA__',
            GDC_FileAccessType.CONTROLLED : 'CA__'
            }

        self.gdcFileAccess = dict()
        

    def recordFileAccessType(self, file_attribute_base_name, access_type):
         
        assert access_type in [GDC_FileAccessType.OPEN, GDC_FileAccessType.CONTROLLED], "unexpected file access type"
        
        if file_attribute_base_name in self.gdcFileAccess:
            assert self.gdcFileAccess[file_attribute_base_name] == access_type, "inconsistent file access type"
        else:
            self.gdcFileAccess[file_attribute_base_name] = access_type
        
    def getAccessTypePrefix(self, file_attribute_base_name):
        assert file_attribute_base_name in self.gdcFileAccess, "file_attribute_base_name not in gdcFileAccess dict"
        return self.accessTypePrefix[self.gdcFileAccess[file_attribute_base_name]]

GDC_FILE_ACCESS = GDC_FileAccess()

# Sample Types
# from https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
class SampleType:
    SAMPLE_TYPES_DESCRIPTION = 0
    SAMPLE_TYPES_LETTER_CODE = 1
    SAMPLE_TYPES_TN = 2
    TUMOR = 'tumor'
    NORMAL = 'normal'
    NA = 'na'
    # this table needs updating - TARGET data had some unknown sample type ids
    SAMPLE_TYPES = {'01': ['Primary Solid Tumor', 'TP', TUMOR],
                '02' : ['Recurrent Solid Tumor', 'TR', TUMOR],
                '03' : ['Primary Blood Derived Cancer - Peripheral Blood', 'TB', TUMOR],
                '04' : ['Recurrent Blood Derived Cancer - Bone Marrow', 'TRBM', TUMOR],
                '05' : ['Additional - New Primary', 'TAP', TUMOR],
                '06' : ['Metastic', 'TM', TUMOR],
                '07' : ['Additional Metastic', 'TAM', TUMOR],
                '08' : ['Human Tumor Original Cells', 'THOC', TUMOR],
                '09' : ['Primary Blood Derived Cancer - Bone Marrow', 'TBM', TUMOR],
                '10' : ['Blood Derived Normal', 'NB', NORMAL],
                '11' : ['Solid Tissue Normal', 'NT', NORMAL],
                '12' : ['Buccal Cell Normal', 'NBC', NORMAL],
                '13' : ['EBV Immortalized Normal', 'NEBV', NORMAL],
                '14' : ['Bone Marrow Normal', 'NBM', NORMAL],
                '15' : ['sample type 15', '15SH', NA],
                '16' : ['sample type 16', '16SH', NA],
                '20' : ['Control Analyte', 'CELLC', NA],
                '40' : ['Recurrent Blood Derived Normal - Peripheral Blood', 'TRB', NORMAL],
                '41' : ['Blood Derived Cancer - Bone Marrow, Post-treatment', 'TBD', TUMOR],
                '42' : ['Blood Derived Cancer - Peripheral Blood, Post-treatement', 'TBD', TUMOR],
                '50' : ['Cell Lines', 'CELL', NA],
                '60' : ['Primary Xenograft Tissue', 'XP', NA],
                '61' : ['Cell Line Derived Xenograft Tissue', 'XCL', NA],
                '99' : ['sample type 99', '99SH', NA]}

    def getTumorNormalClassification(self, sample_type_id):
        return self.SAMPLE_TYPES[sample_type_id][2]

    def getLetterCode(self, sample_type_id):
        return self.SAMPLE_TYPES[sample_type_id][1]

SAMPLE_TYPE = SampleType()

    
class Platform(DataSource):

    def __init__(self, abbreviations):
        DataSource.__init__(self, abbreviations)
        

PLATFORM_ABBREVIATIONS = {
        'Affymetrix SNP 6.0' : 'AffySNP6',
        'Illumina' : 'Illum',
        'Illumina Human Methylation 450' : 'IllumHuMeth450',
        'Illumina Human Methylation 27' : 'IllumHuMeth27'}

PLATFORM = Platform(PLATFORM_ABBREVIATIONS)

class MetadataRetriever():
    def __init__(self, gdc_api_root, fields):
        self.gdc_api_root = gdc_api_root
        self.fields = fields

    def get_metadata(self, file_uuid):
        url = "{0}/files/{1}?fields={2}".format(self.gdc_api_root, file_uuid, self.fields)
        response = requests.get(url, headers=None)
        responseDict = response.json()
        return responseDict['data']

class CaseMetadataRetriever(MetadataRetriever):
    def __init__(self, gdc_api_root):
        fields = "cases.case_id,cases.submitter_id,cases.project.project_id"
        MetadataRetriever.__init__(self, gdc_api_root, fields)
    
class CaseSampleMetadataRetriever(MetadataRetriever):
    def __init__(self, gdc_api_root):
        fields = "cases.case_id,cases.submitter_id,cases.project.project_id"
        fields = fields + ",cases.samples.sample_id,cases.samples.submitter_id,cases.samples.sample_type_id"
        MetadataRetriever.__init__(self, gdc_api_root, fields)

class FileMetadataRetriever(MetadataRetriever):
    def __init__(self, gdc_api_root):
        if gdc_api_root == GDC_API_ROOT:
            fields = "data_category,data_type,access, data_format,experimental_strategy,analysis.workflow_type,submitter_id"
        else:
            fields = "data_category,data_type,access,experimental_strategy"
        MetadataRetriever.__init__(self, gdc_api_root, fields)

SEPARATOR = '/'
UUID_ATTRIBUTE_SUFFIX = "__uuid_and_filename"
URL_ATTRIBUTE_SUFFIX = "__url"
    
def _read_manifestFile(manifestFile):
    
    manifestFileList = []

    with open(manifestFile, 'r') as fp:
        reader = csv.DictReader(fp, delimiter='\t')
        print(reader.fieldnames)
        for row in reader:
            manifestFileList.append(row)
    
    return manifestFileList

def _add_to_knowncases(case_metadata, known_cases):
    case_id = case_metadata['case_id']
    if case_id not in known_cases:
        submitter_id = case_metadata['submitter_id']
        project_id = case_metadata['project']['project_id']
        new_case = {'submitter_id' : submitter_id, 'project_id' : project_id}
        known_cases[case_id] = new_case
    return case_id

def _add_to_knownsamples(sample_metadata, case_id, known_samples):
    sample_id = sample_metadata['sample_id']
    sample_type_id = sample_metadata['sample_type_id']
    if sample_id not in known_samples:
        sample_submitter_id = sample_metadata['submitter_id']
        new_sample = {'submitter_id' : sample_submitter_id, 
                      'sample_type_id' : sample_type_id,
                      'case_id' : case_id}
        known_samples[sample_id] = new_sample
    return sample_id, SAMPLE_TYPE.getTumorNormalClassification(sample_type_id)

def _add_to_knownpairs(tumor_sample_id, normal_sample_id, known_pairs):
    pair_id = "{0}_{1}".format(tumor_sample_id, normal_sample_id)
    if pair_id not in known_pairs:
        known_pairs[pair_id] = {'tumor': tumor_sample_id, 'normal': normal_sample_id}
    return pair_id

def _constructAttributeName_base(experimental_strategy, workflow_type, data_category, data_type):

    if experimental_strategy is not None:
        experimental_strategy_abbrev = EXP_STRATEGY.getAbbreviation(experimental_strategy) + '__'
    else:
        experimental_strategy_abbrev = ''

    if workflow_type is not None:
        workflow_type_abbrev = WORKFLOW.getAbbreviation(workflow_type) + '__'
    else:
        workflow_type_abbrev = ''

    data_type_lc = data_type.lower().replace(" ", "_")

    attribute_name_base = experimental_strategy_abbrev + workflow_type_abbrev + data_type_lc

    return (attribute_name_base)

def _getImageCodeAndPortionFromImageFilename(filename):
    image_code = filename.split('.')[0].split('-')[-1]
    portion = int(filename.split('.')[0].split('-')[-2])
    return image_code, portion

def _constructImageAttributeName_base(experimental_strategy, workflow_type, data_category, data_type, filename=None):

    if experimental_strategy is not None:
        experimental_strategy_abbrev = EXP_STRATEGY.getAbbreviation(experimental_strategy) + '__'
    else:
        experimental_strategy_abbrev = ''

    if workflow_type is not None:
        workflow_type_abbrev = WORKFLOW.getAbbreviation(workflow_type) + '__'
    else:
        workflow_type_abbrev = ''

    data_type_lc = data_type.lower().replace(" ", "_") + '_'

    # see https://wiki.nci.nih.gov/display/TCGA/TCGA+barcode# for interpretation of TCGA bar code
    # that is incorporated into image filename
    image_code, portion = _getImageCodeAndPortionFromImageFilename(filename)
    if 'BS' in image_code:
        image_abbrev = 'bottom'
    elif 'TS' in image_code:
        image_abbrev = 'top'
    elif 'MS' in image_code:
        image_abbrev = 'middle'
    else:
        image_abbrev=image_code[0:2]

    
    attribute_name_base = experimental_strategy_abbrev + workflow_type_abbrev + data_type_lc + image_abbrev

    return attribute_name_base, portion


def _pick_tcga_submitter(a, b):
    '''Comparator function for barcodes, using the rules described in the GDAC
    FAQ entry for replicate samples: https://confluence.broadinstitute.org/display/GDAC/FAQ
    '''
    # Get the analytes and plates
    # TCGA-BL-A0C8-01A-11<Analyte>-<plate>-01
    analyte1 = a[19]
    analyte2 = b[19]
    plate1   = a[21:25]
    plate2   = b[21:25]

    # Equals case
    if a == b:
        return a
    elif analyte1 == analyte2:
        # Prefer the aliquot with the highest lexicographical sort value
        return a if a >= b else b
    elif analyte1 == "H":
        # Prefer H over R and T
        return a
    elif analyte1 == "R":
        # Prefer R over T
        return a if analyte2 == "T" else b
    elif analyte1 == "T":
        # Prefer H and R over T
        return b
    elif analyte1 == "D":
        # Prefer D over G,W,X, unless plate number is higher
        return a if plate2 <= plate1 else b
    elif analyte2 == "D":
        return b if plate1 <= plate2 else a
    else:
        # Default back to highest lexicographical sort value
        return a if a >= b else b

def _pick_target_submitter(a,b):
    # TODO: implement this after implementation method is decided.
    return a


def _resolve_collision(gdc_api_root, data_category, data_type, uuid1, name1, uuid2, name2):

    # NOTE: we chose not to employ the created_datetime or updated_datetime fields in 
    # our decision logic.  From what we can tell, neither should be used to make a selection between 
    # two files.


    # get aliquot submitter IDs
    aliquot_submitter_id1 = None
    aliquot_submitter_id2 = None

    # SNV files associated with two samples: tumor and normal. The choice of the right file is made according to the Tumor sample.
    # note: conditional logic below to accommodate legacy archive
    if data_category.lower() == GDC_DataCategory.SNV.lower() and data_type not in set([GDC_DataType.AGGREGATED_SOMATIC_MUTATION, GDC_DataType.MASKED_SOMATIC_MUTATION]):

        file_fields = "cases.project.program.name,cases.samples.sample_type,cases.samples.portions.analytes.aliquots.submitter_id,cases.samples.sample_type_id"
        meta_retriever = MetadataRetriever(gdc_api_root, file_fields)

        data1 = meta_retriever.get_metadata(uuid1)
        samples_list1 = data1['cases'][0]['samples']
        for s in samples_list1:
            sample_type = SAMPLE_TYPE.getTumorNormalClassification(s['sample_type_id'])
            if sample_type == SAMPLE_TYPE.TUMOR:
                aliquot_submitter_id1 = s['portions'][0]['analytes'][0]['aliquots'][0]['submitter_id']
        assert aliquot_submitter_id1 is not None, "no tumor sample in pair"

        data2 = meta_retriever.get_metadata(uuid2)
        samples_list2 = data2['cases'][0]['samples']
        for s in samples_list2:
            sample_type = SAMPLE_TYPE.getTumorNormalClassification(s['sample_type_id'])
            if sample_type == SAMPLE_TYPE.TUMOR:
                aliquot_submitter_id2 = s['portions'][0]['analytes'][0]['aliquots'][0]['submitter_id']
        assert (aliquot_submitter_id2 is not None), "no tumor sample in pair"
    
    # Here we handle other file types that are associated with single sample, or multiple cases. 
    # Files that are associated with multiple cases won't use information encoded in aliquot barcode
    else:
        file_fields = "cases.project.program.name,cases.samples.portions.analytes.aliquots.submitter_id"
        meta_retriever = MetadataRetriever(gdc_api_root, file_fields)

        data1 = meta_retriever.get_metadata(uuid1)
        data2 = meta_retriever.get_metadata(uuid2)

        num_cases1 = len(data1['cases'])
        num_cases2 = len(data2['cases'])

        # Take care of the case that at least one of the files is associated with more than one case.
        if num_cases1 > 1 or num_cases2 > 1:
            print("number of cases...")
            print("{0}: {1}, {2}: {3}".format(uuid1,num_cases1,uuid2,num_cases2))
            # If one of the two files has more cases associated with it, we assume it's the correct file to pick.
            if num_cases1 >= num_cases2:
                return uuid1, name1
            else:
                return uuid2, name2
        
        # Both files are associated with only one case and one sample
        else:
            assert len(data1['cases'][0]['samples']) == 1, "more than one sample associated with file"
            assert len(data2['cases'][0]['samples']) == 1, "more than one sample associated with file"
            aliquot_submitter_id1 = data1['cases'][0]['samples'][0]['portions'][0]['analytes'][0]['aliquots'][0]['submitter_id']
            aliquot_submitter_id2 = data2['cases'][0]['samples'][0]['portions'][0]['analytes'][0]['aliquots'][0]['submitter_id']

    # At this point we've already taken care of files that were associated with multiple cases and are now left with single-case files, 
    # single-sample files
    print("new aliquot submitter ID is: {0}".format(aliquot_submitter_id1))
    print("old aliquot submitter ID is: {0}".format(aliquot_submitter_id2))

    program1 = data1['cases'][0]['project']['program']['name']

    if program1 == "TCGA":
        chosen = _pick_tcga_submitter(aliquot_submitter_id1, aliquot_submitter_id2)

    elif program1 == "TARGET":
        chosen = _pick_target_submitter(aliquot_submitter_id1, aliquot_submitter_id2)
    else:
        chosen = aliquot_submitter_id1

    if chosen == aliquot_submitter_id1:
        return uuid1, name1
    else:
        return uuid2, name2


def _add_file_attribute(gdc_api_root, entity, file_uuid, filename, file_url,
                        data_category, data_type, experimental_strategy, workflow_type, access):
    # I needed to insert some special-case processing for image data files
    # this probably isn't the cleanest way to handle it, but good enough for now
    if data_type in set([GDC_DataType.TISSUE_SLIDE_IMAGE, GDC_DataType.DIAGNOSTIC_IMAGE]):
        basename, portion = _constructImageAttributeName_base(experimental_strategy, workflow_type,
                                                              data_category, data_type, filename)
        attribute_name = basename + UUID_ATTRIBUTE_SUFFIX

        if attribute_name in entity:
            print("multiple portions!")
            filename_present = entity[attribute_name]
            _, portion_present = _getImageCodeAndPortionFromImageFilename(filename_present)
            if portion > portion_present:
                GDC_FILE_ACCESS.recordFileAccessType(basename, access)
                entity[basename + UUID_ATTRIBUTE_SUFFIX] = file_uuid + SEPARATOR + filename
                entity[basename + URL_ATTRIBUTE_SUFFIX] = file_url
        else:
            GDC_FILE_ACCESS.recordFileAccessType(basename, access)
            entity[basename + UUID_ATTRIBUTE_SUFFIX] = file_uuid + SEPARATOR + filename
            entity[basename + URL_ATTRIBUTE_SUFFIX] = file_url            
    else:
        basename = _constructAttributeName_base(experimental_strategy, workflow_type,
                                                data_category, data_type)
        attribute_name = basename + UUID_ATTRIBUTE_SUFFIX
        
        # see if attribute already defined for entity
        if attribute_name in entity:
            existing_file = entity[attribute_name]
            existing_uuid = existing_file.split("/")[0]
            existing_filename = existing_file.split("/")[1]

            print("multiple files for same attribute!") 
            print("attribute name is {0}".format(attribute_name))
            print("new file: {0}/{1}".format(file_uuid, filename))
            print("existing file: {0}".format(entity[attribute_name]))
            
            chosen_uuid, chosen_filename = _resolve_collision(gdc_api_root, data_category, data_type, 
                                                              file_uuid, filename, existing_uuid, existing_filename)
            print("chosen file is: {0}/{1}".format(chosen_uuid, chosen_filename))

            if chosen_uuid == file_uuid:
                GDC_FILE_ACCESS.recordFileAccessType(basename, access)
                entity[basename + UUID_ATTRIBUTE_SUFFIX] = file_uuid + SEPARATOR + filename
                entity[basename + URL_ATTRIBUTE_SUFFIX] = file_url
            else:
                return
        else:
            GDC_FILE_ACCESS.recordFileAccessType(basename, access)
            entity[basename + UUID_ATTRIBUTE_SUFFIX] = file_uuid + SEPARATOR + filename
            entity[basename + URL_ATTRIBUTE_SUFFIX] = file_url

def get_file_metadata(gdc_api_root, file_uuid, filename, file_url, known_cases, known_samples, known_pairs, deferred_file_uuids):
    
    # get from GDC the data file's category, type, access type, format, experimental strategy,
    # analysis workflow type
    fileMetadataRetriever = FileMetadataRetriever(gdc_api_root)
    responseDict = fileMetadataRetriever.get_metadata(file_uuid)
    
    try:
        data_category = responseDict['data_category']
        data_type = responseDict['data_type']
        access = responseDict['access']
    except KeyError as x:
        # we expect all files to have at least a data_category, data_type and access type assigned them"
        print("KeyError = ", x)
        print("SKIPPING FILE: file uuid = {0}, file name = {1}".format(file_uuid, filename))
        return
    
    if 'data_format' in responseDict:
        data_format = responseDict['data_format']
    else:
        data_format = None
    if 'experimental_strategy' in responseDict:
        experimental_strategy = responseDict['experimental_strategy']
    else: 
        experimental_strategy = None
    if 'analysis' in responseDict and 'workflow_type' in responseDict['analysis']:
        workflow_type = responseDict['analysis']['workflow_type']
    else:
        workflow_type = None

    if data_category in set([GDC_DataCategory.CLINICAL, GDC_DataCategory.BIOSPECIMEN]): 
        metadataRetriever = CaseMetadataRetriever(gdc_api_root)
    else:
        metadataRetriever = CaseSampleMetadataRetriever(gdc_api_root)

    # quick fix for legacy image data - will clean up
    if data_type in set([GDC_DataType.TISSUE_SLIDE_IMAGE, GDC_DataType.DIAGNOSTIC_IMAGE]):
        metadataRetriever = CaseSampleMetadataRetriever(gdc_api_root)

    metadata = metadataRetriever.get_metadata(file_uuid)

    cases = metadata['cases']
    num_associated_cases = len(cases)
    assert num_associated_cases > 0, file_uuid

    num_associated_samples = 0
    samples = None
    if num_associated_cases == 1 and 'samples' in cases[0]:
        samples = cases[0]['samples']
        num_associated_samples = len(samples)

    # first consider files that are associated with a single case
    if num_associated_cases == 1:
        if num_associated_samples == 0:
            case_id = _add_to_knowncases(cases[0], known_cases)
            _add_file_attribute(gdc_api_root, known_cases[case_id], file_uuid, filename, file_url,
                                data_category, data_type, experimental_strategy, workflow_type, access)
        elif num_associated_samples == 1:
            case_id = _add_to_knowncases(cases[0], known_cases)
            sample_id, _ = _add_to_knownsamples(samples[0], case_id, known_samples)
            _add_file_attribute(gdc_api_root, known_samples[sample_id], file_uuid, filename, file_url, 
                                data_category, data_type, experimental_strategy, workflow_type, access)
        elif num_associated_samples == 2:
            case_id = _add_to_knowncases(cases[0], known_cases)
            sample1_id, sample1_type_tn = _add_to_knownsamples(samples[0], case_id, known_samples)
            sample2_id, sample2_type_tn = _add_to_knownsamples(samples[1], case_id, known_samples)
            assert sample1_type_tn != sample2_type_tn, "not tumor/normal pair"
            if sample1_type_tn == SAMPLE_TYPE.TUMOR:
                tumor_sample_id = sample1_id
                normal_sample_id = sample2_id
            else:
                tumor_sample_id = sample2_id
                normal_sample_id = sample1_id

            pair_id = _add_to_knownpairs(tumor_sample_id, normal_sample_id, known_pairs)
            _add_file_attribute(gdc_api_root, known_pairs[pair_id], file_uuid, filename, file_url,
                                data_category, data_type, experimental_strategy, workflow_type, access)
        else:
            # file associated with more than two samples from a single case
            # not sure how to process this...don't believe there are any such files in GDC
            raise ValueError(file_uuid, filename, data_category)

    else:
        # file associated with multiple cases
        # we will record file_uuid and deal with later
        deferred_file_uuids.append([file_uuid, filename])

# may eventually drop this and incorporate into get_file_metadata.  Wasn't sure what to do with files
# associated with multiple cases or files associated with samples across multiple cases.
# For now we only include files in data model if they are associated with cases and samples that were
# pulled in from other files in manifest that were associated with single cases.
def process_deferred_file_uuid(gdc_api_root, file_uuid, filename, file_url, known_cases, known_samples):
    
    # get data file's name, category, type, access, format experimental strategy, workflow type
    fileMetadataRetriever = FileMetadataRetriever(gdc_api_root)
    responseDict = fileMetadataRetriever.get_metadata(file_uuid)

    data_category = responseDict['data_category']
    data_type = responseDict['data_type']
    access = responseDict['access']

    if 'data_format' in responseDict:
        data_format = responseDict['data_format']
    else:
        data_format = None
    
    if 'experimental_strategy' in responseDict:
        experimental_strategy = responseDict['experimental_strategy']
    else: 
        experimental_strategy = None
    if 'analysis' in responseDict and 'workflow_type' in responseDict['analysis']:
        workflow_type = responseDict['analysis']['workflow_type']
    else:
        workflow_type = None
        
    if data_category == GDC_DataCategory.CLINICAL or data_category == GDC_DataCategory.BIOSPECIMEN:
        metadataRetriever = CaseMetadataRetriever(gdc_api_root)
    else:
        metadataRetriever = CaseSampleMetadataRetriever(gdc_api_root)

    metadata = metadataRetriever.get_metadata(file_uuid)

    cases = metadata['cases']
    num_associated_cases = len(cases)
    assert num_associated_cases > 1, file_uuid

    for case in cases:
        case_id = case['case_id']
        if case_id in known_cases:
            case_id = _add_to_knowncases(case, known_cases)
            if 'samples' in case:
                # associated with samples
                samples = case['samples']
                for sample in samples:
                    sample_id = sample['sample_id']
                    if sample_id in known_samples:
                        _add_file_attribute(gdc_api_root, known_samples[sample_id], file_uuid, filename, file_url,
                                            data_category, data_type, experimental_strategy, workflow_type, access)
            else:
                # associated with multiple cases only
                _add_file_attribute(gdc_api_root, known_cases[case_id], file_uuid, filename, file_url,
                                    data_category, data_type, experimental_strategy, workflow_type, access)

def create_participants_file(cases, manifestFileBasename):
    
    attribute_names = []
    for case_id, case in cases.items():
        for attribute_name in case:
            if attribute_name not in attribute_names:
                attribute_names.append(attribute_name)
    
    participants_filename = manifestFileBasename + '_participants.txt'
    participant_sets_membership_filename = manifestFileBasename + '_participant_sets_membership.txt'
    
    with open(participants_filename, 'w') as participantsFile, open(participant_sets_membership_filename, 'w') as membershipFile:
        fieldnames = ['entity:participant_id'] + attribute_names
        participants_writer = csv.DictWriter(participantsFile, fieldnames=fieldnames, delimiter='\t')
        participants_writer.writeheader()

        fieldnames = ['membership:participant_set_id', 'participant_id']
        membership_writer = csv.DictWriter(membershipFile, fieldnames=fieldnames, delimiter='\t')
        membership_writer.writeheader()
        
        for case_id, case in cases.items():
            entity_row = {'entity:participant_id': case_id}
            for attribute_name in attribute_names:
                if attribute_name in case:
                    entity_row[attribute_name] = case[attribute_name]
                    if attribute_name.endswith(UUID_ATTRIBUTE_SUFFIX):
                        basename = attribute_name[0:-len(UUID_ATTRIBUTE_SUFFIX)]
                        membership_row = {'membership:participant_set_id' : GDC_FILE_ACCESS.getAccessTypePrefix(basename) + basename,
                                          'participant_id' : case_id}
                        membership_writer.writerow(membership_row)
                else:
                    entity_row[attribute_name] = '__DELETE__'
                entity_row[attribute_name] = case[attribute_name] if attribute_name in case else '__DELETE__'
            participants_writer.writerow(entity_row)

            membership_row = {'membership:participant_set_id' : 'ALL',
                              'participant_id' : case_id}
            membership_writer.writerow(membership_row)            

def create_samples_file(samples, manifestFileBasename):
    attribute_names = []
    for sample_id, sample in samples.items():
        for attribute_name in sample:
            if attribute_name not in {'submitter_id', 'case_id', 'sample_type_id'} and attribute_name not in attribute_names:
                attribute_names.append(attribute_name)

    samples_filename = manifestFileBasename + '_samples.txt'
    sample_sets_membership_filename = manifestFileBasename + '_sample_sets_membership.txt'
    with open(samples_filename, 'w') as samplesFile, open(sample_sets_membership_filename, 'w') as membershipFile:
        
        fieldnames = ['entity:sample_id', 'participant_id', 'submitter_id', 'sample_type'] + attribute_names
        sample_writer = csv.DictWriter(samplesFile, fieldnames=fieldnames, delimiter='\t')
        sample_writer.writeheader()

        fieldnames = ['membership:sample_set_id', 'sample_id']
        membership_writer = csv.DictWriter(membershipFile, fieldnames=fieldnames, delimiter='\t')
        membership_writer.writeheader()
        
        for sample_id, sample in samples.items():
            entity_row = {'entity:sample_id' : sample_id, 'participant_id': sample['case_id'],
                          'submitter_id' : sample['submitter_id'],
                          'sample_type' : SAMPLE_TYPE.getLetterCode(sample['sample_type_id'])}
            for attribute_name in attribute_names:
                if attribute_name in sample:
                    entity_row[attribute_name] = sample[attribute_name]
                    if attribute_name.endswith(UUID_ATTRIBUTE_SUFFIX):
                        basename = attribute_name[0:-len(UUID_ATTRIBUTE_SUFFIX)]
                        membership_row = {'membership:sample_set_id' : GDC_FILE_ACCESS.getAccessTypePrefix(basename) + basename,
                                          'sample_id' : sample_id}
                        membership_writer.writerow(membership_row)
                else:
                    entity_row[attribute_name] = '__DELETE__'
            sample_writer.writerow(entity_row)

            membership_row = {'membership:sample_set_id' : 'ALL',
                              'sample_id': sample_id}
            membership_writer.writerow(membership_row)
                        
def create_pairs_file(pairs, samples, manifestFileBasename):
    attribute_names = []
    for pair_id, pair in pairs.items():
        for attribute_name in pair:
            if attribute_name not in {'tumor', 'normal'} and attribute_name not in attribute_names:
                attribute_names.append(attribute_name)

    pairs_filename = manifestFileBasename + '_pairs.txt'
    pair_sets_membership_filename = manifestFileBasename + '_pair_sets_membership.txt'
    with open(pairs_filename, 'w') as pairsFile, open(pair_sets_membership_filename, 'w') as membershipFile:
        fieldnames = ['entity:pair_id', 'participant_id', 'case_sample_id', 'control_sample_id',
                    'tumor_submitter_id', 'normal_submitter_id',
                    'tumor_type', 'normal_type'] + attribute_names
        pairs_writer = csv.DictWriter(pairsFile, fieldnames=fieldnames, delimiter='\t')
        pairs_writer.writeheader()

        fieldnames = ['membership:pair_set_id', 'pair_id']
        membership_writer = csv.DictWriter(membershipFile, fieldnames=fieldnames, delimiter='\t')
        membership_writer.writeheader()
        
        for pair_id, pair in pairs.items():

            tumor_submitter_id = samples[pair['tumor']]['submitter_id']
            normal_submitter_id = samples[pair['normal']]['submitter_id']
            entity_row = {'entity:pair_id' : pair_id,
                          'participant_id' : samples[pair['tumor']]['case_id'],
                          'case_sample_id' : pair['tumor'],
                          'control_sample_id' : pair['normal'],
                          'tumor_submitter_id' : tumor_submitter_id,
                          'normal_submitter_id' : normal_submitter_id,
                          'tumor_type' : SAMPLE_TYPE.getLetterCode(samples[pair['tumor']]['sample_type_id']),
                          'normal_type' : SAMPLE_TYPE.getLetterCode(samples[pair['normal']]['sample_type_id'])}
            for attribute_name in attribute_names:
                if attribute_name in pair:
                    entity_row[attribute_name] = pair[attribute_name]
                    if attribute_name.endswith(UUID_ATTRIBUTE_SUFFIX):
                        basename = attribute_name[0:-len(UUID_ATTRIBUTE_SUFFIX)]
                        membership_row = {'membership:pair_set_id' : GDC_FILE_ACCESS.getAccessTypePrefix(basename) + basename,
                                          'pair_id': pair_id}
                        membership_writer.writerow(membership_row)
                else:
                    entity_row[attribute_name] = '__DELETE__'
            pairs_writer.writerow(entity_row)

            row = {'membership:pair_set_id' : 'ALL',
                   'pair_id': pair_id}
            membership_writer.writerow(row)


def create_workspace_attributes_file(manifestFileBasename, is_legacy):
    #This part is hardcoded due to the small number of attributes we need to specify.
    #Please feel free to change this specification according to your needs.
    legacy_flag="false"
    if is_legacy:
        legacy_flag="true"

    #Due to a somewhat weird bug in FireCloud, please keep the workspace-colunm-defaults attribute as the last one in the list.
    #Any new attributes should be added before workspace-colunm-defaults
    with open(manifestFileBasename + "_workspace_attributes.txt", 'w') as workspaceColumnOrderFile:
        workspaceColumnOrderFile.write("workspace:legacy_flag\tworkspace-column-defaults\n")
        workspaceColumnOrderFile.write(legacy_flag + "\t" + "{\"participant\": {\"shown\": [\"submitter_id\", \"project_id\", \"participant_id\"]}, \"sample\":{\"shown\":[\"submitter_id\", \"sample_id\", \"participant\", \"sample_type\"]}, \"pair\":{\"shown\":[\"tumor_submitter_id\", \"normal_submitter_id\", \"pair_id\"]}}")

def main():
    parser = argparse.ArgumentParser(description='create FireCloud workspace load files from GDC manifest')
    parser.add_argument("manifest", help="manifest file from the GDC Data Portal")
    parser.add_argument("-r", "--resolve_uuids", help="TSV file mapping GDC UUIDs to URLs")
    parser.add_argument("-l", "--legacy", help="point to GDC Legacy Archive", action="store_true")
    args = parser.parse_args()

    print("manifestFile = {0}".format(args.manifest))
    print("tsvFile = {0}".format(args.resolve_uuids))

    manifestFile = args.manifest
    uuidResolver = None
    if args.resolve_uuids is not None:
        uuidResolver = gdc_uuidresolver.UuidResolver(args.resolve_uuids, '__DELETE__')

    pp = pprint.PrettyPrinter()

    cases = dict()
    samples = dict()
    pairs = dict()
    deferred_file_uuids = []

    manifestFileList = _read_manifestFile(manifestFile)

    gdc_api_root = GDC_LEGACY_API_ROOT if args.legacy else GDC_API_ROOT

    for i, item in enumerate(manifestFileList):

        file_uuid = item['id']
        filename = item['filename']
        file_url = uuidResolver.getURL(file_uuid) if uuidResolver is not None else "__DELETE__"
    
        print('{0} of {1}: {2}, {3}'.format(i+1, len(manifestFileList), file_uuid, filename))

        for attempt in range(5):
            try:
                get_file_metadata(gdc_api_root, file_uuid, filename, file_url, cases, samples, 
                                    pairs, deferred_file_uuids)
                break
            except (KeyboardInterrupt, SystemExit):
                raise
            except Exception as x:
                print("Exception=", x)
                print("attempt=", attempt, 'file uuid = ', file_uuid)
                time.sleep(5)
        else:
            #failed all attempts
            # - just move on
            print("failed 5 attempts! SKIPPING FILE: file uuid = ", file_uuid)

    for uuid_and_filename in deferred_file_uuids:
        file_uuid = uuid_and_filename[0]
        filename = uuid_and_filename[1]
        file_url =  file_url = uuidResolver.getURL(file_uuid) if uuidResolver is not None else "__DELETE__"

        for attempt in range(5):
            try:
                process_deferred_file_uuid(gdc_api_root, file_uuid, filename, file_url, cases, samples)
            except (KeyboardInterrupt, SystemExit):
                raise
            except Exception as x:
                print("Exception=", x)
                print("attempt=", attempt, 'file uuid = ', file_uuid)
                time.sleep(5)
                raise
            else:
                break
        else:
            #failed all attempts
            # - just move on
            print("failed 5 attempts! SKIPPING FILE: file uuid = ", file_uuid)
            continue

    manifestFileBasename = os.path.splitext(os.path.basename(manifestFile))[0]

    create_participants_file(cases, manifestFileBasename)
    create_samples_file(samples, manifestFileBasename)
    if len(pairs) != 0:
        create_pairs_file(pairs, samples, manifestFileBasename)

    #This part creates a file that specifies the workspace attributes. 
    #The attributes are:
    # 1.Default order of columns when shown in the workspace.
    # 2.Whether the workspace is meant to deal with data fom the legacy site or not.
    create_workspace_attributes_file(manifestFileBasename, args.legacy)
    

if __name__ == '__main__':
    main()
