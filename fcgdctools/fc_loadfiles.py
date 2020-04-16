import csv
import requests
import argparse
import pprint
import os.path
import sys
import time
import traceback

UUID_TO_FILENAME = dict()

DEFERRED_FILE_NUM_OF_CASES = dict()

GDC_API_ROOT = "https://api.gdc.cancer.gov"
GDC_LEGACY_API_ROOT = "https://api.gdc.cancer.gov/legacy"

#program
class GDC_ProgramName:
    TARGET = 'TARGET'
    TCGA = 'TCGA'
    FM = 'FM'
    CPTAC = 'CPTAC'

#data categories                                                                                                                   
class GDC_DataCategory:
    SNV = "Simple Nucleotide Variation"
    TRANSCRIPTOME_PROFILING = "Transcriptome Profiling"
    BIOSPECIMEN = "Biospecimen"
    RAW_SEQUENCING_DATA = "Raw Sequencing Data"    
    COPY_NUMBER_VARIATION = "Copy Number Variation"    
    CLINICAL = "Clinical"
    DNA_METHYLATION = "DNA Methylation"
    COMBINED_NUCLEOTIDE_VARIATION = "Combined Nucleotide Variation"
    SEQUENCING_READS = "Sequencing Reads"

    LEGACY_SNV = "Simple nucleotide variation"

#data types                                                                                                                        
class GDC_DataType:
    
    #associated with Simple Nucleotide Variation data category                                                                     
    RAW_SIMPLE_SOMATIC_MUTATION = "Raw Simple Somatic Mutation"
    ANNOTATED_SOMATIC_MUTATION = "Annotated Somatic Mutation"
    AGGREGATED_SOMATIC_MUTATION = "Aggregated Somatic Mutation"
    MASKED_SOMATIC_MUTATION = "Masked Somatic Mutation"  
    
    #associated with Transcriptome Profiling data category                                                                         
    MIRNA_EXPRESSION_QUANTIFICATION = "miRNA Expression Quantification"
    ISOFORM_EXPRESSION_QUANTIFICATION = "Isoform Expression Quantification"
    GENE_EXPRESSION_QUANTIFICATION = "Gene Expression Quantification"

    #associated with Biospecimen data category                                                                                     
    BIOSPECIMEN_SUPPLEMENT = "Biospecimen Supplement"
    SLIDE_IMAGE = "Slide Image"

    #associated with RAW Sequencing Data data category                                                                             
    ALIGNED_READS = "Aligned Reads"
    
    #associated with Copy Number Variation data category                                                                           
    COPY_NUMBER_SEGMENT = 'Copy Number Segment'
    MASKED_COPY_NUMBER_SEGMENT = 'Masked Copy Number Segment'

    #data types associated with Clinical data category                                                                             
    CLINICAL_SUPPLEMENT = "Clinical Supplement"

    LEGACY_TISSUE_SLIDE_IMAGE = "Tissue slide image"
    LEGACY_DIAGNOSTIC_IMAGE = "Diagnostic image"
    LEGACY_PATHOLOGY_REPORT = "Pathology report"
    LEGACY_CLINICAL_DATA = "Clinical data"
    LEGACY_BIOSPECIMEN_DATA = "Biospecimen data"
    
    #associated with DNA Methylation data category                                                                                 
    METHYLATION_BETA_VALUE = "Methylation Beta Value"
    
    #associated with Combined Nucleotide Variation data category
    RAW_CGI_VARIANT = "Raw CGI Variant"

    LEGACY_SIMPLE_NUCLEOTIDE_VARIATON = "Simple nucleotide variation"

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
        'Genotyping Array' : 'GeneArray',
        'Targeted Sequencing' : 'TargetedSeq',
        'miRNA-Seq' : 'miRNAseq',
        'Tissue Slide': 'TissueSlide',
        'Methylation Array' : 'MethArray',
        'Diagnostic Slide' : 'DiagSlide',
        'WGS': 'WGS'} 

EXP_STRATEGY = DataSource(EXP_STRATEGY_ABBREVIATIONS)

WORKFLOW_ABBREVIATIONS = {
        'DNACopy' : 'DNACopy',        
        'BCGSC miRNA Profiling' : 'BCGSCmiRNA',
        'BWA with Mark Duplicates and Cocleaning' : 'BWAMDupCoClean',
        'FM Simple Somatic Mutation' : 'FMSimpleSomaticMutation',
        'FoundationOne Annotation' : 'F1Annotation',
        'Liftover' : 'Lift',
        'STAR 2-Pass' : 'STAR2Pass',        
        'HTSeq - Counts' : 'HTSeqCounts',
        'HTSeq - FPKM' : 'HTSeqFPKM',
        'HTSeq - FPKM-UQ' : 'HTSEQFPKMUQ',
        'BWA-aln' : 'BWAaln',      
        'SomaticSniper': 'SomSnip',
        'SomaticSniper Annotation' : 'SomSnipAnnot',
        'MuTect2' : 'MuTect2',
        'MuTect2 Annotation' : 'MuTect2Annot',       
        'VarScan2' : 'VarScan2',
        'VarScan2 Annotation' : 'VarScan2Annot',    
        'MuSE' : 'MuSE',
        'MuSE Annotation' : 'MuSEAnnot',
        'VCF LiftOver' : 'VCFLift',
        'MuSE Variant Aggregation and Masking' : 'MuSEAggrMask',
        'MuTect2 Variant Aggregation and Masking' : 'MuTect2AggrMask',
        'SomaticSniper Variant Aggregation and Masking' : 'SomSnipAggrMask',
        'VarScan2 Variant Aggregation and Masking' : 'VarScan2AggrMask',
        'FoundationOne Variant Aggregation and Masking' : 'F1AggrMask'}


WORKFLOW = DataSource(WORKFLOW_ABBREVIATIONS)

PLATFORM_ABBREVIATIONS = {
        'Affymetrix SNP 6.0' : 'AffySNP6',
        'Illumina' : 'Illum',
        'Illumina Human Methylation 450' : 'IllumHuMeth450',
        'Illumina Human Methylation 27' : 'IllumHuMeth27'}

PLATFORM = DataSource(PLATFORM_ABBREVIATIONS)

# Sample Types                                                                                                                     
# from https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes                                              
class SampleType:
    TUMOR = 'tumor'
    NORMAL = 'normal'
    OTHER = 'other'
    NA = 'na'
    # this table needs updating - TARGET data had some unknown sample type ids                                                     
    GDC_SAMPLE_TYPE_IDS = {'01': ['Primary Solid Tumor', 'TP', TUMOR],
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

    GDC_TISSUE_TYPES = {'Tumor': TUMOR,
                        'Normal': NORMAL,
                        'Abnormal': OTHER,
                        'Peritumoral': OTHER,
                        'Unknown': NA,
                        'Not Reported': NA}

    def getTumorNormalClassification(self, gdc_tissue_type, gdc_sample_type, gdc_sample_type_id):
        if gdc_sample_type_id is not None:
            return self.GDC_SAMPLE_TYPE_IDS[gdc_sample_type_id][2]
        elif gdc_tissue_type is not None:
            return self.GDC_TISSUE_TYPES[gdc_tissue_type]
        elif gdc_sample_type is not None:
            if 'Normal' in gdc_sample_type:
                return self.NORMAL
            else:
                return self.TUMOR
        else:
            return self.NA

    def getLetterCode(self, sample_type_id):
        if sample_type_id is not None:
            return self.GDC_SAMPLE_TYPE_IDS[sample_type_id][1]
        else:
            return None

SAMPLE_TYPE = SampleType()

class MetadataRetriever():
    def __init__(self, api_endpoint, fields, gdc_api_root):
        self.gdc_api_root = gdc_api_root
        self.fields = fields
        self.api_endpoint = api_endpoint

    def get_metadata(self, uuid):
        url = "{0}/{1}/{2}?fields={3}".format(self.gdc_api_root, self.api_endpoint, uuid, self.fields)
        #debug
        #print('url: {0}'.format(url))
        response = requests.get(url, headers=None, timeout=5)
        responseDict = response.json()
        return responseDict['data']

class FileCaseMetadataRetriever(MetadataRetriever):
    def __init__(self, gdc_api_root):
        fields = "cases.case_id,cases.submitter_id,cases.project.project_id,cases.tissue_source_site"
        MetadataRetriever.__init__(self, 'files', fields, gdc_api_root)

class FileCaseSampleMetadataRetriever(MetadataRetriever):
    def __init__(self, gdc_api_root):
        fields = "cases.case_id,cases.submitter_id,cases.project.project_id"
        fields = fields + ",cases.samples.sample_id,cases.samples.submitter_id,cases.samples.sample_type_id,cases.samples.sample_type,cases.samples.tissue_type"
        fields = fields + ",cases.samples.portions.analytes.aliquots.submitter_id"
        MetadataRetriever.__init__(self, 'files', fields, gdc_api_root)

class FileMetadataRetriever(MetadataRetriever):
    def __init__(self, gdc_api_root):
        fields = "data_category,data_type,data_format,access,experimental_strategy,analysis.workflow_type,cases.project.program.name"
        MetadataRetriever.__init__(self, 'files', fields, gdc_api_root)

class CaseMetadataRetriever(MetadataRetriever):
    def __init__(self, gdc_api_root):
        fields = "primary_site"
        MetadataRetriever.__init__(self, 'cases', fields, gdc_api_root)
        
class IndexFileUuidRetriever():
    def __init__(self, gdc_api_root):
        self.gdc_api_root = gdc_api_root

    def get_index_uuid(self, bam_uuid):
        url = "{0}/files/{1}?expand=index_files".format(self.gdc_api_root, bam_uuid)
        response = requests.get(url, headers=None, timeout=5)
        responseDict = response.json()
        indexFilesList = responseDict['data']['index_files']
        assert(len(indexFilesList) == 1)
        return(indexFilesList[0]['file_id'])
       
SEPARATOR = '/'
DRS_URL_ATTRIBUTE_SUFFIX = "drs_url"

def _read_manifestFile(manifestFile):

    manifestFileList = []

    with open(manifestFile, 'r') as fp:
        reader = csv.DictReader(fp, delimiter='\t')
        for row in reader:
            manifestFileList.append(row)

    return manifestFileList

def _add_to_knowncases(case_metadata, known_cases, gdc_api_root):
    case_id = case_metadata['case_id']
    if case_id not in known_cases:
        submitter_id = case_metadata['submitter_id']
        project_id = case_metadata['project']['project_id']

        caseMetadataRetriever = CaseMetadataRetriever(gdc_api_root)
        caseMetadata = caseMetadataRetriever.get_metadata(case_id)
        primary_site = caseMetadata.get('primary_site')

        new_case = {'submitter_id' : submitter_id, 'project_id' : project_id, 'primary_site' : primary_site}
        known_cases[case_id] = new_case
    return case_id

def _add_to_knownsamples(sample_metadata, case_id, known_samples):
    sample_id = sample_metadata['sample_id']
    tissue_type = sample_metadata['tissue_type']
    sample_type = sample_metadata['sample_type']
    sample_type_id = sample_metadata.get('sample_type_id')

    if sample_id not in known_samples:
        sample_submitter_id = sample_metadata['submitter_id']
        new_sample = {'submitter_id' : sample_submitter_id,
                      'tissue_type': tissue_type,
                      'sample_type': sample_type,
                      'sample_type_id': sample_type_id,
                      'case_id' : case_id}
        known_samples[sample_id] = new_sample
    return sample_id, SAMPLE_TYPE.getTumorNormalClassification(tissue_type, sample_type, sample_type_id)

def _get_sample_type(sample_metadata):
    tissue_type = sample_metadata['tissue_type']
    sample_type = sample_metadata['sample_type']
    sample_type_id = sample_metadata.get('sample_type_id')
    return SAMPLE_TYPE.getTumorNormalClassification(tissue_type, sample_type, sample_type_id)

def _add_pooled_sample_to_knownsamples(file_uuid, filename, samples_metadata, case_id, known_samples):
    sample_ids = [sample_metadata['sample_id'] for sample_metadata in samples_metadata]
    submitter_ids = [sample_metadata['submitter_id'] for sample_metadata in samples_metadata]
    tissue_types = [sample_metadata.get('tissue_type') for sample_metadata in samples_metadata]
    sample_types = [sample_metadata.get('sample_type') for sample_metadata in samples_metadata]
    sample_type_ids = [sample_metadata.get('sample_type_id') for sample_metadata in samples_metadata]

    if len(set(tissue_types)) != 1:
        print('inconsistent tissue types in polled sample: {0}'.format(tissue_types))
        raise ValueError(file_uuid, filename)
    if len(set(sample_types)) != 1:
        print('inconsistent sample types in polled sample: {0}'.format(sample_types))
        raise ValueError(file_uuid, filename)
    if len(set(sample_type_ids)) != 1:
        print('inconsistent sample type ids in polled sample: {0}'.format(sample_type_ids))
        raise ValueError(file_uuid, filename)

    pooled_sample_id = '__'.join(sorted(sample_ids))
    pooled_sample_submitter_id = '__'.join(sorted(submitter_ids))

    if pooled_sample_id not in known_samples:
        new_sample = {'submitter_id' : pooled_sample_submitter_id,
                      'tissue_type': tissue_types[0],
                      'sample_type': sample_types[0],
                      'sample_type_id': sample_type_ids[0],
                      'case_id' : case_id}
        known_samples[pooled_sample_id] = new_sample
    return pooled_sample_id, SAMPLE_TYPE.getTumorNormalClassification(tissue_types[0], sample_types[0], sample_type_ids[0])
        

def _add_to_knownpairs(tumor_sample_id, normal_sample_id, known_pairs):
    pair_id = "{0}_{1}".format(tumor_sample_id, normal_sample_id)
    if pair_id not in known_pairs:
        known_pairs[pair_id] = {'tumor': tumor_sample_id, 'normal': normal_sample_id}
    return pair_id

def _constructAttributeName_base(experimental_strategy, workflow_type, data_category, data_type, data_format):

    if experimental_strategy is not None:
        experimental_strategy_abbrev = EXP_STRATEGY.getAbbreviation(experimental_strategy) + '__'
    else:
        experimental_strategy_abbrev = ''

    if workflow_type is not None:
        workflow_type_abbrev = WORKFLOW.getAbbreviation(workflow_type) + '__'
    else:
        workflow_type_abbrev = ''

    data_type_lc = data_type.lower().replace(" ", "_") + '__'
    
    data_format_lc = data_format.lower().replace(" ", "_") + '__'

    attribute_name_base = experimental_strategy_abbrev + workflow_type_abbrev + data_type_lc + data_format_lc

    return (attribute_name_base)

def _getImageCodeAndPortionFromImageFilename(filename):
    image_code = filename.split('.')[0].split('-')[-1]
    portion = int(filename.split('.')[0].split('-')[-2])
    return image_code, portion

def _constructImageAttributeName_base(experimental_strategy, workflow_type, data_category, data_type, data_format, filename=None):

    if experimental_strategy is not None:
        experimental_strategy_abbrev = EXP_STRATEGY.getAbbreviation(experimental_strategy) + '__'
    else:
        experimental_strategy_abbrev = ''

    if workflow_type is not None:
        workflow_type_abbrev = WORKFLOW.getAbbreviation(workflow_type) + '__'
    else:
        workflow_type_abbrev = ''

    data_type_lc = data_type.lower().replace(" ", "_") + '__' 
    data_format_lc = data_format.lower().replace(" ", "_") + '__'
    
    # see https://wiki.nci.nih.gov/display/TCGA/TCGA+barcode# for interpretation of TCGA bar code
    # that is incorporated into image filename
    image_code, portion = _getImageCodeAndPortionFromImageFilename(filename)
    image_code_lc = image_code.lower() + '__'

    attribute_name_base = experimental_strategy_abbrev + workflow_type_abbrev + image_code_lc + data_type_lc

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
    
# this seemed clearer to me (Chet), but felt it better to use the same code
# used by GDAC
def _pick_tcga_aliquot(a, b):
    """Comparator function for aliquot barcodes.

    Uses rules described in the GDAC FAQ entry for replicate samples: FAQ entry for replicate samples: 
    https://confluence.broadinstitute.org/display/GDAC/FAQ                                        
    """
    # Get the analytes and plates - 
    # _see https://docs.gdc.cancer.gov/Encyclopedia/pages/images/TCGA-TCGAbarcode-080518-1750-4378.pdf

    # TCGA-XX-XXXX-XXV-XXA-XXXX-XX
    # 012345678901234567890123456789012345678901234567890123456789
    #           1         2         3         4         5
    # tss = 5:7
    # participant = 8:12
    # sample = 13:15
    # vial = 15:17
    # portion = 17:19
    # analyte = 19:20
    # plate = 12:25
    # center = 26:28

    # see https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/portion-analyte-codes for analyte codes

    analyte1 = a[19]
    analyte2 = b[19]
    plate1   = a[21:25]
    plate2   = b[21:25]
    rna_analytes = set(['H', 'R', 'T'])
    dna_analytes = set(['D', 'W', 'X']) 

    if a == b:
        return a

    elif analyte1 == analyte2:
        # Prefer the aliquot with the highest lexicographical sort value
        return a if a >= b else b

    elif analyte1 in rna_analytes and analyte2 in rna_analytes:
        # Prefer H over R and T
        if analyte1 == 'H' or analyte2 == 'H':
            return a if analyte1 == 'H' else b
        # Prefer R over T    
        assert analyte == 'R' or analyte2 == 'R'
        return a if analyte1 == 'R' else b

    elif analyte1 in dna_analytes and analyte2 in dna_analytes:
        # Prefer D over G,W,X, unless plate number is higher
        if analyte1 == 'D' or analyte2 == 'D':
            if plate1 == plate2:
                return a if analyte1 == 'D' else b
            else:
                return a if plate2 < plate1 else b
    else:
        # Default back to highest lexicographical sort value     
        return a if a >= b else b


def _pick_tcga_aliquot_pair(aliquot_pair_1, aliquot_pair_2):
    
    tumor_aliquot_1 = aliquot_pair_1['tumor']
    tumor_aliquot_2 = aliquot_pair_2['tumor']
    normal_aliquot_1 = aliquot_pair_1['normal']
    normal_aliquot_2 = aliquot_pair_2['normal']

    if tumor_aliquot_1 != tumor_aliquot_2:
        tumor_aliquot_choice = _pick_tcga_submitter(tumor_aliquot_1, tumor_aliquot_2)
        return aliquot_pair_1 if tumor_aliquot_choice == tumor_aliquot_1 else aliquot_pair_2
    elif normal_aliquot_1 != normal_aliquot_2:
        normal_aliquot_choice = _pick_tcga_submitter(tumor_aliquot_1, tumor_aliquot_2)
        return aliquot_pair_1 if normal_aliquot_choice == normal_aliquot_1 else aliquot_pair_2
    else:
        print("WARNING: aliquot ids are identical, unable to make rational choice")
        return aliquot_pair_1

def _pick_target_submitter(a,b):
    # Get the analytes and plates                                                                                                  
    # TARGET-##-TSS-ABCDEF-TS.TP.N-<portion><analyte>                                                                              

    analyte1 = a[-1]
    analyte2 = b[-1]
    portion1 = a[-3:-1]
    portion2 = b[-3:-1]

    # Equals case                                                                                                                  
    if a == b:
        return a
    elif analyte1 == analyte2:
        # Prefer the aliquot with the higher portion value                                                                         
        return a if portion1 >= portion2 else b
    #If RNA then analyte codes can be only R or S                                                                                  
    elif analyte1 == "S":
        #Then analyte2 has to be R (because they can't be equal)                                                                   
        return b
    elif analyte1 == "R":
        #Then analyte2 has to be S (because they can't be equal)                                                                   
        return a
    elif analyte1 == "D":
        #prefer D to E,W,X,Y                                                                                                       
        return a
    elif analyte1 == "E":
        return b if analyte2 == "D" else a
    elif analyte1 == "W":
        # "W" is the lowest priority                                                                                               
        return b
    elif analyte1 == "X":
        return a if analyte2 == "W" else b
    else:
        # analyte1 is Y                                                                                                            
        return a if analyte2 == "X" or analyte2 == "W" else b
    
def _pick_target_aliquot_pair(aliquot_pair_1, aliquot_pair_2):
    
    tumor_aliquot_1 = aliquot_pair_1['tumor']
    tumor_aliquot_2 = aliquot_pair_2['tumor']
    normal_aliquot_1 = aliquot_pair_1['normal']
    normal_aliquot_2 = aliquot_pair_2['normal']

    if tumor_aliquot_1 != tumor_aliquot_2:
        tumor_aliquot_choice = _pick_target_submitter(tumor_aliquot_1, tumor_aliquot_2)
        return aliquot_pair_1 if tumor_aliquot_choice == tumor_aliquot_1 else aliquot_pair_2
    elif normal_aliquot_1 != normal_aliquot_2:
        normal_aliquot_choice = _pick_target_submitter(tumor_aliquot_1, tumor_aliquot_2)
        return aliquot_pair_1 if normal_aliquot_choice == normal_aliquot_1 else aliquot_pair_2
    else:
        print("WARNING: aliquot ids are identical, unable to make rational choice")
        return aliquot_pair_1

def _resolve_collision(data_category, data_type, program, uuid1, name1, uuid2, name2, gdc_api_root):

    # NOTE: we chose not to employ the created_datetime or updated_datetime fields in 
    # our decision logic.  From what we can tell, neither should be used to make a selection between 
    # two files.


    # Files that are associated with multiple cases won't use information encoded in aliquot barcode, 
    if uuid1 in DEFERRED_FILE_NUM_OF_CASES or uuid2 in DEFERRED_FILE_NUM_OF_CASES:
        if program == GDC_ProgramName.TARGET and data_category in [GDC_DataCategory.CLINICAL, GDC_DataCategory.BIOSPECIMEN]:
            # special-case logic to deal with TARGET clinical and biospecimin files
            # select file associated with Discovery cohort over file associated with Validation cohort
            # where cohort association is encoded in the filename
            DISCOVERY = "Discovery"
            VALIDATION = "Validation"
            if DISCOVERY in name1 and VALIDATION in name2:
                return uuid1, name1
            if VALIDATION in name1 and DISCOVERY in name2:
                return uuid2, name2

            print("No criteria for selection - skip new file1")
            return uuid2, name2

        if uuid1 in DEFERRED_FILE_NUM_OF_CASES:
            num_cases1 = DEFERRED_FILE_NUM_OF_CASES[uuid1]
        else:
            num_cases1 = 1

        if uuid2 in DEFERRED_FILE_NUM_OF_CASES:
            num_cases2 = DEFERRED_FILE_NUM_OF_CASES[uuid2]
        else:
            num_cases2 = 1

        print("number of cases...")
        print("{0}: {1}, {2}: {3}".format(uuid1,num_cases1,uuid2,num_cases2))
        # If one of the two files has more cases associated with it, we assume it's the correct file to pick.
        if num_cases1 >= num_cases2:
            return uuid1, name1
        else:
            return uuid2, name2

    # Now we are left to deal only with files that are associated with one case
    # get aliquot submitter IDs
    tumor_aliquot_submitter_id1 = None
    tumor_aliquot_submitter_id2 = None
    normal_aliquot_submitter_id1 = None
    normal_aliquot_submitter_id2 = None

    # SNV and Combined Nucleotide Variation (TARGET only) files are associated with two samples: tumor and normal. 
    if ((data_category in GDC_DataCategory.SNV and 
         data_type not in set([GDC_DataType.AGGREGATED_SOMATIC_MUTATION, GDC_DataType.MASKED_SOMATIC_MUTATION])) or
        (data_category in GDC_DataCategory.COMBINED_NUCLEOTIDE_VARIATION) or
        (data_category in GDC_DataCategory.LEGACY_SNV and
         data_type in GDC_DataType.LEGACY_SIMPLE_NUCLEOTIDE_VARIATION)):

        file_fields = "cases.samples.sample_type,cases.samples.portions.analytes.aliquots.submitter_id,cases.samples.sample_type_id"
        meta_retriever = MetadataRetriever('files', file_fields, gdc_api_root)

        data1 = meta_retriever.get_metadata(uuid1)
        samples_list1 = data1['cases'][0]['samples']
        assert len(samples_list1) == 2
        for s in samples_list1:
            sample_type_tn = SAMPLE_TYPE.getTumorNormalClassification(s.get('tissue_type'), s.get('sample_type'), s.get('sample_type_id'))
            if sample_type_tn == SAMPLE_TYPE.TUMOR:
                tumor_aliquot_submitter_id1 = s['portions'][0]['analytes'][0]['aliquots'][0]['submitter_id']
            else:
                assert sample_type_tn == SAMPLE_TYPE.NORMAL, "expected normal sample type"
                normal_aliquot_submitter_id1 = s['portions'][0]['analytes'][0]['aliquots'][0]['submitter_id']

        data2 = meta_retriever.get_metadata(uuid2)
        samples_list2 = data2['cases'][0]['samples']
        assert len(samples_list2) == 2
        for s in samples_list2:
            sample_type_tn = SAMPLE_TYPE.getTumorNormalClassification(s.get('tissue_type_'), s.get('sample_type'), s.get('sample_type_id'))
            if sample_type_tn == SAMPLE_TYPE.TUMOR:
                tumor_aliquot_submitter_id2 = s['portions'][0]['analytes'][0]['aliquots'][0]['submitter_id']
            else:
                assert sample_type_tn == SAMPLE_TYPE.NORMAL, "expected normal sample type"
                normal_aliquot_submitter_id2 = s['portions'][0]['analytes'][0]['aliquots'][0]['submitter_id']

        aliquot_pair_1 = {'tumor' : tumor_aliquot_submitter_id1, 'normal' : normal_aliquot_submitter_id1}
        aliquot_pair_2 = {'tumor' : tumor_aliquot_submitter_id2, 'normal' : normal_aliquot_submitter_id2}
        
        print('aliquot pair name for {0}: {1} / {2}'.format(uuid1, aliquot_pair_1['tumor'], aliquot_pair_1['normal']))
        print('aliquot pair pair for {0}: {1} / {2}'.format(uuid2, aliquot_pair_2['tumor'], aliquot_pair_2['normal']))

        if program == GDC_ProgramName.TCGA:
            chosen_aliquot_pair = _pick_tcga_aliquot_pair(aliquot_pair_1, aliquot_pair_2)
        elif program == GDC_ProgramName.TARGET:
            chosen_aliquot_pair = _pick_target_aliquot_pair(aliquot_pair_1, aliquot_pair_2)
        else:
            # no known structure of metadata encoded in aliquot name; just choose 1 arbitrarily
            print('WARNING: no known structure of metadata encoded in aliquote name; choice is arbitrary!')
            chosen_aliquot_pair = aliquot_pair_1
        
        if chosen_aliquot_pair == aliquot_pair_1:
            return uuid1, name1
        else:
            return uuid2, name2

    # Here we handle other file types that are associated with single sample.
    else:
        file_fields = "cases.project.program.name,cases.samples.portions.analytes.aliquots.submitter_id"
        meta_retriever = MetadataRetriever('files', file_fields, gdc_api_root)

        data1 = meta_retriever.get_metadata(uuid1)
        data2 = meta_retriever.get_metadata(uuid2)

        assert len(data1['cases'][0]['samples']) == 1, "more than one sample associated with file"
        assert len(data2['cases'][0]['samples']) == 1, "more than one sample associated with file"

        aliquot_submitter_id1 = data1['cases'][0]['samples'][0]['portions'][0]['analytes'][0]['aliquots'][0]['submitter_id']
        aliquot_submitter_id2 = data2['cases'][0]['samples'][0]['portions'][0]['analytes'][0]['aliquots'][0]['submitter_id']

        print('aliquot name for {0}: {1}'.format(uuid1, aliquot_submitter_id1))
        print('aliquot name for {0}: {1}'.format(uuid2, aliquot_submitter_id2))

        if aliquot_submitter_id1 != aliquot_submitter_id2:
            if program == GDC_ProgramName.TCGA:
                chosen = _pick_tcga_submitter(aliquot_submitter_id1, aliquot_submitter_id2)

            elif program == GDC_ProgramName.TARGET:
                chosen = _pick_target_submitter(aliquot_submitter_id1, aliquot_submitter_id2)
            else:
                print('WARNING: no known structure of metadata encoded in aliquote name; choice is arbitrary!')
                chosen = aliquot_submitter_id1
        else:
            print('WARNING: aliquot ids are identical; unable to make rational choice!')
            chosen = aliquot_submitter_id1

        if chosen == aliquot_submitter_id1:
            return uuid1, name1
        else:
            return uuid2, name2


def _create_drs_url(file_uuid):
    return 'drs://dataguids.org/{0}'.format(file_uuid)

def _get_file_uuid_from_drs_url(drs_url):
    return drs_url.partition('drs://dataguids.org/')[2]

def _add_file_attribute(entity_id, entity, file_uuid, filename,
                        data_category, data_type, data_format, experimental_strategy, workflow_type, access, program, gdc_api_root):
    # I needed to insert some special-case processing for image data files
    # this probably isn't the cleanest way to handle it, but good enough for now
    if data_type in set([GDC_DataType.SLIDE_IMAGE]):
        basename, portion = _constructImageAttributeName_base(experimental_strategy, workflow_type,
                                                              data_category, data_type, data_format, filename)
        attribute_name = basename + DRS_URL_ATTRIBUTE_SUFFIX

        if attribute_name in entity:
            existing_drs_url = entity[attribute_name]
            existing_uuid = _get_file_uuid_from_drs_url(existing_drs_url)

            existing_filename = UUID_TO_FILENAME[existing_uuid]

            print("multiple files for same attribute!") 
            print("entity id: {0}, attribute name: {1}".format(entity_id, attribute_name))
            print("new file: {0}/{1}".format(file_uuid, filename))
            print("existing file: {0}/{1}".format(existing_uuid, existing_filename))

            _, portion_present = _getImageCodeAndPortionFromImageFilename(existing_filename)
            if portion > portion_present:
                print("newer file has larger portion ID; use newer file")
                entity[basename + DRS_URL_ATTRIBUTE_SUFFIX] = _create_drs_url(file_uuid)
            elif portion < portion_present:
                print("newer file has smaller portion ID; retain existing file")
            else:
                print("Both files have samer portion ID: retain existing file")

        else:
            entity[basename + DRS_URL_ATTRIBUTE_SUFFIX] = _create_drs_url(file_uuid)
    else:
        basename = _constructAttributeName_base(experimental_strategy, workflow_type,
                                                data_category, data_type, data_format)
        attribute_name = basename + DRS_URL_ATTRIBUTE_SUFFIX
        
        # see if attribute already defined for entity
        if attribute_name in entity:
            existing_drs_url = entity[attribute_name]
            existing_uuid = _get_file_uuid_from_drs_url(existing_drs_url)
            existing_filename = UUID_TO_FILENAME[existing_uuid]

            print("multiple files for same attribute!")
            print("entity id: {0}, attribute name: {1}".format(entity_id, attribute_name))
            print("new file: {0}/{1}".format(file_uuid, filename))
            print("existing file: {0}/{1}".format(existing_uuid, existing_filename))
            
            chosen_uuid, chosen_filename = _resolve_collision(data_category, data_type, program,
                                                              file_uuid, filename, existing_uuid, existing_filename, gdc_api_root)
            print("chosen file is: {0}/{1}".format(chosen_uuid, chosen_filename))


            if chosen_uuid == file_uuid:
                entity[basename + DRS_URL_ATTRIBUTE_SUFFIX] = _create_drs_url(file_uuid)
                print("experimental strategy: {0}".format(experimental_strategy))
                # GDC does not provide index files for RNA-Seq BAMs
                if data_format == 'BAM'and experimental_strategy != 'RNA-Seq':
                    bai_basename = basename.replace('__bam__', '__bai__')
                    indexFileUuidRetriever = IndexFileUuidRetriever(gdc_api_root)
                    bai_uuid = indexFileUuidRetriever.get_index_uuid(file_uuid)
                    entity[bai_basename + DRS_URL_ATTRIBUTE_SUFFIX] = _create_drs_url(bai_uuid)

            else:
                return
        else:
            entity[basename + DRS_URL_ATTRIBUTE_SUFFIX] = _create_drs_url(file_uuid)

            # GDC does not provide index files for RNA-Seq BAMs
            if data_format == 'BAM'and experimental_strategy != 'RNA-Seq':
                bai_basename = basename.replace('__bam__', '__bai__')
                indexFileUuidRetriever = IndexFileUuidRetriever(gdc_api_root)
                bai_uuid = indexFileUuidRetriever.get_index_uuid(file_uuid)
                entity[bai_basename + DRS_URL_ATTRIBUTE_SUFFIX] = _create_drs_url(bai_uuid)

def get_file_metadata(file_uuid, filename, known_cases, known_samples, known_pairs, deferred_file_uuids, gdc_api_root):
    
    pp = pprint.PrettyPrinter()

    # get from GDC the data file's category, type, access type, format, experimental strategy,
    # analysis workflow type
    fileMetadataRetriever = FileMetadataRetriever(gdc_api_root)
    responseDict = fileMetadataRetriever.get_metadata(file_uuid)
    
    try:
        data_category = responseDict['data_category']
        data_type = responseDict['data_type']
        data_format = responseDict['data_format']
        access = responseDict['access']
        program = responseDict['cases'][0]['project']['program']['name']
    except KeyError as x:
        # we expect all files to have at least a data_category, data_type, access type and program assigned to them
        print("KeyError = ", x)
        print("SKIPPING FILE: file uuid = {0}, file name = {1}".format(file_uuid, filename))
        return
    
    if 'experimental_strategy' in responseDict:
        experimental_strategy = responseDict['experimental_strategy']
    else: 
        experimental_strategy = None
    if 'analysis' in responseDict and 'workflow_type' in responseDict['analysis']:
        workflow_type = responseDict['analysis']['workflow_type']
    else:
        workflow_type = None
    

    if data_category in set([GDC_DataCategory.CLINICAL, GDC_DataCategory.BIOSPECIMEN]): 
        if data_type == GDC_DataType.SLIDE_IMAGE:
            fileMetadataRetriever = FileCaseSampleMetadataRetriever(gdc_api_root)            
        else:
            fileMetadataRetriever = FileCaseMetadataRetriever(gdc_api_root)
    else:
        fileMetadataRetriever = FileCaseSampleMetadataRetriever(gdc_api_root)

    # quick fix for legacy image data - will clean up
    if data_type in set([GDC_DataType.LEGACY_TISSUE_SLIDE_IMAGE, GDC_DataType.LEGACY_DIAGNOSTIC_IMAGE]):
        fileMetadataRetriever = FileCaseSampleMetadataRetriever(gdc_api_root)


    fileMetadata = fileMetadataRetriever.get_metadata(file_uuid)

    #debug
    #print('metadata:')
    #pp.pprint(fileMetadata)

    cases = fileMetadata['cases']

    num_associated_cases = len(cases)
    assert num_associated_cases > 0, file_uuid

    num_associated_samples = 0
    samples = None
    if num_associated_cases == 1 and 'samples' in cases[0]:
        samples = cases[0]['samples']

        num_associated_samples = len(samples)
        #print(samples)

    # first consider files that are associated with a single case
    if num_associated_cases == 1:
        if num_associated_samples == 0:
            case_id = _add_to_knowncases(cases[0], known_cases, gdc_api_root)
            _add_file_attribute(case_id, known_cases[case_id], file_uuid, filename,
                                data_category, data_type, data_format, experimental_strategy, workflow_type, access, program, gdc_api_root)
        elif num_associated_samples == 1:
            case_id = _add_to_knowncases(cases[0], known_cases, gdc_api_root)
            sample_id, _ = _add_to_knownsamples(samples[0], case_id, known_samples)
            _add_file_attribute(sample_id, known_samples[sample_id], file_uuid, filename,
                                data_category, data_type, data_format, experimental_strategy, 
                                workflow_type, access, program, gdc_api_root)

        elif num_associated_samples == 2:
            if data_category == GDC_DataCategory.SNV or data_category == GDC_DataCategory.COMBINED_NUCLEOTIDE_VARIATION:
                case_id = _add_to_knowncases(cases[0], known_cases, gdc_api_root)
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
                    _add_file_attribute(pair_id, known_pairs[pair_id], file_uuid, filename,
                                        data_category, data_type, data_format, experimental_strategy, 
                                        workflow_type, access, program, gdc_api_root)
            else:
                    # Within some programs (e.g., CPTAC), if a given specimen did not have enough material for genomics analysis
                    # multiple specimens or cores from a given patient are being combined to get 
                    # sufficient material. 
                    case_id = _add_to_knowncases(cases[0], known_cases, gdc_api_root)
                    sample_id, _pwd = _add_pooled_sample_to_knownsamples(file_uuid, filename, samples, case_id, known_samples)
                    _add_file_attribute(sample_id, known_samples[sample_id], file_uuid, filename,
                                        data_category, data_type, data_format, experimental_strategy, 
                                        workflow_type, access, program, gdc_api_root)

        elif num_associated_samples > 2: 
            if program == GDC_ProgramName.CPTAC and (data_category == GDC_DataCategory.SNV or data_category == GDC_DataCategory.COMBINED_NUCLEOTIDE_VARIATION):
                # Within the CPTAC program, if a given specimen did not have enough material for proteomics or 
                # genomics, multiple specimens or cores from a given patient are being combined to get 
                # sufficient material. 
                case_id = _add_to_knowncases(cases[0], known_cases, gdc_api_root)
                tumor_samples = []
                normal_samples = []
                na_samples = []
                for sample in samples:
                    sample_type = _get_sample_type(sample)
                    if sample_type == SAMPLE_TYPE.TUMOR:
                        tumor_samples.append(sample)
                    elif sample_type == SAMPLE_TYPE.NORMAL:
                        normal_samples.append(sample)
                    else:
                        na_samples.append(sample)
                if len(tumor_samples) >= 1:
                    if len(tumor_samples) == 1:
                        tumor_sample_id, _pwd = _add_to_knownsamples(tumor_samples[0], case_id, known_samples)
                    else:
                        tumor_sample_id, _pwd = _add_pooled_sample_to_knownsamples(file_uuid, filename, tumor_samples, case_id, known_samples)
                else:
                    print('issue with tumor samples, number of tumor samples = 0')

                if len(normal_samples) >= 1:
                    if len(normal_samples) == 1:
                        normal_sample_id, _pwd = _add_to_knownsamples(normal_samples[0], case_id, known_samples)
                    else:
                        normal_sample_id, _pwd = _add_pooled_sample_to_knownsamples(file_uuid, filename, normal_samples, case_id, known_samples)
                else:
                    print('issue with normal samples, number of normal samples = 0')
                    raise ValueError(file_uuid, filename, data_category)

                pair_id = _add_to_knownpairs(tumor_sample_id, normal_sample_id, known_pairs)
                _add_file_attribute(pair_id, known_pairs[pair_id], file_uuid, filename,
                                    data_category, data_type, data_format, experimental_strategy, 
                                    workflow_type, access, program, gdc_api_root)

            else:
                case_id = _add_to_knowncases(cases[0], known_cases, gdc_api_root)
                sample_id, _pwd = _add_pooled_sample_to_knownsamples(file_uuid, filename, samples, case_id, known_samples)
                _add_file_attribute(sample_id, known_samples[sample_id], file_uuid, filename,
                                    data_category, data_type, data_format, experimental_strategy, 
                                    workflow_type, access, program, gdc_api_root)

    else:
        # file associated with multiple cases
        # we will record file_uuid and deal with later
        DEFERRED_FILE_NUM_OF_CASES[file_uuid] = num_associated_cases
        deferred_file_uuids.append([file_uuid, filename])

# may eventually drop this and incorporate into get_file_metadata.  Wasn't sure what to do with files
# associated with multiple cases or files associated with samples across multiple cases.
# For now we only include files in data model if they are associated with cases and samples that were
# pulled in from other files in manifest that were associated with single cases; this behavior, however, 
# can be overridden by setting all_cases to true, in which case a paricipant entity will be created for each
# case a file is associated with.

def process_deferred_file_uuid(file_uuid, filename, known_cases, known_samples, all_cases, gdc_api_root):
    
    # get data file's name, category, type, access, format experimental strategy, workflow type
    fileMetadataRetriever = FileMetadataRetriever(gdc_api_root)
    responseDict = fileMetadataRetriever.get_metadata(file_uuid)

    data_category = responseDict['data_category']
    data_type = responseDict['data_type']
    data_format = responseDict['data_format']
    access = responseDict['access']
    program = responseDict['cases'][0]['project']['program']['name']

    # I have decided to ignore (i.e., not incorporate into workspace) Clinical and Biospecimen files of data 
    # format "BCR Biotab"; these files are typically associated with multple cases, and there can be
    # multiple files that would map to the same attribute where all of the files are relevant; i.e., one doesn't 
    # replace another.  This doesn't fit into our data model.
    if data_category == GDC_DataCategory.BIOSPECIMEN and data_type == GDC_DataType.BIOSPECIMEN_SUPPLEMENT:
        if data_format in ['BCR Biotab']:
            print('skipping {0} file {1}'.format(data_format, file_uuid))
            return
                  
    if data_category == GDC_DataCategory.CLINICAL and data_type == GDC_DataType.CLINICAL_SUPPLEMENT:
        if data_format in ['BCR Biotab']:
            print('skipping {0} file {1}'.format(data_format, file_uuid))
            return

    if 'experimental_strategy' in responseDict:
        experimental_strategy = responseDict['experimental_strategy']
    else: 
        experimental_strategy = None
    if 'analysis' in responseDict and 'workflow_type' in responseDict['analysis']:
        workflow_type = responseDict['analysis']['workflow_type']
    else:
        workflow_type = None
        
    if data_category == GDC_DataCategory.CLINICAL or data_category == GDC_DataCategory.BIOSPECIMEN:
        fileMetadataRetriever = FileCaseMetadataRetriever(gdc_api_root)
    else:
        fileMetadataRetriever = FileCaseSampleMetadataRetriever(gdc_api_root)

    filedMmetadata = fileMetadataRetriever.get_metadata(file_uuid)

    cases = fileMetadata['cases']
    num_associated_cases = len(cases)
    assert num_associated_cases > 1, file_uuid

    for case in cases:
        case_id = case['case_id']
        if all_cases or case_id in known_cases:
            case_id = _add_to_knowncases(case, known_cases, gdc_api_root)
            if 'samples' in case:
                # associated with samples
                samples = case['samples']
                for sample in samples:
                    sample_id = sample['sample_id']
                    if sample_id in known_samples:
                        _add_file_attribute(sample_id, known_samples[sample_id], file_uuid, filename,
                                            data_category, data_type, data_format, experimental_strategy, 
                                            workflow_type, access, program, gdc_api_root)
            else:
                # associated with multiple cases only
                _add_file_attribute(case_id, known_cases[case_id], file_uuid, filename,
                                    data_category, data_type, data_format,experimental_strategy, 
                                    workflow_type, access, program, gdc_api_root)


def create_participants_file(cases, manifestFileBasename):
    
    attribute_names = []
    for case_id, case in cases.items():
        for attribute_name in case:
            if attribute_name not in attribute_names:
                attribute_names.append(attribute_name)
    
    participants_filename = manifestFileBasename + '_participants.tsv'
    participant_sets_membership_filename = manifestFileBasename + '_participant_sets_membership.tsv'
    
    with open(participants_filename, 'w') as participantsFile, open(participant_sets_membership_filename, 'w') as membershipFile:
        fieldnames = ['entity:participant_id'] + attribute_names
        participants_writer = csv.DictWriter(participantsFile, fieldnames=fieldnames, delimiter='\t')
        participants_writer.writeheader()

        fieldnames = ['membership:participant_set_id', 'participant']
        membership_writer = csv.DictWriter(membershipFile, fieldnames=fieldnames, delimiter='\t')
        membership_writer.writeheader()
        
        for case_id, case in cases.items():
            entity_row = {'entity:participant_id': case_id}
            for attribute_name in attribute_names:
                if attribute_name in case:
                    entity_row[attribute_name] = case[attribute_name]
                else:
                    entity_row[attribute_name] = '__DELETE__'
                entity_row[attribute_name] = case[attribute_name] if attribute_name in case else '__DELETE__'
            participants_writer.writerow(entity_row)

            membership_row = {'membership:participant_set_id' : 'ALL',
                              'participant' : case_id}
            membership_writer.writerow(membership_row)            

def create_samples_file(samples, manifestFileBasename):
    attribute_names = []
    for sample_id, sample in samples.items():
        for attribute_name in sample:
            if attribute_name not in {'submitter_id', 'case_id', 'sample_type_id'} and attribute_name not in attribute_names:
                attribute_names.append(attribute_name)

    samples_filename = manifestFileBasename + '_samples.tsv'
    sample_sets_membership_filename = manifestFileBasename + '_sample_sets_membership.tsv'
    with open(samples_filename, 'w') as samplesFile, open(sample_sets_membership_filename, 'w') as membershipFile:
        
        fieldnames = ['entity:sample_id', 'participant', 'submitter_id', 'sample_type_code', 'sample_type', 'tissue_type' ] + attribute_names
        sample_writer = csv.DictWriter(samplesFile, fieldnames=fieldnames, delimiter='\t')
        sample_writer.writeheader()

        fieldnames = ['membership:sample_set_id', 'sample']
        membership_writer = csv.DictWriter(membershipFile, fieldnames=fieldnames, delimiter='\t')
        membership_writer.writeheader()
        
        for sample_id, sample in samples.items():
            entity_row = {'entity:sample_id' : sample_id, 'participant': sample['case_id'],
                          'submitter_id' : sample['submitter_id'],
                          'sample_type_code' : SAMPLE_TYPE.getLetterCode(sample['sample_type_id']) if sample['sample_type_id'] is not None else '__DELETE__',
                          'sample_type' : sample['sample_type'] if sample['sample_type'] is not None else '__DELETE__',
                          'tissue_type' : sample['tissue_type'] if sample['tissue_type'] is not None else '__DELETE__'}
            for attribute_name in attribute_names:
                if attribute_name in sample:
                    entity_row[attribute_name] = sample[attribute_name]
                else:
                    entity_row[attribute_name] = '__DELETE__'
            sample_writer.writerow(entity_row)

            membership_row = {'membership:sample_set_id' : 'ALL',
                              'sample': sample_id}
            membership_writer.writerow(membership_row)
                        
def create_pairs_file(pairs, samples, manifestFileBasename):
    attribute_names = []
    for pair_id, pair in pairs.items():
        for attribute_name in pair:
            if attribute_name not in {'tumor', 'normal'} and attribute_name not in attribute_names:
                attribute_names.append(attribute_name)

    pairs_filename = manifestFileBasename + '_pairs.tsv'
    pair_sets_membership_filename = manifestFileBasename + '_pair_sets_membership.tsv'
    with open(pairs_filename, 'w') as pairsFile, open(pair_sets_membership_filename, 'w') as membershipFile:
        fieldnames = ['entity:pair_id', 'participant', 'case_sample', 'control_sample',
                    'tumor_submitter_id', 'normal_submitter_id',
                    'tumor_type', 'normal_type'] + attribute_names
        pairs_writer = csv.DictWriter(pairsFile, fieldnames=fieldnames, delimiter='\t')
        pairs_writer.writeheader()

        fieldnames = ['membership:pair_set_id', 'pair']
        membership_writer = csv.DictWriter(membershipFile, fieldnames=fieldnames, delimiter='\t')
        membership_writer.writeheader()
        
        for pair_id, pair in pairs.items():

            tumor_submitter_id = samples[pair['tumor']]['submitter_id']
            normal_submitter_id = samples[pair['normal']]['submitter_id']
            entity_row = {'entity:pair_id' : pair_id,
                          'participant' : samples[pair['tumor']]['case_id'],
                          'case_sample' : pair['tumor'],
                          'control_sample' : pair['normal'],
                          'tumor_submitter_id' : tumor_submitter_id,
                          'normal_submitter_id' : normal_submitter_id,
                          'tumor_type' : SAMPLE_TYPE.getLetterCode(samples[pair['tumor']]['sample_type_id']),
                          'normal_type' : SAMPLE_TYPE.getLetterCode(samples[pair['normal']]['sample_type_id'])}
            for attribute_name in attribute_names:
                if attribute_name in pair:
                    entity_row[attribute_name] = pair[attribute_name]
                else:
                    entity_row[attribute_name] = '__DELETE__'
            pairs_writer.writerow(entity_row)

            row = {'membership:pair_set_id' : 'ALL',
                   'pair': pair_id}
            membership_writer.writerow(row)


def create_workspace_attributes_file(manifestFileBasename, is_legacy):
    #This part is hardcoded due to the small number of attributes we need to specify.
    #Please feel free to change this specification according to your needs.
    legacy_flag="false"
    if is_legacy:
        legacy_flag="true"

    #Due to a somewhat weird bug in FireCloud, please keep the workspace-colunm-defaults attribute as the last one in the list.
    #Any new attributes should be added before workspace-column-defaults
    with open(manifestFileBasename + "_workspace_attributes.tsv", 'w') as workspaceColumnOrderFile:
        workspaceColumnOrderFile.write("workspace:legacy_flag\tworkspace-column-defaults\n")
        workspaceColumnOrderFile.write(legacy_flag + "\t" + "{\"participant\": {\"shown\": [\"submitter_id\", \"project_id\", \"participant_id\"]}, \"sample\":{\"shown\":[\"submitter_id\", \"sample_id\", \"participant\", \"sample_type\"]}, \"pair\":{\"shown\":[\"tumor_submitter_id\", \"normal_submitter_id\", \"pair_id\"]}}")


def main():
    parser = argparse.ArgumentParser(description='create FireCloud workspace load files from GDC manifest')
    parser.add_argument("manifest", help="manifest file from the GDC Data Portal")
    parser.add_argument("-l", "--legacy", help="point to GDC Legacy Archive", action="store_true")
    parser.add_argument("-c", "--all_cases", help="create participant entities for all referenced cases", action="store_true")
    args = parser.parse_args()

    gdc_api_root = GDC_LEGACY_API_ROOT if args.legacy else GDC_API_ROOT

    print("manifestFile = {0}".format(args.manifest))

    manifestFile = args.manifest

    pp = pprint.PrettyPrinter()

    cases = dict()
    samples = dict()
    pairs = dict()
    deferred_file_uuids = []

    manifestFileList = _read_manifestFile(manifestFile)

    for i, item in enumerate(manifestFileList):

        file_uuid = item['id']
        filename = item['filename']
    
        UUID_TO_FILENAME[file_uuid] = filename

        print('{0} of {1}: {2}, {3}'.format(i+1, len(manifestFileList), file_uuid, filename))

        for attempt in range(5):
            try:
                get_file_metadata(file_uuid, filename, cases, samples, 
                                  pairs, deferred_file_uuids, gdc_api_root)
                break
            except (KeyboardInterrupt, SystemExit):
                raise
            except ValueError as x:
                print('Value Error, skip: {0}'.format(x))
                break
            except Exception as x:
                print(''.join(traceback.format_exception(etype=type(x), value=x, tb=x.__traceback__)))
                print("attempt=", attempt, 'file uuid = ', file_uuid)
                time.sleep((attempt+1)**2)
        else:
            #failed all attempts
            # - just move on
            print("failed 5 attempts! SKIPPING FILE: file uuid = ", file_uuid)


    print("Processing deferred files...")
    for uuid_and_filename in deferred_file_uuids:
        file_uuid = uuid_and_filename[0]
        filename = uuid_and_filename[1]
        print("{0}, {1} ".format(file_uuid, filename))

        for attempt in range(5):
            try:
                process_deferred_file_uuid(file_uuid, filename, cases, samples, args.all_cases, gdc_api_root)
            except (KeyboardInterrupt, SystemExit):
                raise
            except Exception as x:
                print(''.join(traceback.format_exception(etype=type(x), value=x, tb=x.__traceback__)))
                print("attempt=", attempt, 'file uuid = ', file_uuid)
                time.sleep((attempt+1)**2)
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
    create_workspace_attributes_file(manifestFileBasename, False)
    

if __name__ == '__main__':
    main()
