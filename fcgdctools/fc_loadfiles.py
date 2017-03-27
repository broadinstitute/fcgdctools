# Version: 6
# Date: March 26, 2017
# Comments: combine uuid and filename into single attributeincorporated uuid-to-url resolution

import csv
import requests
import argparse
import pprint
import os.path
import sys

import gdc_uuidresolver 

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
    
#data types
class GDC_DataType:
    #associated with Clinical data category
    CLINICAL_SUPPLEMENT = "Clinical Supplement"

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

class FC_DataModelAttributes:

    def __init__(self):
        
        self.attributeNames = {
            GDC_DataCategory.CLINICAL : {GDC_DataType.CLINICAL_SUPPLEMENT : 'clinical_supplement'},
            GDC_DataCategory.BIOSPECIMEN : {GDC_DataType.BIOSPECIMEN_SUPPLEMENT : 'biospecimen_supplement'},
            GDC_DataCategory.DNA_METHYLATION : {GDC_DataType.METHYLATION_BETA_VALUE : 'methylation_beta_value'},
            GDC_DataCategory.COPY_NUMBER_VARIATION : {GDC_DataType.COPY_NUMBER_SEGMENT : 'copy_number_segment',
                                                     GDC_DataType.MASKED_COPY_NUMBER_SEGMENT : 'masked_copy_number_segment'},
            GDC_DataCategory.RAW_SEQUENCING_DATA : {GDC_DataType.ALIGNED_READS : 'aligned_reads'},
            GDC_DataCategory.TRANSCRIPTOME_PROFILING : {GDC_DataType.MIRNA_EXPRESSION_QUANTIFICATION : 'mirna_expression_quantification',
                                                       GDC_DataType.ISOFORM_EXPRESSION_QUANTIFICATION : 'isoform_expression_quantification',
                                                       GDC_DataType.GENE_EXPRESSION_QUANTIFICATION : 'gene_expression_quantification'},
            GDC_DataCategory.SNV : {GDC_DataType.RAW_SIMPLE_SOMATIC_MUTATION : 'raw_simple_somatic_mutation',
                                   GDC_DataType.ANNOTATED_SOMATIC_MUTATION : 'annotated_somatic_mutation',
                                   GDC_DataType.AGGREGATED_SOMATIC_MUTATION : 'aggregated_somatic_mutation',
                                   GDC_DataType.MASKED_SOMATIC_MUTATION : 'masked_somatic_mutation'}
            }

    def __category_and_type_recognized(self, data_category, data_type):
        return data_category in self.attributeNames and data_type in self.attributeNames[data_category]


    def constructAttributeName_base(self, experimental_strategy, workflow_type, data_category, data_type):

        assert self.__category_and_type_recognized(data_category, data_type), "category={0}, type={1}".format(data_category, data_type)

        if experimental_strategy is not None:
            if experimental_strategy in EXP_STRATEGY.abbreviations:            
                experimental_strategy_abbrev = EXP_STRATEGY.getAbbreviation(experimental_strategy) + '__'
            else:
                experimental_strategy_abbrev = experimental_strategy.translate(EXP_STRATEGY.ABBREV_TRANSLATE_TABLE) + '__'
        else:
            experimental_strategy_abbrev = ''
    
        if workflow_type is not None:
            if workflow_type in WORKFLOW.abbreviations:
                workflow_type_abbrev = WORKFLOW.getAbbreviation(workflow_type) + '__'
            else:
                workflow_type_abbrev = workflow_type.translate(WORKFLOW.ABBREV_TRANSLATE_TABLE) + '__'
        else:
            workflow_type_abbrev = ''

        attribute_name_base = experimental_strategy_abbrev + workflow_type_abbrev + self.attributeNames[data_category][data_type]

        return (attribute_name_base)

    
GDC_TO_FC = FC_DataModelAttributes()

class GDC_FileAccess:

    def __init__(self):

        self.accessTypePrefix = {
            GDC_FileAccessType.OPEN : 'OA__',
            GDC_FileAccessType.CONTROLLED : 'CA__'
            }

        self.gdcFileAccess = dict()
        

    def recordFileAccessType(self, file_attribute_base_name, access_type):
         
        assert access_type in [GDC_FileAccessType.OPEN, GDC_FileAccessType.CONTROLLED]
        
        if file_attribute_base_name in self.gdcFileAccess:
            assert self.gdcFileAccess[file_attribute_base_name] == access_type
        else:
            self.gdcFileAccess[file_attribute_base_name] = access_type
        
    def getAccessTypePrefix(self, file_attribute_base_name):
        assert file_attribute_base_name in self.gdcFileAccess, file_attribute_base_name
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

class DataSource:
    
    def __init__(self, abbreviations):
        self.abbreviations = abbreviations

    def getAbbreviation(self, type):
        return self.abbreviations[type]
        
class ExperimentalStrategy(DataSource):

    def __init__(self, abbreviations):
        DataSource.__init__(self, abbreviations)

EXP_STRATEGY_ABBREVIATIONS = {
        'WXS' : 'WXS',
        'RNA-Seq' : 'RNAseq',
        'miRNA-Seq' : 'miRNAseq',
        'Genotyping Array' : 'GeneArray',
        'Methylation Array' : 'MethArray'}

EXP_STRATEGY = ExperimentalStrategy(EXP_STRATEGY_ABBREVIATIONS)
    
class Workflow(DataSource):

    ABBREV_TRANSLATE_TABLE = ''.maketrans({'.' : '', '-' : '', ' ' : '', '_' :''})

    def __init__(self, abbreviations):
        DataSource.__init__(self, abbreviations)

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
        'MuTect2 Variant Aggregation and Masking' : 'MuSEAggrMask',
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

WORKFLOW = Workflow(WORKFLOW_ABBREVIATIONS)
        
class Platform(DataSource):

    def __init__(self, abbreviations):
        DataSource.__init__(self, abbreviations)
        

PLATFORM_ABBREVIATIONS = {
        'Affymetrix SNP 6.0' : 'AffySNP6',
        'Illumina' : 'Illum',
        'Illumina Human Methylation 450' : 'IllumHuMeth450',
        'Illumina Human Methylation 27' : 'IllumHuMeth27'}

PLATFORM = Platform(PLATFORM_ABBREVIATIONS)

GDC_API_ROOT = "https://gdc-api.nci.nih.gov"
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


def _get_case_level_metadata(file_uuid):
    #get case_id and (case) submitter_id
    fields = "cases.case_id,cases.submitter_id"
    url = "{0}/files/{1}?fields={2}".format(GDC_API_ROOT, file_uuid, fields)
    #print("get_case_level_metadata:", url)

    response = requests.get(url, headers=None)
    #print(response)
    responseDict = response.json()
    #print("responseDict:")
    #pp = pprint.PrettyPrinter()
    #pp.pprint(responseDict)

    return responseDict['data']
    
def _get_case_and_sample_level_metadata(file_uuid):
    
    fields = "cases.case_id,cases.submitter_id"
    fields = fields + ",cases.samples.sample_id,cases.samples.submitter_id,cases.samples.sample_type_id"
    url = "{0}/files/{1}?fields={2}".format(GDC_API_ROOT, file_uuid, fields)
    #print("get_case_and sample_level_metadata:", url)

    response = requests.get(url, headers=None)
    responseDict = response.json()

    #pp = pprint.PrettyPrinter()
    #print("responseDict:")
    #pp.pprint(responseDict)
    return responseDict['data']


def _get_clinical_biospecimen_supplement_metadata(file_uuid):
    
    file_metadata = _get_case_level_metadata(file_uuid)
    return file_metadata

def _get_methylation_beta_value_metadata(file_uuid):
    
    file_metadata = _get_case_and_sample_level_metadata(file_uuid)
    return file_metadata

def _get_copy_number_segment_metadata(file_uuid):
    
    file_metadata = _get_case_and_sample_level_metadata(file_uuid)
    return file_metadata

def _get_aligned_reads_metadata(file_uuid):
    
    file_metadata = _get_case_and_sample_level_metadata(file_uuid)
    return file_metadata

def _get_expression_quantification_metadata(file_uuid):
    
    file_metadata = _get_case_and_sample_level_metadata(file_uuid)
    return file_metadata

def _get_paired_SNV_metadata(file_uuid):
    
    file_metadata = _get_case_and_sample_level_metadata(file_uuid)
    return file_metadata

def _get_aggregated_SNV_metadata(file_uuid):
    file_metadata = _get_case_and_sample_level_metadata(file_uuid)
    return file_metadata


def get_file_metadata(file_uuid, filename, file_url, known_cases, known_samples, known_pairs, deferred_file_uuids):
    
    # get from GDC the data file's category, type, access type, format, experimental strategy,
    # analysis workflow type
    fields = "data_category,data_type,access, data_format,experimental_strategy,analysis.workflow_type"
    url = "{0}/files/{1}?fields={2}".format(GDC_API_ROOT, file_uuid, fields)
    response = requests.get(url, headers=None)
    responseDict = response.json()
    #pp.pprint(responseDict)
    
    data_category = responseDict['data']['data_category']
    data_type = responseDict['data']['data_type']
    access = responseDict['data']['access']

    if 'data_format' in responseDict['data']:
        data_format = responseDict['data']['data_format']
    else:
        data_format = None
    if 'experimental_strategy' in responseDict['data']:
        experimental_strategy = responseDict['data']['experimental_strategy']
    else: 
        experimental_strategy = None
    if 'analysis' in responseDict['data'] and 'workflow_type' in responseDict['data']['analysis']:
        workflow_type = responseDict['data']['analysis']['workflow_type']
    else:
        workflow_type = None

    # further analysis based on data category and data type

    if data_category == GDC_DataCategory.CLINICAL:
        if data_type == GDC_DataType.CLINICAL_SUPPLEMENT:
            metadata = _get_clinical_biospecimen_supplement_metadata(file_uuid)
            cases = metadata['cases']
            assert len(cases) > 0, file_uuid

            if len(cases) == 1:
                case_id = cases[0]['case_id']
                submitter_id = cases[0]['submitter_id']
                if case_id not in known_cases:
                    new_case = {'submitter_id' : submitter_id}
                    known_cases[case_id] = new_case
                
                basename = GDC_TO_FC.constructAttributeName_base(experimental_strategy,
                                                                 workflow_type,
                                                                 data_category,
                                                                 data_type)
                GDC_FILE_ACCESS.recordFileAccessType(basename, access)
                known_cases[case_id][basename + UUID_ATTRIBUTE_SUFFIX] = file_uuid + SEPARATOR + filename
                known_cases[case_id][basename + URL_ATTRIBUTE_SUFFIX] = file_url

            elif len(cases) > 1:
                # we have a clinical supplement file that is attached to multiple cases.
                # we will just record the file_uuid and deal with it after having processed all of the
                # files in the manifest
                deferred_file_uuids.append(file_uuid)
        else:
            raise ValueError(file_uuid, filename, data_category, data_type)
    
    elif data_category == GDC_DataCategory.BIOSPECIMEN:
        if data_type == GDC_DataType.BIOSPECIMEN_SUPPLEMENT:
            metadata = _get_clinical_biospecimen_supplement_metadata(file_uuid)
            cases = metadata['cases']
            assert len(cases) > 0, file_uuid
            if len(cases) == 1:
                case_id = cases[0]['case_id']
                submitter_id = cases[0]['submitter_id']
                if case_id not in known_cases:
                    new_case = {'submitter_id' : submitter_id}
                    known_cases[case_id] = new_case
                basename = GDC_TO_FC.constructAttributeName_base(experimental_strategy,
                                                                 workflow_type,
                                                                 data_category,
                                                                 data_type)
                GDC_FILE_ACCESS.recordFileAccessType(basename, access)
                known_cases[case_id][basename + UUID_ATTRIBUTE_SUFFIX] = file_uuid + SEPARATOR + filename
                known_cases[case_id][basename + URL_ATTRIBUTE_SUFFIX] = file_url

            elif len(cases) > 1:
                # we have a biospecimen supplement file that is attached to multiple cases.
                # we will just record the file_uuid and deal with it after having processed all of the
                # files in the manifest
                deferred_file_uuids.append(file_uuid)
        else:
            raise ValueError(file_uuid, filename, data_category, data_type)

    # no methylation data for TARGET Program        
    elif data_category == GDC_DataCategory.DNA_METHYLATION:

        if data_type == GDC_DataType.METHYLATION_BETA_VALUE:
            metadata = _get_methylation_beta_value_metadata(file_uuid)        
            
            cases = metadata['cases']
            assert len(cases) == 1, file_uuid
            case_id = cases[0]['case_id']
            case_submitter_id = cases[0]['submitter_id']
            samples = metadata['cases'][0]['samples']
            assert len(samples) == 1, file_uuid
            sample_id = samples[0]['sample_id']
            sample_submitter_id = samples[0]['submitter_id']
            sample_type_id = samples[0]['sample_type_id']

            if case_id not in known_cases:
                new_case ={'submitter_id' : case_submitter_id}
                known_cases[case_id] = new_case   
            if sample_id not in known_samples:
                new_sample = {'submitter_id' : sample_submitter_id, 
                              'sample_type_id' : sample_type_id,
                              'case_id' : case_id}
                known_samples[sample_id] = new_sample

            basename = GDC_TO_FC.constructAttributeName_base(experimental_strategy,
                                                             workflow_type,
                                                             data_category,
                                                             data_type)
            GDC_FILE_ACCESS.recordFileAccessType(basename, access)
            known_samples[sample_id][basename + UUID_ATTRIBUTE_SUFFIX] = file_uuid + SEPARATOR + filename
            known_samples[sample_id][basename + URL_ATTRIBUTE_SUFFIX] = file_url

        else:
            raise ValueError(file_uuid, filename, data_category, data_type)
    
    elif data_category == GDC_DataCategory.COPY_NUMBER_VARIATION:

        if data_type == GDC_DataType.COPY_NUMBER_SEGMENT:
            metadata = _get_copy_number_segment_metadata(file_uuid)
            #pp.pprint(metadata)
            
            cases = metadata['cases']
            assert len(cases) == 1, file_uuid
            case_id = cases[0]['case_id']
            case_submitter_id = cases[0]['submitter_id']
            samples = metadata['cases'][0]['samples']
            assert len(samples) == 1, file_uuid
            sample_id = samples[0]['sample_id']
            sample_submitter_id = samples[0]['submitter_id']
            sample_type_id = samples[0]['sample_type_id']            
            
            if case_id not in known_cases:
                new_case ={'submitter_id' : case_submitter_id}
                known_cases[case_id] = new_case   
            if sample_id not in known_samples:
                new_sample = {'submitter_id' : sample_submitter_id, 
                              'sample_type_id' : sample_type_id,
                              'case_id' : case_id}
                known_samples[sample_id] = new_sample

            basename = GDC_TO_FC.constructAttributeName_base(experimental_strategy,
                                                             workflow_type,
                                                             data_category,
                                                             data_type)
            GDC_FILE_ACCESS.recordFileAccessType(basename, access)

            known_samples[sample_id][basename + UUID_ATTRIBUTE_SUFFIX] = file_uuid + SEPARATOR + filename
            known_samples[sample_id][basename + URL_ATTRIBUTE_SUFFIX] = file_url

        elif data_type == GDC_DataType.MASKED_COPY_NUMBER_SEGMENT:
            metadata = _get_copy_number_segment_metadata(file_uuid)
            #pp.pprint(metadata)
            
            cases = metadata['cases']
            assert len(cases) == 1, file_uuid
            case_id = cases[0]['case_id']
            case_submitter_id = cases[0]['submitter_id']
            samples = metadata['cases'][0]['samples']
            assert len(samples) == 1, file_uuid
            sample_id = samples[0]['sample_id']
            sample_submitter_id = samples[0]['submitter_id']
            sample_type_id = samples[0]['sample_type_id']            
            
            if case_id not in known_cases:
                new_case ={'submitter_id' : case_submitter_id}
                known_cases[case_id] = new_case   
            if sample_id not in known_samples:
                new_sample = {'submitter_id' : sample_submitter_id, 
                              'sample_type_id' : sample_type_id,
                              'case_id' : case_id}
                known_samples[sample_id] = new_sample

            basename = GDC_TO_FC.constructAttributeName_base(experimental_strategy,
                                                             workflow_type,
                                                             data_category,
                                                             data_type)
            GDC_FILE_ACCESS.recordFileAccessType(basename, access)

            known_samples[sample_id][basename + UUID_ATTRIBUTE_SUFFIX] = file_uuid + SEPARATOR + filename
            known_samples[sample_id][basename + URL_ATTRIBUTE_SUFFIX] = file_url

        else:
            raise ValueError(file_uuid, filename, data_category, data_type)
    
    elif data_category == GDC_DataCategory.RAW_SEQUENCING_DATA:

        if data_type == GDC_DataType.ALIGNED_READS:
            metadata = _get_aligned_reads_metadata(file_uuid)
            
            cases = metadata['cases']
            assert len(cases) == 1, file_uuid
            case_id = cases[0]['case_id']
            case_submitter_id = cases[0]['submitter_id']
            samples = metadata['cases'][0]['samples']
            assert len(samples) == 1, file_uuid
            sample_id = samples[0]['sample_id']
            sample_submitter_id = samples[0]['submitter_id']
            sample_type_id = samples[0]['sample_type_id']            
            
            if case_id not in known_cases:
                new_case ={'submitter_id' : case_submitter_id}
                known_cases[case_id] = new_case   
            if sample_id not in known_samples:
                new_sample = {'submitter_id' : sample_submitter_id, 
                              'sample_type_id' : sample_type_id,
                              'case_id' : case_id}
                known_samples[sample_id] = new_sample            
            
            basename = GDC_TO_FC.constructAttributeName_base(experimental_strategy,
                                                             workflow_type,
                                                             data_category,
                                                             data_type)
            GDC_FILE_ACCESS.recordFileAccessType(basename, access)
            
            known_samples[sample_id][basename + UUID_ATTRIBUTE_SUFFIX] = file_uuid + SEPARATOR + filename
            known_samples[sample_id][basename + URL_ATTRIBUTE_SUFFIX] = file_url

        else:
            raise ValueError(file_uuid, filename, data_category, data_type)

    elif data_category == GDC_DataCategory.TRANSCRIPTOME_PROFILING:

        if data_type == GDC_DataType.MIRNA_EXPRESSION_QUANTIFICATION:
            metadata = _get_expression_quantification_metadata(file_uuid)
            
            cases = metadata['cases']
            assert len(cases) == 1, file_uuid
            case_id = cases[0]['case_id']
            case_submitter_id = cases[0]['submitter_id']
            samples = metadata['cases'][0]['samples']
            assert len(samples) == 1, file_uuid
            sample_id = samples[0]['sample_id']
            sample_submitter_id = samples[0]['submitter_id']
            sample_type_id = samples[0]['sample_type_id']                      
            
            
            if case_id not in known_cases:
                known_cases[case_id] = {"submitter_id" : case_submitter_id}
            if sample_id not in known_samples:
                known_samples[sample_id] = {"submitter_id" : sample_submitter_id,
                                            "sample_type_id" : sample_type_id,
                                            "case_id" : case_id}

            basename = GDC_TO_FC.constructAttributeName_base(experimental_strategy,
                                                             workflow_type,
                                                             data_category,
                                                             data_type)
            GDC_FILE_ACCESS.recordFileAccessType(basename, access)
            
            known_samples[sample_id][basename + UUID_ATTRIBUTE_SUFFIX] = file_uuid + SEPARATOR + filename
            known_samples[sample_id][basename + URL_ATTRIBUTE_SUFFIX] = file_url

        elif data_type == GDC_DataType.ISOFORM_EXPRESSION_QUANTIFICATION:
            metadata = _get_expression_quantification_metadata(file_uuid)
            
            cases = metadata['cases']
            assert len(cases) == 1, file_uuid
            case_id = cases[0]['case_id']
            case_submitter_id = cases[0]['submitter_id']
            samples = metadata['cases'][0]['samples']
            assert len(samples) == 1, file_uuid
            sample_id = samples[0]['sample_id']
            sample_submitter_id = samples[0]['submitter_id']
            sample_type_id = samples[0]['sample_type_id']                      
            
            
            if case_id not in known_cases:
                known_cases[case_id] = {"submitter_id" : case_submitter_id}
            if sample_id not in known_samples:
                known_samples[sample_id] = {"submitter_id" : sample_submitter_id,
                                            "sample_type_id" : sample_type_id,
                                            "case_id" : case_id}

            basename = GDC_TO_FC.constructAttributeName_base(experimental_strategy,
                                                             workflow_type,
                                                             data_category,
                                                             data_type)
            GDC_FILE_ACCESS.recordFileAccessType(basename, access)
            
            known_samples[sample_id][basename + UUID_ATTRIBUTE_SUFFIX] = file_uuid + SEPARATOR + filename
            known_samples[sample_id][basename + URL_ATTRIBUTE_SUFFIX] = file_url

        elif data_type == GDC_DataType.GENE_EXPRESSION_QUANTIFICATION:
            metadata = _get_expression_quantification_metadata(file_uuid)
            
            cases = metadata['cases']
            assert len(cases) == 1, file_uuid
            case_id = cases[0]['case_id']
            case_submitter_id = cases[0]['submitter_id']
            samples = metadata['cases'][0]['samples']
            assert len(samples) == 1, file_uuid
            sample_id = samples[0]['sample_id']
            sample_submitter_id = samples[0]['submitter_id']
            sample_type_id = samples[0]['sample_type_id']                      
            
            
            if case_id not in known_cases:
                known_cases[case_id] = {"submitter_id" : case_submitter_id}
            if sample_id not in known_samples:
                known_samples[sample_id] = {"submitter_id" : sample_submitter_id,
                                            "sample_type_id" : sample_type_id,
                                            "case_id" : case_id}

            basename = GDC_TO_FC.constructAttributeName_base(experimental_strategy,
                                                             workflow_type,
                                                             data_category,
                                                             data_type)
            GDC_FILE_ACCESS.recordFileAccessType(basename, access)
            
            known_samples[sample_id][basename + UUID_ATTRIBUTE_SUFFIX] = file_uuid + SEPARATOR + filename
            known_samples[sample_id][basename + URL_ATTRIBUTE_SUFFIX] = file_url

        else:
            raise ValueError(file_uuid, data_category, data_type)
        
    elif data_category == GDC_DataCategory.SNV:

        if data_type == GDC_DataType.RAW_SIMPLE_SOMATIC_MUTATION:
            metadata = _get_paired_SNV_metadata(file_uuid)
            
            cases = metadata['cases']
            assert len(cases) == 1, file_uuid
            case_id = cases[0]['case_id']
            case_submitter_id = cases[0]['submitter_id']
            samples = metadata['cases'][0]['samples']
            assert len(samples) == 2, file_uuid
            

            pair = dict()
            for sample in samples:
        
                sample_id = sample['sample_id']
                sample_submitter_id = sample['submitter_id']  
                sample_type_id = sample['sample_type_id']
                sample_type_tn = SAMPLE_TYPE.getTumorNormalClassification(sample_type_id)
                assert sample_type_tn not in pair
                pair[sample_type_tn] = {'uuid' : sample_id, 'submitter_id' : sample_submitter_id,
                                        'type_id': sample_type_id}
    
            assert SampleType.TUMOR in pair and SampleType.NORMAL in pair
            
            if case_id not in known_cases:
                known_cases[case_id] = {'submitter_id' : case_submitter_id}
            tumor_sample = pair[SampleType.TUMOR] 
            if tumor_sample['uuid'] not in known_samples:
                sample_uuid = tumor_sample['uuid']
                known_samples[sample_uuid]= {'submitter_id' : tumor_sample['submitter_id'],
                                             'sample_type_id' : tumor_sample['type_id'],
                                             'case_id' : case_id}
            normal_sample = pair[SampleType.NORMAL]
            if normal_sample['uuid'] not in known_samples:
                sample_uuid = normal_sample['uuid']
                known_samples[sample_uuid]= {'submitter_id' : normal_sample['submitter_id'],
                                             'sample_type_id' : normal_sample['type_id'],
                                             'case_id' : case_id}

            pair_id = "{0}_{1}".format(tumor_sample['uuid'], normal_sample['uuid'])
            if pair_id not in known_pairs:
                known_pairs[pair_id] = {'tumor': tumor_sample['uuid'], 'normal': normal_sample['uuid']}

            basename = GDC_TO_FC.constructAttributeName_base(experimental_strategy,
                                                             workflow_type,
                                                             data_category,
                                                             data_type)
            GDC_FILE_ACCESS.recordFileAccessType(basename, access)

            known_pairs[pair_id][basename + UUID_ATTRIBUTE_SUFFIX] = file_uuid + SEPARATOR + filename
            known_pairs[pair_id][basename + URL_ATTRIBUTE_SUFFIX] = file_url

        elif data_type == GDC_DataType.ANNOTATED_SOMATIC_MUTATION:
            metadata = _get_paired_SNV_metadata(file_uuid)
            
            cases = metadata['cases']
            assert len(cases) == 1, file_uuid
            case_id = cases[0]['case_id']
            case_submitter_id = cases[0]['submitter_id']
            samples = metadata['cases'][0]['samples']
            assert len(samples) == 2, file_uuid

            pair = dict()
            for sample in samples:
        
                sample_id = sample['sample_id']
                sample_submitter_id = sample['submitter_id']  
                sample_type_id = sample['sample_type_id']
                sample_type_tn = SAMPLE_TYPE.getTumorNormalClassification(sample_type_id)
                assert sample_type_tn not in pair
                pair[sample_type_tn] = {'uuid' : sample_id, 'submitter_id' : sample_submitter_id,
                                        'type_id': sample_type_id}
    
            assert SampleType.TUMOR in pair and SampleType.NORMAL in pair
            
            if case_id not in known_cases:
                known_cases[case_id] = {'submitter_id' : case_submitter_id}
            tumor_sample = pair[SampleType.TUMOR] 
            if tumor_sample['uuid'] not in known_samples:
                sample_uuid = tumor_sample['uuid']
                known_samples[sample_uuid]= {'submitter_id' : tumor_sample['submitter_id'],
                                             'sample_type_id' : tumor_sample['type_id'],
                                             'case_id' : case_id}
            normal_sample = pair[SampleType.NORMAL]
            if normal_sample['uuid'] not in known_samples:
                sample_uuid = normal_sample['uuid']
                known_samples[sample_uuid]= {'submitter_id' : normal_sample['submitter_id'],
                                             'sample_type_id' : normal_sample['type_id'],
                                             'case_id' : case_id}

            pair_id = "{0}_{1}".format(tumor_sample['uuid'], normal_sample['uuid'])
            if pair_id not in known_pairs:
                known_pairs[pair_id] = {'tumor': tumor_sample['uuid'], 'normal': normal_sample['uuid']}

            basename = GDC_TO_FC.constructAttributeName_base(experimental_strategy,
                                                             workflow_type,
                                                             data_category,
                                                             data_type)
            GDC_FILE_ACCESS.recordFileAccessType(basename, access)

            known_pairs[pair_id][basename + UUID_ATTRIBUTE_SUFFIX] = file_uuid + SEPARATOR + filename
            known_pairs[pair_id][basename + URL_ATTRIBUTE_SUFFIX] = file_url

        elif data_type == GDC_DataType.AGGREGATED_SOMATIC_MUTATION:
            metadata = _get_aggregated_SNV_metadata(file_uuid)
            cases = metadata['cases']
            assert len(cases) > 0, file_uuid
            
            # don't bother creating pair sets; will just attach to all 
            # associated tumor samples that non-aggregated files are linked to
 
            deferred_file_uuids.append(file_uuid)           
                
        elif data_type == GDC_DataType.MASKED_SOMATIC_MUTATION:
            metadata = _get_aggregated_SNV_metadata(file_uuid)
            cases = metadata['cases']
            assert len(cases) > 0, file_uuid

            # don't bother creating pair sets; will just attach to all 
            # associated tumor samples that non-aggregated files are linked to
        
            deferred_file_uuids.append(file_uuid)            
            

            
    else:
        raise ValueError(file_uuid, filename, data_category)
        
    return

def process_deferred_file_uuid(file_uuid, file_url, known_cases, known_samples):
    
    # get data file's name, category, type, access, format experimental strategy, workflow type

    fields = "file_name,data_category,data_type,access, data_format,experimental_strategy,analysis.workflow_type"
    url = "{0}/files/{1}?fields={2}".format(GDC_API_ROOT, file_uuid, fields)
    response = requests.get(url, headers=None)
    responseDict = response.json()
    #pp.pprint(responseDict)
    filename = responseDict['data']['file_name']
    data_category = responseDict['data']['data_category']
    data_type = responseDict['data']['data_type']
    access = responseDict['data']['access']
    if 'data_format' in responseDict['data']:
        data_format = responseDict['data']['data_format']
    else:
        data_format = None
    
    if 'experimental_strategy' in responseDict['data']:
        experimental_strategy = responseDict['data']['experimental_strategy']
    else: 
        experimental_strategy = None
    if 'analysis' in responseDict['data'] and 'workflow_type' in responseDict['data']['analysis']:
        workflow_type = responseDict['data']['analysis']['workflow_type']
    else:
        workflow_type = None
        


    if data_category == GDC_DataCategory.CLINICAL:
        if data_type == GDC_DataType.CLINICAL_SUPPLEMENT:
            metadata = _get_clinical_biospecimen_supplement_metadata(file_uuid)
            cases = metadata['cases']
            assert len(cases) > 1, file_uuid
            for case in cases:
                case_id = case['case_id']
                if case_id in known_cases:
                    basename = GDC_TO_FC.constructAttributeName_base(experimental_strategy,
                                                                     workflow_type,
                                                                     data_category,
                                                                     data_type)
                    GDC_FILE_ACCESS.recordFileAccessType(basename, access)

                    known_cases[case_id][basename + UUID_ATTRIBUTE_SUFFIX] = file_uuid + SEPARATOR + filename
                    known_cases[case_id][basename + URL_ATTRIBUTE_SUFFIX] = file_url

        else:
            raise ValueError(file_uuid, filename, data_category, data_type)
    
    elif data_category == GDC_DataCategory.BIOSPECIMEN:
        if data_type == GDC_DataType.BIOSPECIMEN_SUPPLEMENT:
            metadata = _get_clinical_biospecimen_supplement_metadata(file_uuid)
            cases = metadata['cases']
            assert len(cases) > 1, file_uuid
            for case in cases:
                case_id = case['case_id']
                if case_id in known_cases:
                    basename = GDC_TO_FC.constructAttributeName_base(experimental_strategy,
                                                                     workflow_type,
                                                                     data_category,
                                                                     data_type)

                    GDC_FILE_ACCESS.recordFileAccessType(basename, access)

                    known_cases[case_id][basename + UUID_ATTRIBUTE_SUFFIX] = file_uuid + SEPARATOR + filename
                    known_cases[case_id][basename + URL_ATTRIBUTE_SUFFIX] = file_url
        else:
            raise ValueError(file_uuid, filename, data_category, data_type)

    elif data_category == GDC_DataCategory.SNV:                 
        if data_type == GDC_DataType.AGGREGATED_SOMATIC_MUTATION:
            metadata = _get_aggregated_SNV_metadata(file_uuid)
            cases = metadata['cases']
            assert len(cases) > 0, file_uuid

            for case in cases:
                samples = case['samples']
                for sample in samples:
                    sample_id = sample['sample_id']
                    if sample_id in known_samples:
                        sample_type_id = sample['sample_type_id']
                        sample_type_tn = SAMPLE_TYPE.getTumorNormalClassification(sample_type_id)
                        if sample_type_tn == SampleType.TUMOR:

                            basename = GDC_TO_FC.constructAttributeName_base(experimental_strategy,
                                                                             workflow_type,
                                                                             data_category,
                                                                             data_type)
                            GDC_FILE_ACCESS.recordFileAccessType(basename, access)

                            known_samples[sample_id][basename + UUID_ATTRIBUTE_SUFFIX] = file_uuid + SEPARATOR + filename
                            known_samples[sample_id][basename + URL_ATTRIBUTE_SUFFIX] = file_url

        elif data_type == GDC_DataType.MASKED_SOMATIC_MUTATION:
            metadata = _get_aggregated_SNV_metadata(file_uuid)
            cases = metadata['cases']
            assert len(cases) > 0, file_uuid

            for case in cases:
                samples = case['samples']            
                for sample in samples:
                    sample_id = sample['sample_id']
                    if sample_id in known_samples:
                        sample_type_id = sample['sample_type_id']
                        sample_type_tn = SAMPLE_TYPE.getTumorNormalClassification(sample_type_id)
            
                        if sample_type_tn == SampleType.TUMOR:
                            basename = GDC_TO_FC.constructAttributeName_base(experimental_strategy,
                                                                             workflow_type,
                                                                             data_category,
                                                                             data_type)
                            GDC_FILE_ACCESS.recordFileAccessType(basename, access)
                            
                            known_samples[sample_id][basename + UUID_ATTRIBUTE_SUFFIX] = file_uuid + SEPARATOR + filename
                            known_samples[sample_id][basename + URL_ATTRIBUTE_SUFFIX] = file_url

def create_participants_file(cases, manifestFileBasename):
    
    attribute_names = []
    for case_id, case in cases.items():
        for attribute_name in case:
            if attribute_name not in {'submitter_id'} and attribute_name not in attribute_names:
                attribute_names.append(attribute_name)

    #pp.pprint(attribute_names)
    participants_filename = manifestFileBasename + '_participants.txt'
    participant_sets_membership_filename = manifestFileBasename + '_participant_sets_membership.txt'
    
    with open(participants_filename, 'w') as participantsFile, open(participant_sets_membership_filename, 'w') as membershipFile:
        fieldnames = ['entity:participant_id', 'submitter_id'] + attribute_names
        participants_writer = csv.DictWriter(participantsFile, fieldnames=fieldnames, delimiter='\t')
        participants_writer.writeheader()

        fieldnames = ['membership:participant_set_id', 'participant_id']
        membership_writer = csv.DictWriter(membershipFile, fieldnames=fieldnames, delimiter='\t')
        membership_writer.writeheader()
        
        for case_id, case in cases.items():
            entity_row = {'entity:participant_id': case_id, 'submitter_id' : case['submitter_id']}
            for attribute_name in attribute_names:
                if attribute_name in case:
                    entity_row[attribute_name] = case[attribute_name]
                    if attribute_name.endswith(UUID_ATTRIBUTE_SUFFIX):
                        basename = attribute_name[0:-len(UUID_ATTRIBUTE_SUFFIX)]
                        membership_row = {'membership:participant_set_id' : GDC_FILE_ACCESS.getAccessTypePrefix(basename) + "GDC_FILE_RETRIEVAL_" + attribute_name,
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
                        membership_row = {'membership:sample_set_id' : GDC_FILE_ACCESS.getAccessTypePrefix(basename) + "GDC_FILE_RETRIEVAL_" + attribute_name,
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
                        membership_row = {'membership:pair_set_id' : GDC_FILE_ACCESS.getAccessTypePrefix(basename) + "GDC_FILE_RETRIEVAL_" + attribute_name,
                                          'pair_id': pair_id}
                        membership_writer.writerow(membership_row)
                else:
                    entity_row[attribute_name] = '__DELETE__'
            pairs_writer.writerow(entity_row)

            row = {'membership:pair_set_id' : 'ALL',
                   'pair_id': pair_id}
            membership_writer.writerow(row)


def main():

    parser = argparse.ArgumentParser(description='create FireCloud workspace load files from GDC manifest')
    parser.add_argument("manifest", help="manifest file downloaded from the GDC Data Portal")
    parser.add_argument("-r", "--resolve_uuids", help="resolve uuids in the manifest to urls")
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

    for i, item in enumerate(manifestFileList):

        file_uuid = item['id']
        filename = item['filename']
        file_url = uuidResolver.getURL(file_uuid) if uuidResolver is not None else "__DELETE__"
    
        print('{0} of {1}: {2}, {3}'.format(i+1, len(manifestFileList), file_uuid, filename))

        for attempt in range(5):
            try:
   
                url = "{0}/files/{1}".format(GDC_API_ROOT, file_uuid)
                response = requests.get(url, headers=None)
                responseDict = response.json()
    
                get_file_metadata(file_uuid, filename, file_url, cases, samples, 
                                  pairs, deferred_file_uuids)

            except (KeyboardInterrupt, SystemExit):
                raise
            except Exception as x:
                print("Exception=", x)
                print("attempt=", attempt, 'file uuid = ', file_uuid)
                raise
            else:
                break
        else:
            #failed all attempts
            # - just move on
            print("failed 5 attempts! SKIPPING FILE: file uuid = ", file_uuid)
            continue

    for file_uuid in deferred_file_uuids:
        file_url =  file_url = uuidResolver.getURL(file_uuid) if uuidResolver is not None else "__DELETE__"
        process_deferred_file_uuid(file_uuid, file_url, cases, samples)

    #print("cases:")
    #pp.pprint(cases)
    #print("samples:")
    #pp.pprint(samples)
    #print("pairs:")
    #pp.pprint(pairs)


    manifestFileBasename = os.path.splitext(os.path.basename(manifestFile))[0]

    create_participants_file(cases, manifestFileBasename)
    create_samples_file(samples, manifestFileBasename)
    if len(pairs) != 0:
        create_pairs_file(pairs, samples, manifestFileBasename)

if __name__ == '__main__':
    main()
