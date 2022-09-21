import sys
import argparse
import os
import firecloud.api as fapi
from fcgdctools.manifest_downloader import build_filter_json, download_manifest
from fcgdctools.fc_loadfiles import GDC_API_ROOT, GDC_AWG_API_ROOT, main as gen_loadfiles

FILE_TYPE_DICT = {
    "default": ["open"],
    "TCGA-dbGaP-Authorized": ["open", "controlled"],
    "TARGET-dbGaP-Authorized": ["open", "controlled"]

}

TCGA_AUTH_DOMAIN_NAME = "TCGA-dbGaP-Authorized"
TARGET_AUTH_DOMAIN_NAME = "TARGET-dbGaP-Authorized"

def prepare_workspace_attribute_list(attr_file, auth_domain):
    fp = open(attr_file, 'r')
    attr_names = fp.readline().strip().split(":")[1].split("\t")
    attr_values = fp.readline().split("\t")
    attrs = dict()
    for i in range(len(attr_names)):
        attrs[attr_names[i]] = attr_values[i]

    if auth_domain is not None:
        attrs["token_file"] = "file_path_for_gdc_token_file"

    return attrs

def list_downloadable_attrs(prefix, entities):
    downloadable_attr_names = []
    for ent in entities:
        filename = prefix + "_" + ent + "s.tsv"
        if os.path.exists(filename):
            fp = open(filename, 'r')
            attribute_list = fp.readline().strip().split(":")[1].split("\t")
            for attr in attribute_list:
                if attr.endswith("gdc_url"):
                    downloadable_attr_names.append((attr, ent))

    return downloadable_attr_names


def create_method_configs(billing_project, ws_name, attr_list, auth_domain, awg_flag):

    method_namespace = "getzlab"
    file_downloader_method = "gdc_api_file_download"
    file_downloader_method_snapshot_id = 4

    for attr in attr_list:
        
        attr_name = attr[0]
        attr_name_base = attr_name.rsplit('__', 1)[0]
        attr_entity = attr[1]
        
        new_config_name = "gdc_api_file_download__" + attr_name_base
        print("Uploading and configuring method config {0}, based on {1}".format(new_config_name, file_downloader_method))
        new_config = fapi.get_config_template(method_namespace,
                                              file_downloader_method,
                                              file_downloader_method_snapshot_id)
        fapi._check_response_code(new_config, 200)
        new_config = new_config.json()
        new_config['name'] = new_config_name
        new_config['namespace'] = method_namespace
        new_config['rootEntityType'] = attr_entity

        inputs = new_config['inputs']
        outputs = new_config['outputs']

        if auth_domain is not None or awg_flag:
            inputs['gdc_api_file_download.download_file.gdc_token'] = "workspace.token_file"
            
        inputs['gdc_api_file_download.download_file.url'] = "this.{0}".format(attr_name)
        inputs['gdc_api_file_download.download_file.filename'] = "this.{0}".format(attr_name_base + "__gdc_filename")
        inputs['gdc_api_file_download.download_file.file_size'] = "this.{0}".format(attr_name_base + "__gdc_filesize")
        
        outputs['gdc_api_file_download.downloaded_file'] = "this.{0}".format(attr_name_base)

        fapi.create_workspace_config(billing_project, ws_name, new_config)
            

def main():

    parser = argparse.ArgumentParser(description='create FireCloud workspace for a specific project + cohort + access type')
    parser.add_argument("project_name", help="the name of the project. e.g: TCGA")
    parser.add_argument("cohort_name", help="name_of_cancer_cohort. e.g: LUAD")
    parser.add_argument("billing_project", help="name of billing project to create the workspace under. e.g: broad-firecloud-tcga")
    parser.add_argument("ws_suffix", help="descriptive suffix to add to the workspace auto-generated name. e.g: ControlledAccess_hg38_V1-0_DATA")
    parser.add_argument("-d", "--auth_domain", help="authorization domain. for dbGaP controlled access the domain name is TCGA-dbGaP-Authorized.")
    parser.add_argument("-a", "--awg", action="store_true", help="use AWG api endpoint. Requires a valid GDC download token.")
    parser.add_argument("-t", "--token", type=argparse.FileType('r'), help="GDC download token to access controlled APIs. Required for AWG API.")
    
    args = parser.parse_args()
    api_root = GDC_API_ROOT
    if args.awg:
        if args.token is None:
            sys.exit('AWG API requires a GDC download token\n' + parser.format_help())
        api_root = GDC_AWG_API_ROOT
        
    token = args.token.read() if args.token is not None else None
    token_path = os.path.abspath(args.token.name) if args.token is not None else None

    #STEP 1:
    #Create new directory for the cohort and switch wd to this directory
    if args.auth_domain is not None:
        new_dir_name = args.project_name + "-" + args.cohort_name + "_" + args.auth_domain
    else:
        new_dir_name = args.project_name + "-" + args.cohort_name
    os.mkdir(new_dir_name)
    print("Created new directory for the {0} cohort".format(args.cohort_name))
    os.chdir(new_dir_name)
    print("Switched working directory to ./{0}".format(new_dir_name))

    #STEP 2:
    #Create criteria for downloading manifest and then download it.    
    #Right now the file types that are selected for a new workspace depend on whether that workspace is to have open/controlled access to the GDC data portal.
    #This code will need to be redesigned, or new keys will have to be added to the dictionary if this assumption ever changes.
    if args.auth_domain is not None:
        file_types = FILE_TYPE_DICT[args.auth_domain] 
    else: 
        file_types = FILE_TYPE_DICT["default"]

    filters = dict()
    filters["cases.project.program.name"] = [args.project_name]
    filters["cases.project.project_id"] = [args.project_name+"-"+args.cohort_name]
    filters["files.access"] = file_types
    #Following directions from the GDC, we were told that controlled access workspaces should not contain BAM files
    if args.auth_domain is not None:
        filters["files.data_format"] = ["BCR XML","CDC JSON","TXT","VCF","TSV","MAF","XLSX","BEDPE","BAM"]
    else:
        filters["files.data_format"] = ["BCR XML","CDC JSON","TXT","VCF","TSV","MAF","XLSX","BEDPE"]

    filt_json = build_filter_json(filters)

    #Download manifest file to the new directory
    manifest_filename = download_manifest(filt_json, api_root, token)
    print("manifest downloaded")
    
    #Step 3:
    #Generate loadfiles from the manifest file
    gen_loadfile_args = []
    if args.project_name == "TARGET":
        gen_loadfile_args += ["-c"]
    if args.token is not None:
        gen_loadfile_args += ["-t", token_path]
    if args.awg:
        gen_loadfile_args += ["-a", "awg"]
    gen_loadfile_args += [manifest_filename]
    
    gen_loadfiles(gen_loadfile_args)

    #Step 4:
    #Prepare attributes to be loaded
    workspace_attribute_filename = manifest_filename.split(".")[0] + "_workspace_attributes.tsv"
    attribute_list = prepare_workspace_attribute_list(workspace_attribute_filename, args.auth_domain)
    
    #Step 5:
    #Create the new workspace on FireCloud
    workspace_name = "{0}_{1}_{2}".format(args.project_name, args.cohort_name, args.ws_suffix)
    print("New workspace name is: {0}\nPreparing to create workspace.".format(workspace_name))
    fapi.create_workspace(args.billing_project, workspace_name, args.auth_domain, attribute_list)

    #Step 6:
    #Upload data model .tsv files to the newly created workspace
    data_model_file_prefix = manifest_filename.split(".")[0]
    data_files = ["participants", "participant_sets_membership", "samples", "sample_sets_membership", "pairs", "pair_sets_membership"]
    for filetype in data_files:
        full_name = data_model_file_prefix + "_" + filetype + ".tsv"
        if os.path.exists(full_name):
            print("Uploading file {0}".format(full_name))
            fapi.upload_entities_tsv(args.billing_project, workspace_name, full_name)

    #Step 7:
    #Create and Upload method configurations for downloading files to the new workspace
    downloadable_attrs = list_downloadable_attrs(data_model_file_prefix, ["participant", "sample", "pair"])
    print("The downloadable attributes are:")
    for attr in downloadable_attrs:
        print(attr[0])
    create_method_configs(args.billing_project, workspace_name, downloadable_attrs, args.auth_domain, args.awg)

if __name__ == '__main__':
    main()


