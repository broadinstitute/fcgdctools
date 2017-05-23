import requests
import json
import sys
import argparse
import os
import datetime
import firecloud.api as api
from manifest_downloader import build_filter_json, download_manifest

FILE_TYPE_DICT = {
	"default": ["open"],
	"dbGapAuthorizedUsers": ["open", "controlled"]

}

def prepare_workspace_attribute_list(attr_file, auth_domain):
	fp = open(attr_file, 'r')
	attr_names = fp.readline().strip().split(":")[1].split("\t")
	attr_values = fp.readline().split("\t")
	attrs = dict()
	for i in range(len(attr_names)):
		attrs[attr_names[i]] = attr_values[i]

	if auth_domain == "dbGapAuthorizedUsers":
		attrs["token_file"] = "file_path_for_gdc_token_file"

	return attrs

def list_downloadable_attrs(prefix, entities):
	downloadable_attr_names = []
	for ent in entities:
		filename = prefix + "_" + ent + "s.txt"
		if os.path.exists(filename):
			fp = open(filename, 'r')
			attribute_list = fp.readline().strip().split(":")[1].split("\t")
			for attr in attribute_list:
				if attr.endswith("uuid_and_filename"):
					downloadable_attr_names.append((attr, ent))

	return downloadable_attr_names


def create_method_configs(project, ws_name, attr_list, auth_domain):

	config_namespace = "broadinstitute_cga"
	file_downloader_name = "gdc_file_downloader__default_cfg"
	bam_downloader_name = "gdc_bam_downloader__default_cfg"
	snapshot_id = 1

	for attr in attr_list:
		
		attr_name = attr[0]
		attr_entity = attr[1]
		if "aligned_reads" in attr_name:
			new_config_name = "gdc_bam_downloader__" + attr_name[:-17] + "cfg"
			print("Uploading and configuring method config {0}, based on {1}".format(new_config_name, bam_downloader_name))
			api.copy_config_from_repo(project, ws_name, config_namespace, bam_downloader_name, snapshot_id, config_namespace, new_config_name)
			
			current_config = api.get_workspace_config(project, ws_name, config_namespace, new_config_name)
			current_config = current_config.json()

			inputs = current_config['inputs']
			outputs = current_config['outputs']

			inputs['gdc_bam_downloader_workflow.uuid_and_filename'] = "this.{0}".format(attr_name)
			
			outputs['gdc_bam_downloader_workflow.gdc_bam_downloader.bam_file'] = "this.{0}__bam_url".format(attr_name)
			outputs['gdc_bam_downloader_workflow.gdc_bam_downloader.bai_file'] = "this.{0}__bai_url".format(attr_name)

			current_config['inputs'] = inputs
			current_config['outputs'] = outputs
			current_config['rootEntityType'] = attr_entity
			
			api.update_workspace_config(project, ws_name, config_namespace, new_config_name, current_config)

		else:
			new_config_name = "gdc_file_downloader__" + attr_name[:-17] + "cfg"
			print("Uploading and configuring method config {0}, based on {1}".format(new_config_name, file_downloader_name))
			api.copy_config_from_repo(project, ws_name, config_namespace, file_downloader_name, snapshot_id, config_namespace, new_config_name)
			
			current_config = api.get_workspace_config(project, ws_name, config_namespace, new_config_name)
			current_config = current_config.json()

			inputs = current_config['inputs']
			outputs = current_config['outputs']

			if auth_domain != "dbGapAuthorizedUsers":
				inputs.pop('gdc_file_downloader_workflow.gdc_file_downloader.gdc_user_token', None)
				
			inputs['gdc_file_downloader_workflow.uuid_and_filename'] = "this.{0}".format(attr_name)
			
			outputs['gdc_file_downloader_workflow.gdc_file_downloader.file'] = "this.{0}__url".format(attr_name)
			
			current_config['inputs'] = inputs
			current_config['outputs'] = outputs
			current_config['rootEntityType'] = attr_entity

			api.update_workspace_config(project, ws_name, config_namespace, new_config_name, current_config)
			

def main():

    parser = argparse.ArgumentParser(description='create FireCloud workspace for a specific project + cohort + access type')
    parser.add_argument("project_name", help="the name of the project. e.g: TCGA")
    parser.add_argument("cohort_name", help="name_of_cancer_cohort. e.g: LUAD")
    parser.add_argument("billing_project", help="name of billing project to create the workspace under. e.g: broad-firecloud-tcga")
    parser.add_argument("ws_suffix", help="descriptive suffix to add to the workspace auto-generated name. e.g: ControlledAccess_hg38_V1-0_DATA")
    parser.add_argument("-a", "--auth_domain", help="authorization domain. for dbGaP controlled access the domain name is dbGapAuthorizedUsers.", default="")
    
    args = parser.parse_args()

    #STEP 1:
    #Create new directory for the cohort and switch wd to this directory
    if args.auth_domain:
    	new_dir_name = args.project_name + "-" + args.cohort_name + "_" + args.auth_domain
    else:
    	new_dir_name = args.project_name + "-" + args.cohort_name
    os.mkdir(new_dir_name)
    print("Created new directory for the {0} cohort".format(args.cohort_name))
    os.chdir(new_dir_name)
    print("Switched working directory to ./{0}".format(new_dir_name))

	#STEP 2:
	#Create creteria for downloading manifest and then download it.    
    #Right now the file types that are selected for a new workspace depend on whether that workspace is to have open/controlled access to the GDC data portal.
    #This code will need to be redesigned, or new keys will have to be added to the dictionary if this assumption ever changes.
    if args.auth_domain:
    	file_types = FILE_TYPE_DICT[args.auth_domain] 
    else: 
    	file_types = FILE_TYPE_DICT["default"]

    filters = dict()
    filters["cases.project.program.name"] = [args.project_name]
    filters["cases.project.project_id"] = [args.project_name+"-"+args.cohort_name]
    filters["files.access"] = file_types

    filt_json = build_filter_json(filters)

    #Download manifest file to the new directory
    manifest_filename = download_manifest(filt_json)
    print("manifest downloaded")
    
    #Step 3:
    #Run fcgdctools on the manifest file
    fcgdctools_command = "python3 ~/Projects/fcgdctools/fcgdctools/fc_loadfiles.py " + manifest_filename + ">genFcWsLoadFiles_output.txt"
    print("Executing command {0}\nPlease check the output file to see progress and check for errors.".format(fcgdctools_command))
    os.system(fcgdctools_command)
    
    #Step 4:
    #Prepare attributes to be loaded
    workspace_attribute_filename = manifest_filename.split(".")[0] + "_workspace_attributes.txt"
    attribute_list = prepare_workspace_attribute_list(workspace_attribute_filename, args.auth_domain)
    
    #Step 5:
    #Create the new workspace on FireCloud
    workspace_name = "{0}_{1}_{2}".format(args.project_name, args.cohort_name, args.ws_suffix)
    print("New workspace name is: {0}\nPreparing to create workspace.".format(workspace_name))
    api.create_workspace(args.billing_project, workspace_name, args.auth_domain, attribute_list)

    #Step 6:
    #Upload data model .tsv files to the newly created workspace
    data_model_file_prefix = manifest_filename.split(".")[0]
    data_files = ["participants", "participant_sets_membership", "samples", "sample_sets_membership", "pairs", "pair_sets_membership"]
    for filetype in data_files:
    	full_name = data_model_file_prefix + "_" + filetype + ".txt"
    	if os.path.exists(full_name):
    		print("Uploading file {0}".format(full_name))
    		api.upload_entities_tsv(args.billing_project, workspace_name, full_name)

    #Step 7:
    #Upload method configurations for downloading files to the new workspace
    downloadable_attrs = list_downloadable_attrs(data_model_file_prefix, ["participant", "sample", "pair"])
    print("The downloadable attributes are:")
    for attr in downloadable_attrs:
    	print(attr[0])

    #Step 8:
    create_method_configs(args.billing_project, workspace_name, downloadable_attrs, args.auth_domain)

if __name__ == '__main__':
    main()


