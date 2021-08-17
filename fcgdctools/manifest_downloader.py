import requests
import json
import sys
import argparse
import os
import datetime
import firecloud.api as api

def build_filter_json(filter_attrs):
    filt = {
        "op":"and",
        "content":[]
    }

    for key in filter_attrs:
        search_attr = {
            "op":"in",
            "content":{
                "field":key,
                "value":filter_attrs[key]
            }
        }

        filt["content"].append(search_attr)

    return filt


def download_manifest(filt_json, api_root, token=None):
    
    #This is the API endpoint for performing a search on the GDC data portal and retrieving file information.
    files_endpt = api_root + '/files'

    #Creating a new name for the manifest file
    timestamp='{:%Y-%m-%d_%H-%M-%S}'.format(datetime.datetime.now())
    manifest_filename="gdc_manifest_"+timestamp+".tsv"
    print("downloading manifest {0}".format(manifest_filename))

    params = {'filters':json.dumps(filt_json), 'from':0, 'size':30000,
              'return_type':'manifest'}
    
    # requests URL-encodes automatically
    headers = {}
    
    if token:
        headers['X-Auth-Token'] = token
    
    #Writing the output to the manifest file
    with open(manifest_filename, 'wb') as handle:
        while True:
            response = requests.get(files_endpt, params=params, headers=headers)
            for ln_ct, line in enumerate(response.iter_lines()):
                # Skip manifest header if not on first page
                if ln_ct == 0 and params['from'] > 0:
                    continue
                handle.write(line + b'\n')
            if ln_ct == params['size']:
                params['from'] += params['size']
            else:
                break
            

    return manifest_filename

