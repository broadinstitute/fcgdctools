#! /usr/bin/env python3
import argparse

from .manifest_downloader import build_filter_json, download_manifest_post
from .fc_loadfiles import GDC_API_ROOT
from pprint import pprint

def main():
    """Generate a manifest from a list of aliquots"""
    parser = argparse.ArgumentParser(description="Generate a manifest from a list of aliquots")
    parser.add_argument('aliquots', help='File containing list of aliquots')
    parser.add_argument('-p', '--program', help='Program Name')
    parser.add_argument('-j', '--project', help='Project ID')
    parser.add_argument('-e', '--experimental_strategy', default=['WGS'], nargs='*',
                        help='Experimental Strategies')
    parser.add_argument('-d', '--data_format', default=['BAM'], nargs='*',
                        help='Data Formats')

    args = parser.parse_args()
    pprint(args)
    filters = {}
    aliquots = []
    with open(args.aliquots, 'r') as aliquot_ids:
        for a_id in aliquot_ids:
            aliquots.append(a_id.strip())

    filters["cases.samples.portions.analytes.aliquots.submitter_id"] = aliquots
    filters["files.data_format"] = args.data_format
    if args.program:
        filters["cases.project.program.name"] = [args.program]
    if args.project:
        filters["cases.project.project_id"] = [args.project]
    if args.experimental_strategy:
        filters["files.experimental_strategy"] = args.experimental_strategy

    filt_json = build_filter_json(filters)
    manifest_filename = download_manifest_post(filt_json, GDC_API_ROOT)

if __name__ == '__main__':
    main()