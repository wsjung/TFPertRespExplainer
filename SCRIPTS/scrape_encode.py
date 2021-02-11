import requests
from bs4 import BeautifulSoup
import re
import os.path
import argparse
import logging
import logging.config
import logging.handlers
import json
import sys
from datetime import datetime
import pandas as pd
import numpy as np


NO_BEDFILES = -1
CONFLICTING_BEDFILES = -2

## default manual review output file
MANUAL_REVIEW_OUTPUT_FILE = 'conflicts_manual_review.warning'

## default output file for logging TFs with no binding locations data on ENCODE
NO_BEDFILE_TF_OUTPUT_FILE = 'tfs_no_bedfiles.lst'


## initialize logger
#logging.config.fileConfig('logging.ini', disable_existing_loggers=False)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

fh = logging.FileHandler('output.log', mode='w', encoding='utf-8')
sh = logging.StreamHandler(sys.stdout)

formatter = logging.Formatter('[%(asctime)s %(levelname)s] %(message)s')

fh.setFormatter(formatter)
sh.setFormatter(formatter)

logger.addHandler(fh)
logger.addHandler(sh)



ENCODE_BASE_URL = 'https://www.encodeproject.org'

ENCODE_TF_URL = ENCODE_BASE_URL + ('/search'
                        '?type=Experiment'
                        '&status=released'
                        '&perturbed=false' 
                        '&assay_title=TF+ChIP-seq'
                        '&biosample_ontology.term_name=K562'
                        '&target.investigated_as=transcription+factor'
                        '&format=json')

ENCODE_FILE_URL = ENCODE_BASE_URL + ('/search'
                            '?type=File'
                            '&format=json'
                            '&frame=object'
                            '&limit=all'
                            '&output_type=conservative%20IDR%20thresholded%20peaks'
                            '&assembly=GRCh38'
                            '&file_format=bed')


def parse_args(argv):
    parser = argparse.ArgumentParser(description='')
    parser.add_argument(
        '-i', '--tfs', required=True,
        help='File containing list of TFs whose binding locations data to download from ENCODE.')
    parser.add_argument(
        '-o', '--output_dir', required=True,
        help='Output directory path where downloaded files will go.')
    parser.add_argument(
        '-n', '--no_beds_file', required=False, default=NO_BEDFILE_TF_OUTPUT_FILE,
        help='Output file to list TFs with no binding locations data on ENCODE.')
    parser.add_argument(
        '-m', '--manual_review_file', required=False, default=MANUAL_REVIEW_OUTPUT_FILE,
        help='Output file for TFs with multiple BED files that need manual review.')
    parsed = parser.parse_args(argv[1:])
    return parsed


def get_tfs_list(tfs_filepath):
    """
    Args:
        tfs_filepath: path to list of tfs whose binding location data to scrape from ENCODE
    """


    logger.info("==> Reading list of TFs to scrape <==")

    logger.info('tf_file: {}'.format(tfs_filepath))

    # list of all tfs in induction dataset intersected with lambert tfs
    df = pd.read_csv(tfs_filepath, 'r', header=None, names=['tf'])
    tf_list = sorted(df['tf'].tolist())

    return tf_list


def query_tf(tf_name, write_json=False):

    logger.info('==> QUERYING TF: ##### {} ##### <=='.format(tf_name))

    # Force return from the server in JSON format
    headers = {'accept': 'application/json'}

    # output file 
    out_file = 'response.{}.out'.format(tf_name)

    # query url with given TF name
    query_url = ENCODE_TF_URL + '&target.label={}'.format(tf_name)

    logger.info(f'query url: {query_url}')

    # GET search result
    response = requests.get(query_url, headers=headers).json()

    if write_json:
        with open(out_file, 'w') as out:
            out.write(json.dumps(response, indent=4))

        logger.info('wrote response to {}'.format(out_file))

    notif = response['notification']
    num_results = response['total']

    logger.info(f'Notif: {notif}:\t({num_results} results found)')

    experiment_list = []

    for experiment in response['@graph']:
        experiment_list.append(experiment['@id'])

    return experiment_list

#    if num_results < 1: 
#        print('### NO RESULTS FOUND: EXITING')
#        return
#
#    if num_results > 1:
#        print('### SEVERAL EXPERIMENTS FOUND: EXITING')

#    experiment_ID = response['@graph'][0]['@id']    # format: "/experiments/ID"
#    
#    print('Result experiment ID: {}'.format(experiment_ID))




def file_query(tf_name, experiment_ID, write_json=False):
    """
        Args:
            tf_name: name of the tf
            experiment_ID: ID of the experiment
            write_json: True to write json file. (dfeault = False)

        Returns:
            bedfile_list: list of bedfile hrefs
    """

    logger.info('==> QUERYING EXPERIMENT ID: ### {} ### <=='.format(experiment_ID))

    # force return from the server in JSON format
    headers = {'accept': 'application/json'}

    # output file
    out_file = 'response.exp.{}.out'.format(tf_name)

    # query url with given experiment id
    query_url = ENCODE_FILE_URL + '&dataset={}'.format(experiment_ID)

    logger.info('query url: {}'.format(query_url))

    # GET search result
    response = requests.get(query_url, headers=headers).json()

    if write_json:
        with open(out_file, 'w') as out:
            out.write(json.dumps(response, indent=4))

        logger.info('wrote file response to {}'.format(out_file))

    notif = response['notification']
    num_results = response['total']

    logger.info(f'Notif: {notif}\t({num_results} results found)')

    bedfile_list = []

    for bed in response['@graph']:
        bedfile_list.append(bed['@id'])

    return bedfile_list

  #  if num_results < 1:
  #      print("### EXPERIMENT HAS NO CONSERVATIVE IDR THRESHOLDED PEAKS")
  #      return

  #  if num_results == 1:
  #      bedfile = response['@graph'][0]['href']

  #  if num_resuts > 1:
  #      print('### SEVERAL FILES FOUND: FILE SELECTION')

  #  download_file(bedfile, tf_name, dl_dir)


def check_beds(bedfiles, tf):

    def all_equal(lst):
        return lst.count(lst[0]) == len(lst)

    if len(bedfiles) < 1:   # no bedfiles in list
        return NO_BEDFILES

    if len(bedfiles) == 1:  # one bedfile in list
        return bedfiles[0]

    # several bedfiles in list
    #   check if they are all from the same lab and pipeline
    labs = []
    pipelines = []
    dates = []

    for bedfile in bedfiles:
        r = requests.get(ENCODE_BASE_URL + bedfile, allow_redirects=True)

        soup = BeautifulSoup(r.text, 'html.parser')

        logger.info('url: {}'.format(ENCODE_BASE_URL + bedfile))

        # find pipeline info
        pipeline = soup.find_all('div', attrs={'data-test':lambda x: x=='pipelines'})[0].dd.span.a.contents[0]
        pipelines.append(pipeline)

        # find lab info
        lab = soup.find_all('div', attrs={'data-test':lambda x: x=='lab'})[0].dd.contents[0]
        labs.append(lab)

        # find dates
        date = soup.find_all('div', attrs={'data-test':lambda x: x=='datecreated'})[0].dd.contents[0]
        dates.append(datetime.strptime(date, '%Y-%m-%d'))

    if all_equal(labs) and all_equal(pipelines):    # if all labs and pipelines are equal --> return the bedfile with the most recent date
        return bedfiles[np.argmax(dates)]

    else:   # otherwise, write all bedfiles for manual checking
        dates_str = [datetime.strftime(date, '%Y-%m-%d') for date in dates]

        logger.warning('### CONFLICTING BEDFILES (MANUAL REVIEW REQUIRED) ###')
        logger.warning('bedfiles:\t{}'.format(bedfiles))
        logger.warning('labs:\t{}'.format(labs))
        logger.warning('pipelines:\t{}'.format(pipelines))
        logger.warning('dates:\t{}'.format(dates_str))

        # write to conflicts file
        with open(MANUAL_REVIEW_OUTPUT_FILE, 'a') as out:
            out.write("TF: {}\n\tbedfiles:\t{}\n\tlabs:\t{}\n\tpipelines:\t{}\n\tdates:\t{}\n\n".format(
                         tf,
                         '\t'.join(bedfiles),
                         '\t'.join(labs),
                         '\t'.join(pipelines),
                         '\t'.join(dates_str)
            )
        )

        return CONFLICTING_BEDFILES



def download_file(bedfile, tf_name, dl_dir):

    def getfilename(cd):
    #returns filename from content disposition

        if not cd:
            return None

        fname = re.findall('filename=(.+)', cd)

        if len(fname) == 0:
            return None

        return fname[0]
    
    dl_url = '{}{}'.format(ENCODE_BASE_URL,bedfile)

    print('download url: {}'.format(dl_url))

    response = requests.get(dl_url, allow_redirects=True)
    filename = getfilename(response.headers.get('content-disposition'))

    dl_path = os.path.join(dl_dir, '{}/{}'.format(tf_name, filename))

    print('downloading {} to {}'.format(filename, dl_path))

    open(dl_path, 'wb').write(response.content)

    print('download complete')


def scrape(tfs):

    logger.info('Number of TFs to scrape: {}'.format(len(tfs)))

    logger.info('==> Begin scrape <==')

    for tf in tfs:
        experiments = query_tf(tf)
        logger.info('experiments: {}'.format(experiments))

        if len(experiments) < 1:     # no results for TF
            logger.warning('NO EXPERIMENTS for {}'.format(tf))
            # write tf to list of tfs with no bedfiles
            with open(NO_BEDFILE_TF_OUTPUT_FILE, 'a') as out:
                out.write('{}\n'.format(tf))

            continue

        if len(experiments) == 1:    # there exists results
            # query file
            bedfile_list = file_query(tf, experiments[0])

            bedfile = check_beds(bedfile_list, tf)

            if bedfile == NO_BEDFILES:   # no files of interest
                logger.warning('NO BEDFILES for {}'.format(tf))
                # write tf to list of tfs with no bedfiles
                with open(NO_BEDFILE_TF_OUTPUT_FILE, 'a') as out:
                    out.write('{}\n'.format(tf))
                continue

            if len(bedfile_list) > 0:   # several matches of interest
                logger.info('Available BED files for {}: {}'.format(tf, bedfile_list))
                continue


def main(argv):
    # parse arguments
    args = parse_args(argv)

    # set optional arguments
    NO_BEDFILE_TF_OUTPUT_FILE = args.no_beds_file
    MANUAL_REVIEW_OUTPUT_FILE = args.manual_review_file

    logger.info('Input arguments: {}'.format(args))

    tf_list_path = args.tfs
    dl_dir = args.output_dir

    # get list of tfs 
    tfs_list = get_tfs_list(tf_list_path)

    # do the scrape
    scrape(tfs_list)



if __name__ == '__main__':

    main(sys.argv)
