#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import requests
import json
import click
import time
import subprocess
from random import random

@click.command()
@click.option('-f','--file',
              required=True,
              help="Pass vcf file to send (normally just one chromosome)")
@click.option('-u','--url',
              default="http://127.0.0.1:8080",
              help="Address and port of a functioning Michigan Imputation Server instance")
@click.option('-t','--token',
              required=True,
              help="Personal token from the imputation server (tokens expire in 30 days)")
@click.option('-p','--password',
              default="password",
              help="Put any possword, and remember it for the extraction of the results")
@click.option('-a','--api',
              default="/api/v2/jobs/submit/imputationserver@1.4.1",
              help="API path in the server")
@click.option('-pt','--patient',
              help="Patient name")
def main(file,
         url,
         token,
         password,
         api,
         patient):

    files = [('files', open(file, 'rb'))]

    headers = {'X-Auth-Token' : token }
    data = {
      'refpanel': f'apps@1000g-phase-3-v5@2.0.0',
      'population': 'eur',
      'build': 'hg38',
      'password': f'{password}'
    }

    string_to_int = f"{patient}_{file}"
    time_to_sleep = int(str(hash(string_to_int))[1:4])/1000
    time.sleep(time_to_sleep) # to prevent weird error of sending all chromosomes and failing
    r = requests.post(url + api, files=files, data=data, headers=headers)

    success = r.json()["success"]
    job_id = r.json()["id"]

    # check when it is done
    is_imputing = True
    while is_imputing:
        status = requests.get(url + f"/api/v2/jobs/{job_id.strip()}/status", headers=headers).json()
        time.sleep(5) # check every 5 seconds
        if status["complete"] == True and status["state"] == 5:
            # it failed, need to re-send it
	        #print("RE-SENT to impute")
            #r = requests.post(url + api, files=files, data=data, headers=headers)
		    #time.sleep(2)
		    #job_id = r.json()["id"]
            is_imputing = False
            exit(1)
        elif status["complete"] == True and status["state"] == 4:
            is_imputing = False
    time.sleep(4) # it will give time for folder structure to be recognized
    # if not nextflow will search for the path and will not exist

    # send chr,job_id to stdout so nextflow to catch it and save it to file
    subprocess.run(f"echo {job_id}", shell=True, check=True, stderr=subprocess.STDOUT)
    #subprocess.run(f"echo {string_to_int}_{time_to_sleep}", shell=True, check=True, stderr=subprocess.STDOUT)

if __name__ == "__main__":
    main()
