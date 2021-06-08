# Download DSSP file #
# The 1st script in workflow

import json
import requests
import time
import os
import argparse

start = time.ctime()

main_url = 'https://www3.cmbi.umcn.nl/xssp/'

def pdb_to_dssp(pdb_file,main_url):

	url_create = f'{main_url}api/create/pdb_file/dssp/'

	response_create = requests.post(url_create,files={'file_':open(pdb_file,'rb')})
	response_create.raise_for_status

	job_id = json.loads(response_create.text)['id']
	#print(f"Job submitted successfully. Id is: '{job_id}'")

	ready = False
	while not ready:

		url_status = f'{main_url}/api/status/pdb_file/dssp/{job_id}/'
		response_status = requests.get(url_status)
		response_status.raise_for_status

		status = json.loads(response_status.text)['status']
		#print(f'Job status is {status}')

		if status == 'SUCCESS':
			ready = True
		elif status in ['FAILURE','REVOKED']:
			print('Error')
			break
		else:
			time.sleep(5)
	else:
		url_result = f'{main_url}/api/result/pdb_file/dssp/{job_id}/'

		response_result = requests.get(url_result)
		response_result.raise_for_status
		result = json.loads(response_result.text)['result']

		print(pdb_file[:4]+'.dssp is downloaded successfully!')
		return result

def main():

	parser = argparse.ArgumentParser()
	parser.add_argument("input", help="The name of your PDB file", type=str)
	args = parser.parse_args()
	id_pdb = args.input.split('.')[0]

	result = pdb_to_dssp(f'{id_pdb}.pdb', main_url)

	with open(f'{id_pdb}.dssp','w') as file:
		file.write(result)

if __name__ == "__main__":
	main()
	end = time.ctime()

	print(f'\n******** {__file__} ********\nStart: {start}\nFinish: {end}\n')
