# DSSP and PDB features #

import shutil
import os
import time
import argparse

start = time.ctime()

main_path = os.getcwd()
scripts_path = main_path+"/scripts"
model_path = main_path+"/fitted_model"

def main():
    parser = argparse.ArgumentParser(description='### Here should be the description! ###')
    parser.add_argument("input", help="The name of your PDB file", type=str)
    parser.add_argument("-v", "--verbosity", help="Add structural features", action='store_true')
    parser.add_argument("-l", "--launch", help="Launch DSSP locally (specify '-l') or via server", action="store_true")
    args = parser.parse_args()
    pdb_file = args.input
    pdb_id = pdb_file.split('.')[0]

    if f'{pdb_id}_result' not in os.listdir(main_path):
        os.mkdir(f'{main_path}/{pdb_id}_result')
    result_path = f'{main_path}/{pdb_id}_result'

    print("######## Step 1 -- Creating of DSSP file ########")
    if args.launch:
        os.system(f"mkdssp -i {pdb_file} -o {pdb_id}.dssp")
    else:
        os.system(f"python {scripts_path}/get_dssp_1.py {pdb_file}")

    print("######## Step 2 -- Selecting of structural features ########")
    os.system(f"python {scripts_path}/extract_bf_2_1.py {pdb_file}")
    os.system(f"python {scripts_path}/form_dataset_2_2.py {pdb_id}.dssp")
    os.remove(f"{pdb_id}_bf.txt")

    print("######## Step 3 -- Predicting of the probabilities ########")
    if args.verbosity:
        os.system(f"python {scripts_path}/make_predictions_3.py {pdb_id} {model_path} -v")
    else:
        os.system(f"python {scripts_path}/make_predictions_3.py {pdb_id} {model_path}")

    os.remove(f"{pdb_id}.dssp")

    if os.path.isfile(f'Helix_{pdb_id}_test.csv'):
        os.remove(f"Helix_{pdb_id}_test.csv")

    if os.path.isfile(f'B-sheet_{pdb_id}_test.csv'):
        os.remove(f"B-sheet_{pdb_id}_test.csv")

    if os.path.isfile(f'Loop_{pdb_id}_test.csv'):
        os.remove(f"Loop_{pdb_id}_test.csv")

    if args.verbosity:
        shutil.move(f"{main_path}/proba_{pdb_id}_f.csv",f"{result_path}/proba_{pdb_id}_f.csv")
    else:
        shutil.move(f"{main_path}/proba_{pdb_id}.csv",f"{result_path}/proba_{pdb_id}.csv")

if __name__ == "__main__":

    main()

    end = time.ctime()
    print(f"\n******** {__file__} ********\nStart: {start}\nFinish: {end}\n")
