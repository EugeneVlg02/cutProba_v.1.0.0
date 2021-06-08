import argparse
import time
import os

import pickle
import pandas as pd

start = time.ctime()

def main(pdb_id):

    result = pd.DataFrame()

    if os.path.isfile(f"Helix_{pdb_id}_test.csv"):
        helix_data = pd.read_csv(f"Helix_{pdb_id}_test.csv")
        helix_X = helix_data[["ACC_N_com","bfac_N_sep"]].values

        helix_model = pickle.load(open(f"{args.model_path}/helix_gnb.sav",'rb'))
        helix_proba = helix_model.predict_proba(helix_X)[:,1]

        if args.verbosity:
            helix_result = pd.concat([helix_data,pd.Series(helix_proba,name='probability')],axis=1)
        else:
            helix_result = pd.concat([helix_data[['num_aa','chain','AA']],pd.Series(helix_proba,name='probability')],axis=1)
        result = result.append(helix_result,ignore_index=True)

    if os.path.isfile(f"B-sheet_{pdb_id}_test.csv"):
        bsheet_data = pd.read_csv(f"B-sheet_{pdb_id}_test.csv")
        bsheet_X = bsheet_data[["ACC_N_com","bfac_N_sep"]].values

        bsheet_model = pickle.load(open(f"{args.model_path}/bsheet_lda.sav",'rb'))
        bsheet_proba = bsheet_model.predict_proba(bsheet_X)[:,1]

        if args.verbosity:
            bsheet_result = pd.concat([bsheet_data,pd.Series(bsheet_proba*0.67,name='probability')],axis=1)
        else:
            bsheet_result = pd.concat([bsheet_data[['num_aa','chain','AA']],pd.Series(bsheet_proba*0.67,name='probability')],axis=1)
        result = result.append(bsheet_result,ignore_index=True)

    if os.path.isfile(f"Loop_{pdb_id}_test.csv"):
        loop_data = pd.read_csv(f"Loop_{pdb_id}_test.csv")
        loop_X = loop_data[["ACC_N_com","bfac_N_sep","dist_loop_N_com","len_loop_N_sep"]].values

        loop_model = pickle.load(open(f"{args.model_path}/loop_lgr.sav",'rb'))
        loop_proba = loop_model.predict_proba(loop_X)[:,1]

        if args.verbosity:
            loop_result = pd.concat([loop_data,pd.Series(loop_proba,name='probability')],axis=1)
        else:
            loop_result = pd.concat([loop_data[['num_aa','chain','AA']],pd.Series(loop_proba,name='probability')],axis=1)
        result = result.append(loop_result,ignore_index=True)

    result = result.sort_values(by=['chain','num_aa'])
    #print(result)
    max_proba = result.loc[result['probability'] == result['probability'].max()]
    min_proba = result.loc[result['probability'] == result['probability'].min()]
    #print(f"\n******** Maximum of probability ********\n{max_proba}")
    #print(f"\n******** Minimum of probability ********\n{min_proba}")
    if args.verbosity:
        result.to_csv(f'proba_{pdb_id}_f.csv',index=False)
    else:
        result.to_csv(f'proba_{pdb_id}.csv',index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Predicting the probabilities of proteolytic events')
    parser.add_argument("input", help="ID of your PDB file", type=str)
    parser.add_argument("model_path", help="Path to model trained", type=str)
    parser.add_argument("-v", "--verbosity", help="Add structural features", action="store_true")
    args = parser.parse_args()

    main(args.input)

    end = time.ctime()
    print(f"\n******** {__file__} ********\nStart: {start}\nFinish {end}\n")
