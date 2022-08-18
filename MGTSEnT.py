import argparse
import numpy as np
import pandas as pd
import time
import scipy.stats as stats


def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-m','--mgtdb',help='.csv file of the whole MGT typing dataset in MGTdb of Salmonella Enteritidis')
    parser.add_argument('-i','--isolatelist',help='.txt file of isolates of interest in MGTdb')
    parser.add_argument('-f', '--flag', help='.csv flag file included in the Github pakage, e.g. /srv/MGTSEnT/10-01-mgtmdrflag.csv')
    parser.add_argument('-o','--outpath_prefix',help='output path + prefix of the outfile, e.g. /srv/outdir/test_')

    args = parser.parse_args()
    return args

# ### Historical input for parseargs
# python3 A_MGTSEnter.py -i /mnt/e/2018/2019-06-14-Australia_SEN/Australia_isolates.txt -o /mnt/e/2018/2019-06-14-Australia_SEN/4_MGTSEnter/2_Austest/Australia_
# isolatelist = '/mnt/e/2018/2019-06-14-Australia_SEN/4_MGTSEnter/1_08-28/allisolist.txt'
# isolatelist = '/mnt/e/2018/2019-06-14-Australia_SEN/Australia_isolates.txt'  ###new_isolates.txt
# # outpath = '/mnt/e/2018/2019-06-14-Australia_SEN/4_MGTSEnter/11_16_all/40K_'
# outpath = '/mnt/e/2018/2019-06-14-Australia_SEN/4_MGTSEnter/2_20221023/Australia_'
# # typelistpath = 'pathofodc_list.txt'  ### if offer the type list in advance

## the inputlist is the same input list for A2021_04_23_time_pwd_odc.py, the outputs of which were used here for tables production.

def main():
    args=parseargs()
    metapath = args.mgtdb
    isolatelist = open(args.isolatelist, 'r').read().splitlines()
    mgt_flags = args.flag
    outpath = args.outpath_prefix

    ######## the following files should be included in the github package
    #metapath = '/mnt/e/2018/2019-06-14-Australia_SEN/0_meta/MGTdataset_corre.csv'
    #mgt_flags = '/mnt/e/2018/2019-06-14-Australia_SEN/2020-03-06-AR/2021-01-40K_output/MDR/10-01-mgtmdrflag.csv'
    feature = 'MDR_MGT-ST'
    ARGs = 'False'
    ### The AMR genes based prediction of MDR, AMR_2,1,0. At least two columns are required: Strain, AMR_group
    argpredfile = '/mnt/e/2018/2019-06-14-Australia_SEN/2020-03-06-AR/2021-01-40K_output/MDR/isolistamrdf.csv'

    # mgt_level = 'ODC5'
    # number_iso_threshold = 1
    t1 = time.time()
    metadf = pd.read_csv(metapath, low_memory=False)
    isodic = {a :'True' for a in isolatelist}
    metadf['new_isolate'] = metadf['Strain'].map(isodic)
    metadf = metadf[metadf['new_isolate']=='True']
    flags = pd.read_csv(mgt_flags)

    ### step 1. add PPS and MDR isolates
    key_list = ['Clades','Lineages',feature]
    subdf = metadf.copy(deep=True) ### deep=Ture will avoid the metadf being modified after the following function
    accflag_dic = flag_reprot(flags, subdf, key_list)
    flagdf = pd.DataFrame.from_dict(accflag_dic)
    flagdf.index.name = 'Strain'
    ppsdf = flagdf.merge(metadf,left_on='Strain',right_on='Strain', how='left')
    ppsdf.to_csv(outpath + 'MGTSEnter.csv',index=False)
    if ARGs == 'True':
        ppsdf = pd.read_csv(outpath + 'MGTSEnter.csv',low_memory=False)
        argpredf = pd.read_csv(argpredfile,low_memory=False)
        compppsdf = argcomp_func (ppsdf,argpredf,feature)
        compppsdf.to_csv(outpath + 'MGTSEnter.csv', index=False)

    ### step 2. report potential outbreak or associated international isolates




### flag_reprot: to reported features for each isolate
def flag_reprot(flag_input, subdf,key_list):
    flag_input = flag_input.set_index('MGT_type')
    flag_dic_file = flag_input.to_dict()
    mgt_levels_list = subdf.columns.tolist()
    mgt_levels_list = [a for a in mgt_levels_list if "MGT" in a]
    #### to combine MGT level with ST/CCs into MGT1_ST3302
    for level in mgt_levels_list:
        if 'CC' not in level:
            subdf[level] = level + '_ST' + subdf[level].astype(str)
        if 'CC' in level:
            subdf[level] = level + subdf[level].astype(str)
    subdflist = subdf.values.tolist()
    accflag_dic = {}
    for k1 in flag_dic_file:
        if k1 in key_list:
            accflag_dic[k1] = {}
            dic_flag = {k:v for k, v in flag_dic_file[k1].items() if "nan" not in str(v)}
            mgtst_list = []
            for line in subdflist:
                for mgtst in line:
                    if mgtst in dic_flag.keys():
                        accflag_dic[k1][line[0]] = dic_flag[mgtst]
    return accflag_dic


### argcomp_func: if the precence of amr genes were screened, compare the args with MGT-MDR prediction
def argcomp_func(ppsdf, argpredf,feature):
    if 'Strain' in argpredf.columns and 'AMR_group' in argpredf.columns:
        # argpredf = pd.DataFrame(argpredi, columns= ['Strain', 'AMR_group'])
        argdic = argpredf.set_index('Strain').to_dict()['AMR_group']
        ppsdf['ARGsAMR_group'] = ppsdf['Strain'].map(argdic)

        colcondition = [(ppsdf[feature].isnull()) & (ppsdf['ARGsAMR_group'] != 'MDR'),
                        (ppsdf[feature].notnull()) & (ppsdf['ARGsAMR_group'] != 'MDR'),
                        (ppsdf[feature].notnull()) & (ppsdf['ARGsAMR_group'] == 'MDR'),
                        (ppsdf[feature].isnull()) & (ppsdf['ARGsAMR_group']=='MDR'),
                        ]
        values = ['True-','False+','True+','False-']
        ppsdf['arg_mgtmdr_comp'] = np.select(colcondition,values)
        arg_mgtmdr_compdic = ppsdf['arg_mgtmdr_comp'].value_counts().to_dict()
        print(arg_mgtmdr_compdic)
        proportionagree = arg_mgtmdr_compdic['True+']/(arg_mgtmdr_compdic['True+'] + arg_mgtmdr_compdic['False-'])
        # oddsratio, pvalue = stats.fisher_exact([[arg_mgtmdr_compdic['True-'],arg_mgtmdr_compdic['True+']],
        #                                         [arg_mgtmdr_compdic['False-'],arg_mgtmdr_compdic['False+']]
        #                                         ])
        print('% of agree possitive: ' + str('{:.1%}'.format(proportionagree)))
        # print('Odds Ratio: ' + str(oddsratio))
        # print('P value: ' + str(pvalue))
    return ppsdf



if __name__ == '__main__':
    main()





