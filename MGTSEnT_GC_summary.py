import glob
import pandas as pd
import scipy.stats as stats
import operator
import numpy as np
from time import sleep as sl
import argparse
from sklearn.metrics import pairwise_distances,pairwise_distances_chunked
from sklearn.cluster import AgglomerativeClustering,DBSCAN
import time
import os
from datetime import timedelta
import sys
from datetime import date


def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-m', '--mgtdb',
                        help='.csv file of the whole MGT typing dataset in MGTdb of Salmonella Enteritidis')
    parser.add_argument('-i','--isolatelist',help='.txt file of isolates of interest in MGTdb')
    parser.add_argument('-f', '--flag', help='.csv flag file included in the Github pakage, e.g. /srv/MGTSEnT/10-01-mgtmdrflag.csv')
    parser.add_argument('-o','--outpath_prefix',help='output path + prefix of the outfile, e.g. /srv/outdir/test_')
    args = parser.parse_args()
    return args


# ### Historical input for parseargs
# python MGTSEnT_GC_summary.py -i /mnt/e/2018/2019-06-14-Australia_SEN/Australia_isolates.txt -o /mnt/e/2018/2019-06-14-Australia_SEN/4_MGTSEnter/2_Austest/Australia_

# isolatelist = '/mnt/e/2018/2019-06-14-Australia_SEN/new_isolates.txt'
# outpath = '/mnt/e/2018/2019-06-14-Australia_SEN/2020-03-17-transmission/4_week_ODC5/6_ODC0_10'
# isolatelist = '/mnt/e/2018/2019-06-14-Australia_SEN/new_isolates.txt'  ## Australia_isolates
# outpath = '/mnt/e/2018/2019-06-14-Australia_SEN/4_MGTSEnter/2_Austest/GC_'


def main():
    args = parseargs()
    metapath = args.mgtdb
    isolatepath = args.isolatelist
    print(isolatepath)
    isolatelist = open(isolatepath, 'r').read().splitlines()
    mgt_flags = args.flag
    outpath = args.outpath_prefix
    ### inputlist
    # metapath = '/mnt/e/2018/2019-06-14-Australia_SEN/0_meta/MGTdataset_corre.csv'
    # metapath = '/mnt/e/2018/2019-06-14-Australia_SEN/2020-03-17-transmission/4_week_ODC5/6_ODC0_10/Australia_MGTSEN.csv'
    # mgt_flags = '/mnt/e/2018/2019-06-14-Australia_SEN/2020-03-03-historical_flag/03-04-info_dic_flag.csv'
    # typelistpath = 'pathofodc_list.txt'  ### if offer the type list in advance

    mgtlevelist = ['MGT9','ODC1','ODC2','ODC5','ODC10']  # ['MGT2','MGT3','MGT4','MGT5','MGT6','MGT7','MGT8','MGT9'] ### ['MGT1'] #'ODC2', 'ODC5',
    for mgt_level in mgtlevelist:
        print("*********" + mgt_level + "********")
        number_iso_threshold = 1
        # mgt_level = 'MGT2'
        datethresh = 28
        daytime = 'True'
        # if not os.path.exists(outpath):
        #     os.mkdir(outpath)

        t1 = time.time()
        metadf = pd.read_csv(metapath, low_memory=False)
        ### to get the odc types with No. of isolates >= number_iso_threshold
        odc_size2, typeofinterest_isodic = typeofinterest(metadf, mgt_level, isolatelist, number_iso_threshold)
        types = odc_size2
        print(
            'There are ' + str(len(odc_size2)) + ' ' + mgt_level + ' types >= 2 isolates each in the selected genomes')
        ### if types were offerred without screening
        # types = open(typelistpath,'r').read().splitlines()
        # types = ['20241']

        # ### to get the epidemiological information for the types of interest
        mgt_threshod, type_dic, isolate_dic, type_isolist_dic = transmi_tracking_flags(types, mgt_level,
                                                                                       number_iso_threshold, metadf,
                                                                                       mgt_flags)

        if typeofinterest_isodic:
            for k, v in type_dic.items():
                type_dic[k]['No. of seclected genomes'] = typeofinterest_isodic[k]

        ### to calculate the daytimespan if necessary
        if daytime == 'True':
            daytimespandic = daytimespan(types, metadf, mgt_level,datethresh)
            for k, v in type_dic.items():
                type_dic[k]['Daytimespan'] = daytimespandic[k]
            print(str((time.time() - t1) / 60) + ' min')

        ### to export the epidemiological information table with or without daytimespan information
        typedf = pd.DataFrame.from_dict(type_dic, orient='index')
        typedf.index.name = mgt_threshod
        typedf.drop(mgt_threshod, inplace=True, axis=1)
        isolate_df = pd.DataFrame.from_dict(isolate_dic, orient='index')
        isolate_df.index.name = 'Strain'

        if isolatelist:
            isolate_df['isolate_of_interest'] = np.where(isolate_df.index.isin(isolatelist), 'True', 'False')
            conditions = [(isolate_df['no_isolates'] < 2),
                          (isolate_df['no_isolates'] >= 2) & (isolate_df['no_isolates'] < 50),
                          (isolate_df['no_isolates'] >= 50)]
            values = ['Sporadic', 'Middle', 'Large']
            isolate_df['cluster_size'] = np.select(conditions, values)
            isolate_df = isolate_df.merge(metadf, left_on='Strain', right_on='Strain', how='left')
        typedflist = typedf.columns.tolist() + typedf.values.tolist()
        typedf.to_csv(outpath +  mgt_level + '_mgttype_transmission_link_' + str(number_iso_threshold) + '.csv')

        #### select columns to report for the tableau figures production
        # columnnames = ['Strain', 'Clades', 'internation_nation', mgt_level + '_>=_' + str(number_iso_threshold), 'no_isolates', 'No. of seclected genomes', 'isolate_of_interest', 'cluster_size']
        # isolate_df = pd.DataFrame(isolate_df, columns= columnnames)
        isolate_df = isolate_df.sort_values(by='Strain')
        isolate_df = isolate_df[isolate_df['isolate_of_interest'] == 'True']
        isolate_df.to_csv(
            outpath +  mgt_level + '_isolate_transmission_link_' + str(number_iso_threshold) + '.csv', index=False)
        print(str((time.time() - t1) / 60) + ' min')



def typeofinterest(metadf, odc, isolatelist,number_iso_threshold):
    # metadf = pd.read_csv(metapath, low_memory=False)
    isodic = {a :'True' for a in isolatelist}
    metadf['new_isolate'] = metadf['Strain'].map(isodic)
    metadf = metadf[metadf['new_isolate']=='True']
    odcdic = metadf.groupby([odc])['Strain'].count().to_dict()
    if 'None' in odcdic:
        odcdic.pop('None')

    if 0 in odcdic.keys():
        odcdic.pop(0)
    odcdicfilter = {a: b for a, b in odcdic.items() if b >= number_iso_threshold }
    odc_size2 = list(odcdicfilter.keys())
    return odc_size2, odcdic

def time_metric(a, b):
    match = 0
    missmatch = 0
    a = [int(x) for x in a]
    b = [int(x) for x in b]
    d0 = date(a[0], a[1], a[2])
    d1 = date(b[0], b[1], b[2])
    dayinterv = abs((d1 - d0).days)
    return dayinterv

def time_pw(metadf,odc,type):
    # metadf = pd.read_csv(metapath,low_memory=False)
    metadf[odc] = metadf[odc].astype(str)
    # metadf = metadf[(metadf[odc]==str(type)) | (metadf[odc]== str(type))]
    metadf = metadf[metadf[odc]==str(type)]
    no_isolates = metadf.shape[0]
    timedf = pd.DataFrame(metadf, columns=['Strain','Year','Month','Day'])
    # timedf = pd.read_csv('E:/2018/2019-06-14-Australia_SEN/test/time.csv')
    # # timedf['d'] = pd.to_datetime(timedf['Date'],format = '%Y/%m/%d')
    timedf['Day']=timedf['Day'].replace(np.nan,15)
    timedf = timedf[timedf['Month'].notnull()]
    # print({'time_input_dfsize': timedf.shape[0]})
    time_input_dfsize = timedf.shape[0]
    datedf=pd.DataFrame(timedf,columns=['Year','Month','Day'])
    acclist = timedf['Strain'].values.tolist()
    start_time = time.time()
    if datedf.shape[0]>=2:
        dayinterv = pairwise_distances(datedf, metric=time_metric, n_jobs=int(5)) #args.no_jobs=5
        # print(dayinterv)
        # pairw_outfrefix = 'E:/2018/2019-06-14-Australia_SEN/test/time_pairwise_'+ type + '_'
        # if  len(dayinterv) >=2 :
            # print("pairwise distance time", (" --- %s seconds ---" % (time.time() - start_time)))
            # np.savetxt(outfrefixout + odc +'_' +type + '_'+"time_pwdistances.txt", dayinterv.astype(int), fmt='%i', header=",".join(acclist), delimiter=",")
    if datedf.shape[0]<= 1:
        dayinterv = 'None'
    return time_input_dfsize, dayinterv, no_isolates

def daytimespan (odc_pairwise_list, metadf, mgt_level,datethresh):
    daytimespandic = {}
    for type in odc_pairwise_list:
        time_input_dfsize, dayinterv, no_isolates = time_pw(metadf, mgt_level, type)
        if time_input_dfsize < 2:
            if no_isolates == 1:
                daytimespandic[type] = 'Singleton'
            if no_isolates > 1:
                daytimespandic[type] = 'Not enough meta'
        if time_input_dfsize >= 2:
            part = dayinterv[np.triu_indices(dayinterv.shape[0], k=1)]
            if min(part) <= datethresh and max(part) > datethresh and max(part) > 0:
                # print(str(type) + '_Partially within ' + str(datethresh) +' days')
                daytimespandic[type] = 'Partially within ' + str(datethresh) +' days'
            if min(part) <= datethresh and max(part) <= datethresh:
                # print(str(type) + '_Within ' + str(datethresh) +' days')
                daytimespandic[type] = 'Within ' + str(datethresh) +' days'
            if min(part) > datethresh:
                # print(str(type) + '_Longer than ' + str(datethresh) +' days')
                daytimespandic[type] = 'Longer than ' + str(datethresh) +' days'
    return daytimespandic

def flag_reprot(flag_input, subdf,key_list):
    flag_input = flag_input.set_index('MGT_type')
    dic_file = flag_input.to_dict()
    mgt_levels_list = subdf.columns.tolist()
    mgt_levels_list = [a for a in mgt_levels_list if "MGT" in a]
    #### to combine MGT level with ST/CCs into MGT1_ST3302
    for level in mgt_levels_list:
        if 'CC' not in level:
            subdf.loc[:,level] = level + '_ST' + subdf[level].astype(str)
        if 'CC' in level:
            subdf.loc[:,level] = level + subdf[level].astype(str)
    subdflist = subdf.values.tolist()
    keyflag_dic = {}
    for k1 in dic_file:
        if k1 in key_list:
            dic_flag = {k:v for k, v in dic_file[k1].items() if "nan" not in str(v)}
            mgtst_list = []
            for line in subdflist:
                for value in line:
                    if value in dic_flag.keys():
                        mgtst = value
                        mgtst_list.append(dic_flag[mgtst])
            mgtst_ser = pd.Series(mgtst_list,dtype='str')
            keyflag_dic[k1] = "/".join(mgtst_ser.unique())
            # print("".join(mgtst_ser.unique()))
    return keyflag_dic

def transmi_tracking_flags(typelist,mgt,number_iso_threshold, metadf,mgt_flags):
    # metadf = pd.read_csv(metapath, low_memory=False)
    flags = pd.read_csv(mgt_flags)
    mgt_threshod = mgt + '_>=_' + str(number_iso_threshold)
    type_dic = {}
    type_isolist_dic = {}
    isolate_dic = {}

    for type in typelist:
        type_isolist_dic[type] = []
        subdf = metadf[metadf[mgt]== type]
        key_list = ['Clades','Lineages','MDR', 'AR2_1', 'Top_MGT4_STs']
        keyflag_dic = flag_reprot(flags, subdf, key_list)
        country_dic = subdf.groupby(['Country'])['Strain'].count().to_dict()
        state_dic = subdf.groupby(['State'])['Strain'].count().to_dict()
        year_dic = subdf.groupby(['Year'])['Strain'].count().to_dict()
        source_dic = subdf.groupby(['Type'])['Strain'].count().to_dict()

        type_dic[type] = keyflag_dic
        type_dic[type]['no_isolates'] = subdf.shape[0]
        type_dic[type]['country_detail'] = country_dic
        if 'None' in country_dic.keys():
            type_dic[type]['no_country'] = len(country_dic) - 1
        else:
            type_dic[type]['no_country'] = len(country_dic)

        if type_dic[type]['no_country'] == 0:
            type_dic[type]['internation_nation'] = 'Not enough meta'
        if type_dic[type]['no_country'] == 1 and int(type_dic[type]['no_isolates']) >= 2:
            type_dic[type]['internation_nation'] = 'National'
        if type_dic[type]['no_country'] == 1 and type_dic[type]['no_isolates'] == 1:
            type_dic[type]['internation_nation'] = 'Singleton'
        if type_dic[type]['no_country'] >= 2:
            type_dic[type]['internation_nation'] = 'International'

        type_dic[type]['state_detail'] = state_dic
        if 'None' in state_dic.keys():
            type_dic[type]['no_state'] = len(state_dic) - 1
        else:
            type_dic[type]['no_state'] = len(state_dic)
        type_dic[type]['year_detail'] = year_dic
        if 'None' in year_dic.keys():
            type_dic[type]['no_year'] = len(year_dic) - 1
        else:
            type_dic[type]['no_year'] = len(year_dic)

        type_dic[type]['source'] = source_dic
        if 'None' in source_dic.keys():
            type_dic[type]['no_source'] = len(source_dic) - 1
        else:
            type_dic[type]['no_source'] = len(source_dic)

        ### to product isolate_dic
        acclist = subdf['Strain'].tolist()
        for acc in acclist:
            type_isolist_dic[type].append(acc)
            isolate_dic[acc] = {}
            isolate_dic[acc] = type_dic[type]
            isolate_dic[acc][mgt_threshod] = type

    return mgt_threshod, type_dic, isolate_dic, type_isolist_dic



if __name__ == "__main__":
    main()
