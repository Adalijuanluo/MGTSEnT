import pandas as pd
import scipy.stats as stats
import operator
import numpy as np
from time import sleep as sl
import argparse
from sklearn.metrics import pairwise_distances,pairwise_distances_chunked
from sklearn.cluster import AgglomerativeClustering,DBSCAN
import time
from datetime import timedelta
import sys
from datetime import date
def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-e','--mgt_epi_data',help = 'csv file with MGT and epidemiological information')
    parser.add_argument('-m','--mgt_level',help='the level of mgt for clustering, e.g. MGT9 ')
    parser.add_argument('-n','--number_iso_threshold',help = 'the number of isolates threshold for each type in the tested MGT level')
    parser.add_argument('-o','--outpath',help='the path of outfiles, e.g. /srv/scratch')
    parser.add_argument('-f','--mgt_flags',help='the csv file of the flag MGT-STs for meta data')
    parser.add_argument('-p','--prefix',help='the prefix for the outfile names')
    parser.add_argument('-c','--country',help='country for filtering isolates for pairwise distance analysis, e.g. -c Australia')
    parser.add_argument('-t','--transtype',help='international or national cluster, or a list of ODC10 STs in a .txt file without heading, e.g. -t international')
    parser.add_argument('-a','--whole_allele_profile',help='Allele profile of MGT9, e.g. /srv/scratch/mgt9_alleprofile.txt')
    parser.add_argument("-d", "--distances",
                        help="file containing distances corresponding to the alleleprofiles file (from previous run of this script if applicable)")
    parser.add_argument("-l", "--dist_limits",
                        help="comma separated list of cluster cutoffs or range or both i.e 1,2,5 or 1-10 or 1,2,5-10, note ODC0/MGT9 were automatically given",default="1,2,5,10")
    parser.add_argument("-j", "--no_jobs",
                        help="num jobs to split distance calc into", default=1)
    args = parser.parse_args()
    return  args
def main():
    t1 = time.time()
    args = parseargs()
    mgt_epi = pd.read_csv(args.mgt_epi_data)
    threshod_dic = {args.mgt_level: int(args.number_iso_threshold)}
    if args.mgt_flags == None:
        mgt_threshod, type_dic, isolate_dic, type_isolist_dic = transmi_tracking(threshod_dic, mgt_epi)
    if args.mgt_flags != None:
        flags = pd.read_csv(args.mgt_flags)
        mgt_threshod, type_dic, isolate_dic, type_isolist_dic = transmi_tracking_flags(threshod_dic, mgt_epi, flags)
    typedf = pd.DataFrame.from_dict(type_dic, orient='index')
    typedf.index.name = mgt_threshod
    isolate_df = pd.DataFrame.from_dict(isolate_dic, orient='index')
    isolate_df.index.name = 'Accession'
    typedflist = typedf.columns.tolist() + typedf.values.tolist()
    isolate_df.to_csv(args.outpath + '/' + args.prefix + args.mgt_level + '_isolate_transmission_link.csv')
    typedf.to_csv(args.outpath + '/' + args.prefix + args.mgt_level + '_mgttype_transmission_link.csv')

    ###### allele profile getting for pairwise calculation,
    if args.country != None:
        odc_pairwise_list = []
        if args.transtype == 'international':
            typedf2 = typedf[typedf['no_country'] >= 2] ### >=2 is international
            type_country_dic = typedf2['country_detail'].to_dict()

            for type, subdic in type_country_dic.items():
                for c in subdic.keys():
                    if args.country == c :
                        no_iso_country = subdic[c]
                        if no_iso_country>= 1 and args.mgt_level == "ODC10":  ####### to set the threshold of >= 1 for each cluster in Australia.
                            if type not in odc_pairwise_list and type != "None":
                                odc_pairwise_list.append(type)
                        if args.mgt_level != "ODC10":
                            if type not in odc_pairwise_list and type != "None":
                                odc_pairwise_list.append(type)
            print({"Total No. of types for pairwise distance calculation" : len(odc_pairwise_list)})
            print(odc_pairwise_list)
        if args.transtype == 'national':
            typedf2 = typedf[typedf['no_country'] == 1] ### >=2 is international; == 1 is national
            type_country_dic = typedf2['country_detail'].to_dict()

            for type, subdic in type_country_dic.items():
                for c in subdic.keys():
                    if args.country == c :
                        no_iso_country = subdic[c]
                        if no_iso_country>= 2 and args.mgt_level == "ODC10":  ####### >=2 isolates for national transmission
                            if type not in odc_pairwise_list and type != "None":
                                odc_pairwise_list.append(type)
                        if args.mgt_level != "ODC10":
                            if type not in odc_pairwise_list and type != "None":
                                odc_pairwise_list.append(type)
            print({"Total No. of types for pairwise distance calculation" : len(odc_pairwise_list)})
            print(odc_pairwise_list)
        if args.transtype == None:
            typedf2 = typedf[typedf['no_country'] >= 1] ### including both international and national
            type_country_dic = typedf2['country_detail'].to_dict()
            for type, subdic in type_country_dic.items():
                for c in subdic.keys():
                    if args.country == c :
                        no_iso_country = subdic[c]
                        if no_iso_country>= 2 and args.mgt_level == "ODC10":  ####### >=2 isolates for national transmission
                            if type not in odc_pairwise_list and type != "None":
                                odc_pairwise_list.append(type)
                        if args.mgt_level != "ODC10":
                            if type not in odc_pairwise_list and type != "None":
                                odc_pairwise_list.append(type)
            print({"Total No. of types for pairwise distance calculation" : len(odc_pairwise_list)})
            print(odc_pairwise_list)
        if args.transtype != None and ".txt" in args.transtype:
            odc_pairwise_list=open(args.transtype,'r').read().splitlines()
        # odc_pairwise_list = ['4969']
        for type in odc_pairwise_list :
            if type in type_isolist_dic:
                print(args.mgt_level + '_' + type)
                # time_pw(args, mgt_epi, args.mgt_level, type, args.outpath)
                ### to save the type correlated acc list
                isolatelistfile = open(args.outpath + '/' + args.mgt_level + '_' + type + '_' + args.country + '_correlated_isolatelist.txt','w')
                isolatelistfile.write(args.mgt_level + '_' + type + '\n')
                for acc in type_isolist_dic[type]:
                    isolatelistfile.write(acc + '\n')
                ### to calculate the pairwise distance of alleles
                if args.whole_allele_profile != "":
                    allele_prof = open(args.whole_allele_profile, "r").read().splitlines()
                    allele_proflist = get_part_alleleprofil(allele_prof, type_isolist_dic[type])
                    allele_prof_outfile = open(args.outpath + '/' + args.mgt_level + '_' + type + '_alleleprof.txt', 'w')
                    allele_prof_outlist=[]
                    for a in allele_proflist:
                        allele_prof_outlist.append(a)
                        allele_prof_outfile.write(a + '\n')
                    profs, id_to_strain = process_profiles(allele_prof_outlist)
                    pairw_outfrefix=args.outpath + '/' + args.mgt_level + '_' + type + '_'
                    pairwise_process(args, profs, id_to_strain, pairw_outfrefix)
                    t2 = timecal(t1)
                    ##  iso_process_profiles() and  iso_pairwise_process() are for pairwise distance of isolates in Australia
                    # iso_profs, iso_id_to_strain = iso_process_profiles(allele_prof_outlist)
                    # print(iso_id_to_strain)
                    # iso_pairwise_process(args, iso_profs, iso_id_to_strain, pairw_outfrefix)
def timecal(uptime):
    timespan = time.time() - uptime
    print(timedelta(seconds=timespan))
    return time.time()
def time_metric(a, b):
    match = 0
    missmatch = 0
    a = [int(x) for x in a]
    b = [int(x) for x in b]
    d0 = date(a[0], a[1], a[2])
    d1 = date(b[0], b[1], b[2])
    dayinterv = abs((d1 - d0).days)
    return dayinterv
###note columns have to include 'Accession','Collection Year','Collection Month','Collection Day'.
def time_pw(args, metadf,odc,type,outfrefixout):
    # metadf = pd.read_csv(metapath)
    metadf[odc] = metadf[odc].astype(str)
    metadf = metadf[(metadf[odc]==type) | (metadf[odc]== str(type))]
    timedf = pd.DataFrame(metadf, columns=['Accession','Collection Year','Collection Month','Collection Day'])
    # timedf = pd.read_csv('E:/2018/2019-06-14-Australia_SEN/test/time.csv')
    # # timedf['d'] = pd.to_datetime(timedf['Date'],format = '%Y/%m/%d')
    timedf['Collection Day']=timedf['Collection Day'].replace(np.nan,15)
    timedf = timedf[timedf['Collection Month'].notnull()]
    print({'time_input_dfsize': timedf.shape[0]})
    datedf=pd.DataFrame(timedf,columns=['Collection Year','Collection Month','Collection Day'])
    acclist = timedf['Accession'].values.tolist()
    start_time = time.time()
    if datedf.shape[0]>1:
        dayinterv = pairwise_distances(datedf, metric=time_metric, n_jobs=int(args.no_jobs))
        # pairw_outfrefix = 'E:/2018/2019-06-14-Australia_SEN/test/time_pairwise_'+ type + '_'
        if  len(dayinterv) >=2 :
            print("pairwise distance time", (" --- %s seconds ---" % (time.time() - start_time)))
            np.savetxt(outfrefixout + odc +'_' +type + '_'+"time_pwdistances.txt", dayinterv.astype(int), fmt='%i', header=",".join(acclist), delimiter=",")
def unneg(a):
    if "-" in a:
        return a.split("_")[0][1:]
    else:
        return a
def mgt_dist_metric(a, b):
    match = 0
    missmatch = 0

    for i in range(len(a)):
        aAllele = a[i]
        bAllele = b[i]
        # print(aAllele,bAllele)
        # sl(0.1)
        if aAllele == 0 or bAllele == 0 or aAllele == bAllele:
            match += 1
        else:
            missmatch += 1
            # print(aAllele,bAllele)

    return missmatch
def process_profiles(inprofiles, s=False):
    profs = {}
    id_to_strain = {}

    for line in inprofiles[1:]:
        col = line.split("\t")
        if s:
            if col[0] in s:
                # print(col[0])

                if col[1] not in profs:
                    noneg = [unneg(x) for x in col[3:]]
                    profs[col[1]] = noneg
                    id_to_strain[col[1]] = [str(col[0])]
                else:
                    id_to_strain[col[1]].append(str(col[0]))
        else:
            # print(col[0])
            if col[1] not in profs:
                noneg = [unneg(x) for x in col[3:]]
                profs[col[1]] = noneg
                id_to_strain[col[1]] = [str(col[0])]
            else:
                id_to_strain[col[1]].append(str(col[0]))

    return profs, id_to_strain
def pairwise_process(args,profs,id_to_strain, pairw_outfrefix):
    idlist = list(profs.keys())
    inprofs = [profs[x] for x in idlist]

    dfo = pd.DataFrame(inprofs)
    # distances only calculated if args.distances not set
    lend = ""
    if args.distances:
        # read in distances previosly calculated
        d = np.loadtxt(args.distances)
        lend = len(d) ### number of MGT9 STs in this cluster
    else:
        start_time = time.time()
        d = pairwise_distances(inprofs, metric=mgt_dist_metric, n_jobs=int(args.no_jobs))
        lend = len(d)
        # if  len(d) >=2 :
        print("pairwise distance time", (" --- %s seconds ---" % (time.time() - start_time)))
        np.savetxt(pairw_outfrefix + "mgt9_distances.txt", d.astype(int), fmt='%i', header=",".join(idlist), delimiter=",")
    # distance cutoffs to calculate
    if lend >=2:
        pairw_outfile = open(pairw_outfrefix + 'iso_odc_recal.txt','w')
        diststring = args.dist_limits
        dists = diststring.split(",")
        distances = []
        for i in dists:
            if "-" in i:
                n = i.split("-")
                nlist = list(range(int(n[0]) + 1, int(n[1]) + 2))
                # distance cutoffs seems to be non inclusive i.e. cutoff of 3 means max distance is 2
                # therefore need to add 1 to all values
            else:
                nlist = [int(i) + 1]
            distances += nlist
        clusterlists = {}
        preference = []
        for id in idlist:
            preference.append(len(id_to_strain[id]))
        start_time = time.time()
        for dist in distances:
            clusters = AgglomerativeClustering(n_clusters=None, distance_threshold=dist, affinity="precomputed",
                                               linkage="single").fit_predict(d)
            clusterls = list(clusters)
            clusterlists[dist] = clusterls
        print("clustering time", (" --- %s seconds ---" % (time.time() - start_time)))
        realdists = ["ODC" + str(x - 1) for x in distances]

        pairw_outfile.write("Strain\tMGT9\t{}\n".format("\t".join(realdists)))
        for i in range(len(idlist)):
            id = idlist[i]
            for strain in id_to_strain[id]:
                pairw_outfile.write(strain + '\t' + str(id))
                for d in distances:
                    clust = clusterlists[d][i]
                    pairw_outfile.write("\t" + str(clust + 1))
                pairw_outfile.write("\n")
        pairw_outfile.close()
    if lend < 2:  ### belong to the same MGT9 ST
        pairw_outfile = open(pairw_outfrefix + 'iso_odc_recal.txt','w')
        pairw_outfile.write('Strain' + '\t' + 'MGT9' + '\n')
        for st, isolist in id_to_strain.items():
            for iso in isolist:
                pairw_outfile.write(str(iso) + '\t' + str(st) + '\n')
        return
        ##### pairwise distance calculation
###  iso_process_profiles() and  iso_pairwise_process() are for pairwise distance of isolates in Australia
def iso_process_profiles(inprofiles, s=False):
    profs = {}
    id_to_strain = {}
    for line in inprofiles[1:]:
        col = line.split("\t")
        if s:
            if col[0] in s:
                # print(col[0])
                if col[0] not in profs:
                    noneg = [unneg(x) for x in col[3:]]
                    profs[col[0]] = noneg
                    id_to_strain[col[0]] = [str(col[1])]
                else:
                    id_to_strain[col[0]].append(str(col[1]))
        else:
            # print(col[0])

            if col[0] not in profs:
                noneg = [unneg(x) for x in col[3:]]
                profs[col[0]] = noneg
                id_to_strain[col[0]] = [str(col[1])]
            else:
                id_to_strain[col[0]].append(str(col[1]))
    return profs, id_to_strain
def iso_pairwise_process(args,profs,id_to_strain, pairw_outfrefix):
    idlist = list(profs.keys())
    # print(idlist)
    inprofs = [profs[x] for x in idlist]
    # distances only calculated if args.distances not set
    lend = ""
    if args.distances:
        # read in distances previosly calculated
        d = np.loadtxt(args.distances)
        lend = len(d)
    else:
        start_time = time.time()
        d = pairwise_distances(inprofs, metric=mgt_dist_metric, n_jobs=int(args.no_jobs))
        lend = len(d)
        if  len(d) >=2 :
            print("pairwise distance time", (" --- %s seconds ---" % (time.time() - start_time)))
            np.savetxt(pairw_outfrefix + "iso_distances.txt", d.astype(int), fmt='%i', header=",".join(idlist), delimiter=",")
def epi_filt(mgt_epi,fi_dic):
    col_name = mgt_epi.columns.values.tolist()
    mgt_epi = mgt_epi.values.tolist()
    for key in fi_dic:
        epi_filter_out = []
        for line in mgt_epi:
            line_index = col_name.index(key)
            line_value = str(line[line_index])
            if line_value in fi_dic[key]:
                epi_filter_out.append(line)
        mgt_epi = epi_filter_out
    mgt_epi = pd.DataFrame(mgt_epi)
    mgt_epi.columns = col_name
    print(mgt_epi.shape)
    return mgt_epi
def flag_reprot(flag_input, test_file,key_list):
    flag_input = flag_input.set_index('MGT_type')
    dic_file = flag_input.to_dict()
    mgt_levels_list = test_file.columns.tolist()
    mgt_levels_list = [a for a in mgt_levels_list if "MGT" in a]
    for level in mgt_levels_list:
        test_file[level] = level + test_file[level].astype(str)
    test_list = test_file.values.tolist()
    keyflag_dic = {}
    for k1 in dic_file:
        if k1 in key_list:
            # outfile = open(outpath + '/' + k1 + '.txt','w')
            dic_file2 = {k:v for k, v in dic_file[k1].items() if "nan" not in str(v)}
            mgtst_list = []
            for line in test_list:
                for value in line:
                    if value in dic_file2.keys():
                        mgtst = value
                        mgtst_list.append(mgtst)
                        strain_name = line [0]
                        predict_types = dic_file2[value]
                        # output = "{}\t{}\t{}".format(strain_name,predict_types,mgtst)
                        # print(output)
                        # outfile.write(output + '\n')
            mgtst_ser = pd.Series(mgtst_list)
            keyflag_dic[k1] = mgtst_ser.value_counts().to_dict()
    return keyflag_dic
def transmi_tracking_flags(threshod_dic, mgt_epi,flags):
    for mgt in threshod_dic:
        gp = mgt_epi.groupby([ mgt])['Strain'].count().fillna(0)
        pass_filter_type_dic = gp [gp>= threshod_dic[mgt]].to_dict()
        if 0 in pass_filter_type_dic.keys():
            pass_filter_type_dic.pop(0)
        mgt_threshod = mgt + '_>=_' + str(threshod_dic[mgt])
        type_dic = {}
        type_isolist_dic = {}
        isolate_dic = {}
        interspread = 0
        limited_year = 0
        large_sclale = 0
        for type in pass_filter_type_dic.keys():
            type_isolist_dic[type] = []
            subdf = mgt_epi[mgt_epi[mgt]== type]
            key_list = ['Population_Structure', 'MDR', 'AR2_1', 'Top_MGT4_STs']
            keyflag_dic = flag_reprot(flags, subdf, key_list)
            country_dic = subdf.groupby(['Country'])['Strain'].count().to_dict()
            year_dic = subdf.groupby(['Collection Year'])['Strain'].count().to_dict()
            source_dic = subdf.groupby(['Source Type'])['Strain'].count().to_dict()
            type_dic[type] = keyflag_dic

            # type_dic[type]={"no_isolates":{}}
            type_dic[type]['no_isolates'] = pass_filter_type_dic[type]
            type_dic[type]['country_detail'] = country_dic
            if 'None' in country_dic.keys():
                type_dic[type]['no_country'] = len(country_dic) - 1
            else:
                type_dic[type]['no_country'] = len(country_dic)

            type_dic[type]['year_detail'] = year_dic
            if 'None' in year_dic.keys():
                type_dic[type]['no_year'] = len(year_dic) - 1
            else:
                type_dic[type]['no_year'] = len(year_dic)
            if len(year_dic) <= 2 and  len(year_dic)> 0:
                limited_year =limited_year+1
            if len(country_dic) >= 2 :
                interspread = interspread + 1
            if len(year_dic) > 1 and len(year_dic) > 0 and len(country_dic) >= 2 and pass_filter_type_dic[type] > 50:
                large_sclale = large_sclale + 1

            type_dic[type]['source'] = source_dic
            if 'None' in source_dic.keys():
                type_dic[type]['no_source'] = len(source_dic) - 1
            else:
                type_dic[type]['no_source'] = len(source_dic)
            ########### to product isolate_dic
            acclist = subdf['Accession'].tolist()
            for acc in acclist:
                type_isolist_dic[type].append(acc)
                isolate_dic[acc] = {}
                isolate_dic[acc] = type_dic[type]
                isolate_dic[acc][mgt_threshod] = type

        print("No. of passed types: " + str(len(pass_filter_type_dic)))
        print('No. of potential international spreading clusters: ' + str(interspread))
        print('No. of potential international spreading clusters within years: ' + str(limited_year))
        print('No. of potential international spreading large clusters within years >50: ' + str(large_sclale))

        return mgt_threshod,type_dic, isolate_dic, type_isolist_dic
def transmi_tracking(threshod_dic, mgt_epi):
    for mgt in threshod_dic:
        gp = mgt_epi.groupby([mgt])['Strain'].count().fillna(0)
        pass_filter_type_dic = gp [gp>= threshod_dic[mgt]].to_dict()
        if 0 in pass_filter_type_dic.keys():
            pass_filter_type_dic.pop(0)
        mgt_threshod = mgt + '_>=_' + str(threshod_dic[mgt])
        type_dic = {}
        isolate_dic = {}
        type_isolist_dic ={}
        interspread = 0
        limited_year = 0
        large_sclale = 0
        for type in pass_filter_type_dic.keys():
            type_isolist_dic[type]=[]
            subdf = mgt_epi[mgt_epi[mgt]== type]
            key_list = ['Population_Structure', 'MDR', 'AR2_1', 'Top_MGT4_STs']
            # keyflag_dic = flag_reprot(flags, subdf, key_list)
            country_dic = subdf.groupby(['Country'])['Strain'].count().to_dict()
            year_dic = subdf.groupby(['Collection Year'])['Strain'].count().to_dict()
            source_dic = subdf.groupby(['Source Type'])['Strain'].count().to_dict()
            type_dic[type] = {}
            type_dic[type]['no_isolates'] = pass_filter_type_dic[type]
            type_dic[type]['country_detail'] = country_dic
            if 'None' in country_dic.keys():
                type_dic[type]['no_country'] = len(country_dic) - 1
            else:
                type_dic[type]['no_country'] = len(country_dic)

            type_dic[type]['year_detail'] = year_dic
            if 'None' in year_dic.keys():
                type_dic[type]['no_year'] = len(year_dic) - 1
            else:
                type_dic[type]['no_year'] = len(year_dic)
            if len(year_dic) <= 2 and  len(year_dic)> 0:
                limited_year =limited_year+1
            if len(country_dic) >= 2 :
                interspread = interspread + 1
            if len(year_dic) > 1 and len(year_dic) > 0 and len(country_dic) >= 2 and pass_filter_type_dic[type] > 50:
                large_sclale = large_sclale + 1
            type_dic[type]['source'] = source_dic
            if 'None' in source_dic.keys():
                type_dic[type]['no_source'] = len(source_dic) - 1
            else:
                type_dic[type]['no_source'] = len(source_dic)

            ########### to product isolate_dic
            acclist = subdf['Accession'].tolist()
            for acc in acclist:
                type_isolist_dic[type].append(acc)
                isolate_dic[acc] = {}
                isolate_dic[acc] = type_dic[type]
                isolate_dic[acc][mgt_threshod] = type
        print("No. of passed types: " + str(len(pass_filter_type_dic)))
        print('No. of potential international spreading clusters: ' + str(interspread))
        print('No. of potential international spreading clusters within years: ' + str(limited_year))
        print('No. of potential international spreading large clusters within years >50: '+ str(large_sclale))
        return mgt_threshod,type_dic, isolate_dic, type_isolist_dic
def get_part_alleleprofil(whol_alleprof, isolist):
    outlist = []
    outlist.append(whol_alleprof[0])
    for line in whol_alleprof:
        col = line.split('\t')
        for acc in isolist:
            if acc == col[0]:           # or acc + '_cladeC.fasta' == col[0]:
                outlist.append(line)
    return outlist
if __name__ == "__main__":
    main()

