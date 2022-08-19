# MGTSEnT_Microreactinput.py can generate the nwk file and the metadata as the input for the microreact to visualise the ODC10.
# It required for the recalled GC0 to 10, based on the A2020_05_06_prog_transmission_type_iso.py

from geopy.exc import GeocoderTimedOut
from geopy.geocoders import Nominatim
import matplotlib
import pandas as pd
from sklearn.metrics import pairwise_distances,pairwise_distances_chunked
import time
import numpy as np
######### dendrogram plotting
from numpy import loadtxt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy
from sklearn.cluster import AgglomerativeClustering, DBSCAN
import argparse


def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-m', '--mgtdb',
                        help='.csv file of the whole MGT typing dataset in MGTdb of Salmonella Enteritidis')
    parser.add_argument('-c','--cluster_type',help='.txt file of GC types, e.g. GC10_1813')
    parser.add_argument('-o','--outpath_prefix',help='output path + prefix of the outfile, e.g. /srv/outdir/test_')
    args = parser.parse_args()
    return args

# ### Historical input for parseargs
# python MGTSEnT_Microreactinput.py -c /mnt/e/2018/2019-06-14-Australia_SEN/4_MGTSEnter/GCtypes.txt -o /mnt/e/2018/2019-06-14-Australia_SEN/4_MGTSEnter/2_Austest/Australia_
# isolatelist = '/mnt/e/2018/2019-06-14-Australia_SEN/new_isolates.txt'
# outpath = '/mnt/e/2018/2019-06-14-Australia_SEN/2020-03-17-transmission/4_week_ODC5/6_ODC0_10'
# isolatelist = '/mnt/e/2018/2019-06-14-Australia_SEN/new_isolates.txt'  ## Australia_isolates
# outpath = '/mnt/e/2018/2019-06-14-Australia_SEN/4_MGTSEnter/2_Austest/GC_'


### inputs
# mgt_epi='/mnt/e/2018/2019-06-14-Australia_SEN/2020-03-17-transmission/05-05-pairwise_transmission/mgt_epi_06_17_replace5.csv'
args = parseargs()
gctypelist = open(args.cluster_type,'r').read().splitlines()
outpath = args.outpath_prefix
metapath = args.mgtdb
#metapath = '/mnt/e/2018/2019-06-14-Australia_SEN/0_meta/MGTdataset_corre.csv'
### ocd recall from GC0 to GC10, if add an extra column of ODC0.0, MGT9 can be separated in the dendrogram
#odcrecallpath = '/mnt/e/2018/2019-06-14-Australia_SEN/2020-03-17-transmission/05-05-pairwise_transmission/06-21_pwd_ODC10'
odcrecallpath = "/srv/scratch/z5167454/MGTSEnT/updating_input"
# outpath = odcrecallpath + '/recall_analy/mciroreact'
# outpath = '/mnt/e/2018/2019-06-14-Australia_SEN/2020-03-17-transmission/05-05-pairwise_transmission/2_ODC10_1813'
# gctypelist = ['ODC10_1813',"ODC10_36"] # '36','89','95','1893','1889'
core_country = 'Australia'
longitude = []
latitude = []

### to run the functions to produce input files for microreact
def main():
    for gctype in gctypelist:
        gclevel = gctype.split('_')[0]
        type = gctype.split('_')[1]
        odcrecall = odcrecallpath + '/' + gctype + '_iso_odc_recal.txt'

        ### step one: get the isolate distance based on ODC0 to ODC10
        allele_prof_outlist = open(odcrecall, 'r').read().splitlines()
        # pairw_outfrefix=args.outpath + args.mgt_level + '_' + type + '_'
        pairw_outfrefix = outpath + gctype + '_'
        profs, id_to_strain = iso_process_profiles(allele_prof_outlist)
        iso_pairwise_process(profs, id_to_strain, pairw_outfrefix)

        ### step two: produce dentrogram and the nwk file as input for the microreact
        odcisodistance = pairw_outfrefix + "iso_distances.txt"
        nwkfile = outpath + gctype + '_microreact2.nwk'
        lines = loadtxt(odcisodistance, comments="#", delimiter=",", unpack=False)
        col = open(odcisodistance, 'r').read().splitlines()[0][2:].split(',')
        dists = squareform(lines)

        # from numpy import savetxt
        # file1 = outpath + gctype + '_dist.csv'
        # savetxt(file1, dists,delimiter=',')
        ### ‘complete’, ‘average’, ‘weighted’ and ‘ward’,
        linkage_matrix = linkage(dists, 'average')
        print(linkage_matrix)
        # file2 = outpath + gctype + '_dist_matrix.csv'
        # savetxt(file2, linkage_matrix, delimiter=',')
        # linkage_matrix = linkage(lines, "single")
        mkdendrogram(linkage_matrix, outpath,col)
        ### export the dendrogram into newick file
        tree = hierarchy.to_tree(linkage_matrix, False)
        tree2 = getNewick(tree, "", tree.dist, col)
        with open(nwkfile, 'w') as outf:
            outf.write(tree2)
            print(tree2)

        ### step three: produce the microreact metadata .csv file including the nodes colour, geographic location, date and ODC0-10
        acclist,odcdf = cluster_epi_des(odcrecall)
        submgtepi, type_dic = epides(acclist, metapath)
        submgtepi = submgtepi.merge(odcdf, left_on='Strain',right_on='Strain', how = 'left')
        submgtepi =pd.DataFrame(submgtepi, columns=['Strain','Date', 'Day', 'Year', 'Month','Source','Continent', 'Country','State','MGT9_y','ODC5_x', 'ODC1_y',
               'ODC2_y', 'ODC3', 'ODC4', 'ODC5_y','ODC6','ODC7','ODC8','ODC9','ODC10_x'])
        submgtepi = submgtepi.rename(columns={'Strain': 'id'}).fillna('None')
        submgtepi = date_epidf(submgtepi)
        submgtepi = region_epidf(submgtepi)
        submgtepi = source_epidf(submgtepi)
        submgtepi.to_csv(outpath + gclevel + '_' + str(type) + '_microreact.csv', index=False)



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
def iso_pairwise_process(profs,id_to_strain,pairw_outfrefix ):
    idlist = list(profs.keys())
    # print(idlist)
    inprofs = [profs[x] for x in idlist]
    # distances only calculated if args.distances not set
    lend = ""
    start_time = time.time()
    d = pairwise_distances(inprofs, metric=mgt_dist_metric, n_jobs=int(1))
    lend = len(d)
    if  len(d) >=2 :
        print("pairwise distance time", (" --- %s seconds ---" % (time.time() - start_time)))
        np.savetxt(pairw_outfrefix + "iso_distances.txt", d.astype(int), fmt='%i', header=",".join(idlist), delimiter=",")
        print(d)
def mkdendrogram(linkage_matrix,outpath,col):
    R = dendrogram(linkage_matrix, orientation='top',labels=col,color_threshold=4.1 , count_sort = True) #color_threshold=15 ,
    outdic = {}
    for k, v in R.items():
        print(k)
    outdic['ivl'] = R['ivl']
    outdic['leaves_color_list'] = R['leaves_color_list']
    redf = pd.DataFrame.from_dict(outdic)
    redf = redf.reset_index()
    # redf.to_csv(outpath + 'denrogram.csv', index=False)
    plt.tick_params(axis='x', which='major', labelsize=5)
    plt.tick_params(axis='y', which='major', labelsize=8)
    plt.title("ODC10_45")
    plt.show()

def getNewick(node, newick, parentdist, leaf_names):
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = getNewick(node.get_left(), newick, node.dist, leaf_names)
        newick = getNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick



### get the latitude and longitude of a region, however, it may fail when the internet is not good
def findGeocode(region, attempt=1, max_attempts=10):
    print(region)
    try:
        geolocator = Nominatim(user_agent="your_app_name")
        return geolocator.geocode(region)
    except GeocoderTimedOut:
        if attempt <= max_attempts:
            return findGeocode(region, attempt=attempt + 1)
        raise


# print(geolocator.geocode(city))

def epides (acclist, mgt_epi):
    mgt_epi = pd.read_csv(mgt_epi, low_memory=False)
    clusdf = pd.DataFrame(acclist, columns=['Strain'])
    submgtepi = clusdf.merge(mgt_epi, left_on='Strain', right_on='Strain', how='left')
    # submgtepi = submgtepi.set_index('Strain')
    type_dic ={}
    type_dic['No. of isolates'] = len(acclist)
    ######### date processing
    timedf = pd.DataFrame(submgtepi)
    timedf['Day'] = timedf['Day'].replace(np.nan, 15)
    timedf = timedf[timedf['Month'].notnull()]
    timedf = pd.DataFrame(timedf, columns=["Year", "Month", "Day"])
    timedf.columns = ['year', 'month', 'day']

    timedf['date'] = pd.to_datetime(timedf[['year', 'month', 'day']])
    datespan = timedf['date'].max() - timedf['date'].min()
    type_dic['No. of isolates with date'] = timedf.shape[0]
    type_dic['Date_span'] = datespan.days
    #### meta processing
    # type_dic = dicmake('Country', submgtepi, type_dic)
    # type_dic = dicmake('Region', submgtepi, type_dic)
    # type_dic = dicmake('Source Type', submgtepi, type_dic)
    # type_dic = dicmake('Year', submgtepi, type_dic)
    return submgtepi, type_dic
def cluster_epi_des (odccall):
    odcdf= pd.read_table(odccall)
    colname = odcdf.columns.tolist()
    acclist = odcdf['Strain'].tolist()
    # colname[0] = 'Strain'
    # odcdf.columns = colname
    # odcdf.index = odcdf['Strain']
    # accdic = {}
    # for group in colname[1:]:
    #     # accdic[group] = {}
    #     var = odcdf[group].unique()
    #     # print(len(var))
    #     for v in var:
    #         odcclus = str(group)  + '_' + str(v)
    #         acc = odcdf[odcdf[group] == v].index.tolist()
    #         accdic[odcclus] = acc
    return acclist,odcdf

### get the latitude and longitude of a region
def regionloc(regionlist):
    regionlocdic = {}
    for region in regionlist:
        regionlocdic[region] = {}
        if findGeocode(region) != None:
            loc = findGeocode(region)
            regionlocdic[region]['latitude'] = loc.latitude
            regionlocdic[region]['longitude'] = loc.longitude
            out = ("{}\t{}\t{}").format(region, loc.latitude, loc.longitude)
        if findGeocode(region) == None:
            print(region + ' not found')
            regionlocdic[region]['latitude'] = 'Unknown'
            regionlocdic[region]['longitude'] = 'Unknown'
    regionlocdic['None'] = {}
    regionlocdic['None']['latitude'] = ""
    regionlocdic['None']['longitude'] = ""
    return regionlocdic

# regionlist = ["United Kingdom","South Africa","United States","Australia","Ireland, North Atlantic","Uganda","Belgium, Africa","Denmark","Germany"]
# regionloc(regionlist)
def colordic(colourtypelist):
    a = 1
    color_dic = {}
    c = 7  ######### decide the difference or color span distance
    hex_color = list(matplotlib.colors.cnames.values())
    for k in colourtypelist:
        a = a + c
        color_dic[k] = hex_color[a + c]  ### to set the color
    color_dic['None'] = ''
    return color_dic

def region_epidf(submgtepi):
    regionlist = pd.unique(submgtepi['Country'])
    if len(regionlist) > 1:
        regionlocdic = regionloc(regionlist)
        countrycolor = colordic(regionlist)
        submgtepi['Latitude'] = submgtepi['Country'].map(lambda x: regionlocdic[x]['latitude'])
        submgtepi['Longitude'] = submgtepi['Country'].map(lambda x: regionlocdic[x]['longitude'])
        submgtepi['Country__Colour'] = submgtepi['Country'].map(countrycolor)
        submgtepi = submgtepi.replace('None', "")

    if len(regionlist) == 1 and core_country in regionlist:
        regionlist = list(submgtepi['State'].value_counts().to_dict().keys())
        regionlocdic = regionloc(regionlist)
        countrycolor = colordic(regionlist)
        submgtepi['Latitude'] = submgtepi['State'].map(lambda x: regionlocdic[x]['latitude'])
        submgtepi['Longitude'] = submgtepi['State'].map(lambda x: regionlocdic[x]['longitude'])
        submgtepi['State__Colour'] = submgtepi['State'].map(countrycolor)
        submgtepi = submgtepi.replace('None', "")
        # submgtepi.to_csv(outpath + gclevel + '_' + str(type) + '_microreact.csv', index=False)
    return submgtepi


def source_epidf(submgtepi):
    sourcelist = pd.unique(submgtepi['Source'])
    invasivetype = []
    invasavelist = ['blood', 'invasive','fluid','CSF','Wound']
    for sourcetype in sourcelist:
        if 'Human' in sourcetype:
            for a in invasavelist:
                if a in sourcetype:
                    invasivetype.append(sourcetype)
    if len(invasivetype) > 0:
        ### to create a column if
        submgtepi["Invasive"] = np.where(submgtepi["Source"].isin(invasivetype), "TRUE", "")
    return submgtepi



def date_epidf(submgtepi):

    submgtepi["Date2"] = submgtepi.apply(lambda row: str(row["Month"]) + '/' + str(row["Year"]) , axis =1)
    submgtepilist  = submgtepi.values.tolist()
    collist = submgtepi.columns.tolist()
    val  = collist.index('Date')
    val2 = collist.index('Date2')
    for line in submgtepilist:
        if line[val] == 'None':
            line[val] = line[val2]
    submgtepi2 = pd.DataFrame(submgtepilist, columns=collist)

    submgtepi2= submgtepi2.drop(['Date2'], axis =1 )

    # submgtepi['Date'] = np.where((submgtepi['Date']=='None'),  row('Date'), else row('Date'), axis =1)
    # lambda row: row.Actual_Price - ((row.Discount_Percentage / 100) * row.Actual_Price), axis = 1
    # submgtepi['Event'] = np.where((df.Event == 'Painting'), 'Art', df.Event)

    return submgtepi2


if __name__ == "__main__":
    main()

