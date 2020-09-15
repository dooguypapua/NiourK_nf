#=====================================================
# -*- coding: utf-8 -*-                              |
# title           : Nk_mongo.py                      |    
# description     : NiourK-DB mongo functions        |
# author          : dooguypapua                      |
# copyright       : CHU Angers                       |
# date            : 20200910                         |
# version         : 0.2                              |
# python_version  : 3.8.2                            |
#=====================================================
import os
import json
import pymongo
import vcfpy
from tqdm import *
from tabulate import tabulate
from Nk_functions import *


#---------------------------------------------------------------#
#---------------------------------------------------------------#
#                   CONNECTION & INITIALIZATION                 #
#---------------------------------------------------------------#
#---------------------------------------------------------------#

#***** CONNECT *****#
def connectMongo(dicoInit):
    try :
        myclient = pymongo.MongoClient("mongodb://"+dicoInit["mongoHost"]+":"+dicoInit["mongoPort"]+"/",serverSelectionTimeoutMS=dicoInit["maxSevSelDelay"])
        myclient.server_info()
        dicoInit["db"] = myclient.NiourK_db
        dicoInit["collect_nk_run"] = dicoInit["db"].nk_run
        dicoInit["collect_nk_sample"] = dicoInit["db"].nk_sample
        dicoInit["collect_nk_var"] = dicoInit["db"].nk_var
        dicoInit["collect_nk_depth"] = dicoInit["db"].nk_depth
    except:
        exit("\nUnable to connect to `"+"mongodb://"+dicoInit["mongoHost"]+":"+dicoInit["mongoPort"]+"/"+"`\n\nAre you sure mongod is running ?\n `sudo mongod --port 27018 --dbpath /media/dooguy/ultima_thule/niourkdb`\n")

#***** INIT whole genome nk_depth *****# (deprecated)
def initNkDepth(dicoInit,dicoChrSize):
    printcolor("\n`nk_depth` collection is incomplete.\n","0",dicoInit['red'],None,dicoInit['colorBool'])
    printcolor(" Initialization can take a long time.\n Do you want to start initialization (y/n) :","0",dicoInit['red'],None,dicoInit['colorBool'])
    choice = input().lower()
    if choice=="y":
        dicoChrSize = {
                   "chr1":248956422,"chr2":242193529,"chr3":198295559,"chr4":190214555,"chr5":181538259,"chr6":170805979, \
                   "chr7":159345973,"chr8":145138636,"chr9":138394717,"chr10":133797422,"chr11":135086622,"chr12":133275309, \
                   "chr13":114364328,"chr14":107043718,"chr15":101991189,"chr16":90338345,"chr17":83257441,"chr18":80373285, \
                   "chr19":58617616,"chr20":64444167,"chr21":46709983,"chr22":50818468,"chrM":16569,"chrX":156040895,"chrY":57227415 \
                  }
        printcolor("\nInit nk_depth collection\n","1",dicoInit['blue1'],None,dicoInit['colorBool'])
        dicoInit["db"].drop_collection("nk_depth")
        for chrom in dicoChrSize:
            printcolor("  Insert "+chrom+"\n","0",dicoInit['blue1'],None,dicoInit['colorBool'])
            lstPos = list(range(1,dicoChrSize[chrom]+1,1))
            lstChunk = [lstPos[i * dicoInit["nbChunk"]:(i + 1) * dicoInit["nbChunk"]] for i in range((len(lstPos) + dicoInit["nbChunk"] - 1) // dicoInit["nbChunk"] )]  
            for i in tqdm(range(len(lstChunk)),ncols=60,bar_format="    {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt}"):
                chunk = lstChunk[i]
                insertList = []
                for pos in chunk: insertList.append({ "_id":chrom+"_"+str(pos) })
                dicoInit["collect_nk_depth"].insert_many(insertList,ordered=False)





#---------------------------------------------------------------#
#---------------------------------------------------------------#
#                        STATS & SUMMARY                        #
#---------------------------------------------------------------#
#---------------------------------------------------------------#

#***** STATS *****#
def dbStats(dicoInit):
    stats = dicoInit["db"].command("dbstats") # Call stats
    printcolor("\nSub-command: stats\n\n","1",dicoInit['blue1'],None,dicoInit['colorBool'])
    printcolor("  DATABASE\n","1",dicoInit['blue2'],None,dicoInit['colorBool'])
    printcolor("    name    : "+stats['db']+"\n","0",dicoInit['blue2'],None,dicoInit['colorBool'])
    printcolor("    objects : "+str(stats['objects'])+"\n","0",dicoInit['blue2'],None,dicoInit['colorBool'])
    printcolor("    size    : "+convertByteSize(stats['dataSize'])+"\n\n","0",dicoInit['blue2'],None,dicoInit['colorBool'])
    printcolor("  COLLECTIONS\n","1",dicoInit['white'],None,dicoInit['colorBool'])
    for collection in dicoInit["db"].collection_names():
        printcolor("    name    : "+collection+"\n","0",dicoInit['white'],None,dicoInit['colorBool'])
        printcolor("    objects : "+str(dicoInit["collect_"+collection].count())+"\n","0",dicoInit['white'],None,dicoInit['colorBool'])
        printcolor("    size    : "+convertByteSize(dicoInit["db"].command("collstats", collection)["size"])+"\n\n","0",dicoInit['white'],None,dicoInit['colorBool'])

#***** LIST RUN *****#
def listRun(dicoInit):
    runKeys = ["instrument","project","ref","target","name","num","date","seqname","seqnum","seqdate","seqstatus","chip","seqkit","libKit"]
    prettyRunKeys = ["instrument","project","ref","name","num","date"]
    table = []
    # Output CSV
    OUT = open(dicoInit["pathOutput"],'w')
    OUT.write("\t".join(runKeys)+"\n")
    for runEntry in dicoInit["collect_nk_run"].find().sort("num",pymongo.ASCENDING):
        ToWrite = ""
        row = []
        for key in runKeys:
            if runEntry["name"]+"_"+str(runEntry["num"]) in ["Auto_user_Proton-455-171010CMTL2panel6poolsHi-Q_658_999","Auto_user_Proton-459-RUN024_ENDOC_628_1007","Auto_user_Proton-460-RUN094_MITO_627_1009","Auto_user_Proton-461-RUN023_NOH_661_1011","Auto_user_Proton-462-RUN019_NUCL_662_1013","Auto_user_Proton-463-171020TSCcustompanel4poolsHi-Q_666_1015","Auto_user_Proton-464-RUN098_MITO_663_1017","Auto_user_Proton-465-RUN007_MITOME_PUCE1_664_1019","Auto_user_Proton-466-RUN007_MITOME_PUCE2_665_1021","Auto_user_Proton-467-171025TSCCUSTOMPANEL4POOLSHI-Q_667_1023","Auto_user_Proton-468-171010CMTL1panel6poolsHi-Q_659_1025","Auto_user_Proton-469-171106TSCcustompanel4poollsHi-Q_668_1027","Auto_user_Proton-470-171019CMTL1panel6poolHi-Q_669_1029","Auto_user_Proton-471-RUN005_VIRGINIE_671_1031","Auto_user_Proton-472-RUN006_VIRGINIE_672_1033","Auto_user_Proton-473-RUN048_ONCO_670_1035","Auto_user_Proton-474-RUN099_MITO_673_1037","Auto_user_Proton-475-RUN020_NUCL_675_1039","Auto_user_Proton-476-RUN025_ENDOC_674_1041","Auto_user_Proton-477-171116TSCCUSTOMPANEL4POOLSHI-Q_679_1043","Auto_user_Proton-478-171019CMTL2panel6poolHi-Q_680_1045","Auto_user_Proton-479-RUN024_NOH_676_1047","Auto_user_Proton-480-RUN100_MITO_681_1049","Auto_user_Proton-481-171127TSCcustompanel4poolsHi-Q_683_1051","Auto_user_Proton-482-RUN025_NOH_682_1053","Auto_user_Proton-483-RUN026_NOH_686_1055","Auto_user_Proton-484-RUN026_ENDOC_684_1057","Auto_user_Proton-485-RUN049_ONCO_687_1059","Auto_user_Proton-486-RUN101_MITO_685_1061","Auto_user_Proton-487-171213TSCcustompanel4poolsHi-Q_695_1063","Auto_user_Proton-488-171207CMTCUSTOMPANEL6POOLSHI-Q_688_1065","Auto_user_Proton-489-RUN102_MITO_696_1067","Auto_user_Proton-490-RUN021_NUCL_694_1069","Auto_user_Proton-491-171207CMT2CUSTOMPANEL6POOLSHI-Q_689_1071","Auto_user_Proton-492-171208CMT1CUSTOMPANEL6POOLSHI-Q_690_1073","Auto_user_Proton-493-RUN027_ENDOC_697_1075","Auto_user_Proton-494-171208CMT2CUSTOMPANEL6POOLSHI-Q_691_1077","Auto_user_Proton-495-171226TSCcustompanel4poolsHi-Q_698_1079","Auto_user_Proton-496-171211CMTL1panel6poolsHi-Q_692_1081","Auto_user_Proton-497-RUN103_MITO_701_1083","Auto_user_Proton-498-RUN028_ENDOC_702_1085","Auto_user_Proton-499-171211CMTL2panel6poolsHi-Q_693_1087","Auto_user_Proton-500-171220CMTL1panel6poolsHi-Q_699_1089","Auto_user_Proton-501-RUN050_ONCO_708_1091","Auto_user_Proton-502-RUN027_NOH_707_1093","Auto_user_Proton-503-180110TSCcustompanel4poolsHi-Q_709_1095","Auto_user_Proton-504-171220CMTL2panel6poolsHi-Q_700_1097","Auto_user_Proton-505-180115TSCCUSTOMPANEL4POOLSHI-Q_710_1099","Auto_user_Proton-506-171228CMT1CUSTOMPANEL6POOLSHI-Q_703_1101","Auto_user_Proton-507-RUN104_MITO_712_1103","Auto_user_Proton-508-RUN105_MITO_711_1105","Auto_user_Proton-509-RUN010BIS_MITOME_PUCE1_713_1107","Auto_user_Proton-510-RUN010BIS_MITOME_PUCE2_714_1109","Auto_user_Proton-511-180119TSCcustompanel4poolsHi-Q_718_1111","Auto_user_Proton-512-180122CMT1CUSTOMPANEL6POOLSHI-Q_716_1113","Auto_user_Proton-513-180126TSCCUSTOMPANEL4POOLSHI-Q_721_1115","Auto_user_Proton-514-180125CMTL2panel6poolsHi-Q_720_1117","Auto_user_Proton-515-171228CMT2CUSTOMPANEL6POOLSHI-Q_704_1119","Auto_user_Proton-516-180122CMT2CUSTOMPANEL6POOLSHI-Q_717_1121","Auto_user_Proton-517-RUN106_MITO_723_1123","Auto_user_Proton-518-RUN051_ONCO_722_1125","Auto_user_Proton-519-RUN029_ENDOC_BIS_715_1127","Auto_user_Proton-520-RUN022_NUCL_724_1129","Auto_user_Proton-521-180125CMTL1panel6poolsHi-Q_719_1131","Auto_user_Proton-522-180201CMTpanel6poolsHi-Q_725_1133","Auto_user_Proton-523-RUN029_NOH_726_1135","Auto_user_Proton-524-RUN028_NOH_728_1137","Auto_user_Proton-525-180207TSCcustompanel4poolsHi-Q_727_1139","Auto_user_Proton-526-180214TSCCUSTOMPANEL4POOLSHI-Q_731_1141","Auto_user_Proton-527-RUN108_MITO_732_1143","Auto_user_Proton-528-RUN107_MITO_730_1145","Auto_user_Proton-529-RUN011_MITOME_PUCE1_733_1147","Auto_user_Proton-530-RUN011_MITOME_PUCE2_734_1149","Auto_user_Proton-531-RUN030_NOH_729_1151","Auto_user_Proton-532-RUN052_ONCO_738_1153","Auto_user_Proton-533-RUN030_ENDOC_737_1155","Auto_user_Proton-534-RUN109_MITO_740_1157","Auto_user_Proton-535-180226TSCCUSTOMPANEL4POOLSHI-Q_739_1159","Auto_user_Proton-536-180108CMT1CUSTOMPANEL6POOLSHI-Q_705_1161","Auto_user_Proton-537-RUN113_MITO_747_1163","Auto_user_Proton-539-180319TSCCUSTOMPANEL4POOLSHI-Q_749_1167","Auto_user_Proton-540-RUN011_MITOME_PUCE1_2EPASSAGE_750_1169","Auto_user_Proton-541-RUN023_NUCL_744_1171","Auto_user_Proton-542-RUN110_MITO_741_1173","Auto_user_Proton-543-RUN111_MITO_745_1175","Auto_user_Proton-544-RUN031_ENDOC_748_1177","Auto_user_Proton-545-RUN012_MITOME_PUCE1_751_1179","Auto_user_Proton-546-RUN012_MITOME_PUCE2_752_1181","Auto_user_Proton-547-RUN112_MITO_754_1183","Auto_user_Proton-548-RUN031_NOH_753_1185","Auto_user_Proton-549-RUN053_ONCO_758_1187","Auto_user_Proton-550-RUN114_MITO_755_1189","Auto_user_Proton-551-180326CMTL1panel6poolsHi-Q_756_1191","Auto_user_Proton-552-180326CMTL2panel6poolsHi-Q_757_1193","Auto_user_Proton-553-180403TSCcustompanel4poolsHi-Q_759_1195","Auto_user_Proton-554-180108CMT2CUSTOMPANEL6POOLSHI-Q_706_1197","Auto_user_Proton-555-RUN053_ONCO_BIS_760_1199","Auto_user_Proton-556-RUN114_MITO_BIS_761_1201","Auto_user_Proton-557-RUN008bis_MITOME_MERIEM_PUCE1_763_1203","Auto_user_Proton-558-RUN008bis_MITOME_MERIEM_PUCE2_764_1205","Auto_user_Proton-559-180416TSCCUSTOMPANEL4POOLSHI-Q_765_1207","Auto_user_Proton-560-180220CMT1PANEL6POOLSHI-Q_742_1209","Auto_user_Proton-561-RUN032_NOH_762_1211","Auto_user_Proton-562-180220CMT2PANEL6POOLSHI-Q_743_1213","Auto_user_Proton-563-180420TSCCUSTOMPANEL4POOLSHI-Q_767_1215","Auto_user_Proton-564-RUN032_ENDOC_766_1217","Auto_user_Proton-565-180426TSCCUSTOMPANEL4POOLSHI-Q_773_1219","Auto_user_Proton-566-180410CMTL1panel6poolsHi-Q_771_1221","Auto_user_Proton-567-RUN115_MITO_769_1223","Auto_user_Proton-568-RUN024_NUCL_768_1225","Auto_user_Proton-569-RUN054_ONCO_774_1227","Auto_user_Proton-570-RUN116_MITO_770_1229","Auto_user_Proton-571-180524TSCCUSTOMPANEL4POOLSHI-Q_777_1231","Auto_user_Proton-572-180410CMTL2panel6poolsHi-Q_772_1233","Auto_user_Proton-573-RUN033_ENDOC_775_1235","Auto_user_Proton-574-RUN117_MITO_776_1237","Auto_user_Proton-575-RUN118_MITO_781_1239","Auto_user_Proton-576-RUN033_NOH_778_1241","Auto_user_Proton-577-RUN055_ONCO_779_1243","Auto_user_Proton-578-RUN034_NOH_780_1245","Auto_user_Proton-579-180605TSCCUSTOMPANEL4POOLSHI-Q_782_1247","Auto_user_Proton-580-180605TSCCUSTOMPANEL4POOLSHI-Qbis_783_1249","Auto_user_Proton-582-180523CMTL1panel6poolsHi-Q_787_1253","Auto_user_Proton-583-RUN013_MITOME_PUCE1_784_1255","Auto_user_Proton-584-RUN013_MITOME_PUCE2_785_1257","Auto_user_Proton-585-RUN119_MITO_791_1259","Auto_user_Proton-586-RUN001_ONCOENDOC_CAP_792_1261","Auto_user_Proton-587-180621TSCCUSTOMPANEL4POOLSHI-Q_794_1263","Auto_user_Proton-588-180523CMTL2panel6poolsHi-Q_788_1265","Auto_user_Proton-589-RUN014_MITOME_PUCE1_796_1267","Auto_user_Proton-590-RUN014_MITOME_PUCE2_797_1269","Auto_user_Proton-591-RUN034_ENDOC_793_1271","Auto_user_Proton-592-RUN035_NOH_795_1273","Auto_user_Proton-593-RUN025_NUCL_786_1275","Auto_user_Proton-594-RUN001_NOH_CAP_805_1277","Auto_user_Proton-595-180702TSCCUSTOMPANEL4POOLSHI-Q_803_1279","Auto_user_Proton-596-180615CMTL1panel6poolsHi-Q_798_1281","Auto_user_Proton-597-RUN120BIS_MITO_809_1283","Auto_user_Proton-598-RUN056BIS_ONCO_810_1285","Auto_user_Proton-599-RUN036BIS_NOH_814_1287","Auto_user_Proton-600-RUN121_MITO_SW_813_1289","Auto_user_Proton-601-180713TSCCUSTOMPANEL4POOLSHI-Q_812_1291","Auto_user_Proton-602-180615CMTL2panel6poolsHi-Q_799_1293","Auto_user_Proton-603-RUN122_MITO_818_1295","Auto_user_Proton-604-RUN035_ENDOC_811_1297","Auto_user_Proton-605-RUN015_MITOME_V2_PUCE1_815_1299","Auto_user_Proton-606-RUN015_MITOME_V2_PUCE2_816_1301","Auto_user_Proton-607-RUN123_MITO_820_1303","Auto_user_Proton-608-RUN016_MITOME_V1_PUCE2_824_1305","Auto_user_Proton-609-RUN016_MITOME_V1_PUCE1_823_1307","Auto_user_Proton-610-RUN057_ONCO_825_1309","Auto_user_Proton-611-RUN036_ENDOC_819_1311","Auto_user_Proton-612-RUN037_NOH_822_1313","Auto_user_Proton-613-180807TSCcustompanel4poolsHi-Q_826_1315","Auto_user_Proton-614-180718TSCCUSTOMPANEL4POOLSHI-Q_817_1317","Auto_user_Proton-615-180807TSCCUSTOMPANEL4POOLSHI-QBIS_827_1319","Auto_user_Proton-616-RUN125_MITO_828_1321","Auto_user_Proton-617-180718TSCCUSTOMPANEL4POOLSHI-QBIS_830_1323","Auto_user_Proton-618-180813TSCCUSTOMPANEL4POOLSHI-Q_829_1325","Auto_user_Proton-619-180822TSCCUSTOMPANEL4POOLSHI-Q_831_1327","Auto_user_Proton-620-180616CMTL1panel6poolsHi-Q_800_1329","Auto_user_Proton-621-RUN124_MITO_821_1331","Auto_user_Proton-622-RUN127_MITO_832_1333","Auto_user_Proton-623-RUN125BIS_MITO_835_1335","Auto_user_Proton-624-RUN038_NOH_834_1337","Auto_user_Proton-625-RUN128_MITO_833_1339","Auto_user_Proton-626-RUN126_MITO_836_1341","Auto_user_Proton-627-180904CMTL1panel6poolsHi-Q_842_1343","Auto_user_Proton-628-180904CMTL2panel6poolsHi-Q_843_1345","Auto_user_Proton-629-RUN039_NOH_837_1347","Auto_user_Proton-630-RUN058_ONCO_841_1349","Auto_user_Proton-631-RUN026_NUCL_847_1351","Auto_user_Proton-632-RUN037_ENDOC_846_1353","Auto_user_Proton-633-180616CMTL2panel6poolsHi-Q_801_1355","Auto_user_Proton-634-RUN129_MITO_838_1357","Auto_user_Proton-635-180925TSCCUSTOMPANEL4POOLSHI-Q_850_1359","Auto_user_Proton-636-180925CMTCUSTOMPANEL6POOLSHI-Q_851_1361","Auto_user_Proton-638-RUN128BIS_MITO_853_1365","Auto_user_Proton-639-RUN017_MITOME_PUCE1_848_1367","Auto_user_Proton-640-RUN017_MITOME_PUCE2_849_1369","Auto_user_Proton-641-180615CMTL2PANEL6POOLHI-QBIS_856_1371","Auto_user_Proton-642-RUN129_MITO_TEST_855_1373","Auto_user_Proton-649-RUN038_ENDOC_870_1388","Auto_user_Proton-650-180925L1panel6poolsHi-Q-BIS_871_1390","Auto_user_Proton-651-RUN040_NOH_867_1392","Auto_user_Proton-652-RUN041_NOH_875_1394","Auto_user_Proton-653-181023TSCCUSTOMPANEL4POOLSHI-Q_877_1397","Auto_user_Proton-654-180910CMT1panel6poolsHi-Q_839_1399","Auto_user_Proton-655-180910CMT2panel6poolsHi-Q_840_1401","Auto_user_Proton-656-RUN131_MITO_874_1403","Auto_user_Proton-657-RUN132_MITO_876_1406","Auto_user_Proton-658-RUN042_NOH_878_1408","Auto_user_Proton-659-181004CMTL1panel6poolsHi-Q_872_1410","Auto_user_Proton-660-181004CMTL2panel6poolsHi-Q_873_1412","Auto_user_Proton-661-180905CMTL1panel6poolsHi-Q_844_1414","Auto_user_Proton-662-180905CMTL2panel6poolsHi-Q_845_1416","Auto_user_Proton-663-RUN060_ONCO_881_1418","Auto_user_Proton-664-RUN133_MITO_882_1420","Auto_user_Proton-665-181114TSCCUSTOMPANEL4POOLSHI-Q_884_1422","Auto_user_Proton-666-181031CMTL1panel6poolsHi-Q_879_1424","Auto_user_Proton-667-RUN134_MITO_885_1426","Auto_user_Proton-668-RUN043_NOH_883_1428","Auto_user_Proton-669-RUN009_MITOME_MERIEM_PUCE_B_866_1430","Auto_user_Proton-670-RUN009_MITOME_MERIEM_PUCE_A_865_1432","Auto_user_Proton-671-RUN027_NUCL_890_1434","Auto_user_Proton-672-RUN039_ENDOC_886_1436","Auto_user_Proton-673-181121TSCcustompanel4poolsHi-Q_889_1438","Auto_user_Proton-674-181031CMTL2panel6poolsHi-Q_880_1440","Auto_user_Proton-675-RUN018_Mitome_v2_PUCE1_892_1442","Auto_user_Proton-676-RUN018_Mitome_v2_PUCE2_893_1444","Auto_user_Proton-677-RUN061_ONCO_896_1446","Auto_user_Proton-678-RUN135_MITO_891_1448","Auto_user_Proton-679-181204TSCCUSTOMPANEL4POOLSHI-Q_897_1450","Auto_user_Proton-680-181113CMTL2panel6poolsHi-Q_888_1452","Auto_user_Proton-681-181213TSCCUSTOMPANEL4POOLSHI-Q_908_1454","Auto_user_Proton-682-181113CMTL1panel6poolsHi-Q_887_1456","Auto_user_Proton-683-RUN002_NOH_CAP_898_1458","Auto_user_Proton-684-RUN040_ENDOC_901_1460","Auto_user_Proton-685-RUN003_NOH_CAP_907_1462","Auto_user_Proton-686-RUN136_MITO_911_1464","Auto_user_Proton-689-181227custompanel4poolsHi-Q_915_1470","Auto_user_Proton-690-181212CMTL1panel6poolsHi-Q_919_1472","Auto_user_Proton-691-RUN004_NOH_CAP_913_1474","Auto_user_Proton-692-RUN137_MITO_912_1476","Auto_user_Proton-693-RUN138_MITO_914_1478","Auto_user_Proton-694-RUN005_NOH_CAP_916_1480","Auto_user_Proton-697-190115TSCCUSTOMPANEL4POOLSHI-Q_927_1486","Auto_user_Proton-698-181212CMTL2panel6poolsHi-Q_920_1488","Auto_user_Proton-699-RUN062_ONCO_921_1490","Auto_user_Proton-700-RUN041_ENDOC_926_1492","Auto_user_Proton-701-RUN139_MITO_922_1494","Auto_user_Proton-702-RUN006_NOH_CAP_923_1496","Auto_user_Proton-703-190123TSCcustompanel4poolsHi-Q_933_1498","Auto_user_Proton-704-190121CMTL1PANEL6POOLSHI-Q_930_1500","Auto_user_Proton-705-RUN140_MITO_928_1502","Auto_user_Proton-706-RUN041BIS_ENDOC_929_1504","Auto_user_Proton-708-RUN007_NOH_CAP_932_1508","Auto_user_Proton-709-RUN141_MITO_936_1510","Auto_user_Proton-710-RUN063_ONCO_938_1512","Auto_user_Proton-711-190201CMTL1panel6poolsHi-Q_940_1514","Auto_user_Proton-712-RUN008_NOH_CAP_935_1516","Auto_user_Proton-713-190206tsccustompanel4poolshi-q_941_1518","Auto_user_Proton-714-190121CMTL2PANEL6POOLSHI-Q_931_1520","Auto_user_Proton-715-190208TSCcustompanel4poolsHi-Q_944_1522","Auto_user_Proton-716-NGS_NOH009_CAP_939_1524","Auto_user_Proton-718-RUN019_MITOME_PUCE_B_943_1528","Auto_user_Proton-719-190225tsccustompanel4poolshi-q_946_1530","Auto_user_Proton-720-190225CMTL1panel6poolsHi-Q_948_1532","Auto_user_Proton-721-RUN143_MITO_950_1534","Auto_user_Proton-722-RUN042_ENDOC_937_1536","Auto_user_Proton-723-RUN010_NOH_CAP_947_1538","Auto_user_Proton-724-RUN011_NOH_CAP_945_1540","Auto_user_Proton-725-190304TSCcustompanel4poolsHi-Q_954_1542","Auto_user_Proton-726-190225CMTL2panel6poolsHi-Q_949_1544","Auto_user_Proton-727-RUN020_MITOME_PUCE1_952_1546","Auto_user_Proton-728-RUN020_MITOME_PUCE2_953_1548","Auto_user_Proton-729-RUN012_NGS_NOH_CAP_955_1550","Auto_user_Proton-730-RUN043_ENDOC_997_1552","Auto_user_Proton-731-190311TSCcustompanel4poolsHi-Q_1003_1554","Auto_user_Proton-732-190315TSCcustompanel4poolsHi-Q_1002_1556","Auto_user_Proton-733-RUN064_ONCO_998_1558","Auto_user_Proton-734-RUN142_MITO_1005_1560","Auto_user_Proton-735-RUN021_MITOME_V2_PUCE1_1000_1562","Auto_user_Proton-736-RUN021_MITOME_V2_PUCE2_1001_1564","Auto_user_Proton-738-RUN015_NOH_CAP_1004_1568","Auto_user_Proton-741-190326TSCbisCUSTOMPANEL4POOLSHI-Q_1016_1576","Auto_user_Proton-742-RUN144_MITO_1009_1578","Auto_user_Proton-743-RUN065_ONCO_1022_1580","Auto_user_Proton-744-RUN014_NOH_CAP_1014_1582","Auto_user_Proton-745-RUN147_MITO_1018_1584","Auto_user_Proton-746-RUN016_NOH_CAP_1010_1586","Auto_user_Proton-747-RUN022_MITOME_V2_PUCE_1_1006_1588","Auto_user_Proton-748-RUN022_MITOME_V2_PUCE_2_1007_1590","Auto_user_Proton-749-190409TSCcustompanel4poolsHi-Q_1025_1592","Auto_user_Proton-750-190327TSCbisCUSTOMPANEL4POOLSHI-Q_1021_1594","Auto_user_Proton-751-RUN145_MITO_1015_1596","Auto_user_Proton-752-RUN017_NOH_CAP_1013_1598","Auto_user_Proton-753-RUN146_MITO_1020_1600","Auto_user_Proton-754-RUN018_NOH_CAP_1017_1602","Auto_user_Proton-755-RUN023_MITOME_V2_PUCE2_1024_1604","Auto_user_Proton-756-RUN023_MITOME_V2_PUCE1_1023_1606","Auto_user_Proton-757-RUN148_MITO_1027_1608","Auto_user_Proton-758-RUN044_ENDOC_1026_1610","Auto_user_Proton-759-RUN019_NOH_CAP_1028_1612","Auto_user_Proton-760-RUN149_MITO_1031_1614","Auto_user_Proton-761-190502TSCcustompanel4poolsHi-Q_1033_1616","Auto_user_Proton-762-190408CMTL1panel6poolsHi-Q_1029_1618","Auto_user_Proton-763-RUN045_ENDOC_1032_1620","Auto_user_Proton-764-RUN150_MITO_1034_1622","Auto_user_Proton-765-RUN020_NOH_CAP_1036_1624","Auto_user_Proton-766-RUN066_ONCO_1035_1626","Auto_user_Proton-771-RUN024_MITOME_PUCE1_1040_1636","Auto_user_Proton-772-RUN024_MITOME_PUCE2_1041_1638","Auto_user_Proton-773-190528TSCcustompanel4poolsHi-Q_1046_1640","Auto_user_Proton-774-190502CMTL1panel6poolsHi-Q_1043_1642","Auto_user_Proton-775-RUN152_MITO_1042_1644","Auto_user_Proton-776-RUN022_NOH_CAP_1045_1646","Auto_user_Proton-777-RUN153_MITO_1048_1648","Auto_user_Proton-778-RUN046_ENDOC_1047_1650","Auto_user_Proton-779-RUN023_NOH_CAP_1049_1652","Auto_user_Proton-780-RUN067_ONCO_1050_1654","Auto_user_Proton-781-190614TSCCUSTOMPANEL4POOLSHI-Q_1053_1656","Auto_user_Proton-782-190502CMTL2panel6poolsHi-Q_1044_1658","Auto_user_Proton-783-RUN148BIS_MITO_1051_1660","Auto_user_Proton-784-RUN024_NOH_CAP_1054_1662","Auto_user_Proton-785-190624TSCcustompanel4poolsHi-Q_1058_1664","Auto_user_Proton-786-190618CMTL1panel6poolsHi-Q_1059_1666","Auto_user_Proton-789-RUN155_MITO_1057_1672","Auto_user_Proton-790-RUN026_NOH_CAP_1055_1674","Auto_user_Proton-791-RUN154_MITO_BIS_1065_1676","Auto_user_Proton-792-RUN025_NOH_CAP_BIS_1066_1678","Auto_user_Proton-793-RUN025_MITOME_PUCE1_1061_1680","Auto_user_Proton-794-RUN025_MITOME_PUCE2_1062_1682","Auto_user_Proton-795-RUN026_MITOME_PUCE1_1071_1684","Auto_user_Proton-796-RUN026_MITOME_PUCE2_1072_1686","Auto_user_Proton-797-190704TSCcustompanel4poolsHi-Q_1075_1688","Auto_user_Proton-798-190618CMTL2panel6poolsHi-Q_1060_1690","Auto_user_Proton-799-RUN156_MITO_1064_1692","Auto_user_Proton-800-RUN068_ONCO_1070_1694","Auto_user_Proton-801-RUN157_MITO_1068_1696","Auto_user_Proton-802-RUN047_ENDOC_1063_1698","Auto_user_Proton-803-RUN027_NOH_CAP_1067_1700","Auto_user_Proton-804-RUN_028_NOH_CAP_1069_1702","Auto_user_Proton-805-190722TSCCUSTOMPANEL4POOLSHI-Q_1077_1705","Auto_user_Proton-806-190701CMTL1panel6poolsHI-Q_1073_1707","Auto_user_Proton-807-RUN158_MITO_1079_1709","Auto_user_Proton-808-RUN029_NOH_CAP_1078_1711","Auto_user_Proton-809-190730CMTL1panel6poolsHi-Q_1084_1713","Auto_user_Proton-813-RUN069_ONCO_1088_1721","Auto_user_Proton-814-RUN030_NOH_CAP_1086_1723","Auto_user_Proton-815-190814TSCCUSTOMPANEL4POOLSHI-Q_1089_1725","Auto_user_Proton-816-190701CMTL2panel6poolsHi-Q_1074_1727","Auto_user_Proton-817-RUN049_ENDOC_1090_1729","Auto_user_Proton-818-RUN031_NOH_CAP_1087_1731","Auto_user_Proton-819-190902TSCcustompanel4poolsHi-Q_1093_1733","Auto_user_Proton-820-190730CMTL2panel6poolsHi-Q_1085_1735","Auto_user_Proton-821-RUN032_NOH_CAP_1092_1737","Auto_user_Proton-822-RUN161_MITO_1091_1739","Auto_user_Proton-823-RUN160_MITO_1096_1741","Auto_user_Proton-824-RUN162_MITO_1094_1743","Auto_user_Proton-825-190913TSCCUSTOMPANEL4POOLSHI-Q_1103_1745","Auto_user_Proton-826-190909CMTL1panel6poolsHi-Q_1100_1747","Auto_user_Proton-827-RUN050_ENDOC_1104_1749","Auto_user_Proton-828-RUN033_NOH_CAP_1097_1751","Auto_user_Proton-830-RUN070_ONCO_1098_1755","Auto_user_Proton-831-RUN034_NOH_CAP_1106_1757","Auto_user_Proton-832-RUN165_MITO_1099_1759","Auto_user_Proton-833-RUN027_MITOME_PUCE1_1105_1761","Auto_user_Proton-834-RUN027_MITOME_PUCE2_1107_1763","Auto_user_Proton-835-190924TSCcustompanel4poolsHi-Q_1111_1765","Auto_user_Proton-837-RUN036_NOH_CAP_1108_1769","Auto_user_Proton-838-RUN164_MITO_1102_1771","Auto_user_Proton-839-191001TSCcustompanel4poolsHi-Q_1115_1773","Auto_user_Proton-840-190909CMTL2panel6poolsHi-Q_1101_1775","Auto_user_Proton-841-RUN035_NOH_CAP_1109_1777","Auto_user_Proton-842-RUN051_ENDOC_1113_1779","Auto_user_Proton-843-RUN037_NOH_CAP_1114_1781","Auto_user_Proton-844-RUN038_NOH_CAP_1116_1783","Auto_user_Proton-845-191008TSCCUSTOMPANEL4POOLSHI-Q_1117_1785","Auto_user_Proton-846-RUN166_MITO_1110_1787","Auto_user_Proton-847-191011TSCcustompanel4poolsHi-Q_1121_1789","Auto_user_Proton-848-RUN167_MITO_1120_1791","Auto_user_Proton-849-RUN039_NOH_CAP_1118_1793","Auto_user_Proton-850-RUN071_ONCO_1122_1795","Auto_user_Proton-851-RUN040_NOH_CAP_1119_1797","Auto_user_Proton-852-RUN052_ENDOC_1125_1799","Auto_user_Proton-853-191010CMTL1panel6poolsHi-Q_1123_1801","Auto_user_Proton-854-RUN168_MITO_1126_1803","Auto_user_Proton-855-RUN042_NOH_CAP_1127_1805","Auto_user_Proton-856-191010CMTL2panel6poolsHi-Q_1124_1807","Auto_user_Proton-857-191021BISTSCcustompanel4poolsHi-Q_1130_1809","Auto_user_Proton-858-RUN041_NOH_CAP_1128_1811","Auto_user_Proton-859-RUN166BIS_MITO_1134_1813","Auto_user_Proton-860-191104TSCCUSTOMPANEL4POOLSHI-Q_1131_1815","Auto_user_Proton-861-RUN072_ONCO_1136_1817","Auto_user_Proton-862-RUN169_MITO_1135_1819","Auto_user_Proton-863-RUN042bis_NOH_CAP_1137_1821","Auto_user_Proton-864-RUN053_ENDOC_1138_1823","Auto_user_Proton-865-191119TSCcustompanel4poolsHi-Q_1142_1825","Auto_user_Proton-866-191022CMTL1panel6poolsHi-Q_1132_1827","Auto_user_Proton-867-RUN028_MITOME_PUCE_A_1139_1829","Auto_user_Proton-868-RUN028_MITOME_PUCE_B_1140_1831","Auto_user_Proton-870-RUN044_NOH_CAP_1144_1835","Auto_user_Proton-871-RUN043_NOH_CAP_1145_1837","Auto_user_Proton-872-RUN054_ENDOC_1146_1839","Auto_user_Proton-873-191128TSCCUSTOMPANEL4POOLSHI-Q_1148_1841","Auto_user_Proton-874-191022CMTL2panel6poolsHi-Q_1133_1843","Auto_user_Proton-875-RUN045_NOH-CAP_1149_1845","Auto_user_Proton-876-RUN171_MITO_1143_1847","Auto_user_Proton-877-RUN046_NOH_CAP_1147_1849","Auto_user_Proton-878-RUN047_NOH_CAP_1152_1851","Auto_user_Proton-879-RUN048_NOH_CAP_1150_1853","Auto_user_Proton-880-RUN172_MITO_1151_1855","Auto_user_Proton-881-191210TSCcustompanel4poolsHi-Q_1156_1857","Auto_user_Proton-882-191202CMTL1panel6poolsHi-Q_1153_1859","Auto_user_Proton-883-RUN173_MITO_1157_1861","Auto_user_Proton-884-RUN073_ONCO_1155_1863","Auto_user_Proton-885-191218TSCCUSTOMPANELPOOLSHI-Q_1161_1865","Auto_user_Proton-886-191202CMTL2panel6poolsHi-Q_1154_1867","Auto_user_Proton-887-RUN055_ENDOC_1158_1869","Auto_user_Proton-888-RUN174_MITO_1162_1871","Auto_user_Proton-889-RUN029_MITOME_PUCE_A_1159_1873","Auto_user_Proton-890-RUN029_MITOME_PUCE_B_1160_1875","Auto_user_Proton-891-RUN056_ENDOC_1164_1877","Auto_user_Proton-892-RUN049_NOH_CAP_1163_1879","Auto_user_Proton-893-RUN175_MITO_1165_1881","Auto_user_Proton-894-RUN003_ONCO_ENDOC_CAP_1166_1883","Auto_user_Proton-895-200110TSCcustompanel4poolsHi-Q_1167_1885","Auto_user_Proton-896-200107CMTL1panel6poolsHi-Q_1168_1887","Auto_user_Proton-897-RUN176_MITO_1172_1889","Auto_user_Proton-898-RUN074_ONCO_1170_1891","Auto_user_Proton-899-200120CMTL1panel6poolsHi-Q_1174_1893","Auto_user_Proton-900-200107CMTL2panel6poolsHi-Q_1169_1895","Auto_user_Proton-901-RUN177_MITO_1173_1897","Auto_user_Proton-902-RUN050_NOH_CAP_1171_1899","Auto_user_Proton-903-RUN178_MITO_1176_1901","Auto_user_Proton-904-RUN051_NOH_CAP_reanalyse_1910","Auto_user_Proton-905-RUN052_NOH_CAP_1177_1906","Auto_user_Proton-906-RUN179_MITO_1178_1908","Auto_user_Proton-907-200205TSCcustompanel4poolsHi-Q_1182_1911","Auto_user_Proton-908-RUN053_NOH_CAP_1179_1913","Auto_user_Proton-909-RUN054_NOH_CAP_1181_1915","Auto_user_Proton-910-RUN180_MITO_1180_1917","Auto_user_Proton-911-RUN181_MITO_1183_1919","Auto_user_Proton-912-RUN057_ENDOC_1184_1921","Auto_user_Proton-913-200214TSCcustompanel4poolsHi-Q_1188_1923","Auto_user_Proton-915-RUN183_MITO_1189_1927","Auto_user_Proton-916-RUN075_ONCO_1186_1929","Auto_user_Proton-917-200226TSCcustompanel4poolsHi-Q_1191_1931","Auto_user_Proton-918-RUN182_MITO_1187_1933","Auto_user_Proton-919-RUN076_panel_ONCO_1195_1935","Auto_user_Proton-920-RUN184_MITO_1193_1937","Auto_user_Proton-921-200317TSCCUSTOMPANEL4POOLSHI-Q_1197_1939","Auto_user_Proton-922-RUN058_ENDOC_1192_1941","Auto_user_Proton-923-RUN056_NOH_CAP_1194_1943","Auto_user_Proton-925-RUN186_MITO_1203_1947","Auto_user_Proton-926-RUN057_NOH_CAP_1200_1949","Auto_user_Proton-930-RUN059_ENDOC_142_004","Auto_user_Proton-931-200305CMTL1panel6poolsHi-Q_141_006","Auto_user_Proton-932-200424TSCCUSTOMPANEL4POOLSHI-Q_146_008","Auto_user_Proton-933-200305CMTL2POOLSHI-Q_147_010","Auto_user_Proton-934-RUN187_MITO_148_012","Auto_user_Proton-935-RUN077_ONCO_149_014","Auto_user_Proton-936-200518TSCCUSTOMPANEL4POOLSHI-Q_151_016","Auto_user_Proton-937-200402CMTL1POOLSHI-Q_144_018","Auto_user_Proton-938-RUN188_MITO_153_020","Auto_user_Proton-939-RUN058_NOH_CAP_150_022","Auto_user_Proton-940-200529TSCCUSTOMPANEL4POOLSHI-Q_155_024","Auto_user_Proton-941-200402CMTL2POOLSHI-Q_145_026","Auto_user_Proton-942-RUN059_NOH_CAP_152_028","Auto_user_Proton-943-RUN060_ENDOC_156_030","Auto_user_Proton-944-200609TSCCUSTOMPANEL4POOLSHI-Q_159_032","Auto_user_Proton-945-200525CMTL1PANEL6POOLSHI-Q_160_034","Auto_user_Proton-946-RUN060_NOH_CAP_154_036","Auto_user_Proton-947-RUN061_NOH_CAP_157_038","Auto_user_Proton-949-RUN062_NOH_CAP_158_042","Auto_user_Proton-950-200626TSCCUSTOMPANEL4POOLSHI-Q_163_044","Auto_user_S5XL-11-PUCE_NOH_RUN08_72_027","Auto_user_S5XL-12-PUCE_NOH_RUN09_73_029","Auto_user_S5XL-187-RUN129_MITO_356_410","Auto_user_S5XL-188-RUN130_MITO_357_412","Auto_user_S5XL-193-RUN017_MITOME_PUCE1_362_422","Auto_user_S5XL-194-RUN017_MITOME_PUCE2_363_424","Auto_user_S5XL-204-mtDNA_HIV_Val_181121_381_444","Auto_user_S5XL-207-RUN_020_mtDNA_MAJIDA_385_451","Auto_user_S5XL-224-mtDNA_413_487","Auto_user_S5XL-225-mtDNA_PUCE3_SHARP_415_489","Auto_user_S5XL-226-mtDNA_PUCE2_SHARP_414_491","Auto_user_S5XL-230-mtDNA_SHARP1_refait_421_499","Auto_user_S5XL-237-RUN021_mtDNA_459_515","Auto_user_S5XL-251-200414TSCCUSTOMPANEL4POOLSHI-Q_485_548","Auto_user_S5XL-255-200424TSCcustompanel4poolsHi-Q_487_560","Auto_user_S5XL-256-200305CMTL2panel6poolsHi-Q_488_562","Auto_user_S5XL-7-PUCE_NOH_RUN07_68_018","Proton-737-RUN013_NOH_CAP_v2_1575","user_Proton-295-160915CMTL2panel6poolsHi-Q_674"]:

                ToWrite+=str(runEntry[key])+"\t"
                if dicoInit["listPretty"]==True and key in prettyRunKeys:
                    text = str(runEntry[key])
                    if len(text)<=dicoInit["truncateWidth"][1]: row.append(text)
                    else: row.append(text[0:dicoInit["truncateWidth"][1]-3]+"...")
        OUT.write(ToWrite[:-1]+"\n")
        table.append(row)
    OUT.close()
    if dicoInit["listPretty"]==True: print("\n"+tabulate(table, prettyRunKeys, tablefmt="fancy_grid"))

#***** LIST SAMPLE *****#
def listSample(dicoInit):
    sampleKeys = ["name","runid","bc"]
    runKeys = ["instrument","project","ref","target","name","num","date","seqname","seqnum","seqdate","seqstatus","chip","seqkit","libKit"]
    prettyRunKeys = ["instrument","project","name","num","date","ref"]
    table = []
    # Output CSV
    OUT = open(dicoInit["pathOutput"],'w')
    OUT.write("\t".join(sampleKeys+runKeys)+"\n")
    for runEntry in dicoInit["collect_nk_sample"].find().sort("name",pymongo.ASCENDING):
        ToWrite = ""
        row = []
        for key in sampleKeys:
            ToWrite+=str(runEntry[key])+"\t"
            if dicoInit["listPretty"]==True:
                text = str(runEntry[key])
                if len(text)<=dicoInit["truncateWidth"][0]: row.append(text)
                else: row.append(text[0:dicoInit["truncateWidth"][0]-3]+"...")
        findRun = dicoInit["collect_nk_run"].find_one({"_id": runEntry["runid"]})
        for key in runKeys:
            ToWrite+=str(findRun[key])+"\t"
            if dicoInit["listPretty"]==True and key in prettyRunKeys:
                text = str(findRun[key])
                if len(text)<=dicoInit["truncateWidth"][0]: row.append(text)
                else: row.append(text[0:dicoInit["truncateWidth"][0]-3]+"...")
        OUT.write(ToWrite[:-1]+"\n")
        table.append(row)
    OUT.close()
    if dicoInit["listPretty"]==True: print("\n"+tabulate(table, sampleKeys+prettyRunKeys, tablefmt="fancy_grid"))
    




#---------------------------------------------------------------#
#---------------------------------------------------------------#
#                    ADD OBJECTS TO DATABASE                    #
#---------------------------------------------------------------#
#---------------------------------------------------------------#

#***** Add Run *****#
def addRun(dicoInit):
    if not dicoInit["quiet"]: printcolor("\nSub-command: add-run\n","1",dicoInit['blue2'],None,dicoInit['colorBool'])
    dataJson,HashExp,HashBarcode = loadTorrentJson(dicoInit["pathInput"])
    dicoRun = {
               "_id": dataJson["runid"],\
               "seqname": dataJson["expName"],\
               "seqnum" : dataJson["experimentAnalysisSettings"]["experiment"],\
               "seqdate" : HashExp["date"].split("T")[0],\
               "seqstatus" : dataJson["experimentAnalysisSettings"]["status"],\
               "instrument" : HashExp["pgmName"],\
               "chip" : dataJson["chiptype"],\
               "seqkit" : HashExp["sequencekitname"],\
               "libKit" : dataJson["experimentAnalysisSettings"]["libraryKitName"],\
               "project" : dataJson["project"],\
               "name" : dataJson["resultsName"],\
               "num" : HashExp["repResult"],\
               "date" : dataJson["experimentAnalysisSettings"]["date"].split("T")[0],\
               "ref" : dataJson["experimentAnalysisSettings"]["reference"],\
               "target" : os.path.basename(dataJson["experimentAnalysisSettings"]["targetRegionBedFile"]) \
              }
    # Insert if asbent
    findRun = dicoInit["collect_nk_run"].find_one({"_id": dataJson["runid"]})
    if findRun==None:
        insertRun = dicoInit["collect_nk_run"].insert_one(dicoRun)
        if not dicoInit["quiet"]: printcolor("  `"+insertRun.inserted_id+"` insert in NiourK-db (nk_run).\n","0",dicoInit['blue1'],None,dicoInit['colorBool'])
    else:
        if not dicoInit["quiet"]: printcolor("  `"+findRun["_id"]+"` already inserted in NiourK-db (nk_run).\n","0",dicoInit['blue1'],None,dicoInit['colorBool'])
    # Insert samples features
    # only for Auto_PRO-30-mtDNA_87942_74_102 => HashBarcode = {"IonXpress_001":{'barcodes':["IonXpress_001"]},"IonXpress_002":{'barcodes':["IonXpress_002"]},"IonXpress_003":{'barcodes':["IonXpress_003"]},"IonXpress_004":{'barcodes':["IonXpress_004"]},"IonXpress_005":{'barcodes':["IonXpress_005"]},"IonXpress_006":{'barcodes':["IonXpress_006"]},"IonXpress_008":{'barcodes':["IonXpress_008"]},"IonXpress_009":{'barcodes':["IonXpress_009"]}}
    for sample in HashBarcode:
        bc = HashBarcode[sample]['barcodes'][0]
        dicoSample = {
                      "_id" : dataJson["runid"]+"_"+sample.replace(" ","_"),\
                      "runid" : dataJson["runid"],\
                      "name" : sample.replace(" ","_"),\
                      "bc" : bc
                     }
        # Insert if asbent
        findSample = dicoInit["collect_nk_sample"].find_one({"_id" : dataJson["runid"]+"_"+sample.replace(" ","_")})
        if findSample==None:
            insertSample = dicoInit["collect_nk_sample"].insert_one(dicoSample)
            if not dicoInit["quiet"]: printcolor("    `"+insertSample.inserted_id+"`".ljust(25-len(insertSample.inserted_id))+" insert in NiourK-db (nk_sample).\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])
        else:
            if not dicoInit["quiet"]: printcolor("    `"+findSample["_id"]+"`".ljust(25-len(findSample["_id"]))+" already inserted in NiourK-db (nk_sample).\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])



#***** Add VCF *****#
def addVcf(dicoInit):
    if not dicoInit["quiet"]: printcolor("\nSub-command: add-vcf\n","1",dicoInit['blue2'],None,dicoInit['colorBool'])
    # Read input VCF file
    if not dicoInit["quiet"]: printcolor("    Load VCF file\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])    
    vcfReader = vcfpy.Reader.from_path(dicoInit["pathInput"])
    vcfHeader = vcfReader.header
    nkVersion = ""
    refGenome = ""
    lstCaller = []
    # Retrieve caller list
    for headerLine in vcfHeader.lines:
        if headerLine.key=="Nk_version": nkVersion = headerLine.value
        if headerLine.key=="reference": refGenome = os.path.basename(headerLine.value).replace(".fasta","")
        if headerLine.key=="Nk_calls": lstCaller.extend(headerLine.value.split("|"))
    # Check well formated header
    if nkVersion=="": mainUsage(dicoInit,"Missing or empty `Nk_version` tag in input vcf header `"+dicoInit["pathInput"]+"`")
    if refGenome=="": mainUsage(dicoInit,"Any reference genome found in input vcf `"+dicoInit["pathInput"]+"`")
    if len(lstCaller)==0: mainUsage(dicoInit,"Missing or empty `Nk_calls` tag in input vcf header `"+dicoInit["pathInput"]+"`")
    # Output GRCh38 liftovered BED
    if refGenome=="GRCh37":
        BED = open(dicoInit["pathDirTmp"]+"/temp_GRCh37.bed",'w')
        dicoLiftOver = {}
    # Browse variants
    dicoCall = {}
    for record in vcfReader:
        dicoVar = { "nkversion":nkVersion , "call":{}, "filter":{} }
        varId = record.CHROM+"_"+str(record.POS)+"_"+record.REF+"_"+str(record.ALT[0].value)
        dicoVar["af"] = round(float(record.calls[0].data.get('AF')[0]),2)
        # Update depth to Nk_depth
        depthId = record.CHROM+"_"+str(record.POS)
        depth = int(record.calls[0].data.get('DP'))
        if refGenome=="GRCh38":
            dicoInit["collect_nk_depth"].update_one({"_id" : depthId},{"$set": { dicoInit["runId"]+"_"+dicoInit["sampleName"]: depth }})
        # Calling results
        for i in range(len(lstCaller)):
            if record.INFO["CALLFILTER"][0].split("|")[i]!=".":
                dicoVar["filter"][lstCaller[i]] = record.INFO["CALLFILTER"][0].split("|")[i]
            else: 
                dicoVar["call"][lstCaller[i]] = record.INFO["CALLQUAL"][0].split("|")[i]
            if record.INFO["CALLAF"][0].split("|")[i]!=".":
                dicoVar["callaf"] = float(record.INFO["CALLAF"][0].split("|")[i])
        # Prepare Insert/update variant to Nk_var
        if refGenome=="GRCh38" or refGenome=="chrM": 
            dicoCall[varId] = dicoVar
        # Write to temp GRCh37 BED (0-based position)
        else:
            bedLine = record.CHROM+"\t"+str(record.POS-1)+"\t"+str(record.POS)+"\t"+record.REF+"#"+record.ALT[0].value+"#"+str(dicoVar["af"])+"#"+str(depth)+"#"+json.dumps(dicoVar).replace(" ","")+"\n"
            BED.write(bedLine)
    # For GRCh38 direct Insert/update variant to Nk_var
    if refGenome=="GRCh38" or refGenome=="chrM":
        if not dicoInit["quiet"]: printcolor("    Insert/Update NiourK-db (nk_var)\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])
        lstVarId = list(dicoCall.keys())
        for i in tqdm(range(len(lstVarId)),ncols=30,leave=False,bar_format="      {percentage:3.0f}%|{bar}|"):
            findVar = dicoInit["collect_nk_var"].find_one({"_id": lstVarId[i]})
            if findVar==None:
                insertRun = dicoInit["collect_nk_var"].insert_one({ "_id":lstVarId[i]})
            dicoInit["collect_nk_var"].update_one({"_id" : lstVarId[i]},{"$set": { dicoInit["runId"]+"_"+dicoInit["sampleName"]: dicoCall[lstVarId[i]] }})
        if not dicoInit["quiet"]: printcolor("      "+str(len(lstVarId))+" variants inserted in NiourK-db (nk_var).\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])
    # Liftover
    else:
        BED.close()
        if not dicoInit["quiet"]: printcolor("    Liftover to GRCh38\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])    
        cmd_liftover = dicoInit["pathLiftover"]+" "+dicoInit["pathDirTmp"]+"/temp_GRCh37.bed "+dicoInit["pathLiftChain"]+" "+dicoInit["pathDirTmp"]+"/temp_GRCh38.bed "+dicoInit["pathDirTmp"]+"/temp_unmapped.bed 2>/dev/null"
        os.system(cmd_liftover)
        # Read GRCh38 variants
        if not dicoInit["quiet"]: printcolor("    Insert/Update NiourK-db (nk_var)\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])    
        BED = open(dicoInit["pathDirTmp"]+"/temp_GRCh38.bed",'r')
        lstLines = BED.read().split("\n")
        BED.close()
        for i in tqdm(range(len(lstLines)),ncols=30,leave=False,bar_format="      {percentage:3.0f}%|{bar}|"):
            if lstLines[i]!="":
                splitLine = lstLines[i].split("\t")
                splitField = splitLine[3].split("#")
                splitLine
                ref = splitField[0]
                alt = splitField[1]
                varId = splitLine[0]+"_"+str(int(splitLine[1])+1)+"_"+ref+"_"+alt
                dicoVar = json.loads(splitField[4])
                dicoVar["af"] = float(splitField[2])
                # Update depth to Nk_depth
                depthId = splitLine[0]+"_"+str(int(splitLine[1])+1)
                depth = int(splitField[3])
                dicoInit["collect_nk_depth"].update_one({"_id" : depthId},{"$set": { dicoInit["runId"]+"_"+dicoInit["sampleName"]: depth }})
                # Insert/update variant to Nk_var
                findVar = dicoInit["collect_nk_var"].find_one({"_id": varId})
                if findVar==None:
                    insertRun = dicoInit["collect_nk_var"].insert_one({ "_id":varId})
                dicoInit["collect_nk_var"].update_one({"_id" : varId},{"$set": { dicoInit["runId"]+"_"+dicoInit["sampleName"]: dicoVar }})
        if not dicoInit["quiet"]: printcolor("      "+str(len(lstLines)-1)+" variants inserted in NiourK-db (nk_var).\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])



#***** Add Nk-Sample *****#
def addNkSample(dicoInit):
    sampleID = dicoInit["runId"]+"_"+dicoInit["sampleName"]
    if not dicoInit["quiet"]:
        printcolor("\nSub-command: add-nksample\n","1",dicoInit['blue1'],None,dicoInit['colorBool'])
        printcolor("  sampleId: "+sampleID+"\n","1",dicoInit['blue2'],None,dicoInit['colorBool'])
    # Load JSON
    if not dicoInit["quiet"]: printcolor("    Load NkSample JSON file\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])    
    JSON = open(dicoInit["pathInput"],'r')
    dataSampleJson = json.load(JSON)
    JSON.close()
    # Output GRCh38 liftovered BED
    BED = open(dicoInit["pathDirTmp"]+"/temp_GRCh37.bed",'w')
    dicoLiftOver = {}
    # Browse sample variants
    if not dicoInit["quiet"]: printcolor("    Read NkSample variants\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])    
    for key in dataSampleJson:
        dicoCall = {"version":"1.7", "call":{}, "filter":{}}
        # Haplo features
        if key=="haplo": # H1b2 (0.64)
            haplo = dataSampleJson[key].split(" ")[0]
            haploscore = float(dataSampleJson[key].split("(")[1].replace(")",""))
            dicoInit["collect_nk_sample"].update_one({"_id" : sampleID},{"$set": { "haplo":haplo, "haploscore":haploscore }})
        elif key!="haplo_mitomap":
            splitKey = key.split("_")
            chrom = splitKey[0]
            pos = int(splitKey[1]) #(1-based position)
            ref = splitKey[2]
            alt = splitKey[3]
            varId = key
            af = dataSampleJson[key]["allele_freq"]
            sb = dataSampleJson[key]["strand_bias"]
            depth = dataSampleJson[key]["dico_cov"]["reads_all"]
            for caller in ["GATKu","LoFreq","platypus","SNVer","TSVC","VarScan"]:
                if caller in dataSampleJson[varId]:
                    dicoCall["call"][caller.lower()] = dataSampleJson[varId][caller]
                elif caller+"_pseudo" in dataSampleJson[varId]:
                    dicoCall["call"][caller.lower()] = dataSampleJson[varId][caller+"_pseudo"]
                elif caller+"filtered" in dataSampleJson[varId]:
                    dicoCall["filter"][caller.lower()] = dataSampleJson[varId][caller+"filtered"]
        # Write to temp GRCh37 BED (0-based position)
        bedLine = chrom+"\t"+str(pos-1)+"\t"+str(pos)+"\t"+ref+"#"+alt+"#"+str(af)+"#"+str(sb)+"#"+str(depth)+"#"+json.dumps(dicoCall).replace(" ","")+"\n"
        BED.write(bedLine)
    BED.close()
    # Liftover
    if not dicoInit["quiet"]: printcolor("    Liftover to GRCh38\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])    
    cmd_liftover = dicoInit["pathLiftover"]+" "+dicoInit["pathDirTmp"]+"/temp_GRCh37.bed "+dicoInit["pathLiftChain"]+" "+dicoInit["pathDirTmp"]+"/temp_GRCh38.bed "+dicoInit["pathDirTmp"]+"/temp_unmapped.bed 2>/dev/null"
    os.system(cmd_liftover)
    # Read GRCh38 variants
    if not dicoInit["quiet"]: printcolor("    Insert/Update NiourK-db (nk_var)\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])    
    BED = open(dicoInit["pathDirTmp"]+"/temp_GRCh38.bed",'r')
    lstLines = BED.read().split("\n")
    BED.close()
    for i in tqdm(range(len(lstLines)),ncols=30,leave=False,bar_format="      {percentage:3.0f}%|{bar}|"):
        if lstLines[i]!="":
            splitLine = lstLines[i].split("\t")
            splitField = splitLine[3].split("#")
            ref = splitField[0]
            alt = splitField[1]
            varId = splitLine[0]+"_"+str(int(splitLine[1])+1)+"_"+ref+"_"+alt
            dicoVar = json.loads(splitField[5])
            dicoVar["af"] = float(splitField[2])
            dicoVar["sb"] = float(splitField[3])
            # Add depth to Nk_depth if absent
            depthId = splitLine[0]+"_"+str(int(splitLine[1])+1)
            depth = int(splitField[4])
            findDepth = dicoInit["collect_nk_depth"].find_one({"_id": depthId})
            if not findDepth:
                insertDepth = dicoInit["collect_nk_depth"].insert_one({"_id" : depthId},{"$set": { sampleID: depth }})
            elif not sampleID in findDepth:
                dicoInit["collect_nk_depth"].update_one({"_id" : depthId},{"$set": { sampleID: depth }})
            # Insert/update variant to Nk_var
            findVar = dicoInit["collect_nk_var"].find_one({"_id": varId})
            if findVar==None:
                insertRun = dicoInit["collect_nk_var"].insert_one({ "_id":varId})
            dicoInit["collect_nk_var"].update_one({"_id" : varId},{"$set": { sampleID: dicoVar }})
    if not dicoInit["quiet"]: printcolor("      "+str(len(lstLines)-1)+" variants inserted in NiourK-db (nk_var).\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])



#***** Add Depth *****#
def addDepth(dicoInit):
    sampleID = dicoInit["runId"]+"_"+dicoInit["sampleName"]
    if not dicoInit["quiet"]:
        printcolor("\nSub-command: add-depth\n","1",dicoInit['blue1'],None,dicoInit['colorBool'])
        printcolor("  sampleId: "+sampleID+"\n","1",dicoInit['blue2'],None,dicoInit['colorBool'])
    # Check if sample already inserted:
    findInsertStatut = dicoInit["collect_nk_depth"].find_one({"_id":"insertsample"})
    boolInsert = True
    if findInsertStatut==None:
        dicoInit["collect_nk_depth"].insert_one({"_id" : "insertsample", "lstrunid": []})
    elif sampleID in findInsertStatut["lstrunid"]:
        if not dicoInit["quiet"]: printcolor("      already inserted in NiourK-db (nk_depth).\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])
        boolInsert = False
    if boolInsert:
        # Insert GRCH38 sample positions
        if not dicoInit["quiet"]: printcolor("    Load depth BED file\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])            
        # Accept bed and bed.gz
        if ".gz" in dicoInit["pathInput"]: BED = gzip.open(dicoInit["pathInput"], 'rb')
        else: BED = open(dicoInit["pathInput"],'r')
        lstLines = BED.read().split("\n")
        BED.close()
        if not dicoInit["quiet"]: printcolor("    Insert/Update NiourK-db (nk_depth)\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])            
        for i in tqdm(range(len(lstLines)),ncols=30,leave=True,bar_format="  {percentage:3.0f}%|{bar}|"):
            try:
                splitLine = lstLines[i].split("\t")
                chrom = splitLine[0]
                start = int(splitLine[1])+1 # 1-based in nk_depth
                end = int(splitLine[2])
                depth = int(splitLine[3]) # zero depth was already filtered from input bed
                for pos in range(start,end+1,1):
                    depthId = chrom+"_"+str(pos)
                    # upsert => creates a new document if no documents match the filter.
                    dicoInit["collect_nk_depth"].update_one({"_id" : depthId},{"$set": { sampleID: depth }},upsert=True)
            except: pass
        if not dicoInit["quiet"]: printcolor("      "+str(len(lstLines))+" positions depth inserted in NiourK-db (nk_depth).\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])
        # add to InsertStatut
        dicoInit["collect_nk_depth"].update_one({"_id" : "insertsample"}, {'$push': {'lstrunid': sampleID}})



#***** Get variant features *****#
def getVariant(dicoInit):
    varID = dicoInit["varInput"]
    if not dicoInit["quiet"]:
        printcolor("\nSub-command: get-variant\n","1",dicoInit['blue1'],None,dicoInit['colorBool'])
        printcolor("  varId    : "+varID+"\n","1",dicoInit['blue2'],None,dicoInit['colorBool'])
    # Search variant in nk_var
    findVar = dicoInit["collect_nk_var"].find_one({"_id":varID})
    if not findVar: printcolor("  find   : 0 occurence in NiourK-db\n","1",dicoInit['red'],None,dicoInit['colorBool'])
    else:
        countSamples = len(findVar)-1
        printcolor("  find     : "+str(countSamples)+" occurences in NiourK-db\n","1",dicoInit['green'],None,dicoInit['colorBool'])
        # Search number of samples with depth at this position in nk_depth
        depthId = "_".join(varID.split("_")[:-2])
        nbOverlapSample = 1
        findDepth = dicoInit["collect_nk_depth"].find_one({"_id":depthId})
        if findDepth:
            for key in findDepth:
                if key!="_id" and findDepth[key]>=dicoInit["mindepth"]:
                    nbOverlapSample+=1
        # Compute variant DB frequency
        varDBfreq = round(((countSamples*100)/nbOverlapSample),1)
        printcolor("  DB depth : "+str(nbOverlapSample)+" samples","1",dicoInit['white'],None,dicoInit['colorBool'])
        printcolor(" (depth>="+str(dicoInit["mindepth"])+")\n","0",dicoInit['white'],None,dicoInit['colorBool'])
        printcolor("  DB freq  : "+str(varDBfreq)+"% in NiourK-db\n","1",dicoInit['white'],None,dicoInit['colorBool'])
        # Display samples summary table
        header = ["sample","af","sb","call","filter","version"]
        table = []
        for key in findVar:
            if key!="_id":
                row = [key,str(findVar[key]['af']),str(findVar[key]['sb'])]
                callCell = " "
                filterCell = " "
                for tool in findVar[key]['call']: callCell+=tool+" ("+str(findVar[key]['call'][tool])+")"+"\n"
                for tool in findVar[key]['filter']: callCell+=tool+" ("+str(findVar[key]['filter'][tool])+")"+"\n"
                row.append(callCell[:-1])
                row.append(filterCell[:-1])
                row.append(findVar[key]['version'])
                table.append(row)
        for line in tabulate(table, header, tablefmt="fancy_grid").split("\n"):
            printcolor("  "+line+"\n","0",dicoInit['grey1'],None,dicoInit['colorBool'])
