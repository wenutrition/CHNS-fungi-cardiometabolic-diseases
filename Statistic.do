
*-----1-dataset describe---------------

*---basic characters--

use "D:\真菌分析\data_pre\fugi_15all_g_diversity.dta",clear 
merge m:m SampleID using  "D:\真菌分析\data_pre\2015_all_phenotype.dta"
keep if _merge==3 // 10695

drop _merge pielou_evenness-g761
gen dataset="data1"
save "D:\真菌分析\data_pre\data1_meta.dta",replace

use  "D:\真菌分析\data_pre\2015_all_phenotype.dta",clear
rename idind Idind
merge m:m Idind using "D:\真菌分析\data_pre\repeatDMM_cluster.dta",force
keep if time==15
drop _merge
save "D:\真菌分析\data_pre\data2_meta.dta",replace
gen dataset="data2"
append using "D:\真菌分析\data_pre\data1_meta.dta"
gen veg=lveg+dveg
qui foreach var of varlist hyper_med xinjigengse zhongfeng cancer smoke alcohol{
     replace `var'=0 if `var'==9
}
table1, by(dataset) vars(age contn \ sex cat \ city cat\ city_score contn\pet cat\education cat\ marrige cat\income contn\BMI contn\ wc contn \ hpc contn \ SBP contn \ DBP contn  \ HbA1c contn \ glucose contn \ ins contn \HDL_C contn \ LDL_C contn\  tc contn \tg contn \t2d cat  \ pret2d cat\ hyp cat \dys cat\ gut_disease cat\ fuxie cat\ xinjigengse cat\ zhongfeng cat\ cancer cat\ hyper_med cat\diabetes_med cat\ antibiotic_current cat\ antibiotic_6month cat\ probiotics cat\ kangyan_med cat\ kangsuan_med cat\ weisuan_med cat\changdao_shoushu cat\  smoke cat\  alcohol cat\ MET contn\ wheat contn \ rice contn \dveg contn \lveg contn \veg contn\saltveg contn \fruit contn \nuts contn \pork contn\poultry contn \milk contn \egg contn \fish contn\cereoth contn \tuber contn \dryleg contn \legpro contn \orgmeat contn \pastes contn\othmeat contn \cake contn \sugar contn \OIL_VEG contn \OIL_ANI contn \salt contn \sauce contn \others contn )format(%8.2f) onecol test pdp(2) saving (D:\真菌分析\data_des.xlsx, replace)

*---sequencing depth--

*-all-

import excel "D:\真菌分析\data_des.xlsx", sheet("all") firstrow clear 
merge 1:1 SampleID using  "D:\真菌分析\data_pre\fugi_15all_g_diversity.dta"
keep if _merge==3
drop _merge
merge m:m SampleID using  "D:\真菌分析\data_pre\2015_all_phenotype.dta"
su count,de 

*-repeat-

import excel "D:\真菌分析\data_des.xlsx", sheet("repeat_sequence") firstrow clear 
merge m:m SampleID using  "D:\真菌分析\data_pre\data2_meta.dta"
keep if _merge==3 //不能直接匹配

*----overlaped genera----

*---2015 all and repeated fungi data mapping--

import excel "D:\真菌分析\rawdata\raw_fungi_mapping.xlsx", sheet("all") firstrow clear 
save  "D:\真菌分析\data_pre\all_fungi_mapping.dta",replace
merge 1:1 g using  "D:\真菌分析\data_pre\raw_fungi_mapping.dta"
keep if _merge==3 // 558 matched (202 from all, 74 from repeated without matched)
drop _merge
save  "D:\真菌分析\data_pre\fungi_mapping_all_repeat.dta",replace

*---all filted---

import excel "D:\真菌分析\data_des.xlsx", sheet("Sheet4") firstrow clear 
merge m:m code_all using "D:\真菌分析\data_pre\all_fungi_mapping.dta"
keep if _merge==3
drop _merge
save "D:\真菌分析\data_pre\all_filtergenus.dta",replace

*---repeat filted---

import excel "D:\真菌分析\data_des.xlsx", sheet("Sheet4") firstrow clear 
merge m:m code_repeat using "D:\真菌分析\data_pre\raw_fungi_mapping.dta"
keep if _merge==3
save "D:\真菌分析\data_pre\repeat_filtergenus.dta",replace

*--overlaped--

use "D:\真菌分析\data_pre\all_filtergenus.dta",clear
append using "D:\真菌分析\data_pre\repeat_filtergenus.dta"
duplicates drop g,force

*--matched the discovery filtered genus (see:D:\真菌分析\rawdata\raw_fungi_mapping)

import excel "D:\真菌分析\rawdata\raw_fungi_mapping.xlsx", sheet("Sheet3") firstrow clear 
merge 1:1 code_all using  "D:\真菌分析\data_pre\fungi_mapping_all_repeat.dta"
keep if _merge==3
drop _merge
rename code_all micro
merge 1:1 micro using "D:\真菌分析\data_pre\all_g_mapping.dta"
keep if _merge==3
drop _merge
save "D:\真菌分析\data_pre\repeat_g_mapping.dta",replace


*----------------2-Repeat datasets compare-----

*--2.1-----DMM clusters analyses------

*--2.1.1--DMM clusters and diversity--

import excel "D:\真菌分析\Final results\Repeated data\DMM_update\All_cluster.xlsx", sheet("All_cluster") firstrow clear 
tab cluster time
drop SampleID
save "D:\真菌分析\data_pre\repeatDMM_cluster.dta",replace

use "D:\真菌分析\data_pre\repeatDMM_cluster.dta",clear
tab time
keep if time==15
merge m:m Idind using "D:\真菌分析\data_pre\fugi_15_diversity.dta",force
keep if _merge==3
sort cluster
by cluster: su shannon observed_features pielou_evenness faith_pd
gr box observed_features,over(cluster)

use "D:\真菌分析\data_pre\repeatDMM_cluster.dta",clear
tab time
keep if time==18
merge m:m Idind using "D:\真菌分析\data_pre\fugi_18_diversity.dta",force
keep if _merge==3
sort cluster
by cluster: su shannon observed_features pielou_evenness faith_pd
gr box observed_features,over(cluster)
table1, by(cluster) vars(shannon contn \ observed_features contn\ pielou_evenness contn \ faith_pd contn  )format(%8.2f) onecol test pdp(2) saving (D:\真菌分析\cluster_phenotype.xlsx, replace)

*--2.1.1--DMM clusters and phenotype-

*--statistic analysis--

*----15----

use  "D:\真菌分析\data_pre\2015_all_phenotype.dta",clear
rename idind Idind
merge m:m Idind using "D:\真菌分析\data_pre\repeatDMM_cluster.dta",force
keep if time==15
*keep SampleID sex
*save "D:\真菌分析\data_pre\repeatsex.dta",replace
tab cluster
table1, by(cluster) vars(age contn \ sex cat \BMI contn\ wc contn \ hpc contn \ district cat\ city cat\education cat\ marrige cat\ smoke cat\  alcohol cat\ income contn\ city_score contn\ SBP contn \ DBP contn \ HDL_C contn \ LDL_C contn \ ins contn \ HbA1c contn \ glucose contn \ DBP contn \tg contn \  tc contn \t2d cat  \  hyper_med cat\ xinjigengse cat\ zhongfeng cat\ cancer cat\ diabetes_med cat\ hyp cat \dys cat\ pret2d cat\ gut_disease cat\ fuxie cat\ antibiotic_current cat\ antibiotic_6month cat\ probiotics cat\ kangyan_med cat\ kangsuan_med cat\ weisuan_med cat\changdao_shoushu cat\ yogurt_drink cat\ pet cat \ famine cat\ wheat contn \ rice contn \cereoth contn \tuber contn \dryleg contn \legpro contn \dveg contn \lveg contn \saltveg contn \fruit contn \nuts contn \pork contn \othmeat contn \orgmeat contn \poultry contn \milk contn \egg contn \fish contn \OIL_VEG contn \OIL_ANI contn \cake contn \sugar contn \salt contn \sauce contn \others contn \pastes contn)format(%8.2f) onecol test pdp(2) saving (D:\真菌分析\cluster_phenotype.xlsx, replace)

*----18----

use  "D:\真菌分析\data_pre\2018_all_phenotype.dta",clear
rename idind Idind
drop _merge
merge m:m Idind using "D:\真菌分析\data_pre\repeatDMM_cluster.dta",force
keep if time==18
tab cluster
table1, by(cluster) vars(age contn \ sex cat \BMI contn\ wc contn \ hpc contn \ district cat\ city cat\education cat\ marrige cat\ smoke cat\  alcohol cat\ income contn\ city_score contn\ SBP contn \ DBP contn \ HDL_C contn \ LDL_C contn \ ins contn \ HbA1c contn \ glucose contn \ DBP contn \tg contn \  tc contn \t2d cat  \  hyper_med cat\ xinjigengse cat\ zhongfeng cat\ cancer cat\ diabetes_med cat\ hyp cat \dys cat\ pret2d cat\ gut_disease cat\ fuxie cat\ antibiotic_current cat\ antibiotic_6month cat\ probiotics cat\ kangyan_med cat\ kangsuan_med cat\ weisuan_med cat\changdao_shoushu cat\ yogurt_drink cat\ pet cat \ famine cat\ wheat contn \ rice contn \cereoth contn \tuber contn \dryleg contn \legpro contn \dveg contn \lveg contn \saltveg contn \fruit contn \nuts contn \pork contn \othmeat contn \orgmeat contn \poultry contn \milk contn \egg contn \fish contn \OIL_VEG contn \OIL_ANI contn \cake contn \sugar contn \salt contn \sauce contn \others contn \pastes contn)format(%8.2f) onecol test pdp(2) saving (D:\真菌分析\cluster_phenotype_18.xlsx, replace)

*--phenotype predict analysis---

*---2015--

use  "D:\真菌分析\data_pre\2015_all_phenotype.dta",clear
rename idind Idind
merge m:m Idind using "D:\真菌分析\data_pre\repeatDMM_cluster.dta",force
keep if time==15
gen veg=dveg+lveg
sort cluster
by cluster:su city_score

*-all data-

keep SampleID age t2d city education marrige veg smoke alcohol hyper_med xinjigengse-cancer MET city_score income hpc-sex BMI-dys rice-antibiotic_current antibiotic_6month probiotics-pet cluster

*--all diet related data--

keep SampleID smoke alcohol MET wheat-pastes veg cluster

preserve
replace cluster=0 if cluster!=1
rename cluster outcome
tab outcome
save "D:\真菌分析\data_pre\15cluster_predict1.dta",replace
restore

preserve
gen outcome=0
replace outcome=1 if cluster==2
drop cluster
tab outcome
save "D:\真菌分析\data_pre\15cluster_predict2.dta",replace
restore

preserve
gen outcome=0
replace outcome=1 if cluster==3
drop cluster
tab outcome
save "D:\真菌分析\data_pre\15cluster_predict3.dta",replace
restore

preserve
gen outcome=0
replace outcome=1 if cluster==4
drop cluster
tab outcome
save "D:\真菌分析\data_pre\15cluster_predict4.dta",replace
restore


preserve
keep if cluster==3| cluster==4
gen outcome=0
replace outcome=1 if cluster==4
drop cluster
tab outcome
save "D:\真菌分析\data_pre\15cluster_predict5.dta",replace // AUC=0.85 (all data)
restore

*---2018--

use  "D:\真菌分析\data_pre\2018_all_phenotype.dta",clear
rename idind Idind
drop _merge
merge m:m Idind using "D:\真菌分析\data_pre\repeatDMM_cluster.dta",force
keep if time==18
rename Idind SampleID
keep SampleID age t2d city education marrige smoke alcohol hyper_med xinjigengse-cancer MET city_score income hpc-sex BMI-dys rice-antibiotic_current antibiotic_6month probiotics-pet cluster

preserve
keep if cluster==3| cluster==4
gen outcome=0
replace outcome=1 if cluster==4
drop cluster
tab outcome
save "D:\真菌分析\data_pre\18cluster_predict5.dta",replace // AUC=0.85
restore

*----2.1.2--DMM clusters transport

use "D:\真菌分析\data_pre\repeatDMM_cluster.dta",clear
keep if time==15
rename cluster cluster_15
save "D:\真菌分析\data_pre\repeatDMM_15cluster.dta",replace

use "D:\真菌分析\data_pre\repeatDMM_cluster.dta",clear
keep if time==18
rename cluster cluster_18
save "D:\真菌分析\data_pre\repeatDMM_18cluster.dta",replace

use "D:\真菌分析\data_pre\repeatDMM_15cluster.dta",clear
merge 1:1 Idind using "D:\真菌分析\data_pre\repeatDMM_18cluster.dta"
keep if _merge==3
drop _merge time
save "D:\真菌分析\data_pre\repeatDMM_indcluster.dta",replace
export excel using "D:\真菌分析\Final results\Repeated data\DMM\ind_repeatcluster.xlsx",  firstrow(variables) sheet("Sheet1") sheetreplace  

use "D:\真菌分析\data_pre\repeatDMM_indcluster.dta",clear
sort cluster_15 cluster_18
egen group=group(cluster_15 cluster_18)
tab group 
save "D:\真菌分析\data_pre\repeatDMM_clusterchange.dta",replace

sort cluster_15
by cluster_15:tab cluster_18

*---transport compare---

*---baseline

use  "D:\真菌分析\data_pre\2015_all_phenotype.dta",clear
rename idind Idind
merge m:m Idind using "D:\真菌分析\data_pre\repeatDMM_cluster.dta",force
keep if time==15
drop _merge
merge 1:1 Idind using "D:\真菌分析\data_pre\repeatDMM_clusterchange.dta",force


keep SampleID group age t2d city education marrige smoke alcohol hyper_med xinjigengse-cancer MET city_score income hpc-sex BMI-dys rice-antibiotic_current antibiotic_6month probiotics-pet cluster
tab group,missing
sort group
by group :su age SBP DBP glucose

preserve
keep if group==3 | group==4
gen outcome=0
replace outcome=1 if group==4
drop cluster group
tab outcome
save "D:\真菌分析\data_pre\15cluster_predict1_.dta",replace
table1, by(outcome) vars(age contn \ sex cat \BMI contn\ wc contn \ hpc contn \  city cat\education cat\ marrige cat\ smoke cat\  alcohol cat\ income contn\ city_score contn\ SBP contn \ DBP contn \ HDL_C contn \ LDL_C contn \ ins contn \ HbA1c contn \ glucose contn \ DBP contn \tg contn \  tc contn \t2d cat  \  hyper_med cat\ xinjigengse cat\ zhongfeng cat\ cancer cat\ diabetes_med cat\ hyp cat \dys cat\ gut_disease cat\ fuxie cat\ antibiotic_current cat\ antibiotic_6month cat\ probiotics cat\ kangyan_med cat\ kangsuan_med cat\ weisuan_med cat\changdao_shoushu cat\ yogurt_drink cat\ pet cat\ wheat contn \ rice contn \cereoth contn \tuber contn \dryleg contn \legpro contn \dveg contn \lveg contn \saltveg contn \fruit contn \nuts contn \pork contn \othmeat contn \orgmeat contn \poultry contn \milk contn \egg contn \fish contn \OIL_VEG contn \OIL_ANI contn \cake contn \sugar contn \salt contn \sauce contn \others contn \pastes contn)format(%8.2f) onecol test pdp(2) saving (D:\真菌分析\cluster_phenotype_1_34_baseline.xlsx, replace)
restore

*---follow

use  "D:\真菌分析\data_pre\2015_all_phenotype.dta",clear
rename idind Idind
merge m:m Idind using "D:\真菌分析\data_pre\repeatDMM_cluster.dta",force
keep if time==18
drop _merge
merge 1:1 Idind using "D:\真菌分析\data_pre\repeatDMM_clusterchange.dta",force
keep SampleID group age t2d city education marrige smoke alcohol hyper_med xinjigengse-cancer MET city_score income hpc-sex BMI-dys rice-antibiotic_current antibiotic_6month probiotics-pet cluster

preserve
keep if group==3 | group==4
gen outcome=0
replace outcome=1 if group==4
drop cluster group
tab outcome
save "D:\真菌分析\data_pre\15cluster_predict1_f.dta",replace
restore
\repeatDMM_cluster.dta",clear
tab time

*-------GEE modle for fungi-diseases associations dataset prepare--

use "D:\真菌分析\data_pre\repeatDMM_cluster.dta",clear
gen ID=time+Idind
save "D:\真菌分析\data_pre\repeatDMM_cluster_.dta",replace

use  "D:\真菌分析\data_pre\2018_all_phenotype.dta",clear
gen time=18
append using "D:\真菌分析\data_pre\2015_all_phenotype.dta"
tab time,missing
replace time=15 if time==.
rename idind Idind
gen ID=time+Idind
merge m:m Idind using "D:\真菌分析\data_pre\repeatDMM_cluster_.dta",force
keep if _merge==3
tab  cluster,missing
spearman BMI city_score
gen veg=lveg+dveg
sort Idind
egen code_id=group(Idind)
replace time=0 if time==15
replace time=1 if time==18
sort time
qui foreach var of varlist age MET income city_score hpc-tc BMI SBP DBP veg rice-pastes collect_month {
     by time: egen `var'_mean= mean(`var')
	 by time: replace `var'=`var'_mean if `var'==.
}

foreach phenotype of varlist t2d dys hyp marrige  smoke alcohol xinjigengse zhongfeng cancer sex hyper_med diabetes_med antibiotic_current antibiotic_6month probiotics kangyan_med kangsuan_med weisuan_med yogurt_drink pet hyper_med city  hyp dys pret2d fuxie gut_disease changdao_shoushu{
	 replace `phenotype'=0 if `phenotype'==9 |`phenotype'==. 
}
misstable sum  t2d dys hyp age sex   smoke alcohol city income fruit  pork veg milk fish antibiotic_current MET
save   "D:\真菌分析\data_pre\cluster_gee_.dta",replace

*----2.1.3-----DMM clusters and metabolism----

*------Predict clusters dataset-----

use "D:\CNHS\data\phenotype data\2015update_20220426\chns_2015_westlake.dta",clear
duplicates drop Idind,force // 13315
merge m:m Idind using  "D:\真菌分析\data_pre\metabolism_15.dta"
keep if _merge==3
drop if age==. // 996 obs
drop _merge
merge m:m Idind using "D:\真菌分析\data_pre\repeatDMM_15cluster.dta"
keep if _merge==3

preserve
keep Idind cluster 
save "D:\真菌分析\data_pre\metabolism_idcluster.dta",replace
restore

tab cluster
sort cluster
by cluster:su MEDN0390

keep SampleID cluster CMPF_N-mwxq04_N

preserve
keep if cluster==4 | cluster==3
gen outcome=0
replace outcome=1 if cluster==4
tab outcome
drop cluster
keep SampleID outcome MEDN1863	MEDN1445	MEDP1693	MEDN2099	MEDN1004	MEDP0052	MEDP1003	MEDP1126	MEDN1231	MEDP1424	MEDN0846	MEDP1404	MEDN0551	MEDP1919	MEDP0401	MEDP1954	MEDN1740	MEDN1753	MEDN1699	MEDN2197	MEDP2266	MEDP1255	MEDN1092	MEDP1670	MEDN1193	MEDN0071	MEDN1029	MEDP2132	MEDP0022	MEDN0729	MEDN1442	MEDN0431	MEDN1854	MEDP0006	MEDP0174	MEDP2325	MEDP0313	MEDN1708	MEDN0065	MEDP0325	MEDN0579	MEDN0754	MEDP2672	MEDN0105	MEDP2494	MEDP2088	MEDP2131	MEDP1321	MEDN0589	MEDN1440	MEDN0261	MEDN1441	MEDN1444	MEDP2633	MEDN0807	MEDP1061	MEDN1279	MEDN0364	MEDN0536	MEDN0340	MEDP0251	MEDN0100	MEDP1778	MEDP0026	MEDN0478	MEDN0682	MEDN1686	MEDN0011	MEDN2103	MEDN0661	MEDN0062	MEDN0072	MEDN1576	MEDP0389	MEDN0826	MEDN0398	MEDN0227	MEDN1221	MEDN0463	MEDP0453	MEDP2431	MEDP2223	MEDN0616	MEDP1045	MEDN1246	MEDP1928	MEDP0598	MEDN1861	MEDP1318	MEDP1319	MEDN0403	MEDP0296	MEDN0115	MEDN1965	MEDN1902	MEDN1126	MEDN1723	MEDN1077	MEDP2176	MEDN2111	MEDN1532	MEDN1943	MEDP1803	MEDN0658	MEDN0769	MEDN1613	MEDN1853	MEDP1885	MEDP1201	MEDP1002	MEDP1242	MEDP1488	MEDN0179	MEDN1264	MEDN1265	MEDN0348	MEDN2142	MEDN0511	MEDN1927	MEDP2334	MEDN1862	MEDP2283	MEDN1156	MEDN1640	MEDN0335	MEDN0201	MEDN1566	MEDP0147	MEDN2048	MEDN1796	MEDN1906	MEDN1448	MEDP1691	MEDP1695	MEDP2143	MEDN0659	MEDN0362	MEDN0758	MEDP0297	MEDN1886	MEDN1415	MEDP2339	MEDN1935	MEDN0745	MEDP1663	MEDN0748	MEDN0856	MEDP0084	MEDP0831	MEDP1657	MEDN1311	MEDN0429	MEDP1753	MEDP1434	MEDP1857	MEDP0336	MEDP1326	MEDP1702	MEDP1171	MEDN0323	MEDN1213	MEDP1391	MEDN2182	MEDN1428	MEDN0128	MEDP1407	MEDP0307	MEDP1251	MEDP1333	MEDP1334	MEDP1389	MEDN1413	MEDP1701	MEDP1339	MEDN1932	MEDP1399	MEDN0294	MEDN1643	MEDN0383	MEDN0752	MEDN0750	MEDN0756	MEDP1659	MEDP1165	MEDP1146	MEDN1690	MEDN0290	MEDP1658	MEDN0375	MEDN1416	MEDP1322	MEDP1412	MEDN0245	MEDP1418	MEDP1692	MEDP0494	MEDN0759	MEDN1597	MEDN1691	MEDN1648	MEDP1518	MEDN1269	MEDN1270	MEDN1278	MEDN1277	MEDN0390	MEDN1069	MEDN1688	MEDN1840	MEDN0378	MEDP1901	MEDP1428
sort outcome
save "D:\真菌分析\data_pre\metabolism_cluster34raw.dta",replace
restore


foreach metabolism of varlist CMPF_N-mwxq04_N {
	gen ln_`metabolism'=log(`metabolism') 
	drop `metabolism'
	rename ln_`metabolism' `metabolism'
}

preserve
keep SampleID cluster MEDN1863	MEDN1445	MEDP1693	MEDN2099	MEDN1004	MEDP0052	MEDP1003	MEDP1126	MEDN1231	MEDP1424	MEDN0846	MEDP1404	MEDN0551	MEDP1919	MEDP0401	MEDP1954	MEDN1740	MEDN1753	MEDN1699	MEDN2197	MEDP2266	MEDP1255	MEDN1092	MEDP1670	MEDN1193	MEDN0071	MEDN1029	MEDP2132	MEDP0022	MEDN0729	MEDN1442	MEDN0431	MEDN1854	MEDP0006	MEDP0174	MEDP2325	MEDP0313	MEDN1708	MEDN0065	MEDP0325	MEDN0579	MEDN0754	MEDP2672	MEDN0105	MEDP2494	MEDP2088	MEDP2131	MEDP1321	MEDN0589	MEDN1440	MEDN0261	MEDN1441	MEDN1444	MEDP2633	MEDN0807	MEDP1061	MEDN1279	MEDN0364	MEDN0536	MEDN0340	MEDP0251	MEDN0100	MEDP1778	MEDP0026	MEDN0478	MEDN0682	MEDN1686	MEDN0011	MEDN2103	MEDN0661	MEDN0062	MEDN0072	MEDN1576	MEDP0389	MEDN0826	MEDN0398	MEDN0227	MEDN1221	MEDN0463	MEDP0453	MEDP2431	MEDP2223	MEDN0616	MEDP1045	MEDN1246	MEDP1928	MEDP0598	MEDN1861	MEDP1318	MEDP1319	MEDN0403	MEDP0296	MEDN0115	MEDN1965	MEDN1902	MEDN1126	MEDN1723	MEDN1077	MEDP2176	MEDN2111	MEDN1532	MEDN1943	MEDP1803	MEDN0658	MEDN0769	MEDN1613	MEDN1853	MEDP1885	MEDP1201	MEDP1002	MEDP1242	MEDP1488	MEDN0179	MEDN1264	MEDN1265	MEDN0348	MEDN2142	MEDN0511	MEDN1927	MEDP2334	MEDN1862	MEDP2283	MEDN1156	MEDN1640	MEDN0335	MEDN0201	MEDN1566	MEDP0147	MEDN2048	MEDN1796	MEDN1906	MEDN1448	MEDP1691	MEDP1695	MEDP2143	MEDN0659	MEDN0362	MEDN0758	MEDP0297	MEDN1886	MEDN1415	MEDP2339	MEDN1935	MEDN0745	MEDP1663	MEDN0748	MEDN0856	MEDP0084	MEDP0831	MEDP1657	MEDN1311	MEDN0429	MEDP1753	MEDP1434	MEDP1857	MEDP0336	MEDP1326	MEDP1702	MEDP1171	MEDN0323	MEDN1213	MEDP1391	MEDN2182	MEDN1428	MEDN0128	MEDP1407	MEDP0307	MEDP1251	MEDP1333	MEDP1334	MEDP1389	MEDN1413	MEDP1701	MEDP1339	MEDN1932	MEDP1399	MEDN0294	MEDN1643	MEDN0383	MEDN0752	MEDN0750	MEDN0756	MEDP1659	MEDP1165	MEDP1146	MEDN1690	MEDN0290	MEDP1658	MEDN0375	MEDN1416	MEDP1322	MEDP1412	MEDN0245	MEDP1418	MEDP1692	MEDP0494	MEDN0759	MEDN1597	MEDN1691	MEDN1648	MEDP1518	MEDN1269	MEDN1270	MEDN1278	MEDN1277	MEDN0390	MEDN1069	MEDN1688	MEDN1840	MEDN0378	MEDP1901	MEDP1428
*keep SampleID cluster MEDP1901	MEDP1428	MEDN0378	MEDN1840	MEDN1688	MEDN1069	MEDN0390	MEDN1691	MEDN1269	MEDN1270	MEDN0294	MEDN1277	MEDN1278	MEDN1416	MEDP1412	MEDP1518	MEDN0375	MEDN0128	MEDN0290	MEDP1322	MEDP1434	MEDN1597	MEDN1648	MEDN1690	MEDP1399	MEDN0759	MEDP0494	MEDP1692	MEDP2176	MEDN0750	MEDN0752	MEDN0756	MEDN0383	MEDP0307	MEDN1428	MEDN1643	MEDN1932	MEDN2048	MEDP1658	MEDN0245	MEDP1333	MEDP1334	MEDP1251	MEDP0147	MEDP1002	MEDP1146	MEDP1165	MEDP1659	MEDN0826	MEDP1389	MEDN1213	MEDN1311	MEDN2182	MEDN0362	MEDN0511	MEDN1886	MEDP1171	MEDP1326	MEDP1702	MEDN0011	MEDN1686	MEDP1753	MEDN0661	MEDP1339	MEDP1701	MEDP1045	MEDP1691	MEDN1193	MEDN0745	MEDP1407	MEDN0729	MEDP0336	MEDP1857	MEDN0659	MEDN1708	MEDP0296	MEDN1413	MEDP2143	MEDP1391	MEDN0201	MEDN0335	MEDP1418	MEDP2088	MEDN1965	MEDN0429	MEDN0748	MEDN0856	MEDP2431	MEDP2633	MEDN1576	MEDN1853	MEDP1928	MEDP1885	MEDN0179	MEDN0323	MEDN1264	MEDN1265	MEDP0598	MEDP1242	MEDP2283	MEDN0062	MEDN0072	MEDN1796	MEDN0348	MEDN1156	MEDN1029	MEDN0769	MEDN1415	MEDN2111	MEDP2334	MEDN0340	MEDP1201	MEDP1061	MEDP1318	MEDP1319	MEDP0084	MEDP0831	MEDN1862	MEDP0026	MEDP1778	MEDN1440	MEDN1699	MEDN1740	MEDN1753	MEDN1613	MEDP0174	MEDN0616	MEDN1246	MEDN1441	MEDN1640	MEDN1927	MEDP0052
rename cluster outcome
save "D:\真菌分析\data_pre\metabolism_predict_cluster1234.dta",replace 
restore 

preserve
keep if cluster==4 | cluster==3
gen outcome=0
replace outcome=1 if cluster==4
tab outcome
drop cluster
keep SampleID outcome MEDN1863	MEDN1445	MEDP1693	MEDN2099	MEDN1004	MEDP0052	MEDP1003	MEDP1126	MEDN1231	MEDP1424	MEDN0846	MEDP1404	MEDN0551	MEDP1919	MEDP0401	MEDP1954	MEDN1740	MEDN1753	MEDN1699	MEDN2197	MEDP2266	MEDP1255	MEDN1092	MEDP1670	MEDN1193	MEDN0071	MEDN1029	MEDP2132	MEDP0022	MEDN0729	MEDN1442	MEDN0431	MEDN1854	MEDP0006	MEDP0174	MEDP2325	MEDP0313	MEDN1708	MEDN0065	MEDP0325	MEDN0579	MEDN0754	MEDP2672	MEDN0105	MEDP2494	MEDP2088	MEDP2131	MEDP1321	MEDN0589	MEDN1440	MEDN0261	MEDN1441	MEDN1444	MEDP2633	MEDN0807	MEDP1061	MEDN1279	MEDN0364	MEDN0536	MEDN0340	MEDP0251	MEDN0100	MEDP1778	MEDP0026	MEDN0478	MEDN0682	MEDN1686	MEDN0011	MEDN2103	MEDN0661	MEDN0062	MEDN0072	MEDN1576	MEDP0389	MEDN0826	MEDN0398	MEDN0227	MEDN1221	MEDN0463	MEDP0453	MEDP2431	MEDP2223	MEDN0616	MEDP1045	MEDN1246	MEDP1928	MEDP0598	MEDN1861	MEDP1318	MEDP1319	MEDN0403	MEDP0296	MEDN0115	MEDN1965	MEDN1902	MEDN1126	MEDN1723	MEDN1077	MEDP2176	MEDN2111	MEDN1532	MEDN1943	MEDP1803	MEDN0658	MEDN0769	MEDN1613	MEDN1853	MEDP1885	MEDP1201	MEDP1002	MEDP1242	MEDP1488	MEDN0179	MEDN1264	MEDN1265	MEDN0348	MEDN2142	MEDN0511	MEDN1927	MEDP2334	MEDN1862	MEDP2283	MEDN1156	MEDN1640	MEDN0335	MEDN0201	MEDN1566	MEDP0147	MEDN2048	MEDN1796	MEDN1906	MEDN1448	MEDP1691	MEDP1695	MEDP2143	MEDN0659	MEDN0362	MEDN0758	MEDP0297	MEDN1886	MEDN1415	MEDP2339	MEDN1935	MEDN0745	MEDP1663	MEDN0748	MEDN0856	MEDP0084	MEDP0831	MEDP1657	MEDN1311	MEDN0429	MEDP1753	MEDP1434	MEDP1857	MEDP0336	MEDP1326	MEDP1702	MEDP1171	MEDN0323	MEDN1213	MEDP1391	MEDN2182	MEDN1428	MEDN0128	MEDP1407	MEDP0307	MEDP1251	MEDP1333	MEDP1334	MEDP1389	MEDN1413	MEDP1701	MEDP1339	MEDN1932	MEDP1399	MEDN0294	MEDN1643	MEDN0383	MEDN0752	MEDN0750	MEDN0756	MEDP1659	MEDP1165	MEDP1146	MEDN1690	MEDN0290	MEDP1658	MEDN0375	MEDN1416	MEDP1322	MEDP1412	MEDN0245	MEDP1418	MEDP1692	MEDP0494	MEDN0759	MEDN1597	MEDN1691	MEDN1648	MEDP1518	MEDN1269	MEDN1270	MEDN1278	MEDN1277	MEDN0390	MEDN1069	MEDN1688	MEDN1840	MEDN0378	MEDP1901	MEDP1428
*keep SampleID outcome MEDP1901	MEDP1428	MEDN0378	MEDN1840	MEDN1688	MEDN1069	MEDN0390	MEDN1691	MEDN1269	MEDN1270	MEDN0294	MEDN1277	MEDN1278	MEDN1416	MEDP1412	MEDP1518	MEDN0375	MEDN0128	MEDN0290	MEDP1322	MEDP1434	MEDN1597	MEDN1648	MEDN1690	MEDP1399	MEDN0759	MEDP0494	MEDP1692	MEDP2176	MEDN0750	MEDN0752	MEDN0756	MEDN0383	MEDP0307	MEDN1428	MEDN1643	MEDN1932	MEDN2048	MEDP1658	MEDN0245	MEDP1333	MEDP1334	MEDP1251	MEDP0147	MEDP1002	MEDP1146	MEDP1165	MEDP1659	MEDN0826	MEDP1389	MEDN1213	MEDN1311	MEDN2182	MEDN0362	MEDN0511	MEDN1886	MEDP1171	MEDP1326	MEDP1702	MEDN0011	MEDN1686	MEDP1753	MEDN0661	MEDP1339	MEDP1701	MEDP1045	MEDP1691	MEDN1193	MEDN0745	MEDP1407	MEDN0729	MEDP0336	MEDP1857	MEDN0659	MEDN1708	MEDP0296	MEDN1413	MEDP2143	MEDP1391	MEDN0201	MEDN0335	MEDP1418	MEDP2088	MEDN1965	MEDN0429	MEDN0748	MEDN0856	MEDP2431	MEDP2633	MEDN1576	MEDN1853	MEDP1928	MEDP1885	MEDN0179	MEDN0323	MEDN1264	MEDN1265	MEDP0598	MEDP1242	MEDP2283	MEDN0062	MEDN0072	MEDN1796	MEDN0348	MEDN1156	MEDN1029	MEDN0769	MEDN1415	MEDN2111	MEDP2334	MEDN0340	MEDP1201	MEDP1061	MEDP1318	MEDP1319	MEDP0084	MEDP0831	MEDN1862	MEDP0026	MEDP1778	MEDN1440	MEDN1699	MEDN1740	MEDN1753	MEDN1613	MEDP0174	MEDN0616	MEDN1246	MEDN1441	MEDN1640	MEDN1927	MEDP0052
sort outcome
save "D:\真菌分析\data_pre\metabolism_predict_cluster34.dta",replace
restore

preserve
gen outcome=0
replace outcome=1 if cluster==4
tab outcome
drop cluster
save "D:\真菌分析\data_pre\metabolism_predict_cluster4others.dta",replace
restore

preserve
keep if cluster==4 | cluster==2
gen outcome=0
replace outcome=1 if cluster==4
tab outcome
drop cluster
save "D:\真菌分析\data_pre\metabolism_predict_cluster24.dta",replace
restore

preserve
keep if cluster==4 | cluster==1
gen outcome=0
replace outcome=1 if cluster==4
tab outcome
drop cluster
save "D:\真菌分析\data_pre\metabolism_predict_cluster14.dta",replace
restore

preserve
keep if cluster==1 | cluster==2
keep SampleID cluster MEDN1863	MEDN1445	MEDP1693	MEDN2099	MEDN1004	MEDP0052	MEDP1003	MEDP1126	MEDN1231	MEDP1424	MEDN0846	MEDP1404	MEDN0551	MEDP1919	MEDP0401	MEDP1954	MEDN1740	MEDN1753	MEDN1699	MEDN2197	MEDP2266	MEDP1255	MEDN1092	MEDP1670	MEDN1193	MEDN0071	MEDN1029	MEDP2132	MEDP0022	MEDN0729	MEDN1442	MEDN0431	MEDN1854	MEDP0006	MEDP0174	MEDP2325	MEDP0313	MEDN1708	MEDN0065	MEDP0325	MEDN0579	MEDN0754	MEDP2672	MEDN0105	MEDP2494	MEDP2088	MEDP2131	MEDP1321	MEDN0589	MEDN1440	MEDN0261	MEDN1441	MEDN1444	MEDP2633	MEDN0807	MEDP1061	MEDN1279	MEDN0364	MEDN0536	MEDN0340	MEDP0251	MEDN0100	MEDP1778	MEDP0026	MEDN0478	MEDN0682	MEDN1686	MEDN0011	MEDN2103	MEDN0661	MEDN0062	MEDN0072	MEDN1576	MEDP0389	MEDN0826	MEDN0398	MEDN0227	MEDN1221	MEDN0463	MEDP0453	MEDP2431	MEDP2223	MEDN0616	MEDP1045	MEDN1246	MEDP1928	MEDP0598	MEDN1861	MEDP1318	MEDP1319	MEDN0403	MEDP0296	MEDN0115	MEDN1965	MEDN1902	MEDN1126	MEDN1723	MEDN1077	MEDP2176	MEDN2111	MEDN1532	MEDN1943	MEDP1803	MEDN0658	MEDN0769	MEDN1613	MEDN1853	MEDP1885	MEDP1201	MEDP1002	MEDP1242	MEDP1488	MEDN0179	MEDN1264	MEDN1265	MEDN0348	MEDN2142	MEDN0511	MEDN1927	MEDP2334	MEDN1862	MEDP2283	MEDN1156	MEDN1640	MEDN0335	MEDN0201	MEDN1566	MEDP0147	MEDN2048	MEDN1796	MEDN1906	MEDN1448	MEDP1691	MEDP1695	MEDP2143	MEDN0659	MEDN0362	MEDN0758	MEDP0297	MEDN1886	MEDN1415	MEDP2339	MEDN1935	MEDN0745	MEDP1663	MEDN0748	MEDN0856	MEDP0084	MEDP0831	MEDP1657	MEDN1311	MEDN0429	MEDP1753	MEDP1434	MEDP1857	MEDP0336	MEDP1326	MEDP1702	MEDP1171	MEDN0323	MEDN1213	MEDP1391	MEDN2182	MEDN1428	MEDN0128	MEDP1407	MEDP0307	MEDP1251	MEDP1333	MEDP1334	MEDP1389	MEDN1413	MEDP1701	MEDP1339	MEDN1932	MEDP1399	MEDN0294	MEDN1643	MEDN0383	MEDN0752	MEDN0750	MEDN0756	MEDP1659	MEDP1165	MEDP1146	MEDN1690	MEDN0290	MEDP1658	MEDN0375	MEDN1416	MEDP1322	MEDP1412	MEDN0245	MEDP1418	MEDP1692	MEDP0494	MEDN0759	MEDN1597	MEDN1691	MEDN1648	MEDP1518	MEDN1269	MEDN1270	MEDN1278	MEDN1277	MEDN0390	MEDN1069	MEDN1688	MEDN1840	MEDN0378	MEDP1901	MEDP1428
save "D:\真菌分析\data_pre\metabolism_cluster12raw.dta",replace
restore

*----------Validated on the 2018 datasets-------------------

use  "D:\真菌分析\data_pre\metabolism_18.dta",clear
merge m:m Idind using "D:\真菌分析\data_pre\repeatDMM_18cluster.dta"
keep if _merge==3
tab cluster
sort cluster
by cluster:su MEDN0390

foreach metabolism of varlist CMPF_N-mwxq04_N {
	gen ln_`metabolism'=log(`metabolism') 
	drop `metabolism'
	rename ln_`metabolism' `metabolism'
}
keep SampleID cluster CMPF_N-mwxq04_N

preserve
keep if cluster==4 | cluster==3
gen outcome=0
replace outcome=1 if cluster==4
tab outcome
drop cluster
keep SampleID outcome MEDN1863	MEDN1445	MEDP1693	MEDN2099	MEDN1004	MEDP0052	MEDP1003	MEDP1126	MEDN1231	MEDP1424	MEDN0846	MEDP1404	MEDN0551	MEDP1919	MEDP0401	MEDP1954	MEDN1740	MEDN1753	MEDN1699	MEDN2197	MEDP2266	MEDP1255	MEDN1092	MEDP1670	MEDN1193	MEDN0071	MEDN1029	MEDP2132	MEDP0022	MEDN0729	MEDN1442	MEDN0431	MEDN1854	MEDP0006	MEDP0174	MEDP2325	MEDP0313	MEDN1708	MEDN0065	MEDP0325	MEDN0579	MEDN0754	MEDP2672	MEDN0105	MEDP2494	MEDP2088	MEDP2131	MEDP1321	MEDN0589	MEDN1440	MEDN0261	MEDN1441	MEDN1444	MEDP2633	MEDN0807	MEDP1061	MEDN1279	MEDN0364	MEDN0536	MEDN0340	MEDP0251	MEDN0100	MEDP1778	MEDP0026	MEDN0478	MEDN0682	MEDN1686	MEDN0011	MEDN2103	MEDN0661	MEDN0062	MEDN0072	MEDN1576	MEDP0389	MEDN0826	MEDN0398	MEDN0227	MEDN1221	MEDN0463	MEDP0453	MEDP2431	MEDP2223	MEDN0616	MEDP1045	MEDN1246	MEDP1928	MEDP0598	MEDN1861	MEDP1318	MEDP1319	MEDN0403	MEDP0296	MEDN0115	MEDN1965	MEDN1902	MEDN1126	MEDN1723	MEDN1077	MEDP2176	MEDN2111	MEDN1532	MEDN1943	MEDP1803	MEDN0658	MEDN0769	MEDN1613	MEDN1853	MEDP1885	MEDP1201	MEDP1002	MEDP1242	MEDP1488	MEDN0179	MEDN1264	MEDN1265	MEDN0348	MEDN2142	MEDN0511	MEDN1927	MEDP2334	MEDN1862	MEDP2283	MEDN1156	MEDN1640	MEDN0335	MEDN0201	MEDN1566	MEDP0147	MEDN2048	MEDN1796	MEDN1906	MEDN1448	MEDP1691	MEDP1695	MEDP2143	MEDN0659	MEDN0362	MEDN0758	MEDP0297	MEDN1886	MEDN1415	MEDP2339	MEDN1935	MEDN0745	MEDP1663	MEDN0748	MEDN0856	MEDP0084	MEDP0831	MEDP1657	MEDN1311	MEDN0429	MEDP1753	MEDP1434	MEDP1857	MEDP0336	MEDP1326	MEDP1702	MEDP1171	MEDN0323	MEDN1213	MEDP1391	MEDN2182	MEDN1428	MEDN0128	MEDP1407	MEDP0307	MEDP1251	MEDP1333	MEDP1334	MEDP1389	MEDN1413	MEDP1701	MEDP1339	MEDN1932	MEDP1399	MEDN0294	MEDN1643	MEDN0383	MEDN0752	MEDN0750	MEDN0756	MEDP1659	MEDP1165	MEDP1146	MEDN1690	MEDN0290	MEDP1658	MEDN0375	MEDN1416	MEDP1322	MEDP1412	MEDN0245	MEDP1418	MEDP1692	MEDP0494	MEDN0759	MEDN1597	MEDN1691	MEDN1648	MEDP1518	MEDN1269	MEDN1270	MEDN1278	MEDN1277	MEDN0390	MEDN1069	MEDN1688	MEDN1840	MEDN0378	MEDP1901	MEDP1428
save "D:\真菌分析\data_pre\metabolism_predict_cluster34_18.dta",replace
restore

*---------Validate on the 2018 datasets after remove the repeat datasets--

use "D:\真菌分析\data_pre\metabolism_predict_cluster34_18.dta",clear
merge 1:1 SampleID using "D:\真菌分析\data_pre\samplelist.dta"
keep if _merge==3
drop _merge
save "D:\真菌分析\data_pre\metabolism_predict_cluster34_18_sens.dta",replace
 

*---AUC results confirm(SampleID)---

*--Discovery--

import excel "D:\真菌分析\Final results\predict_update\metabolis_fungi_discovery_auc.xlsx", sheet("Sheet1") firstrow clear 
merge 1:1 SampleID using "D:\真菌分析\data_pre\metabolism_predict_cluster34.dta"
tab outcome label // OK

*--Validation-all--

import excel "D:\真菌分析\Final results\predict_update\metabolis_fungi_validation_select_auc.xlsx", sheet("Sheet1") firstrow clear 
merge 1:1 SampleID using "D:\真菌分析\data_pre\metabolism_predict_cluster34_18.dta"
tab outcome label // OK
keep if _merge==3
drop _merge
merge 1:1 SampleID using "D:\真菌分析\data_pre\metabolism_cluster.dta"
keep if _merge==3
tab cluster_15 cluster_18 // only 4 participants trans from 3 to 4,and same number participants for 4 to 3.
drop if cluster_15==3 | cluster_15==4 
sort cluster_15 cluster_18
egen group=group(cluster_15 cluster_18)
sort group
by group:su predict
keep cluster_15 group predict
export excel using "D:\真菌分析\Plot\data\trans_predict.xlsx",  firstrow(variables) sheet("2018") sheetreplace  
 
*----------cluster 3 vs cluster 4 metabolism VIP values------

import excel "D:\CNHS\data\代谢组数据\Repeatdata\1.Data_Assess\all_group\ALL_sample_data_leveled.xlsx", sheet("Sheet1") firstrow clear 
keep Index-物质二级分类 Level Formula HMDB PubChemCID kegg_map
save "D:\真菌分析\data_pre\metabolism_inf.dta",replace

import delimited "D:\真菌分析\Final results\All data\metabolism\kw_test.csv",  clear 
rename id Index
save "D:\真菌分析\data_pre\metabolism_kw.dta",replace
 
 
import delimited "D:\真菌分析\Final results\All data\metabolism\metabolim_vip.csv",  clear 
rename vip Index
rename metabolism VIP
merge 1:1 Index using  "D:\真菌分析\data_pre\metabolism_inf.dta"
keep if _merge==3
drop _merge
merge 1:1 Index using "D:\真菌分析\data_pre\metabolism_kw.dta"
keep if _merge==3
drop _merge
sort VIP
save "D:\真菌分析\data_pre\metabolism_vipkw.dta",replace 
export excel using "D:\真菌分析\Final results\All data\metabolism\metabolism_vip.xlsx",  firstrow(variables) sheet("2015") sheetreplace  

*----------metabolism fold change calculate and results present----

import delimited "D:\真菌分析\Final results\All data\metabolism\metabolismfoldchange.csv",  clear 
rename index Index
merge 1:1 Index using "D:\真菌分析\data_pre\metabolism_vipkw.dta"
drop _merge
save "D:\真菌分析\data_pre\metabolism_vipkwfold.dta",replace
gen change=.
replace change=0 if fold<1 & fdr<0.05
replace change=1 if fold>1 & fdr<0.05
tab change
sort  change // 48 decreased and 84 increased
tab ClassI
*duplicates drop VIP, force // 31
export excel using "D:\真菌分析\Final results\All data\metabolism\metabolism_results_all.xlsx",  firstrow(variables) sheet("2015") sheetreplace  
replace ClassI="others" if ClassI!="Amino acid and Its metabolites" & ClassI!="Bile acids" & ClassI!="Benzene and substituted derivatives" & ClassI!="FA" & ClassI!="GP" & ClassI!="Organic acid And Its derivatives " 
save "D:\真菌分析\data_pre\metabolism_vipkwfold.dta",replace 
 
*---------metabolism module and phenotype analyses----------based all data

*---merge the modules' mapping--

import delimited "D:\真菌分析\Final results\All data\metabolism\metabolism_module.csv",  clear 
rename metabolism Index
merge 1:1 Index using "D:\真菌分析\data_pre\metabolism_vipkwfold.dta"
tab color
order Index HMDB color VIP fold fdr
sort color VIP // metabolism_vip

*--continue vars--

import delimited "D:\真菌分析\Final results\All data\metabolism\metabolism_module_降维.csv",  clear 
rename sampleid SampleID
merge 1:1 SampleID using  "D:\真菌分析\data_pre\2015_metabolism_phenotype.dta"

qui foreach var of varlist mepink-megrey glucose HbA1c ins HDL_C LDL_C tc tg SBP DBP{
	egen d_`var'=std(`var')
	drop `var'
	rename d_`var' `var'
}


tempname coef
 tempfile res
 postfile `coef' str200(outcome metabolism)  float(n rr lul uul p) using "`res'", replace
 foreach outcome of varlist glucose HbA1c ins HDL_C LDL_C tc tg SBP DBP {
foreach var of varlist mepink-megrey {
regress `outcome' `var' age i.sex BMI
 test `var'
 post `coef'   ("`outcome'") ("`var'")  (e(N)) (_b[`var'])  (_b[`var']-1.96*_se[`var'])  (_b[`var']+1.96*_se[`var'])  (chi2tail(1,(_b[`var']/_se[`var'])^2)) 
}
}
 postclose `coef'
 preserve
use "`res'", clear

export excel using "D:\真菌分析\Final results\All data\metabolism\metabolism_module_analyses.xlsx",  firstrow(variables) sheet("module_phenotype") sheetreplace  
restore

*--binary vars--

tempname coef
tempfile res
 postfile `coef' str200(outcome metabolism)  float(n rr lul uul p) using "`res'", replace
 foreach outcome of varlist t2d hyp dys {
foreach var of varlist mepink-megrey {
 logit `outcome' `var' age i.sex BMI
test `var'
 post `coef'   ("`outcome'") ("`var'")  (e(N)) (exp(_b[`var']))  (exp(_b[`var']-1.96*_se[`var']))  (exp(_b[`var']+1.96*_se[`var']))  (chi2tail(1,(_b[`var']/_se[`var'])^2))  
}
}
  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\All data\metabolism\metabolism_module_analyses.xlsx",  firstrow(variables) sheet("module_diseases") sheetreplace  
restore

*--fungi cluster--

tab cluster 
keep if cluster==3 | cluster==4
gen outcome=0
replace outcome=1 if cluster==4

 tempname coef
tempfile res
 postfile `coef' str200(outcome metabolism)  float(n rr lul uul p) using "`res'", replace
 foreach outcome of varlist outcome{
foreach var of varlist mepink-megrey {
 logit `outcome' `var' age i.sex BMI
test `var'
 post `coef'   ("`outcome'") ("`var'")  (e(N)) (exp(_b[`var']))  (exp(_b[`var']-1.96*_se[`var']))  (exp(_b[`var']+1.96*_se[`var']))  (chi2tail(1,(_b[`var']/_se[`var'])^2))  
}
}
  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\All data\metabolism\metabolism_module_analyses.xlsx",  firstrow(variables) sheet("module_cluster") sheetreplace  
restore

*------metabolism module and fungi associations----

import delimited "D:\真菌分析\Final results\All data\metabolism\metabolism_module_降维.csv",  clear 
rename sampleid SampleID
merge 1:1 SampleID using  "D:\真菌分析\data_pre\2015_metabolism_phenotype.dta"
drop _merge
merge 1:1 SampleID using "D:\真菌分析\data_pre\2015_repeat_fungi_clr.dta"
qui foreach var of varlist mepink-megrey g2-g633{
	egen d_`var'=std(`var')
	drop `var'
	rename d_`var' `var'
}

tempname coef
 tempfile res
 postfile `coef' str200(module micro)  float(n rr lul uul p) using "`res'", replace
 foreach module of varlist  mepink-megrey {
foreach var of varlist g2-g633{
regress `module' `var' age i.sex BMI
 test `var'
 post `coef'   ("`module'") ("`var'")  (e(N)) (_b[`var'])  (_b[`var']-1.96*_se[`var'])  (_b[`var']+1.96*_se[`var'])  (chi2tail(1,(_b[`var']/_se[`var'])^2)) 
}
}
 postclose `coef'
 preserve
use "`res'", clear

export excel using "D:\真菌分析\Final results\All data\metabolism\metabolism_module_analyses.xlsx",  firstrow(variables) sheet("micro_module") sheetreplace  
restore

*-----results mapping-------

import excel "D:\真菌分析\Final results\All data\metabolism\metabolism_module_analyses.xlsx", sheet("micro_module") firstrow clear 
rename micro code_repeat
merge m:m code_repeat using   "D:\真菌分析\data_pre\fungi_mapping_all_repeat.dta"
keep if _merge==3
keep if Code_All=="G139" | Code_All=="G286" | Code_All=="G142" | Code_All=="G247" | Code_All=="G245" | Code_All=="G280" | Code_All=="G21" | Code_All=="G746" | Code_All=="G186" 
sort p
keep if p<0.05
tab module

*--------------mediation analysis datasets----------

use "D:\真菌分析\data_pre\repeatDMM_15cluster.dta",clear
merge m:m Idind using  "D:\真菌分析\data_pre\metabolism_15.dta"
keep if _merge==3
rename cluster cluster_15 
rename SampleID SampleID_15
keep Idind cluster_15 SampleID_15 
save  "D:\真菌分析\data_pre\metabolism_cluster_15.dta",replace

use  "D:\真菌分析\data_pre\metabolism_18.dta",clear
merge m:m Idind using "D:\真菌分析\data_pre\repeatDMM_18cluster.dta"
keep if _merge==3
drop _merge
rename SampleID SampleID_18
keep Idind cluster_18 SampleID_18
merge 1:1 Idind using "D:\真菌分析\data_pre\metabolism_cluster_15.dta"
keep SampleID_15 SampleID_18 cluster_18 cluster_15 Idind
rename SampleID_18 SampleID
save  "D:\真菌分析\data_pre\metabolism_cluster.dta",replace // cluster 3 or 4, and trans from 1 or 2
 
tab cluster_18 cluster_15 // 43个一直保持cluster3， 15个一直保持cluster4
tab cluster_15 // 3-85, 4-80
tab cluster_18 // 3-114, 4-82
drop if cluster_15==3 & cluster_18==3
drop if cluster_15==4 & cluster_18==4
drop if cluster_18==1 | cluster_18==2
tab cluster_18 // 71 vs 67
keep SampleID Idind cluster_15 cluster_18
save  "D:\真菌分析\data_pre\samplelist.dta",replace // cluster 3 or 4, and trans from 1 or 2

*--metabolism modules and diseases---

import delimited "D:\真菌分析\Final results\All data\metabolism\metabolism_module_降维.csv",  clear 
rename sampleid SampleID
merge 1:1 SampleID using  "D:\真菌分析\data_pre\2015_metabolism_phenotype.dta"
drop _merge
merge 1:1 SampleID using "D:\真菌分析\data_pre\2015_repeat_fungi_clr.dta"
qui foreach var of varlist mepink-megrey g2-g633{
	egen d_`var'=std(`var')
	drop `var'
	rename d_`var' `var'
}
tab cluster
keep if cluster==3 | cluster==4
gen outcome=0
replace outcome=1 if cluster==4
qui foreach var of varlist mepink-megrey g2-g633{
	egen d_`var'=std(`var')
	drop `var'
	rename d_`var' `var'
}

preserve
keep SampleID mepink-megrey
save "D:\真菌分析\data_pre\metabolism_data.dta",replace 
restore

preserve
keep SampleID outcome
save "D:\真菌分析\data_pre\outcome_data.dta",replace 
restore

preserve
keep SampleID hyp t2d
save "D:\真菌分析\data_pre\diseases_data.dta",replace 
restore

preserve
keep SampleID age sex BMI
save "D:\真菌分析\data_pre\cov_data.dta",replace 
restore

tempname coef
tempfile res
 postfile `coef' str200(outcome metabolism)  float(n rr lul uul p) using "`res'", replace
 foreach outcome of varlist t2d hyp dys {
foreach var of varlist mepink-megrey {
 logit `outcome' `var' age i.sex BMI
test `var'
 post `coef'   ("`outcome'") ("`var'")  (e(N)) (exp(_b[`var']))  (exp(_b[`var']-1.96*_se[`var']))  (exp(_b[`var']+1.96*_se[`var']))  (chi2tail(1,(_b[`var']/_se[`var'])^2))  
}
}
  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\All data\metabolism\mediation_input.xlsx",  firstrow(variables) sheet("module_diseases") sheetreplace  
restore

*-----metabolism module and phenotype--

qui foreach var of varlist glucose HbA1c ins HDL_C LDL_C tc tg SBP DBP{
	egen d_`var'=std(`var')
	drop `var'
	rename d_`var' `var'
}

tempname coef
 tempfile res
 postfile `coef' str200(outcome metabolism)  float(n rr lul uul p) using "`res'", replace
 foreach outcome of varlist glucose HbA1c ins HDL_C LDL_C tc tg SBP DBP {
foreach var of varlist mepink-megrey {
regress `outcome' `var' age i.sex BMI
 test `var'
 post `coef'   ("`outcome'") ("`var'")  (e(N)) (_b[`var'])  (_b[`var']-1.96*_se[`var'])  (_b[`var']+1.96*_se[`var'])  (chi2tail(1,(_b[`var']/_se[`var'])^2)) 
}
}
 postclose `coef'
 preserve
use "`res'", clear

export excel using "D:\真菌分析\Final results\All data\metabolism\mediation_input.xlsx",  firstrow(variables) sheet("module_phenotype") sheetreplace  
restore

*---cluster and metabolism--------

tempname coef
 tempfile res
 postfile `coef' str200(metabolism)  float(n rr lul uul p) using "`res'", replace
foreach var of varlist mepink-megrey {
regress  `var' i.outcome age i.sex BMI
 test 1.outcome
 post `coef'   ("`var'")  (e(N)) (_b[1.outcome])  (_b[1.outcome]-1.96*_se[1.outcome])  (_b[1.outcome]+1.96*_se[1.outcome])  (chi2tail(1,(_b[1.outcome]/_se[1.outcome])^2)) 
}

 postclose `coef'
 preserve
use "`res'", clear

export excel using "D:\真菌分析\Final results\All data\metabolism\mediation_input.xlsx",  firstrow(variables) sheet("fungi_module") sheetreplace  
restore

*--sig metabolism metabolism and diseases---

import delimited "D:\真菌分析\Final results\All data\metabolism\metabolism_module.csv",  clear 
rename metabolism Index
merge 1:1 Index using "D:\真菌分析\data_pre\metabolism_vipkwfold.dta"
tab color
keep if color=="brown" | color=="purple" | color=="yellow" 
keep color Index
save "D:\真菌分析\data_pre\module_metabolism_list.dta",replace 
 
use "D:\真菌分析\data_pre\metabolism_predict_cluster34.dta",clear
merge 1:1 SampleID using  "D:\真菌分析\data_pre\2015_metabolism_phenotype.dta"
keep if _merge==3
drop _merge
keep SampleID outcome-cluster MEDN0011	MEDN0062	MEDN0072	MEDN0378	MEDN0383	MEDN0390	MEDN0807	MEDN1069	MEDN1092	MEDN1246	MEDN1269	MEDN1270	MEDN1277	MEDN1278	MEDN1640	MEDN1686	MEDN1840	MEDN1965	MEDN2103	MEDP0006	MEDP0022	MEDP0026	MEDP0052	MEDP0174	MEDP0325	MEDP0336	MEDP0494	MEDP1003	MEDP1061	MEDP1126	MEDP1691	MEDP1692	MEDP1693	MEDP1778	MEDP1885	MEDP1928	MEDP2132
*merge 1:1 SampleID using "D:\真菌分析\data_pre\2015_repeat_fungi_clr.dta"
qui foreach var of varlist MEDN0011-MEDP2132{
	egen d_`var'=std(`var')
	drop `var'
	rename d_`var' `var'
}

preserve
keep SampleID MEDN0011-MEDP2132
save "D:\真菌分析\data_pre\meta_data.dta",replace 
restore

preserve
rename MEDP1692 LPC_15
rename MEDN1270 LPE_22
rename MEDN0378 FFA_183
rename MEDP1691 LPC_17
rename MEDP0026 L_Valine
rename MEDP1778 L_Norvaline
rename MEDN1840 Pinolenic_acid
rename MEDP0336 LPC_14
rename MEDP0006 L_Glycine
keep SampleID t2d FFA_183 LPE_22 Pinolenic_acid L_Glycine L_Valine LPC_14 LPC_17 LPC_15 L_Norvaline
restore


tempname coef
tempfile res
 postfile `coef' str200(outcome metabolism)  float(n rr lul uul p) using "`res'", replace
 foreach outcome of varlist t2d hyp  {
foreach var of varlist MEDN0011-MEDP2132 {
 logit `outcome' `var' age i.sex BMI
test `var'
 post `coef'   ("`outcome'") ("`var'")  (e(N)) (exp(_b[`var']))  (exp(_b[`var']-1.96*_se[`var']))  (exp(_b[`var']+1.96*_se[`var']))  (chi2tail(1,(_b[`var']/_se[`var'])^2))  
}
}
  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\All data\metabolism\mediation_input.xlsx",  firstrow(variables) sheet("metabolism_diseases") sheetreplace  
restore

*--mapping the metabolism modules--

import excel "D:\真菌分析\Final results\All data\metabolism\mediation_input.xlsx", sheet("metabolism_diseases") firstrow clear 
rename metabolism Index
merge m:m Index using  "D:\真菌分析\data_pre\module_metabolism_list.dta" 
sort p
gen mark=color+outcome
export excel using "D:\真菌分析\Final results\All data\metabolism\mediation_input.xlsx",  firstrow(variables) sheet("metabolism_diseases_mpping") sheetreplace  

*---driven fungi and metabolism--------

use "D:\真菌分析\data_pre\metabolism_predict_cluster34.dta",clear
merge 1:1 SampleID using  "D:\真菌分析\data_pre\2015_metabolism_phenotype.dta"
keep if _merge==3
drop _merge
merge 1:1 SampleID using "D:\真菌分析\data_pre\2015_repeat_fungi_clr.dta"
keep SampleID outcome-cluster g114 g116 g146 g186 g188 g19 g220 g225 g623 MEDN0011	MEDN0062	MEDN0072	MEDN0378	MEDN0383	MEDN0390	MEDN0807	MEDN1069	MEDN1092	MEDN1246	MEDN1269	MEDN1270	MEDN1277	MEDN1278	MEDN1640	MEDN1686	MEDN1840	MEDN1965	MEDN2103	MEDP0006	MEDP0022	MEDP0026	MEDP0052	MEDP0174	MEDP0325	MEDP0336	MEDP0494	MEDP1003	MEDP1061	MEDP1126	MEDP1691	MEDP1692	MEDP1693	MEDP1778	MEDP1885	MEDP1928	MEDP2132
drop MEDP0174	MEDP1061	MEDP1003	MEDP1126	MEDN0062	MEDN0072	MEDN1640	MEDN1246	MEDN1965	MEDN1686	MEDN0011	MEDN1092	MEDP1885	MEDP1928	MEDN2103

tab outcome,missing
drop if outcome==.
qui foreach var of varlist MEDN0378-MEDP2132 g19-g623{
	egen d_`var'=std(`var')
	drop `var'
	rename d_`var' `var'
}
preserve
keep SampleID g19-g623
save "D:\真菌分析\data_pre\fungi_data.dta",replace 
restore

tempname coef
 tempfile res
 postfile `coef' str200(module micro)  float(n rr lul uul p) using "`res'", replace
 foreach module of varlist  MEDN0378-MEDP2132 {
foreach var of varlist g19-g623{
regress `module' `var' age i.sex BMI
 test `var'
 post `coef'   ("`module'") ("`var'")  (e(N)) (_b[`var'])  (_b[`var']-1.96*_se[`var'])  (_b[`var']+1.96*_se[`var'])  (chi2tail(1,(_b[`var']/_se[`var'])^2)) 
}
}
 postclose `coef'
 preserve
use "`res'", clear

export excel using "D:\真菌分析\Final results\All data\metabolism\mediation_input.xlsx",  firstrow(variables) sheet("fungi_metabolism") sheetreplace  
restore

*----driven fungi and diseases---

tempname coef
tempfile res
 postfile `coef' str200(outcome metabolism)  float(n rr lul uul p) using "`res'", replace
 foreach outcome of varlist t2d hyp dys {
foreach var of varlist g19-g623 {
 logit `outcome' `var' age i.sex BMI
test `var'
 post `coef'   ("`outcome'") ("`var'")  (e(N)) (exp(_b[`var']))  (exp(_b[`var']-1.96*_se[`var']))  (exp(_b[`var']+1.96*_se[`var']))  (chi2tail(1,(_b[`var']/_se[`var'])^2))  
}
}
  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\All data\metabolism\mediation_input.xlsx",  firstrow(variables) sheet("fungi_diseases") sheetreplace  
restore

*-------基于18年数据test肠道真菌与代谢物的中介效应----

use  "D:\真菌分析\data_pre\2018_all_phenotype.dta",clear
rename idind Idind
merge m:m Idind using  "D:\真菌分析\data_pre\metabolism_18.dta"
keep if _merge==3
drop _merge
merge 1:1 SampleID using "D:\真菌分析\data_pre\2018_repeat_fungi_clr.dta"
keep if _merge==3
drop _merge
merge 1:1 Idind using "D:\真菌分析\data_pre\repeatDMM_18cluster.dta"
keep if _merge==3
keep if cluster==3 | cluster==4
gen outcome=0
replace outcome=1 if cluster==4
drop _merge

*--model 2 drop the repeat data--

merge 1:1 SampleID using "D:\真菌分析\data_pre\samplelist.dta"
keep if _merge==3
drop _merge

*--model 1: no drop--

save  "D:\真菌分析\data_pre\2018_mediation.dta",replace

logit t2d g220 age i.sex BMI,or
logit t2d MEDN0378 age i.sex BMI,or
*logit incident_t2d g220 age i.sex BMI,or // 原始所有数据集里显著
keep SampleID outcome age sex BMI t2d hyp MEDN0378 MEDN0390 MEDN1269 MEDN1270 MEDN1840 MEDP0026 MEDP1778
qui foreach var of varlist MEDN0378-MEDP1778{
	egen d_`var'=std(`var')
	drop `var'
	rename d_`var' `var'
}

preserve
keep SampleID MEDN0378-MEDP1778
save "D:\真菌分析\data_pre\meta18_data.dta",replace 
restore


tempname coef
tempfile res
 postfile `coef' str200(outcome metabolism)  float(n rr lul uul p) using "`res'", replace
 foreach outcome of varlist t2d hyp  {
foreach var of varlist MEDN0378-MEDP1778 {
 logit `outcome' `var' age i.sex BMI
test `var'
 post `coef'   ("`outcome'") ("`var'")  (e(N)) (exp(_b[`var']))  (exp(_b[`var']-1.96*_se[`var']))  (exp(_b[`var']+1.96*_se[`var']))  (chi2tail(1,(_b[`var']/_se[`var'])^2))  
}
}
  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\All data\metabolism\mediation_input.xlsx",  firstrow(variables) sheet("metabolism_diseases_18") sheetreplace  
restore

*--mapping the metabolism modules--

import excel "D:\真菌分析\Final results\All data\metabolism\mediation_input.xlsx", sheet("metabolism_diseases_18") firstrow clear 
rename metabolism Index
merge m:m Index using  "D:\真菌分析\data_pre\module_metabolism_list.dta" 
sort p
gen mark=color+outcome
export excel using "D:\真菌分析\Final results\All data\metabolism\mediation_input.xlsx",  firstrow(variables) sheet("metabolism_diseases_18_mpping") sheetreplace  

*---driven fungi and metabolism--------

use  "D:\真菌分析\data_pre\2018_mediation.dta",clear
keep SampleID  outcome age sex BMI t2d hyp g146  g220 g225 MEDN0378 MEDN0390 MEDN1269 MEDN1270 MEDN1840 MEDP0026 MEDP1778
tab outcome,missing
drop if outcome==.
qui foreach var of varlist MEDN0378-MEDP1778 g146-g225{
	egen d_`var'=std(`var')
	drop `var'
	rename d_`var' `var'
}

preserve
keep SampleID hyp t2d
save "D:\真菌分析\data_pre\diseases18_data.dta",replace 
restore

preserve
keep SampleID g146-g225
save "D:\真菌分析\data_pre\fungi18_data.dta",replace 
restore

preserve
keep SampleID age sex BMI
save "D:\真菌分析\data_pre\cov18_data.dta",replace 
restore

tempname coef
 tempfile res
 postfile `coef' str200(module micro)  float(n rr lul uul p) using "`res'", replace
 foreach module of varlist  MEDN0378-MEDP1778 {
foreach var of varlist g146-g225{
regress `module' `var' age i.sex BMI
 test `var'
 post `coef'   ("`module'") ("`var'")  (e(N)) (_b[`var'])  (_b[`var']-1.96*_se[`var'])  (_b[`var']+1.96*_se[`var'])  (chi2tail(1,(_b[`var']/_se[`var'])^2)) 
}
}
 postclose `coef'
 preserve
use "`res'", clear

export excel using "D:\真菌分析\Final results\All data\metabolism\mediation_input.xlsx",  firstrow(variables) sheet("fungi_metabolism_18") sheetreplace  
restore

*----driven fungi and diseases---

tempname coef
tempfile res
 postfile `coef' str200(outcome metabolism)  float(n rr lul uul p) using "`res'", replace
 foreach outcome of varlist t2d hyp  {
foreach var of varlist g146-g225 {
 logit `outcome' `var' age i.sex BMI
test `var'
 post `coef'   ("`outcome'") ("`var'")  (e(N)) (exp(_b[`var']))  (exp(_b[`var']-1.96*_se[`var']))  (exp(_b[`var']+1.96*_se[`var']))  (chi2tail(1,(_b[`var']/_se[`var'])^2))  
}
}
  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\All data\metabolism\mediation_input.xlsx",  firstrow(variables) sheet("fungi_diseases_18") sheetreplace  
restore

*------------Results combined--------------

*--module results---

use "D:\真菌分析\data_pre\metabolism_vipkwfold.dta",clear
rename Index metabolism
drop v3 Formula PubChemCID v1 pvalue fdr
save  "D:\真菌分析\data_pre\metabolism_names.dta",replace

*---discovery cohort--

*---metabolism and disease--

import excel "D:\真菌分析\Plot\data\代谢组分析.xlsx", sheet("metabolism_disease") firstrow clear 
drop _merge
merge m:m metabolism using  "D:\真菌分析\data_pre\metabolism_names.dta"
keep if _merge==3
keep outcome Compounds rr lul uul p p_adjust color ClassI ClassII

*---fungi and metabolism--

import excel "D:\真菌分析\Plot\data\代谢组分析.xlsx", sheet("fungi_metabolism") firstrow clear 
merge m:m metabolism using  "D:\真菌分析\data_pre\metabolism_names.dta"
keep if _merge==3
drop _merge
rename fungi code_repeat
merge m:m code_repeat using "D:\真菌分析\data_pre\fungi_mapping_all_repeat.dta"
keep if _merge==3
keep outcome Compounds rr lul uul p p_adjust color ClassI ClassII

*--中介分析结果--

import delimited "D:\真菌分析\Final results\All data\metabolism\代谢物中介分析结果.csv",  clear 
order pval_mediate qval_mediate 
sort pval_mediate
keep if qval_mediate<0.05 // drop 2
keep if qval_mediate_inverse>0.05 // drop 0
drop  coefci_mediate pval_direct coefci_direct pval_total coefci_total pval_ratio coefci_ratio pval_direct_inverse
merge m:m metabolism using  "D:\真菌分析\data_pre\metabolism_names.dta"
keep if _merge==3
drop _merge
rename fungi code_repeat
merge m:m code_repeat using "D:\真菌分析\data_pre\fungi_mapping_all_repeat.dta"
keep if _merge==3
keep g Compounds 物质 pval_mediate outcome coef_mediate rr rr_metabolismdisease coef_ratio
order g Compounds 物质 pval_mediate outcome coef_mediate rr rr_metabolismdisease coef_ratio
gen mark=g+Compounds
rename coef_mediate mediabe_d
rename rr rr_d
rename rr_metabolismdisease rr_metadis_d
rename coef_ratio coef_ratio_d
rename outcome outcome_d
rename pval_mediate pval_mediate_d
save  "D:\真菌分析\data_pre\mediate_sig_discovery.dta",replace

*--validation cohort

import delimited "D:\真菌分析\Final results\All data\metabolism\代谢物中介分析结果_18_select.csv",  clear 
keep if pval_mediate<0.05 // drop 2
drop  coefci_mediate pval_direct coefci_direct pval_total coefci_total pval_ratio coefci_ratio pval_direct_inverse
merge m:m metabolism using  "D:\真菌分析\data_pre\metabolism_names.dta"
keep if _merge==3
drop _merge
rename fungi code_repeat
merge m:m code_repeat using "D:\真菌分析\data_pre\fungi_mapping_all_repeat.dta"
keep if _merge==3
keep g Compounds 物质 pval_mediate outcome coef_mediate rr rr_metabolismdisease coef_ratio
order g Compounds 物质 pval_mediate outcome coef_mediate rr rr_metabolismdisease coef_ratio
gen mark=g+Compounds
rename coef_mediate mediabe_v
rename rr rr_v
rename rr_metabolismdisease rr_metadis_v
rename coef_ratio coef_ratio_v
rename outcome outcome_v
rename pval_mediate pval_mediate_v
merge 1:1 mark using  "D:\真菌分析\data_pre\mediate_sig_discovery.dta"
keep if _merge==3 //2
order g Compounds 物质 outcome_d outcome_v pval_mediate_d pval_mediate_v coef_ratio_d coef_ratio_v rr_d rr_v rr_metadis_d rr_metadis_v
save  "D:\真菌分析\data_pre\mediate_sig_all.dta",replace

*----sensitive analysis--

import delimited "D:\真菌分析\Final results\All data\metabolism\代谢物中介分析结果_18_select.csv",  clear 
keep if pval_mediate<0.1 // drop 6
drop  coefci_mediate pval_direct coefci_direct pval_total coefci_total pval_ratio coefci_ratio pval_direct_inverse
merge m:m metabolism using  "D:\真菌分析\data_pre\metabolism_names.dta"
keep if _merge==3
drop _merge
rename fungi code_repeat
merge m:m code_repeat using "D:\真菌分析\data_pre\fungi_mapping_all_repeat.dta"
keep if _merge==3
keep g Compounds 物质 pval_mediate outcome coef_mediate rr rr_metabolismdisease coef_ratio
order g Compounds 物质 pval_mediate outcome coef_mediate rr rr_metabolismdisease coef_ratio
gen mark=g+Compounds
rename coef_mediate mediabe_v
rename rr rr_v
rename rr_metabolismdisease rr_metadis_v
rename coef_ratio coef_ratio_v
rename outcome outcome_v
rename pval_mediate pval_mediate_v
merge 1:1 mark using  "D:\真菌分析\data_pre\mediate_sig_discovery.dta"
keep if _merge==3 //2
order g Compounds 物质 outcome_d outcome_v pval_mediate_d pval_mediate_v coef_ratio_d coef_ratio_v rr_d rr_v rr_metadis_d rr_metadis_v
save  "D:\真菌分析\data_pre\mediate_sig_all_sensi.dta",replace

*----------dataset for the mediation results plot---

import delimited "D:\真菌分析\Final results\All data\metabolism\代谢物中介分析结果.csv",  clear 
keep if qval_mediate<0.05 // drop 1
keep if qval_mediate_inverse>0.05 // drop 0
drop  coefci_mediate pval_direct coefci_direct pval_total coefci_total pval_ratio  pval_direct_inverse
merge m:m metabolism using  "D:\真菌分析\data_pre\metabolism_names.dta"
keep if _merge==3
drop _merge
rename fungi code_repeat
merge m:m code_repeat using "D:\真菌分析\data_pre\fungi_mapping_all_repeat.dta"
keep if _merge==3
split (g), parse(;) 
drop if coef_ratio>1

preserve
keep g6  Compounds rr coef_ratio
rename g6 X
rename Compounds M
save  "D:\真菌分析\data_pre\media1.dta",replace 
restore

preserve
keep outcome  Compounds rr_metabolismdisease  
rename Compounds M
rename outcome Y
rename rr_metabolismdisease rr
save  "D:\真菌分析\data_pre\media2.dta",replace 
restore


use  "D:\真菌分析\data_pre\media1.dta",clear
rename rr rr_xm
merge m:m M  using "D:\真菌分析\data_pre\media2.dta"
order X M Y rr_xm rr 
save "D:\真菌分析\data_pre\mediation_input.dta",replace

*----------validate the results on the independent 2018 dataset--

import delimited "D:\真菌分析\Final results\All data\metabolism\代谢物中介分析结果_18_select.csv",  clear 
keep if pval_mediate<0.05 // drop 2
keep if pval_mediate_inverse>0.05 // drop 0
merge m:m metabolism using  "D:\真菌分析\data_pre\metabolism_names.dta"
keep if _merge==3
drop _merge
rename fungi code_repeat
merge m:m code_repeat using "D:\真菌分析\data_pre\fungi_mapping_all_repeat.dta"
keep if _merge==3
split (g), parse(;) 
drop if coef_ratio>1

preserve
keep g6  Compounds rr coef_ratio
rename g6 X
rename Compounds M
save  "D:\真菌分析\data_pre\media1_validate.dta",replace 
restore

preserve
keep outcome  Compounds rr_metabolismdisease  
rename Compounds M
rename outcome Y
rename rr_metabolismdisease rr
save  "D:\真菌分析\data_pre\media2_validate.dta",replace 
restore

use  "D:\真菌分析\data_pre\media1_validate.dta",clear
rename rr rr_xm
merge m:m M  using "D:\真菌分析\data_pre\media2_validate.dta"
order X M Y rr_xm rr 
save "D:\真菌分析\data_pre\mediation_input_validate.dta",replace

*-------Validate the results in the all 2018 data

import delimited "D:\真菌分析\Final results\All data\metabolism\代谢物中介分析结果_18_all.csv",  clear 
keep if pval_mediate<0.05 // drop 6
keep if pval_mediate_inverse>0.05 // drop 2
merge m:m metabolism using  "D:\真菌分析\data_pre\metabolism_names.dta"
keep if _merge==3
drop _merge
rename fungi code_repeat
merge m:m code_repeat using "D:\真菌分析\data_pre\fungi_mapping_all_repeat.dta"
keep if _merge==3
split (g), parse(;) 
drop if coef_ratio>1

preserve
keep g6  Compounds rr coef_ratio
rename g6 X
rename Compounds M
save  "D:\真菌分析\data_pre\media1_validate_2.dta",replace 
restore

preserve
keep outcome  Compounds rr_metabolismdisease  
rename Compounds M
rename outcome Y
rename rr_metabolismdisease rr
save  "D:\真菌分析\data_pre\media2_validate2.dta",replace 
restore


use  "D:\真菌分析\data_pre\media1_validate_2.dta",clear
rename rr rr_xm
merge m:m M  using "D:\真菌分析\data_pre\media2_validate2.dta"
order X M Y rr_xm rr 
save "D:\真菌分析\data_pre\mediation_input_validate2.dta",replace


*--------------3---All data analysis---------

*-----3.1-----PCOA analysis dataset----------

use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
tab region
keep SampleID region g1-g751
save "D:\真菌分析\data_pre\2015_all_fungi_forpca.dta",replace

import delimited "D:\真菌分析\Final results\All data\fungi_PCA.csv",  clear 
rename sampleid SampleID
save "D:\真菌分析\data_pre\2015_all_fungi_pc.dta",replace

*------3.8.2--mixed model---------

*----fungal features-categries--

tempname coef
tempfile res
postfile `coef' str200( micro var) float(n rr lul uul p)   using "`res'", replace

gen rural=0
replace rural=1 if city==1 // 1-city

foreach outcome of varlist pielou_evenness-shannon_entropy g1-g751  {
foreach var of varlist  rural smoke alcohol { 
    mixed `outcome' i.`var' age i.sex BMI ||region:, covariance(unstructured) 
test 1.`var'
 post `coef'   ("`outcome'") ("`var'")  (e(N)) (_b[1.`var'])  (_b[1.`var']-1.96*_se[1.`var'])  (_b[1.`var']+1.96*_se[1.`var'])  (chi2tail(1,(_b[1.`var']/_se[1.`var'])^2))  
 }
 }

  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\All data\associaiton analysis.xlsx",  firstrow(variables) sheet("diet_micro_categray") sheetreplace  
restore

*-----fungal features to phenotype------
tempname coef
tempfile res
postfile `coef' str200( var micro ) float(n rr lul uul p) str200(cov)  using "`res'", replace


foreach outcome of varlist  HDL_C LDL_C HbA1c glucose tg tc SBP DBP {
foreach var of varlist pielou_evenness-shannon_entropy g1-g751 { 
    mixed `outcome' `var' age i.sex BMI ||region:, covariance(unstructured) 
test `var'
 post `coef'   ("`outcome'") ("`var'")  (e(N)) (_b[`var'])  (_b[`var']-1.96*_se[`var'])  (_b[`var']+1.96*_se[`var'])  (chi2tail(1,(_b[`var']/_se[`var'])^2))  ("${cov`m'}")
 }
 }

  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\All data\associaiton analysis.xlsx",  firstrow(variables) sheet("fungi_phenotype") sheetreplace  
restore

*----fungi to diseases-------
tempname coef
tempfile res
postfile `coef' str200( micro var) float(n rr lul uul p) str200(cov)  using "`res'", replace


foreach outcome of varlist   t2d hyp dys {
foreach var of varlist pielou_evenness-shannon_entropy g1-g751 { 
    melogit `outcome' `var' age i.sex BMI ||region:, covariance(unstructured) 
test `var'
 post `coef'   ("`outcome'") ("`var'")  (e(N)) (exp(_b[`var']))  (exp(_b[`var']-1.96*_se[`var']))  (exp(_b[`var']+1.96*_se[`var']))  (chi2tail(1,(_b[`var']/_se[`var'])^2))  ("${cov`m'}")
 }
 }

  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\All data\associaiton analysis.xlsx",  firstrow(variables) sheet("fungi_diseases_model2") sheetreplace  
restore

*--------subgroup analyses based on diabetes risk score----后续在其它队列验证交互显著性

use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
count if g280>0 // 8830
drop g1-g751
merge 1:1 SampleID using "D:\真菌分析\data_pre\fungi_clr.dta"

su age
gen age_score=0
replace age_score=4 if age>=25 & age<=34 
replace age_score=8 if age>=35 & age<=39 
replace age_score=11 if age>=40 & age<=44 
replace age_score=12 if age>=45 & age<=49 
replace age_score=13 if age>=50 & age<=54 
replace age_score=15 if age>=55 & age<=59 
replace age_score=16 if age>=60 & age<=64 
replace age_score=18 if age>=65 
tab age_score t2d

gen BMI_score=0
replace BMI_score=1 if BMI>=22 & BMI<24 // 6 obs
replace BMI_score=3 if BMI>=24 & BMI<30 // 1468 obs
replace BMI_score=5 if BMI>=30 // 183 obs
tab BMI_score

gen wc_score=0
replace wc_score=3 if (wc>=75 & wc<80 & sex==1) | (wc>=70 & wc<75 & sex==2)
replace wc_score=5 if (wc>=80 & wc<85 & sex==1) | (wc>=75 & wc<80 & sex==2)
replace wc_score=7 if (wc>=85 & wc<90 & sex==1) | (wc>=80 & wc<85 & sex==2)
replace wc_score=8 if (wc>=90 & wc<95 & sex==1) | (wc>=85 & wc<90 & sex==2)
replace wc_score=10 if (wc>=95 & sex==1) | (wc>=90 & sex==2)
tab wc_score

gen sys_score=0
replace sys_score=1 if SBP>=110 & SBP<119
replace sys_score=3 if SBP>=120 & SBP<129
replace sys_score=6 if SBP>=130 & SBP<139
replace sys_score=7 if SBP>=140 & SBP<149
replace sys_score=8 if SBP>=150 & SBP<159
replace sys_score=10 if SBP>=160

gen sex_score=0
replace sex_score=2 if sex==1
gen risk=age_score+BMI_score+wc_score+sys_score+sex_score

spearman g280 risk
spearman g280 BMI_score
spearman g280 wc_score

preserve
keep SampleID age_score-risk g280 t2d
rename t2d outcome
save "D:\真菌分析\data_pre\riskscore_g280prediciton.dta",replace
restore

preserve
keep SampleID age_score-risk  t2d
rename t2d outcome
save "D:\真菌分析\data_pre\riskscore_prediciton.dta",replace
restore

egen std_g280=std(g280)
egen std_g702=std(g702)
egen std_risk=std(risk)
egen std_observed_features=std(observed_features)
melogit t2d c.g280##c.risk   ||region:, covariance(unstructured) or // 交互项显著, p=0.012
melogit t2d c.risk##c.observed_features   ||region:, covariance(unstructured) or // 交互项不显著
melogit t2d c.std_risk##c.std_g702   ||region:, covariance(unstructured) or // 交互项不显著
