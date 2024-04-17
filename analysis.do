
*-----0-dataset describe---------------

*---basic characters--

use "D:\真菌分析\data_pre\fugi_15all_g_diversity.dta",clear // 10833
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

import excel "D:\真菌分析\rawdata\raw_fungi_mapping.xlsx", sheet("repeat") firstrow clear 
save  "D:\真菌分析\data_pre\raw_fungi_mapping.dta",replace

use "D:\真菌分析\data_pre\fugi_15all_g_diversity.dta",clear // 10833
merge m:m SampleID using  "D:\真菌分析\data_pre\2015_all_phenotype.dta"
keep if _merge==3 // 10695
keep SampleID sex
append using "D:\真菌分析\data_pre\repeatsex.dta"
tab  sex

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



*---------1-进化树绘制-----analysis_1.R (all dataset,top 150)--

*---------2--分布描述------analysis_1.R (all dataset, all genera)

*---core microbiota--

*---change analysis--

*---注意研究肠道真菌与膳食对疾病结局的交互作用

*----------------2-Repeat datasets compare-----

*--PCA-(Miv pipline,10% prevelance)
*--门、核心菌群的变化-(microbiome pipline,all genera)

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


*--2.1.1--DMM clusters and phenotype--后续补充18年的比较

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

*---baseline(类似于前瞻性分析)

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

*---prospective analysis-----

*---incident cases

use "D:\真菌分析\data_pre\repeatDMM_cluster.dta",clear
tab time
keep if time==15
save "D:\真菌分析\data_pre\repeatDMM_15cluster.dta",replace


use  "D:\真菌分析\data_pre\2018_all_phenotype.dta",clear
order idind hhid commid 
renvars age-fecal_time, pref(follow_)
merge 1:1 idind using   "D:\真菌分析\data_pre\2015_all_phenotype.dta"
drop _merge
rename idind Idind
merge m:m Idind using "D:\真菌分析\data_pre\repeatDMM_15cluster.dta",force
keep if _merge==3
tab t2d,missing // 235
tab follow_t2d,missing // 285
gen duration=follow_age-age
su duration,de

preserve
drop if t2d==1
gen incident_t2d=0
replace incident_t2d=1 if follow_t2d==1
tab incident_t2d // 108
tab cluster incident_t2d
logit incident_t2d ib3.cluster,or
restore

*---repeat measurement analysis-----重要

*----incidnet cases baseline and follow-up microbial composition comparation---

use  "D:\真菌分析\data_pre\2018_all_phenotype.dta",clear
order idind hhid commid 
renvars age-fecal_time, pref(follow_)
merge 1:1 idind using   "D:\真菌分析\data_pre\2015_all_phenotype.dta"
drop _merge
rename idind Idind
merge m:m Idind using "D:\真菌分析\data_pre\repeatDMM_15cluster.dta",force
keep if _merge==3


gen incident_t2d=999
replace  incident_t2d=0 if t2d==0 & follow_t2d==0
replace incident_t2d=1 if t2d==0 & follow_t2d==1
tab incident_t2d // 108

gen incident_dys=999
replace  incident_dys=0 if dys==0 & follow_dys==0
replace incident_dys=1 if dys==0 & follow_dys==1
tab incident_dys // 171

gen incident_hyp=999
replace  incident_hyp=0 if hyp==0 & follow_hyp==0
replace incident_hyp=1 if hyp==0 & follow_hyp==1
tab incident_hyp // 347

gen incident_case=999
replace  incident_case=0 if incident_t2d==0 & incident_dys==0 &incident_hyp==0
replace incident_case=1 if  incident_t2d==1 | incident_dys==1 | incident_hyp==1
tab incident_case // 562

keep Idind incident_t2d incident_dys incident_hyp incident_case
save "D:\真菌分析\data_pre\incident_diseases.dta",replace

*------2018 microbiome and incident diseases---

use "D:\真菌分析\data_pre\fugi_repeat.dta",clear
keep if time==18
merge 1:1 Idind using "D:\真菌分析\data_pre\incident_diseases.dta"
drop _merge time
sort incident_t2d
by incident_t2d:su g220
save "D:\真菌分析\data_pre\2018_micro_incidentdiseases.dta",replace

*------2015 microbiome and incident diseases---

use "D:\真菌分析\data_pre\fugi_repeat.dta",clear
keep if time==15
merge 1:1 Idind using "D:\真菌分析\data_pre\incident_diseases.dta"
drop _merge time
sort incident_t2d
by incident_t2d:su g220
save "D:\真菌分析\data_pre\2015_micro_incidentdiseases.dta",replace


*ref-GEE-https://www.stata.com/features/overview/generalized-estimating-equations/ In this example, we consider a probit model in which we wish to model whether a worker belongs to the union based on the person's age and whether they are living outside of an SMSA. The people in the study appear multiple times in the dataset (this type of panel dataset is commonly referred to as a longitudinal dataset), and we assume that the observations on a given person are more correlated than those between different persons.

*-ref：中介vs 交互-https://www.iikx.com/news/statistics/1737.html

*--交互作用应尽量检验-毕竟我们想得到不依赖于变量差异的暴露-结局关联
*--中介效应在只关注暴露和结局的关系时可忽略，但在讨论可能机制时要考虑

use  "D:\真菌分析\data_pre\2018_all_phenotype.dta",clear
gen time=18
append using "D:\真菌分析\data_pre\2015_all_phenotype.dta"
tab time,missing
replace time=15 if time==.
rename idind Idind
merge m:m Idind using "D:\真菌分析\data_pre\repeatDMM_cluster.dta",force
keep if _merge==3
order time Time
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
save   "D:\真菌分析\data_pre\cluster_gee.dta",replace


xtset code_id time
xtgee  t2d  age i.sex  , family(binomial) link(logit) corr(ind) vce(robust) ef
xtgee  t2d ib3.cluster age i.sex , family(binomial) link(logit) corr(ind) vce(robust) ef

xtgee  t2d ib3.cluster age i.sex BMI , family(binomial) link(logit) corr(ind) vce(robust) ef


inspect age sex BMI smoke alcohol city income fruit  pork veg milk fish antibiotic_current MET
global cov  age i.sex   i.smoke i.alcohol i.city income fruit  pork veg milk fish i.antibiotic_current MET
order code_id Idind
xtset code_id time
global cov  age i.sex   i.smoke i.alcohol i.city income fruit  pork veg milk fish i.antibiotic_current MET
xtset code_id time
xtgee  t2d ib3.cluster ${cov}, family(binomial) link(logit) corr(ind) robust  ef

xtgee t2d ib3.cluster##c.BMI  age i.sex   i.smoke i.alcohol i.city income fruit  pork veg milk fish i.antibiotic_current MET , family(binomial) link(logit) corr(ind) robust  ef


xtgee BMI ib3.cluster  age i.sex   i.smoke i.alcohol i.city income fruit  pork veg milk fish i.antibiotic_current MET, family(gaussian) link(identity) corr(ind)   robust
xtgee BMI ib4.cluster  age i.sex   i.smoke i.alcohol i.city income fruit  pork veg milk fish i.antibiotic_current MET , family(gaussian) link(identity) corr(ind)   robust

xtgee t2d ib3.cluster  age i.sex   i.smoke i.alcohol i.city income fruit  pork veg milk fish i.antibiotic_current MET , family(binomial) link(logit) corr(ind)  robust ef
xtgee t2d ib4.cluster  age i.sex   i.smoke i.alcohol i.city income fruit  pork veg milk fish i.antibiotic_current MET , family(binomial) link(logit) corr(ind)  robust ef

xtgee t2d ib3.cluster  age i.sex  BMI i.smoke i.alcohol i.city income fruit  pork veg milk fish i.antibiotic_current MET , family(binomial) link(logit) corr(ind)  robust ef
xtgee t2d ib4.cluster  age i.sex  BMI i.smoke i.alcohol i.city income fruit  pork veg milk fish i.antibiotic_current MET , family(binomial) link(logit) corr(ind)  robust ef


xtgee hyp ib3.cluster  age i.sex   i.smoke i.alcohol i.city income fruit dveg milk fish i.antibiotic_current , family(binomial) link(logit) corr(exchangeable) robust  ef
xtgee hyp ib4.cluster  age i.sex   i.smoke i.alcohol i.city income fruit dveg milk fish i.antibiotic_current , family(binomial) link(logit) corr(exchangeable) robust  ef
xtgee dys ib3.cluster  age i.sex   i.smoke i.alcohol i.city income fruit dveg milk fish i.antibiotic_current , family(binomial) link(logit) corr(exchangeable) robust  ef
xtgee dys ib4.cluster  age i.sex   i.smoke i.alcohol i.city income fruit dveg milk fish i.antibiotic_current , family(binomial) link(logit) corr(exchangeable) robust  ef

*---cluster3及cluster4与糖尿病及肥胖同时相关，肥胖可能是中介作用

*-----排除疾病的影响-----

*--随访期间血糖都正常的人群中看fungi与treat的关联--

*-------GEE modle for fungi-diseases associations---更新

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

*----Validate the fungi and disease associaitons based on the GEE analysis---

*可能引入新的繁琐：饮食为什么不验证，这里可以在验证队列在后续部分只报告sacc的GEE分析结果


*----2.1.3-----DMM clusters and metabolism----

*---pca analysis--

use "D:\真菌分析\data_pre\repeatDMM_15cluster.dta",clear
merge m:m Idind using  "D:\真菌分析\data_pre\metabolism_15.dta"
keep if _merge==3
rename cluster_15 cluster
keep SampleID-mwxq04_N cluster
save "D:\真菌分析\data_pre\15_metabolism_cluster.dta",replace

*--predict age--

use "D:\CNHS\data\phenotype data\2015update_20220426\chns_2015_westlake.dta",clear
duplicates drop Idind,force // 13315
merge m:m Idind using  "D:\真菌分析\data_pre\metabolism_15.dta"
keep if _merge==3
drop if age==. // 996 obs
drop _merge
merge m:m Idind using "D:\真菌分析\data_pre\repeatDMM_15cluster.dta"
keep if _merge==3
tab cluster

preserve
keep if cluster==1
keep SampleID-mwxq04_N age 
save "D:\真菌分析\data_pre\metabolism_15_predict_cluster1.dta",replace
restore

preserve
keep if cluster==2
keep SampleID-mwxq04_N age 
save "D:\真菌分析\data_pre\metabolism_15_predict_cluster2.dta",replace
restore

preserve
keep if cluster==3
keep SampleID-mwxq04_N age 
save "D:\真菌分析\data_pre\metabolism_15_predict_cluster3.dta",replace
restore

preserve
keep if cluster==4
keep SampleID-mwxq04_N age 
save "D:\真菌分析\data_pre\metabolism_15_predict_cluster4.dta",replace
restore

*------Phenotype dataset---

use "D:\真菌分析\data_pre\2015_all_phenotype.dta",clear
rename idind Idind
duplicates drop Idind,force // 13315
merge m:m Idind using  "D:\真菌分析\data_pre\metabolism_idcluster.dta"
keep if _merge==3
keep SampleID age BMI t2d district city smoke alcohol MET city_score hpc wc HDL_C LDL_C ins HbA1c glucose tg tc sex SBP DBP hyp dys pret2d cluster hyp dys
order SampleID t2d
save "D:\真菌分析\data_pre\2015_metabolism_phenotype.dta",replace


*------Predict clusters-----

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

*-对比基线相同group trans到cluster3 或cluster4时代谢组的评分差异

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
 

*-对比trans到cluster3 或cluster4组别的人在基线cluster（1或者2）的代谢组评分差异





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

*--------------中介分析数据准备----------可以把pre-t2d引入，呈现三个方向整体的结果（网络图或者trans图）--说明真菌和代谢图谱共变

*--只有cluster3和cluster4-dataset 1：fungi cluster-metabolism module-diseases--(only based on the cluster 3 and cluster 4)

*---15 和18有多少人重复出现在cluster3或cluster4中？

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

*-----只用cluster3和cluster4---dataset2: individual fungi (that driven the clusters)-metabolisms(sig results)--diseases--

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


*========真菌标记物与代谢物的关联结果mapping===

import excel "D:\真菌分析\Plot\data\代谢组分析.xlsx", sheet("fungi_metabolism") firstrow clear 
merge m:m metabolism using  "D:\真菌分析\data_pre\metabolism_names.dta"
keep if _merge==3
drop _merge
rename fungi code_repeat
merge m:m code_repeat using "D:\真菌分析\data_pre\fungi_mapping_all_repeat.dta"
keep if _merge==3
split (g), parse(;) 
order g6 Compounds rr p fold VIP
sort p
keep if p<0.05

*========代谢物与疾病的关联结果mapping===

import excel "D:\真菌分析\Plot\data\代谢组分析.xlsx", sheet("metabolism_disease") firstrow clear 
drop _merge
merge m:m metabolism using  "D:\真菌分析\data_pre\metabolism_names.dta"
keep if _merge==3
sort p

*-----------什么因素决定了cluster间的tranform？-------------




*----------data transform for heatmap plot------

use "D:\真菌分析\data_pre\metabolism_predict_cluster34.dta",clear
foreach v of varlist MEDN0011-MEDP2672 {
           egen s_`v' = std(`v')
		   drop `v'
		   rename s_`v' `v'
         }
save "D:\真菌分析\data_pre\metabolism_predict_cluster34_.dta",replace



*----2.1.4--DMM clusters and fungi/bacteria diversity----

*--fungi diversity

use "D:\真菌分析\data_pre\repeatDMM_cluster.dta",clear
keep if time==15
merge m:m Idind using "D:\真菌分析\data_pre\fugi_15_diversity.dta",force
keep if _merge==3
save "D:\真菌分析\data_pre\cluster_fungi_15diversity.dta",replace

use "D:\真菌分析\data_pre\repeatDMM_cluster.dta",clear
keep if time==18
merge m:m Idind using "D:\真菌分析\data_pre\fugi_18_diversity.dta",force
keep if _merge==3
save "D:\真菌分析\data_pre\cluster_fungi_18diversity.dta",replace

*--bacteria diversity

use "D:\真菌分析\data_pre\repeatDMM_cluster.dta",clear
keep if time==15
merge m:m Idind using "D:\真菌分析\data_pre\bac_15_diversity.dta",force
keep if _merge==3
save "D:\真菌分析\data_pre\cluster_bac_15diversity.dta",replace

use "D:\真菌分析\data_pre\repeatDMM_cluster.dta",clear
keep if time==18
merge m:m Idind using "D:\真菌分析\data_pre\bac_18_diversity_raw.dta",force
keep if _merge==3
save "D:\真菌分析\data_pre\cluster_bac_18diversity.dta",replace

*--DMM change pattern-fungi diversity

use "D:\真菌分析\data_pre\repeatDMM_indcluster.dta",clear
tab cluster_15 cluster_18
egen group=group(cluster_15 cluster_18)
merge m:m Idind using "D:\真菌分析\data_pre\fugi_18_diversity.dta",force
keep if _merge==3
tab group,missing
gen Group=.
replace Group=1 if group<5
replace Group=2 if group>4 & group<9
replace Group=3 if group>8 & group<13
replace Group=4 if group>12 & group<17
tab Group
save "D:\真菌分析\data_pre\clusterchange_fungidiversity.dta",replace
 
use "D:\真菌分析\data_pre\repeatDMM_indcluster.dta",clear
tab cluster_15 cluster_18
order cluster_18 cluster_15
egen group=group(cluster_15 cluster_18)
merge m:m Idind using "D:\真菌分析\data_pre\fugi_18_diversity.dta",force
keep if _merge==3
tab group,missing
tab cluster_18 cluster_15
gen Group=""
replace Group="a" if group<5
replace Group="b" if group>4 & group<9
replace Group="c" if group>8 & group<13
replace Group="d" if group>12 & group<17
tab Group
save "D:\真菌分析\data_pre\clusterchange_fungidiversity1.dta",replace
 
*--DMM change pattern-phenotype---

use "D:\CNHS\data\phenotype data\2015update_20220426\chns_2015_westlake.dta",clear
duplicates drop Idind,force // 13315
merge m:m Idind using  "D:\真菌分析\data_pre\repeatDMM_indcluster.dta"
keep if _merge==3
tab cluster_15 cluster_18
egen group=group(cluster_15 cluster_18)
sort group
by group:su age dbp BMI1 
tab group U24

 
*----2.1.5----fungi cluster and bilogical age---

*--predict data prepare--

use "D:\CNHS\data\phenotype data\2015update_20220426\chns_2015_westlake.dta",clear
duplicates drop Idind,force // 13315
merge m:m Idind using  "D:\真菌分析\data_pre\metabolism_15.dta"
keep if _merge==3
drop if age==. // 996 obs
keep SampleID-mwxq04_N age 
save "D:\真菌分析\data_pre\metabolism_15_predict.dta",replace
drop if age>80
save "D:\真菌分析\data_pre\metabolism_15_predict1.dta",replace

*---fungi cluster and bilogical age--

use "D:\真菌分析\data_pre\repeatDMM_cluster.dta",clear
keep if time==15
merge m:m Idind using "D:\真菌分析\data_pre\metabolism_15.dta",force
keep if _merge==3
keep SampleID cluster
save "D:\真菌分析\data_pre\metabolism_15_cluster.dta",replace

import excel "D:\真菌分析\Final results\microbiota metabolism\metabolism_predictage.xlsx", sheet("Sheet1") firstrow clear 
merge 1:1 SampleID  using "D:\真菌分析\data_pre\metabolism_15_cluster.dta",force
keep if _merge==3
sort cluster
by cluster:su predict actual
spearman predict actual // 0.658

preserve
keep if cluster==1
spearman predict actual // 0.679
restore

preserve
keep if cluster==2
spearman predict actual // 0.641
restore

preserve
keep if cluster==3
spearman predict actual // 0.535
restore

preserve
keep if cluster==4
spearman predict actual // 0.648
restore

*---fungi cluster and bac-fungi interaction--

use "D:\真菌分析\data_pre\repeatDMM_cluster.dta",clear
keep if time==15
merge 1:1 Idind using "D:\真菌分析\data_pre\fugibac_15_g_filterdata.dta"
keep if _merge==3
tab cluster,missing
drop Idind time
save "D:\真菌分析\data_pre\clustr_bacfungi_interaction.dta",replace

*--2.2-----Contribution of microbial instability (baseline phenotype) ------

*---后续可能需要进一步过滤差异菌种后再做分析--这点挺好，同一家庭的人可以对比看不同年龄的波动情况

use "D:\CNHS\data\phenotype data\2015update_20220426\chns_2015_westlake.dta",clear

use "D:\真菌分析\data_pre\2015_all_phenotype.dta",clear
*duplicates drop Idind,force // 
rename idind Idind
merge m:m Idind using "D:\真菌分析\data_pre\repeat_dis1.dta"
keep if _merge==3
drop _merge
merge m:m Idind using "D:\真菌分析\data_pre\fugi_15_diversity.dta",force
keep if _merge==3
spearman  age dis // -0.058 (0.01)
spearman  MET dis // 0.07 (0.0017)
spearman  BMI dis // -0.0534 (0.0206)
spearman faith_pd dis // 没有关系，此前细菌的研究报道有关系

su age,de
xtile q_age = age, nq(4)
sort q_age
by q_age: su dis

preserve
keep age q_age sex dis
save "D:\真菌分析\data_pre\age_dis.dta",replace
restore

sort sex
by sex: su dis

sort pret2d
by pret2d:su dis //sig 

sort famine
by famine:su dis

tempname coef
 tempfile res
 postfile `coef' str200(phenotype dis) float( n corr  p)  using "`res'", replace
qui foreach g of varlist age t2d city education marrige smoke alcohol xinjigengse zhongfeng cancer MET city_score income hpc wc HDL_C LDL_C ins HbA1c glucose tg tc sex BMI SBP DBP hyp dys pret2d famine rice-pastes gut_disease fuxie antibiotic_current antibiotic_6month probiotics kangyan_med kangsuan_med weisuan_med yogurt_drink changdao_shoushu pet{
spearman `g' dis
 post `coef'   ("`g'")  ("dis") (r(N)) (r(rho))  (r(p))
}


 postclose `coef'
 preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\Repeated data\基线因素与真菌波动性的关联.xlsx",  firstrow(variables) sheet("spearman") sheetreplace  
restore

 
*merge m:m Idind using "D:\真菌分析\data_pre\fugi_15_diversity.dta"
*keep if _merge==3
*spearman  shannon dis // p=0.45 与baseline的真菌多样性无关

*---更新！基于重复测量的均值计算相关性------

*--表型均值与真菌波动的关联---

use  "D:\真菌分析\data_pre\2015_all_phenotype.dta",clear
rename idind Idind
renvars age t2d smoke alcohol MET hpc-tc BMI-pret2d rice-pastes , postfix(_15)  // 批量增加后缀
save  "D:\真菌分析\data_pre\2015_all_phenotype_recode.dta",replace
 
use  "D:\真菌分析\data_pre\2018_all_phenotype.dta",clear
rename idind Idind
*renvars age t2d smoke alcohol MET hpc-tc BMI-pret2d rice-pastes , postfix(_18)  // 批量增加后缀
*save  "D:\真菌分析\data_pre\2018_all_phenotype_recode.dta",replace
merge 1:1 Idind using  "D:\真菌分析\data_pre\2015_all_phenotype_recode.dta"
keep if _merge==3
drop _merge
merge 1:1 Idind using "D:\真菌分析\data_pre\repeat_dis1.dta"
keep if _merge==3
drop _merge

 foreach var of varlist age  MET hpc-tc BMI SBP-dys rice-pastes {
  egen `var'_mean = rmean(`var' `var'_15 )
}

foreach var of varlist   t2d smoke alcohol pret2d hyp dys  {
  gen `var'_change=0
  replace `var'_change=1 if (`var'!=`var'_15)&(`var'!=9 |`var'!=.| `var'_15!=9 | `var'_15!=.)
}
  
tempname coef
 tempfile res
 postfile `coef' str200(phenotype dis) float( n corr  p)  using "`res'", replace
qui foreach g of varlist age_mean-dys_change sex city city_score income{
spearman `g' dis
 post `coef'   ("`g'")  ("dis") (r(N)) (r(rho))  (r(p))
}


 postclose `coef'
 preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\Repeated data\基线因素与真菌波动性的关联.xlsx",  firstrow(variables) sheet("表型均值") sheetreplace  
restore

preserve
keep age_mean dis
save  "D:\真菌分析\data_pre\age_dis.dta",replace
restore

twoway (scatter dis age_mean) (lfit dis age) ///
       , xtitle("age") ytitle("dis") legend(off)
su age_mean,de
drop age
rename age_mean age
gen age_group=0
replace age_group=1 if age>30 & age<=40
replace age_group=2 if age>40 & age<=50
replace age_group=3 if age>50 & age<=60
replace age_group=4 if age>60 & age<=70
replace age_group=5 if age>70 & age<=80
replace age_group=6 if age>80
tab age_group
sort age_group
by age_group:su dis

su age,de
xtile q_age = age, nq(4)
sort q_age
by q_age: su dis

gen age_group1=0
replace age_group1=1 if age>=65
sort age_group1
by age_group1: su dis



*--真菌均值与真菌波动的关联---

use  "D:\真菌分析\data_pre\2015_repeat_datasets.dta",clear
renvars pielou_evenness-shannon g2-g633 , postfix(_15)  // 批量增加后缀
merge 1:1 Idind using "D:\真菌分析\data_pre\2018_repeat_datasets.dta"
keep if _merge==3
drop _merge
merge 1:1 Idind using "D:\真菌分析\data_pre\repeat_dis1.dta"
keep if _merge==3
drop g74 g104 g223 g292 g398 g559 g573

 foreach var of varlist pielou_evenness-shannon g2-g633 {
  egen `var'_mean = rmean(`var' `var'_15 )
}

tempname coef
 tempfile res
 postfile `coef' str200(phenotype dis) float( n corr  p)  using "`res'", replace
qui foreach g of varlist g2_mean-g633_mean{
spearman `g' dis
 post `coef'   ("`g'")  ("dis") (r(N)) (r(rho))  (r(p))
}


 postclose `coef'
 preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\Repeated data\基线因素与真菌波动性的关联.xlsx",  firstrow(variables) sheet("真菌均值") sheetreplace  
restore

*-------生活方式改变及基本表型与肠道真菌波动性的关联-----

use  "D:\真菌分析\data_pre\2018_all_phenotype.dta",clear
rename idind Idind
merge 1:1 Idind using  "D:\真菌分析\data_pre\2015_all_phenotype_recode.dta"
keep if _merge==3
drop _merge
merge 1:1 Idind using "D:\真菌分析\data_pre\repeat_dis1.dta"
keep if _merge==3
drop _merge 

 foreach var of varlist age  MET hpc-tc BMI SBP DBP rice-pastes {
  gen `var'_change = `var' -`var'_15 
}

foreach var of varlist   t2d smoke alcohol pret2d hyp dys  {
  gen `var'_change=0
  replace `var'_change=1 if (`var'!=`var'_15)&(`var'!=9 |`var'!=.| `var'_15!=9 | `var'_15!=.)
}
  
tempname coef
 tempfile res
 postfile `coef' str200(phenotype dis) float( n corr  p)  using "`res'", replace
qui foreach g of varlist age_change-dys_change age sex education city city_score income{
spearman `g' dis
 post `coef'   ("`g'")  ("dis") (r(N)) (r(rho))  (r(p))
}


 postclose `coef'
 preserve
use "`res'", clear
export excel using "D:\真菌分析\revised_results\基线因素波动与真菌波动性的关联.xlsx",  firstrow(variables) sheet("表型均值") sheetreplace  
restore


*------重要-bad trans vs good trans---基线相同的cluster看trans的差异

use "D:\真菌分析\data_pre\repeatDMM_indcluster.dta",clear
tab cluster_15 cluster_18
egen group=group(cluster_15 cluster_18)
tab cluster_15 group

use  "D:\真菌分析\data_pre\2018_all_phenotype.dta",clear
rename idind Idind
merge 1:1 Idind using  "D:\真菌分析\data_pre\2015_all_phenotype_recode.dta"
keep if _merge==3
drop _merge
merge 1:1 Idind using "D:\真菌分析\data_pre\repeatDMM_indcluster.dta"
keep if _merge==3
drop _merge
tab cluster_15 cluster_18
egen group=group(cluster_15 cluster_18)

 foreach var of varlist age  MET hpc-tc BMI SBP-dys rice-pastes {
  egen `var'_mean = rmean(`var' `var'_15 )
}

foreach var of varlist   t2d smoke alcohol pret2d hyp dys  {
  gen `var'_change=0
  replace `var'_change=1 if (`var'!=`var'_15)&(`var'!=9 |`var'!=.| `var'_15!=9 | `var'_15!=.)
}

preserve
keep if cluster_15==3
table1, by(cluster_18) vars(age_15 contn \ sex cat\ MET_15 contn \ fruit_15 contn  )format(%8.2f) onecol test pdp(2) saving (D:\真菌分析\trans_phenotype.xlsx, replace)
restore

*------------fungi and fungi diversity----

use  "D:\真菌分析\data_pre\2015_repeat_datasets.dta",clear

tempname coef
 tempfile res
 postfile `coef' str200(diversity genus) float( n corr  p)  using "`res'", replace
qui foreach diversity of varlist pielou_evenness-shannon {
foreach phenotype of varlist g2-g633 {
 spearman `diversity'  `phenotype' 
 post `coef'    ("`diversity'")  ("`phenotype'") (r(N)) (r(rho))  (r(p))
}
}
 postclose `coef'
 preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\All data\真菌_真菌多样性.xlsx",  firstrow(variables) sheet("repeat_raw") sheetreplace  
restore

*--name mappping

import excel "D:\真菌分析\Final results\All data\真菌_真菌多样性.xlsx", sheet("raw") firstrow clear 
rename genus micro
merge m:m micro using  "D:\真菌分析\data_pre\all_g_mapping.dta"
keep if _merge==3
sort diversity corr
export excel using "D:\真菌分析\Plot\data\真菌_真菌多样性.xlsx",  firstrow(variables) sheet("all") sheetreplace  


*--2.4---microbiota-fungi-metabolism paire-wise associations ------

*---microbiota-fungi--

use "D:\真菌分析\data_pre\bac_dis.dta",clear
merge 1:1 Idind using "D:\真菌分析\data_pre\repeat_dis1.dta"
keep if _merge==3
spearman dis dis_bac // no associaiton

*---microbiota-metabolism--

use "D:\真菌分析\data_pre\repeat_metabolismbc.dta",clear
merge 1:1 Idind using "D:\真菌分析\data_pre\bac_dis.dta"
keep if _merge==3
spearman dis dis_bac // no associaiton

*---fugi-metabolism--

use "D:\真菌分析\data_pre\repeat_metabolismbc.dta",clear
rename dis dis_metabolism
merge 1:1 Idind using "D:\真菌分析\data_pre\repeat_dis1.dta"
keep if _merge==3
spearman dis_metabolism dis // no associaiton

*-----2.5-------fungi change calculate---------------

import excel "D:\真菌分析\rawdata\fungirepeatdata.xlsx", sheet("g_raw") firstrow clear 
keep SampleID time Idind g2	g4	g18	g19	g44	g63	g64	g73	g87	g106	g112	g113	g114	g115	g116	g123	g126	g146	g186	g187	g188	g190	g194	g200	g201	g207	g212	g216	g217	g220	g221	g223	g225	g226	g227	g230	g237	g243	g265	g267	g278	g279	g280	g283	g284	g302	g303	g317	g318	g327	g328	g331	g397	g398	g399	g461	g496	g498	g505	g507	g528	g530	g538	g564	g567	g582	g583	g584	g585	g597	g609	g612	g617	g622	g623	g626	g633
 foreach micro of varlist g2-g633{
   replace `micro'=`micro'+1
   } 
egen g_all=rowtotal(g2-g633)
 foreach micro of varlist g2-g633{
   replace `micro'=`micro'/g_all
   } 
export excel using "D:\真菌分析\rawdata\fungirepeatdata.xlsx",  firstrow(variables) sheet("g_rela") sheetreplace  
   
order Idind time 
sort Idind time
 foreach micro of varlist g2-g633 {
	gen d_`micro'=(`micro'[_n]-`micro'[_n-1])/`micro'[_n-1]*100
}

keep if time==18

 foreach micro of varlist d_g2-d_g633 {
	egen meanchange_`micro'=median(`micro')
}
keep meanchange_d_g2-meanchange_d_g633  


import delimited "D:\真菌分析\Final results\Repeated data\KW_relative.csv",  clear 
save "D:\真菌分析\data_pre\repeat_kw_test.dta",replace

import excel "D:\真菌分析\Final results\Repeated data\repeat_差异分析结果.xlsx", sheet("Sheet2") firstrow clear 
merge 1:1 micro using "D:\真菌分析\data_pre\repeat_kw_test.dta"
drop _merge
rename micro code_repeat
merge m:m code_repeat using "D:\真菌分析\data_pre\fungi_mapping_all_repeat.dta"
keep if _merge==3
rename code_all micro
drop _merge
merge m:m micro using "D:\真菌分析\data_pre\all_g_mapping.dta"
keep if _merge==3

*----------2.6----Association of core fungi and metabolism ------

*---2015 linear model---


*----------2.7---------个性化营养分析-------

*---plant-based diet 



*--------------3---All data analysis---------

*--总的数据和重复的数据有多少重叠
*--真菌和细菌有多少重合
*--缺失值填补--

*-------fungi diversity and age----

use "D:\真菌分析\data_pre\fugi_15all_g_diversity.dta",clear
keep SampleID pielou_evenness-shannon_entropy
save "D:\真菌分析\data_pre\fugi_15all_g_diversity.dta"


*--真菌&细菌重复数据--

use "D:\真菌分析\data_pre\fugi_15all_g_diversity.dta",clear
merge 1:1 SampleID using  "D:\真菌分析\data_pre\bac_15all_g_diversity.dta" // 10778 keep (55 drop for fungi,320 for bac)

*-2015 repeated data及all data 重复数据 (剔除重复测量数据)--

*-fungi-

use "D:\真菌分析\data_pre\fugi_15all_g_diversity.dta",clear // 10833
merge m:m SampleID using  "D:\真菌分析\data_pre\2015_all_phenotype.dta"
keep if _merge==3 // 10695
spearman age faith_pd
rename idind Idind
drop _merge
merge m:m Idind using "D:\真菌分析\data_pre\repeatDMM_cluster.dta"
keep if _merge==3
drop _merge
keep if time==15 // only 3 overlaped, 后续分析删除
keep Idind
gen idind=Idind
save  "D:\真菌分析\data_pre\delete_fungi_id.dta",replace

*-bacteria-同fungi

*-----3.1----data pre----

*---fungi---

use "D:\真菌分析\data_pre\fugi_15all_g_diversity.dta",clear // 10833
merge m:m SampleID using  "D:\真菌分析\data_pre\2015_all_phenotype.dta"
keep if _merge==3 // 10695
save "D:\CNHS\data\Data share\for CDC\2015_all_fungi_dataset",replace

tab region,missing 

 foreach micro of varlist g2-g761 {
	count if `micro'==0 
	if r(N)>9626 {	
	drop  `micro' 
	}
}
tab famine,missing
logit t2d g280 age i.sex ,or // 
tab t2d,missing
tab diabetes_med,missing
*drop if diabetes_med==1
keep if diabetes_med==1 | t2d==0
logit t2d g280 age i.sex ,or // 
tab t2d,missing
save "D:\真菌分析\data_pre\2015_all_dataset.dta",replace

*--前瞻性分析---

use "D:\真菌分析\data_pre\2015_all_dataset.dta",clear
tab t2d
rename t2d t2d_base
drop _merge
merge 1:1 idind using  "D:\真菌分析\data_pre\2018_all_phenotype.dta" 
keep if _merge==3
drop if t2d_base==1 // 974
tab t2d // 413 cases 
logit t2d shannon_entropy age i.sex ,or // nosig
logit t2d observed_features age i.sex ,or // nosig
logit t2d  age i.sex ,or // nosig


*---bacteria--

use "D:\真菌分析\data_pre\bac_15all_g_diversity.dta",clear
merge m:m SampleID using  "D:\真菌分析\data_pre\2015_all_phenotype.dta"
keep if _merge==3 // 10958
save "D:\CNHS\data\Data share\for CDC\2015_all_bacteria_dataset",replace

replace region="K" if region=="k"
replace region="L" if region=="l"
tab region
replace region="Beijing" if region=="J"
replace region="Liaoning" if region=="A"
replace region="Heilongjiang" if region=="B"
replace region="Shanghei" if region=="K"
replace region="Jiangsu" if region=="C"
replace region="Zhengjiang" if region=="M"
replace region="Shandong" if region=="D"
replace region="Henan" if region=="E"
replace region="Hubei" if region=="F"
replace region="Hunan" if region=="G"
replace region="Guangxi" if region=="H"
replace region="Guizhou" if region=="I"
replace region="Yunnan" if region=="N"
replace region="Chongqin" if region=="L"
replace region="Shanxi" if region=="O"
gen District="South"
replace District="North" if region=="Beijing" | region=="Liaoning" | region=="Heilongjiang"  | region=="Shandong" | region=="Henan"  | region=="Shanxi" //1093
tab District // north 4445, south 6420
tab region
drop _merge

 foreach micro of varlist g1-g1409{
	count if `micro'==0 
	if r(N)>9862{	
	drop  `micro' 
	}
}
tab famine,missing
sort famine
by famine:su observed_features
regress observed_features i.famine age i.sex BMI // p=0.062
regress pielou_evenness i.famine age i.sex BMI // p=0.062

*drop if famine==9
regress observed_features i.famine age i.sex BMI // p=0.036
spearman age observed_features
save "D:\真菌分析\data_pre\2015_all_bac_dataset.dta",replace

*----datset update for the describe analyses-(此前描述基于最大样本量，这里更新为主数据)---

set excelxlsxlargefile on
clear mata 
set maxvar 15000

use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
keep SampleID
save "D:\真菌分析\data_pre\2015_all_dataset_ID.dta",replace

import excel "D:\真菌分析\rawdata\rawdata.xlsx", sheet("metadata") firstrow clear 
rename OTUID SampleID
merge m:m SampleID using  "D:\真菌分析\data_pre\2015_all_dataset_fill.dta"
keep if _merge==3
rename SampleID OTUID
tab region
xtile q_city = city_score, nq(4)
drop pielou_evenness-g751
export excel using "D:\真菌分析\rawdata\rawdata_.xlsx",  firstrow(variables) sheet("metadata_") sheetreplace  

import excel "D:\真菌分析\rawdata\rawdata_.xlsx", sheet("Sheet1") firstrow clear 
merge m:m SampleID using "D:\真菌分析\data_pre\2015_all_dataset_ID.dta"
keep if _merge==3
drop _merge
import excel "D:\真菌分析\rawdata\rawdata.xlsx", sheet("taxonomy") firstrow clear 

*-----3.2----PCOA analysis---缺失值处理需考虑

*--fungi--

use "D:\真菌分析\data_pre\2015_all_dataset.dta",clear
tab hhid  //这里可以看同一家庭的真菌构成与非同一家庭真菌构成的差异
tab commid //社区ID
*drop if age==. | BMI==. | fruit==. | tg==. | city_score==.
sort region

foreach phenotype of varlist t2d  marrige  smoke alcohol xinjigengse zhongfeng cancer sex hyper_med diabetes_med antibiotic_current antibiotic_6month probiotics kangyan_med kangsuan_med weisuan_med yogurt_drink pet hyper_med city  hyp dys pret2d fuxie gut_disease changdao_shoushu{
	 replace `phenotype'=. if `phenotype'==9 
}
misstable sum age MET income city_score hpc-tc BMI SBP DBP t2d  marrige  smoke alcohol xinjigengse zhongfeng cancer sex hyper_med diabetes_med antibiotic_current antibiotic_6month probiotics kangyan_med kangsuan_med weisuan_med yogurt_drink pet hyper_med city  hyp dys pret2d fuxie gut_disease changdao_shoushu rice-pastes collect_month 


qui foreach var of varlist age MET income city_score hpc-tc BMI SBP DBP rice-pastes collect_month {
     by region: egen `var'_mean= mean(`var')
	 by region: replace `var'=`var'_mean if `var'==.
}

foreach phenotype of varlist t2d  marrige  smoke alcohol xinjigengse zhongfeng cancer sex hyper_med diabetes_med antibiotic_current antibiotic_6month probiotics kangyan_med kangsuan_med weisuan_med yogurt_drink pet hyper_med city  hyp dys pret2d fuxie gut_disease changdao_shoushu{
	 replace `phenotype'=0 if `phenotype'==. 
}

keep SampleID idind pielou_evenness-shannon_entropy g1-g751 region-age t2d city education-alcohol xinjigengse zhongfeng cancer hyper_med MET city_score income  hpc-sex BMI-pret2d famine-antibiotic_current antibiotic_6month probiotics-pet  hhid commid collect_month
replace famine=9 if famine==.
save "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",replace

use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
tab famine
regress pielou_evenness i.famine age i.sex BMI // 

egen commid_group=group(commid)
tab commid_group
egen hhid_group=group(hhid)

preserve
keep  g1-g751
save "D:\真菌分析\data_pre\2015_all_fungi_pcoa.dta",replace 
restore

preserve
replace rice=0.001 if rice==0
gen ratio=wheat/rice
su ratio
keep  region-ratio
drop idind occupation commid hhid
save "D:\真菌分析\data_pre\2015_fungi_phenotype_pcoa.dta",replace 
restore

*----bacteria---

use "D:\真菌分析\data_pre\2015_all_bac_dataset.dta",clear
tab hhid  //这里可以看同一家庭的真菌构成与非同一家庭真菌构成的差异
tab commid //社区ID
*drop if age==. | BMI==. | fruit==. | tg==. | city_score==.
sort region
qui foreach var of varlist age MET income city_score hpc-tc BMI SBP DBP rice-pastes collect_month {
     by region: egen `var'_mean= mean(`var')
	 by region: replace `var'=`var'_mean if `var'==.
}

foreach phenotype of varlist t2d  marrige  smoke alcohol xinjigengse zhongfeng cancer sex hyper_med diabetes_med antibiotic_current antibiotic_6month probiotics kangyan_med kangsuan_med weisuan_med yogurt_drink pet hyper_med city  hyp dys pret2d fuxie gut_disease changdao_shoushu{
	 replace `phenotype'=0 if `phenotype'==9 |`phenotype'==. 
}

keep SampleID hhid commid idind pielou_evenness-shannon_entropy g51-g1400 region-age t2d city education-alcohol xinjigengse zhongfeng cancer hyper_med MET city_score income  hpc-sex BMI-pret2d famine-antibiotic_current antibiotic_6month probiotics-pet  
replace famine=9 if famine==.
save "D:\真菌分析\data_pre\2015_all_bac_dataset_fill.dta",replace

egen commid_group=group(commid)
tab commid_group
egen hhid_group=group(hhid)

preserve
keep idind g51-g1400
save "D:\真菌分析\data_pre\2015_all_bac_pcoa.dta",replace 
restore

preserve
keep  region-hhid_group
drop idind occupation commid hhid
save "D:\真菌分析\data_pre\2015_bac_phenotype_pcoa.dta",replace 
restore


*-----3.3----phenotype-Alpha diversity analysis----

*-特别关注独居、宠物饲养等对真菌构成的影响,predict analysis

use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
tab region
sort region
by region:spearman observed_features fish
preserve
drop if observed_features>500
keep pielou_evenness-shannon_entropy age
drop if age==. | shannon_entropy==.
save "D:\真菌分析\data_pre\fungi_agediversity.dta",replace
restore

tab famine
sort famine
by famine:su observed_features

*------fugi-连续变量--

tempname coef
 tempfile res
 postfile `coef' str200(diversity genus) float( n corr  p)  using "`res'", replace
qui foreach diversity of varlist pielou_evenness-shannon_entropy {
foreach phenotype of varlist age t2d city education marrige smoke-dys rice-pet {
 spearman `diversity'  `phenotype' 
 post `coef'    ("`diversity'")  ("`phenotype'") (r(N)) (r(rho))  (r(p))
}
}
 postclose `coef'
 preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\All data\表型_真菌多样性.xlsx",  firstrow(variables) sheet("连续变量") sheetreplace  
restore


tempname coef
 tempfile res
  postfile `coef' str200(outcome phenotype) float( n1 n2  p)  using "`res'", replace
  
qui foreach outcome of varlist pielou_evenness-shannon_entropy {
foreach phenotype of varlist t2d  marrige  smoke alcohol xinjigengse zhongfeng cancer sex   antibiotic_current antibiotic_6month probiotics kangyan_med kangsuan_med weisuan_med yogurt_drink pet hyper_med city  hyp dys {
	 replace `phenotype'=0 if `phenotype'==9 |`phenotype'==. 
     ranksum `outcome',by(`phenotype')
     post `coef'  ("`outcome'") ("`phenotype'")  (r(N_1)) (r(N_2)) (2*normprob(-abs(r(z))))
}
}
 postclose `coef'
 preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\All data\表型_真菌多样性.xlsx",  firstrow(variables) sheet("类别变量") sheetreplace  
restore



*--bacteria--

use "D:\真菌分析\data_pre\2015_all_bac_dataset_fill.dta",clear

tempname coef
 tempfile res
 postfile `coef' str200(diversity genus) float( n corr  p)  using "`res'", replace
qui foreach diversity of varlist pielou_evenness-shannon_entropy {
foreach phenotype of varlist age t2d city education marrige smoke-dys rice-pet {
 spearman `diversity'  `phenotype' 
 post `coef'    ("`diversity'")  ("`phenotype'") (r(N)) (r(rho))  (r(p))
}
}
 postclose `coef'
 preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\All data\表型_细菌多样性.xlsx",  firstrow(variables) sheet("连续变量") sheetreplace  
restore


*----3.4----Association of fungi diversity and bacterial diversity---

use "D:\真菌分析\data_pre\bac_15all_diversity.dta",clear
renvars pielou_evenness-shannon_entropy,prefix(bac_) 
merge 1:1 SampleID using "D:\真菌分析\data_pre\fugi_15all_diversity.dta"

keep if _merge==3
spearman bac_pielou_evenness pielou_evenness // 0.047
spearman bac_observed_features observed_features // 0.23
spearman bac_faith_pd faith_pd // 0.2349
spearman bac_shannon_entropy shannon_entropy // 0.0998

*----3.5----PCOA analysis----------

use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
tab region
keep SampleID region g1-g751
save "D:\真菌分析\data_pre\2015_all_fungi_forpca.dta",replace

import delimited "D:\真菌分析\Final results\All data\fungi_PCA.csv",  clear 
rename sampleid SampleID
save "D:\真菌分析\data_pre\2015_all_fungi_pc.dta",replace

*-----3.6----PCA  analysis----

use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
merge 1:1 SampleID using "D:\真菌分析\data_pre\2015_all_fungi_pc.dta"
keep if _merge==3
spearman dim1 age // -0.08
spearman dim2 age // 0.07
spearman dim3 age // 0.049
spearman dim4 age // 0.03
spearman dim5 age // -0.037

logit t2d dim2 age sex BMI,or // 1.05 (p=0.018)
logit t2d dim1 dim2 dim3 dim4 dim5 age sex BMI,or // 3负相关，1,2,5正相关
logit dys dim1 dim2 dim3 dim4 dim5 age sex BMI,or // 1,3正相关，4负相关
logit hyp dim1 dim2 dim3 dim4 dim5 age sex BMI,or // 4正相关，5 负相关


rename dim1 outcome
keep SampleID age t2d city-marrige smoke-pet outcome
save "D:\真菌分析\data_pre\2015_all_dataset_pc1_predict.dta",replace

*------3.7------同一社区或者家庭的人、采样时间肠道真菌及细菌分析 （PCOA等等）---

use "D:\真菌分析\data_pre\2015_all_dataset.dta",clear
egen hhid_group=group(hhid)
tab hhid_group
egen commid_group=group(commid)
tab commid_group
keep SampleID age sex BMI region hhid_group commid_group t2d dys hyp
save "D:\真菌分析\data_pre\fungi_geography.dta",replace // 患病了会不会离的更远


use "D:\真菌分析\Final results\All data\bray_dis.dta",  clear 
rename Var2 S1
rename value SampleID 
rename Col4 dis
merge m:m SampleID using "D:\真菌分析\data_pre\fungi_geography.dta"
drop _merge
renvars SampleID region-commid_group,prefix(person1_) 
rename S1 SampleID 
merge m:m SampleID using "D:\真菌分析\data_pre\fungi_geography.dta"
drop if dis==. 

gen G1=.
replace G1=1 if hhid_group==person1_hhid_group
replace G1=2 if (hhid_group!=person1_hhid_group) & (person1_commid_group==commid_group)
replace G1=3 if (hhid_group!=person1_hhid_group) & (person1_commid_group!=commid_group)& (region==person1_region)
replace G1=4 if (hhid_group!=person1_hhid_group) & (person1_commid_group!=commid_group)& (region!=person1_region)
tab G1,missing
sort G1
by G1:su dis

gen disease=0
replace disease=1 if t2d==1 | dys==1 | hyp==1
gen person_disease=0
replace person_disease=1 if person1_t2d==1 | person1_dys==1 | person1_hyp==1

gen G2=.
replace G2=1 if G1==1 & disease==0 & person_disease==0
replace G2=2 if G1==1 & disease==1 & person_disease==1
replace G2=3 if G1==1 & disease!=person_disease
replace G2=4 if G1==2 & disease==0 & person_disease==0
replace G2=5 if G1==2 & disease==1 & person_disease==1
replace G2=6 if G1==2 & disease!=person_disease
replace G2=7 if G1==3 & disease==0 & person_disease==0
replace G2=8 if G1==3 & disease==1 & person_disease==1
replace G2=9 if G1==3 & disease!=person_disease
replace G2=10 if G1==4 & disease==0 & person_disease==0
replace G2=11 if G1==4 & disease==1 & person_disease==1
replace G2=12 if G1==4 & disease!=person_disease
sort G2
by G2:su dis

gen age_group=.
replace age_group=1 if age<=20
replace age_group=2 if age>20 & age<40
replace age_group=3 if age>=40 & age<60
replace age_group=4 if age>=60 & age<80
replace age_group=5 if age>=80

gen age_group1=.
replace age_group1=1 if person1_age<=20
replace age_group1=2 if person1_age>20 & person1_age<40
replace age_group1=3 if person1_age>=40 & person1_age<60
replace age_group1=4 if person1_age>=60 & person1_age<80
replace age_group1=5 if person1_age>=80

save "D:\真菌分析\data_pre\fungi_geography_dis.dta",replace

preserve 
keep if G1==1
sort age_group age_group1
egen G3=group(age_group age_group1)
tab age_group G3
replace G3=2 if G3==6 
replace G3=3 if G3==11  
replace G3=4 if G3==16 
replace G3=5 if G3==21 
replace G3=8 if G3==12 
replace G3=9 if G3==17 
replace G3=10 if G3==22  
replace G3=14 if G3==18 
replace G3=15 if G3==23 
replace G3=20 if G3==24 
tab G3
sort G3
by G3:su dis
restore

*------3.8 Associaitons analysis---------注意结果在N OF 1 里重复（运动的数据可以整合分析）

import delimited "D:\真菌分析\Final results\All data\fungi_clr.csv",  clear 
rename sampleid SampleID 
save "D:\真菌分析\data_pre\fungi_clr.dta",replae


*-----fungi and diversity---

*--clr-based--

use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
drop g1-g751
merge 1:1 SampleID using "D:\真菌分析\data_pre\fungi_clr.dta"
logit t2d g280,or
sort region
by region:logit t2d g280,or
by region:logit t2d observed_features,or

sort District
by District:logit t2d g280,or
by District:logit t2d observed_features,or

tempname coef
 tempfile res
 postfile `coef' str200(diversity genus) float( n corr  p)  using "`res'", replace
qui foreach diversity of varlist pielou_evenness-shannon_entropy {
foreach phenotype of varlist g1-g751 {
 spearman `diversity'  `phenotype' 
 post `coef'    ("`diversity'")  ("`phenotype'") (r(N)) (r(rho))  (r(p))
}
}
 postclose `coef'
 preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\All data\真菌_真菌多样性.xlsx",  firstrow(variables) sheet("clr") sheetreplace  
restore

*--raw data----注意基于该结果深入分析，plot

use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear

tempname coef
 tempfile res
 postfile `coef' str200(diversity genus) float( n corr  p)  using "`res'", replace
qui foreach diversity of varlist pielou_evenness-shannon_entropy {
foreach phenotype of varlist g1-g751 {
 spearman `diversity'  `phenotype' 
 post `coef'    ("`diversity'")  ("`phenotype'") (r(N)) (r(rho))  (r(p))
}
}
 postclose `coef'
 preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\All data\真菌_真菌多样性.xlsx",  firstrow(variables) sheet("raw") sheetreplace  
restore

*-----fungi and phenotype----------

use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
hist city_score
gen city_score_=city_score/10
logit t2d city_score_ age i.sex BMI,or
logit dys city_score_ age i.sex BMI,or
logit hyp city_score_ age i.sex BMI,or
logit hyp i.city age i.sex BMI ,or
logit t2d i.city age i.sex BMI ,or

preserve
keep SampleID g1-g751 city_score
save "D:\真菌分析\data_pre\fungi_cityscore.dta",replace
restore

drop  g1-g751
merge 1:1 SampleID using "D:\真菌分析\data_pre\fungi_clr.dta"
tab commid
sort commid
by commid: su city_score
order region commid city_score
sort region commid

preserve
keep if city==1
keep SampleID city city_score g280 pielou_evenness faith_pd observed_features shannon_entropy
save "D:\真菌分析\data_pre\fungi_cityscore1.dta",replace
restore

preserve
keep if city==2
keep SampleID city city_score g280 pielou_evenness faith_pd observed_features shannon_entropy
save "D:\真菌分析\data_pre\fungi_cityscore2.dta",replace
restore


*--3.8.1---------Urbanisation analysis-----

*--distributin--

sort city
by city:su city_score g280 // 4281 vs 6414 (80.83 vs 65.84)

preserve
keep if city==1
hist city_score 
restore
preserve
keep if city==2
hist city_score 
restore

*--associaiton with diet or lifestyle

gen vegetable=lveg+dveg
tempname coef
 tempfile res
 postfile `coef' str200(city_score die_lifestyle) float( n corr  p)  using "`res'", replace
foreach phenotype of varlist smoke alcohol MET rice-tuber dveg-othmeat poultry-sauce pastes vegetable {
 spearman city_score  `phenotype' 
 post `coef'    ("city_score")  ("`phenotype'") (r(N)) (r(rho))  (r(p))
}

 postclose `coef'
 preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\All data\associaiton analysis.xlsx",  firstrow(variables) sheet("cityscore_diet") sheetreplace  
restore

use "D:\真菌分析\data_pre\fungi_cityscore1.dta",clear
append using "D:\真菌分析\data_pre\fungi_cityscore1.dta"
save "D:\真菌分析\data_pre\cityscore_distribution.dta",replace

*---associaiton with healthy index of diet / diet diversity---

**----------prediction analysis-----

use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
keep SampleID-g751 city_score
save "D:\真菌分析\data_pre\2015_all_dataset_cityscore_predict.dta",replace

*---3.8.2--mixed model---------

*---diet、lifestyle、urbanisation to fungi(continue vars)--

use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
preserve
drop  if observed_features>500
keep t2d dys hyp pielou_evenness-shannon_entropy 
save "D:\真菌分析\data_pre\diversity_disease.dta",replace 
restore	 
drop g1-g751
merge 1:1 SampleID using "D:\真菌分析\data_pre\fungi_clr.dta"
keep if _merge==3
gen vegetable=lveg+dveg
gen fungi_index=g247+g186+g286

mixed SBP fungi_index  age i.sex BMI ||region:, covariance(unstructured) // positive
mixed tg fungi_index  age i.sex BMI ||region:, covariance(unstructured) // positive
melogit hyp fungi_index  age i.sex  ||region:, covariance(unstructured) // 0.004
melogit dys fungi_index  age i.sex  ||region:, covariance(unstructured) // 0.007
drop if observed_features>500

save "D:\真菌分析\data_pre\fungi_disease_analysis.dta",replace 
		 
preserve
keep t2d dys hyp g332	g246	g280	g286	g149	g249	g2	g702 t2d hyp dys
save "D:\真菌分析\data_pre\fungi_disease.dta",replace 
restore	 

*---diet、lifestyle、urbanisation to fungi(continue vars)--		 

*--hub index--效果不好

tempname coef
tempfile res
postfile `coef' str200( micro var) float(n rr lul uul p)   using "`res'", replace


foreach var of varlist   MET city_score  rice-tuber dveg-othmeat poultry-sauce pastes vegetable{ 
    mixed fungi_index `var' age i.sex BMI ||region:, covariance(unstructured) 
test `var'
 post `coef'   ("fungi_index") ("`var'")  (e(N)) (_b[`var'])  (_b[`var']-1.96*_se[`var'])  (_b[`var']+1.96*_se[`var'])  (chi2tail(1,(_b[`var']/_se[`var'])^2))  
 }
 

  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\All data\associaiton analysis.xlsx",  firstrow(variables) sheet("diet_fungiindex_continues") sheetreplace  
restore

*--fungi features--
		 
tempname coef
tempfile res
postfile `coef' str200( micro var) float(n rr lul uul p)   using "`res'", replace


foreach outcome of varlist pielou_evenness-shannon_entropy g1-g751  {
foreach var of varlist   MET city_score  rice-tuber dveg-othmeat poultry-sauce pastes vegetable{ 
    mixed `outcome' `var' age i.sex BMI ||region:, covariance(unstructured) 
test `var'
 post `coef'   ("`outcome'") ("`var'")  (e(N)) (_b[`var'])  (_b[`var']-1.96*_se[`var'])  (_b[`var']+1.96*_se[`var'])  (chi2tail(1,(_b[`var']/_se[`var'])^2))  
 }
 }

  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\All data\associaiton analysis.xlsx",  firstrow(variables) sheet("diet_micro_continues") sheetreplace  
restore

*---diet、lifestyle、urbanisation to fungi(catgary vars)-

*--fungi index--不太显著

tempname coef
tempfile res
postfile `coef' str200( micro var) float(n rr lul uul p)   using "`res'", replace

gen rural=0
replace rural=1 if city==1 // 1-city

foreach var of varlist  rural smoke alcohol { 
    mixed fungi_index i.`var' age i.sex BMI ||region:, covariance(unstructured) 
test 1.`var'
 post `coef'   ("fungi_index") ("`var'")  (e(N)) (_b[1.`var'])  (_b[1.`var']-1.96*_se[1.`var'])  (_b[1.`var']+1.96*_se[1.`var'])  (chi2tail(1,(_b[1.`var']/_se[1.`var'])^2))  
 }
 

  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\All data\associaiton analysis.xlsx",  firstrow(variables) sheet("diet_fungiindex_categray") sheetreplace  
restore

*--fungi features--

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


*-----fungi to phenotype------

*--features----

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

*--index----

tempname coef
tempfile res
postfile `coef' str200( var micro ) float(n rr lul uul p) str200(cov)  using "`res'", replace

foreach outcome of varlist  HDL_C LDL_C HbA1c glucose tg tc SBP DBP {
foreach var of varlist fungi_index { 
    mixed `outcome' `var' age i.sex BMI ||region:, covariance(unstructured) 
test `var'
 post `coef'   ("`outcome'") ("`var'")  (e(N)) (_b[`var'])  (_b[`var']-1.96*_se[`var'])  (_b[`var']+1.96*_se[`var'])  (chi2tail(1,(_b[`var']/_se[`var'])^2))  ("${cov`m'}")
 }
 }

  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\All data\associaiton analysis.xlsx",  firstrow(variables) sheet("fungi_index_phenotype") sheetreplace  
restore


*----fungi to diseases-------注意探究交互作用

*---model 1: adjust age+sex

tempname coef
tempfile res
postfile `coef' str200( micro var) float(n rr lul uul p) str200(cov)  using "`res'", replace


foreach outcome of varlist   t2d hyp dys {
foreach var of varlist pielou_evenness-shannon_entropy g1-g751 { 
    melogit `outcome' `var' age i.sex  ||region:, covariance(unstructured) 
test `var'
 post `coef'   ("`outcome'") ("`var'")  (e(N)) (exp(_b[`var']))  (exp(_b[`var']-1.96*_se[`var']))  (exp(_b[`var']+1.96*_se[`var']))  (chi2tail(1,(_b[`var']/_se[`var'])^2))  ("${cov`m'}")
 }
 }

  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\All data\associaiton analysis.xlsx",  firstrow(variables) sheet("fungi_diseases_model1") sheetreplace  
restore

*---model 1: adjust age+sex+BMI

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

*---model 3: adjust age+sex+BMI+cityscore+diet and lifestyle

tempname coef
tempfile res
postfile `coef' str200( micro var) float(n rr lul uul p) str200(cov)  using "`res'", replace


foreach outcome of varlist   t2d hyp dys {
foreach var of varlist pielou_evenness-shannon_entropy g1-g751 { 
    melogit `outcome' `var' age i.sex BMI i.smoke i.alcohol MET city_score fish milk fruit vegetable ||region:, covariance(unstructured) 
test `var'
 post `coef'   ("`outcome'") ("`var'")  (e(N)) (exp(_b[`var']))  (exp(_b[`var']-1.96*_se[`var']))  (exp(_b[`var']+1.96*_se[`var']))  (chi2tail(1,(_b[`var']/_se[`var'])^2))  ("${cov`m'}")
 }
 }

  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\All data\associaiton analysis.xlsx",  firstrow(variables) sheet("fungi_diseases_model3") sheetreplace  
restore

*---fungi_index: adjust age+sex+BMI

tempname coef
tempfile res
postfile `coef' str200( micro var) float(n rr lul uul p) str200(cov)  using "`res'", replace


foreach outcome of varlist   t2d hyp dys {
foreach var of varlist fungi_index { 
    melogit `outcome' `var' age i.sex  ||region:, covariance(unstructured) 
test `var'
 post `coef'   ("`outcome'") ("`var'")  (e(N)) (exp(_b[`var']))  (exp(_b[`var']-1.96*_se[`var']))  (exp(_b[`var']+1.96*_se[`var']))  (chi2tail(1,(_b[`var']/_se[`var'])^2))  ("${cov`m'}")
 }
 }

  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\All data\associaiton analysis.xlsx",  firstrow(variables) sheet("fungi_index_diseases_model2") sheetreplace  
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

melogit t2d std_risk   ||region:, covariance(unstructured) or // 2.38
melogit t2d std_g280   ||region:, covariance(unstructured) or // 1.10

melogit t2d c.g280##c.age_score   ||region:, covariance(unstructured) or // 交互项边缘显著, p=0.064
melogit t2d c.g280##c.BMI_score   ||region:, covariance(unstructured) or // 交互项不显著
melogit t2d c.g280##c.wc_score   ||region:, covariance(unstructured) or // 交互项不显著
melogit t2d c.g280##c.sys_score   ||region:, covariance(unstructured) or // 交互项不显著

melogit t2d c.g280##c.risk   ||region:, covariance(unstructured) or // 交互项显著, p=0.012
melogit t2d c.risk##c.observed_features   ||region:, covariance(unstructured) or // 交互项不显著
melogit t2d c.std_risk##c.std_g702   ||region:, covariance(unstructured) or // 交互项不显著

xtile q_risk = risk, nq(2)
xtile q_g280 = g280, nq(2)
tab q_risk q_g280
egen group=group(q_risk q_g280)
tab group t2d

gen carry=0 // clr转换后，这样不对
replace carry=1 if g280>0

melogit t2d i.q_risk   ||region:, covariance(unstructured) or // 3.47
melogit t2d i.q_g280   ||region:, covariance(unstructured) or // 1.158
melogit t2d i.q_risk#i.q_g280   ||region:, covariance(unstructured) or // 效果不错

preserve
keep if q_risk==1
melogit t2d std_g280   ||region:, covariance(unstructured) or // 1.22,0.001
melogit t2d i.carry   ||region:, covariance(unstructured) or // 1.42,0.04
melogit t2d std_observed_features   ||region:, covariance(unstructured) or // 不显著
restore

preserve
keep if q_risk==2
melogit t2d std_g280   ||region:, covariance(unstructured) or // 显著-0.016 （1.09）
melogit t2d i.carry   ||region:, covariance(unstructured) or // 不显著
melogit t2d std_observed_features   ||region:, covariance(unstructured) or // 边缘显著
restore

gen risk_group=0
replace risk_group=1 if risk>=25

melogit t2d i.risk_group   ||region:, covariance(unstructured) or // 3.47

preserve
keep if risk_group==0
melogit t2d std_g280   ||region:, covariance(unstructured) or // 1.22
melogit t2d std_observed_features   ||region:, covariance(unstructured) or // 边缘显著
restore

preserve
keep if risk_group==1
melogit t2d std_g280   ||region:, covariance(unstructured) or // 1.09
melogit t2d std_observed_features   ||region:, covariance(unstructured) or // 显著
restore

preserve
keep if q_g280==1
melogit t2d std_risk   ||region:, covariance(unstructured) or // 2.59
restore

preserve
keep if q_g280==2
melogit t2d std_risk   ||region:, covariance(unstructured) or // 2.25
restore


*----name mapping-----生成2*4=8套标签分别做FDR

import excel "D:\真菌分析\Plot\data\permonova.xlsx", sheet("fungi") firstrow clear
rename Factors var 
keep Var var Categrary
save  "D:\真菌分析\data_pre\phenotype_mapping.dta",replace
 
import excel "D:\真菌分析\rawdata\g_level_mapping.xlsx", sheet("g_level_mapping") firstrow clear
rename A micro
save  "D:\真菌分析\data_pre\all_g_mapping.dta",replace

*----continue diet vars--

import excel "D:\真菌分析\Final results\All data\associaiton analysis.xlsx", sheet("diet_micro_continues") firstrow clear
merge m:m micro using "D:\真菌分析\data_pre\all_g_mapping.dta"
order Phylum Class Order Family Genus micro var p
sort p
drop if p==.
tab var,missing
drop Kingdom _merge
merge m:m var using "D:\真菌分析\data_pre\phenotype_mapping.dta"
drop if p==.
replace Var="Vegetable" if var=="vegetable"
gen micro_var="genus"
replace micro_var="diversity" if Genus==""
replace Genus=micro if Genus==""
keep Phylum Class Order Family Genus Var n rr lul uul p micro_var
order Phylum Class Order Family Genus Var n rr lul uul p micro_var
sort p
save  "D:\真菌分析\data_pre\diet_fungi_continue.dta",replace

*----categary diet vars--

import excel "D:\真菌分析\Final results\All data\associaiton analysis.xlsx", sheet("diet_micro_categray") firstrow clear
merge m:m micro using "D:\真菌分析\data_pre\all_g_mapping.dta"
order Phylum Class Order Family Genus micro var p
sort p
drop if p==.
tab var,missing
drop Kingdom _merge
replace var="city" if var=="rural"
merge m:m var using "D:\真菌分析\data_pre\phenotype_mapping.dta"
drop if p==.
gen micro_var="genus"
replace micro_var="diversity" if Genus==""
replace Genus=micro if Genus==""
keep Phylum Class Order Family Genus Var n rr lul uul p micro_var
order Phylum Class Order Family Genus Var n rr lul uul p micro_var
sort p
save  "D:\真菌分析\data_pre\diet_fungi_cat.dta",replace
 
*---diet vars append-----

use  "D:\真菌分析\data_pre\diet_fungi_continue.dta",clear
append using   "D:\真菌分析\data_pre\diet_fungi_cat.dta"
export excel using "D:\真菌分析\Plot\data\mixed_associaitons.xlsx",  firstrow(variables) sheet("diet_city_fungi") sheetreplace  

*-----fungi-phenotype-

import excel "D:\真菌分析\Final results\All data\associaiton analysis.xlsx", sheet("fungi_phenotype") firstrow clear
merge m:m micro using "D:\真菌分析\data_pre\all_g_mapping.dta"
order Phylum Class Order Family Genus micro var p
sort p
drop if p==.
tab var,missing
drop Kingdom _merge
merge m:m var using "D:\真菌分析\data_pre\phenotype_mapping.dta"
drop if p==.
gen micro_var="genus"
replace micro_var="diversity" if Genus==""
replace Genus=micro if Genus==""
keep Phylum Class Order Family Genus Var n rr lul uul p micro_var
order Phylum Class Order Family Genus Var n rr lul uul p micro_var
sort p
export excel using "D:\真菌分析\Plot\data\mixed_associaitons.xlsx",  firstrow(variables) sheet("fungi_phenotype") sheetreplace  

*-----fungi-disease-----

import excel "D:\真菌分析\Final results\All data\associaiton analysis.xlsx", sheet("fungi_diseases_model1") firstrow clear
merge m:m micro using "D:\真菌分析\data_pre\all_g_mapping.dta"
order Phylum Class Order Family Genus micro var p
sort p
drop if p==.
tab var,missing
drop Kingdom _merge
merge m:m var using "D:\真菌分析\data_pre\phenotype_mapping.dta"
drop if p==.
gen micro_var="genus"
replace micro_var="diversity" if Genus==""
replace Genus=micro if Genus==""
keep Phylum Class Order Family Genus Var n rr lul uul p micro_var
order Phylum Class Order Family Genus Var n rr lul uul p micro_var
sort p
export excel using "D:\真菌分析\Plot\data\mixed_associaitons.xlsx",  firstrow(variables) sheet("fungi_diseases_model1") sheetreplace  

*----------------PDI UPDI analysis-------------------

use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
drop g1-g751
merge 1:1 SampleID using "D:\真菌分析\data_pre\fungi_clr.dta"

gen wholegrain=cereoth
gen vegs=dveg+lveg
gen legumes = dryleg+legpro
gen refinegrain=rice+wheat+pastes
gen beverage=alcohol
gen sweeter=cake+sugar
gen meat=pork+othmeat+orgmeat+poultry

 foreach var of varlist wholegrain fruit vegs nuts legumes OIL_VEG  tuber refinegrain beverage sweeter OIL_ANI milk egg fish meat {
xtile `var'_score = `var', nq(5)
}
 

** plant
foreach var of varlist wholegrain fruit vegs nuts legumes OIL_VEG  tuber refinegrain beverage sweeter{
gen `var'_PDI = .
replace `var'_PDI = 1 if `var'_score == 1 
replace `var'_PDI = 2 if `var'_score == 2
replace `var'_PDI = 3 if `var'_score == 3
replace `var'_PDI = 4 if `var'_score == 4
replace `var'_PDI = 5 if `var'_score == 5
}

 **animal
foreach var of varlist OIL_ANI milk egg fish meat {
gen `var'_PDI = .
replace `var'_PDI = 5 if `var'_score == 1 
replace `var'_PDI = 4 if `var'_score == 2
replace `var'_PDI = 3 if `var'_score == 3
replace `var'_PDI = 2 if `var'_score == 4
replace `var'_PDI = 1 if `var'_score == 5
}

 ** hPDI healthy plant food  groups
foreach var of varlist wholegrain fruit vegs nuts legumes OIL_VEG  tuber{
gen `var'_hPDI = .
replace `var'_hPDI = 1 if `var'_score == 1 
replace `var'_hPDI = 2 if `var'_score == 2
replace `var'_hPDI = 3 if `var'_score == 3
replace `var'_hPDI = 4 if `var'_score == 4
replace `var'_hPDI = 5 if `var'_score == 5
}

 **hPDI unhealthy food  groups
foreach var of varlist refinegrain beverage sweeter OIL_ANI milk egg fish meat{
gen `var'_hPDI = .
replace `var'_hPDI = 1 if `var'_score == 5 
replace `var'_hPDI = 2 if `var'_score == 4
replace `var'_hPDI = 3 if `var'_score == 3
replace `var'_hPDI = 4 if `var'_score == 2
replace `var'_hPDI = 5 if `var'_score == 1
}


** uPDI healthy plant food  groups
foreach var of varlist wholegrain fruit vegs nuts legumes OIL_VEG  tuber{
gen `var'_uPDI = .
replace `var'_uPDI = 1 if `var'_score == 5 
replace `var'_uPDI = 2 if `var'_score == 4
replace `var'_uPDI = 3 if `var'_score == 3
replace `var'_uPDI = 4 if `var'_score == 2
replace `var'_uPDI = 5 if `var'_score == 1
}

 **uPDI unhealthy plant food  groups
foreach var of varlist refinegrain beverage sweeter{
gen `var'_uPDI = .
replace `var'_uPDI = 1 if `var'_score == 1 
replace `var'_uPDI = 2 if `var'_score == 2
replace `var'_uPDI = 3 if `var'_score == 3
replace `var'_uPDI = 4 if `var'_score == 4
replace `var'_uPDI = 5 if `var'_score == 5
}

**uPDI animal food groups
foreach var of varlist  OIL_ANI milk egg fish meat  {
gen `var'_uPDI = .
replace `var'_uPDI = 5 if `var'_score == 1 
replace `var'_uPDI = 4 if `var'_score == 2
replace `var'_uPDI = 3 if `var'_score == 3
replace `var'_uPDI = 2 if `var'_score == 4
replace `var'_uPDI = 1 if `var'_score == 5
}

egen PDI = rowtotal (**_PDI)
egen hPDI = rowtotal (**_hPDI)
egen uPDI = rowtotal (**_uPDI)
egen pdi = rowtotal (wholegrain_PDI-tuber_PDI)
egen unhealthy_pdi = rowtotal (refinegrain_PDI beverage_PDI sweeter_PDI)
egen animal_pdi = rowtotal (OIL_ANI_PDI-meat_PDI)

spearman city_score PDI // -0.125
spearman city_score hPDI // -0.14
spearman city_score uPDI // -0.23
spearman city_score pdi // -0.08
spearman city_score unhealthy_pdi // 0.0207
spearman city_score animal_pdi // -0.31
spearman city_score refinegrain // -0.268

spearman observed_features PDI // 0.023
spearman observed_features hPDI // 0.05
spearman observed_features uPDI // 不显著

spearman g280 PDI // 0.023
spearman g280 hPDI // 0.021



*------------------City score anslysis---------------

*--diversity analysis-(alpha and beta-to demonstrate beta more individual)

use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
drop g1-g751
merge 1:1 SampleID using "D:\真菌分析\data_pre\fungi_clr.dta"
gen legumes=dryleg+legpro
spearman city_score pielou_evenness // -0.048
spearman city_score observed_features // -0.065
spearman city_score faith_pd // -0.14
spearman city_score shannon_entropy // -0.0576
spearman cereoth observed_features // 0.023
spearman cereoth city_score // nosig
spearman legumes city_score // 0.066
spearman legumes observed_features // -0.023

*=======假设与城市化水平正相关的菌有害，与城市化负相关的菌有益，构建index链接疾病


*--------fungi-cityscore and fungi-phenotype association consistent--12/19

import excel "D:\真菌分析\Plot\Results\mixed_models.xlsx", sheet("combine_") firstrow clear
keep if Var=="Urbanisation index"
keep if p_adjust<=0.05
keep Genus rr
rename rr city_rr
save  "D:\真菌分析\data_pre\cityscore_fungi_rr.dta",replace

import excel "D:\真菌分析\Plot\Results\mixed_models.xlsx", sheet("combine_") firstrow clear
keep if dataset=="fungi_phenotype"
keep if p_adjust<=0.05
replace rr=-rr if Var=="HDL-C"
keep Genus rr
rename rr phenotype_rr
save  "D:\真菌分析\data_pre\phenotype_fungi_rr.dta",replace

use  "D:\真菌分析\data_pre\cityscore_fungi_rr.dta",clear
merge m:m Genus using  "D:\真菌分析\data_pre\phenotype_fungi_rr.dta"
keep if _merge==3
spearman city_rr phenotype_rr

*--------CLR based----以此为准

use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
drop g1-g751
merge 1:1 SampleID using "D:\真菌分析\data_pre\fungi_clr.dta"

*-------Relative abundance based---

use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
egen g_all=rowtotal(g1-g751)
 foreach micro of varlist g1-g751{
   replace `micro'=`micro'/g_all
   } 
sort t2d
by t2d: su g280 g702

egen city_neg=rowtotal(g78	g332	g247	g149	g411	g717	g54	g349	g305	g347	g746	g21	g79	g407	g138	g137	g139	g348	g252)
egen city_pos=rowtotal(g266	g287	g354	g299	g751	g20	g291	g246	g701	g245	g658	g618	g127	g259	g564	g702	g700	g740	g643	g277	g605	g2	g1	g186	g280)
gen index=city_pos-city_neg

qui foreach micro of varlist city_neg city_pos index {
	egen std_`micro'=std(`micro')
	drop `micro'
	rename std_`micro' `micro'
}

*----diseases---

tempname coef
tempfile res
postfile `coef' str200( micro var) float(n rr lul uul p) str200(cov)  using "`res'", replace

foreach outcome of varlist   t2d hyp dys {
foreach var of varlist  city_neg city_pos index{ 
    melogit `outcome' `var' age i.sex  BMI ||region:, covariance(unstructured) 
test `var'
 post `coef'   ("`outcome'") ("`var'")  (e(N)) (exp(_b[`var']))  (exp(_b[`var']-1.96*_se[`var']))  (exp(_b[`var']+1.96*_se[`var']))  (chi2tail(1,(_b[`var']/_se[`var'])^2))  ("${cov`m'}")
 }
 }

  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\All data\city_score_fungi.xlsx",  firstrow(variables) sheet("cityscore_index_diseases_rela") sheetreplace  
restore

keep if city==1 // 4281
tempname coef
tempfile res
postfile `coef' str200( micro var) float(n rr lul uul p) str200(cov)  using "`res'", replace

foreach outcome of varlist   t2d hyp dys {
foreach var of varlist  city_neg city_pos { 
    melogit `outcome' `var' age i.sex  BMI ||region:, covariance(unstructured) 
test `var'
 post `coef'   ("`outcome'") ("`var'")  (e(N)) (exp(_b[`var']))  (exp(_b[`var']-1.96*_se[`var']))  (exp(_b[`var']+1.96*_se[`var']))  (chi2tail(1,(_b[`var']/_se[`var'])^2))  ("${cov`m'}")
 }
 }

  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\All data\city_score_fungi.xlsx",  firstrow(variables) sheet("city1_clr") sheetreplace  
restore

keep if city==2 // 6414
tempname coef
tempfile res
postfile `coef' str200( micro var) float(n rr lul uul p) str200(cov)  using "`res'", replace

foreach outcome of varlist   t2d hyp dys {
foreach var of varlist  city_neg city_pos { 
    melogit `outcome' `var' age i.sex  BMI ||region:, covariance(unstructured) 
test `var'
 post `coef'   ("`outcome'") ("`var'")  (e(N)) (exp(_b[`var']))  (exp(_b[`var']-1.96*_se[`var']))  (exp(_b[`var']+1.96*_se[`var']))  (chi2tail(1,(_b[`var']/_se[`var'])^2))  ("${cov`m'}")
 }
 }

  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\All data\city_score_fungi.xlsx",  firstrow(variables) sheet("city2_clr") sheetreplace  
restore

*----phenotypes---

tempname coef
tempfile res
postfile `coef' str200( var micro ) float(n rr lul uul p) str200(cov)  using "`res'", replace

foreach outcome of varlist  HDL_C LDL_C HbA1c glucose tg tc SBP DBP {
foreach var of varlist  city_neg city_pos index{ 
    mixed `outcome' `var' age i.sex BMI ||region:, covariance(unstructured) 
test `var'
 post `coef'   ("`outcome'") ("`var'")  (e(N)) (_b[`var'])  (_b[`var']-1.96*_se[`var'])  (_b[`var']+1.96*_se[`var'])  (chi2tail(1,(_b[`var']/_se[`var'])^2))  ("${cov`m'}")
 }
 }

  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\All data\city_score_fungi.xlsx",  firstrow(variables) sheet("cityscore_index_phynotypes_rela") sheetreplace  
restore

*---diet and lifestyle---

gen vegetable=lveg+dveg
foreach v of varlist age MET rice-pastes  vegetable  {
           egen s_`v' = std(`v')
		   drop `v'
		   rename s_`v' `v'
         }

*---diet、lifestyl to urban index(continue vars)--不好解释 
		 
tempname coef
tempfile res
postfile `coef' str200( micro var) float(n rr lul uul p)   using "`res'", replace


foreach outcome of varlist city_neg city_pos index  {
foreach var of varlist   MET   rice-tuber dveg-othmeat poultry-sauce pastes vegetable{ 
    mixed `outcome' `var' age i.sex BMI ||region:, covariance(unstructured) 
test `var'
 post `coef'   ("`outcome'") ("`var'")  (e(N)) (_b[`var'])  (_b[`var']-1.96*_se[`var'])  (_b[`var']+1.96*_se[`var'])  (chi2tail(1,(_b[`var']/_se[`var'])^2))  
 }
 }

  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\All data\city_score_fungi.xlsx",  firstrow(variables) sheet("cityscore_index_diet_continue_rela") sheetreplace  
restore

*---data for plot----

use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
drop g1-g751
merge 1:1 SampleID using "D:\真菌分析\data_pre\fungi_clr.dta"
egen city_neg=rowtotal(g78	g332	g247	g149	g411	g717	g54	g349	g305	g347	g746	g21	g79	g407	g138	g137	g139	g348	g252)
egen city_pos=rowtotal(g266	g287	g354	g299	g751	g20	g291	g246	g701	g245	g658	g618	g127	g259	g564	g702	g700	g740	g643	g277	g605	g2	g1	g186	g280)
 
preserve 
keep  city_pos t2d 
rename t2d outcome 
rename city_pos value
gen disease="T2D"
save "D:\真菌分析\data_pre\city_index_t2d_pos.dta",replace
restore 

preserve 
keep  city_pos dys 
rename dys outcome 
rename city_pos value
gen disease="DYS"
save "D:\真菌分析\data_pre\city_index_dys_pos.dta",replace
restore 

preserve 
keep  city_pos hyp 
rename hyp outcome 
rename city_pos value
gen disease="HYP"
save "D:\真菌分析\data_pre\city_index_hyp_pos.dta",replace
restore 


preserve 
keep  city_neg t2d 
rename t2d outcome 
rename city_neg value
gen disease="T2D"
save "D:\真菌分析\data_pre\city_index_t2d_neg.dta",replace
restore 

preserve 
keep  city_neg dys 
rename dys outcome 
rename city_neg value
gen disease="DYS"
save "D:\真菌分析\data_pre\city_index_dys_neg.dta",replace
restore 

preserve 
keep  city_neg hyp 
rename hyp outcome 
rename city_neg value
gen disease="HYP"
save "D:\真菌分析\data_pre\city_index_hyp_neg.dta",replace
restore 

use "D:\真菌分析\data_pre\city_index_t2d_pos.dta",clear
append using "D:\真菌分析\data_pre\city_index_dys_pos.dta"
append using "D:\真菌分析\data_pre\city_index_hyp_pos.dta"
save "D:\真菌分析\data_pre\city_index_disease_pos.dta",replace

use "D:\真菌分析\data_pre\city_index_t2d_neg.dta",clear
append using "D:\真菌分析\data_pre\city_index_dys_neg.dta"
append using "D:\真菌分析\data_pre\city_index_hyp_neg.dta"
save "D:\真菌分析\data_pre\city_index_disease_neg.dta",replace


*=======菌群预测的城市化水平与实际的mismatch index与疾病结局的关联

import excel "D:\真菌分析\Final results\predict\fungi_cityscore_predict.xlsx", sheet("Sheet1") firstrow clear 
merge 1:1 SampleID using "D:\真菌分析\data_pre\2015_all_dataset_fill.dta"
spearman predict label // 0.35
spearman predict observed_features // -0.13
spearman cereoth predict // 0.025
logit t2d predict,or // p=0.01
logit t2d label,or // 更显著
gen rediual=label-predict 
logit t2d rediual,or // 显著有害
gen rediual1=predict-label
logit t2d rediual1,or // 显著有益

xtile q_predict=predict,nq(2)
xtile q_label=label,nq(2)
tab q_predict q_label

gen match=0
replace match=1 if q_predict==q_label
logit t2d i.match,or // nosig

*----------------Network analysis----------------

import excel "D:\真菌分析\data_pre\networkanalysis.xlsx", sheet("Sheet2") firstrow clear
rename SampleID micro
merge 1:1 micro using   "D:\真菌分析\data_pre\all_g_mapping.dta"
keep if _merge==3
split (micro), parse(g) 
destring micro2,replace
sort micro2

*-------------------DMM cluster and phenotype-------------有点不好解释

import delimited "D:\真菌分析\Final results\All data\DMM\All_cluster.csv",  clear //这个好
import delimited "D:\真菌分析\results1\DMM results\All_cluster.csv",  clear 
rename sampleid SampleID 
merge 1:1 SampleID using  "D:\真菌分析\data_pre\2015_all_dataset_fill.dta"
sort cluster
by cluster:su  pielou_evenness faith_pd observed_features shannon_entropy
logit t2d i.cluster ,or
logit dys i.cluster ,or
logit hyp i.cluster ,or


*--------------------Healthy lifestyle score and fungi------

*----------------All data based---------------------

*---score 1 (ref:PMID: 35103793)-缺少睡眠数据

use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
drop g1-g751
merge 1:1 SampleID using "D:\真菌分析\data_pre\fungi_clr.dta"
gen veg=dveg+lveg
 foreach var of varlist fruit veg  pork poultry milk sugar fish OIL_VEG MET nuts{
	xtile q_`var'=`var',nq(2)
	}
	
gen bmi_score=0
replace bmi_score=2 if BMI>=18.5 & BMI<=24.9
gen fruit_score=0
replace fruit_score=1 if q_fruit==2
gen veg_score=0
replace veg_score=1 if q_veg==2
gen meat_score=0
replace meat_score=1 if q_pork==1
gen total_score=veg_score+fruit_score+meat_score
gen diet_score=0
replace diet_score=2 if total_score==3
replace diet_score=1 if total_score==2
gen smoke_score=0
replace smoke_score=2 if smoke==0
gen MET_score=0
replace MET_score=2 if q_MET==2
gen alc_score=0
replace alc_score=2 if alcohol==0 // 酒精加上效果不好
gen lifestyle_score=bmi_score+diet_score+smoke_score+MET_score // 参考ref1构建
tab lifestyle_score
gen lifestyle_score1=bmi_score+diet_score+smoke_score+MET_score+alc_score // 改进

spearman lifestyle_score observed_features // 0.07
spearman lifestyle_score1 observed_features // 0.067
spearman lifestyle_score shannon_entropy // 边缘显著 p=0.06
spearman lifestyle_score pielou_evenness // 不显著
spearman lifestyle_score faith_pd // 0.059

spearman lifestyle_score city_score // -0.055
spearman lifestyle_score1 city_score // -0.056
spearman lifestyle_score g280 // 不显著
spearman lifestyle_score1 g280 // 不显著

*----继续改进-------该结果较理想

gen sugar_score=0
replace sugar_score=1 if sugar==0
gen total_score1=veg_score+fruit_score+meat_score+sugar_score
tab total_score1
gen diet_score1=0
replace diet_score1=2 if total_score1==4
replace diet_score1=1 if total_score1==3 | total_score1==2 
gen lifestyle_score2=bmi_score+diet_score1+smoke_score+MET_score // 改进
spearman lifestyle_score2 observed_features // 0.077
spearman lifestyle_score2 g280 // -0.03
spearman lifestyle_score2 city_score  //  -0.11

*------再改进--效果不好

gen fish_score=0
replace fish_score=1 if q_fish==2
gen pou_score=0
replace pou_score=1 if q_poultry==2
gen nut_score=0
replace nut_score=1 if q_nut==2
gen total_score2=veg_score+fruit_score+meat_score+sugar_score+fish_score //这里不能加milk（类别不明）
tab total_score2

gen diet_score2=0
replace diet_score2=2 if total_score2==5
replace diet_score2=1 if  total_score2>=1 & total_score2<5
gen lifestyle_score3=bmi_score+diet_score2+smoke_score+MET_score
spearman lifestyle_score3 observed_features // 0.057
spearman lifestyle_score3 g280 // -0.033 
spearman lifestyle_score3 city_score  //  -0.088


*--------------------repeated data based---------------------

*---2015--

use "D:\真菌分析\data_pre\2015_repeat_datasets.dta",clear
drop g2-g633
merge 1:1 SampleID using "D:\真菌分析\data_pre\2015_repeat_fungi_clr.dta"
keep if _merge==3
gen veg=dveg+lveg
 foreach var of varlist fruit veg  pork poultry milk sugar fish OIL_VEG MET nuts{
	xtile q_`var'=`var',nq(2)
	}
	
gen bmi_score=0
replace bmi_score=2 if BMI>=18.5 & BMI<=24.9
gen fruit_score=0
replace fruit_score=1 if q_fruit==2
gen veg_score=0
replace veg_score=1 if q_veg==2
gen meat_score=0
replace meat_score=1 if q_pork==1
gen total_score=veg_score+fruit_score+meat_score
gen diet_score=0
replace diet_score=2 if total_score==3
replace diet_score=1 if total_score==2
gen smoke_score=0
replace smoke_score=2 if smoke==0
gen MET_score=0
replace MET_score=2 if q_MET==2
gen sugar_score=0
replace sugar_score=1 if sugar==0
gen total_score1=veg_score+fruit_score+meat_score+sugar_score
tab total_score1
gen diet_score1=0
replace diet_score1=2 if total_score1==4
replace diet_score1=1 if total_score1==3 | total_score1==2 
gen lifestyle_score2=bmi_score+diet_score1+smoke_score+MET_score // 改进
spearman lifestyle_score2 observed_features // 0.08
spearman lifestyle_score2 g220 // -0.08
spearman lifestyle_score2 city_score  //  -0.2
save "D:\真菌分析\data_pre\2015_repeat_datasets_lifescore.dta",replace


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

logit t2d c.g220##c.risk 

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

melogit t2d std_risk   ||region:, covariance(unstructured) or // 2.38
melogit t2d std_g280   ||region:, covariance(unstructured) or // 1.10

melogit t2d c.g280##c.age_score   ||region:, covariance(unstructured) or // 交互项边缘显著, p=0.064
melogit t2d c.g280##c.BMI_score   ||region:, covariance(unstructured) or // 交互项不显著
melogit t2d c.g280##c.wc_score   ||region:, covariance(unstructured) or // 交互项不显著
melogit t2d c.g280##c.sys_score   ||region:, covariance(unstructured) or // 交互项不显著


*---2018--

use "D:\真菌分析\data_pre\2018_repeat_datasets.dta",clear
drop g2-g633
rename SampleID_18 SampleID
merge 1:1 SampleID using "D:\真菌分析\data_pre\2018_repeat_fungi_clr.dta"
keep if _merge==3
gen veg=dveg+lveg
 foreach var of varlist fruit veg  pork poultry milk sugar fish OIL_VEG MET nuts{
	xtile q_`var'=`var',nq(2)
	}
	
gen bmi_score=0
replace bmi_score=2 if BMI>=18.5 & BMI<=24.9
gen fruit_score=0
replace fruit_score=1 if q_fruit==2
gen veg_score=0
replace veg_score=1 if q_veg==2
gen meat_score=0
replace meat_score=1 if q_pork==1
gen total_score=veg_score+fruit_score+meat_score
gen diet_score=0
replace diet_score=2 if total_score==3
replace diet_score=1 if total_score==2
gen smoke_score=0
replace smoke_score=2 if smoke==0
gen MET_score=0
replace MET_score=2 if q_MET==2
gen sugar_score=0
replace sugar_score=1 if sugar==0
gen total_score1=veg_score+fruit_score+meat_score+sugar_score
tab total_score1
gen diet_score1=0
replace diet_score1=2 if total_score1==4
replace diet_score1=1 if total_score1==3 | total_score1==2 
gen lifestyle_score2=bmi_score+diet_score1+smoke_score+MET_score // 改进
spearman lifestyle_score2 observed_features // 0.06
spearman lifestyle_score2 g220 // -0.058
spearman lifestyle_score2 city_score  //  -0.15
drop _merge
save "D:\真菌分析\data_pre\2018_repeat_datasets_lifescore.dta",replace



*--2015 TO 2018 transport and lifestylescore---

*--2015 lifescore and cluster---

use "D:\真菌分析\data_pre\2015_repeat_datasets_lifescore.dta",clear
drop _merge
merge 1:1 Idind using  "D:\真菌分析\data_pre\repeatDMM_clusterchange.dta"
keep if _merge==3
sort cluster_15
xtile q_score=lifestyle_score2,nq(3)
tab cluster_15 q_score
by cluster_15:su lifestyle_score2 // 4.34 4.58 4.78 4.25
tab cluster_15 lifestyle_score2
gr box lifestyle_score2,over(cluster_15)
hist lifestyle_score2 if cluster_15==1
hist lifestyle_score2 if cluster_15==4
regress  glucose lifestyle_score2
regress  glucose lifestyle_score2 if cluster_15==1
regress  glucose lifestyle_score2 if cluster_15==2
regress  glucose lifestyle_score2 if cluster_15==3
regress  glucose lifestyle_score2 if cluster_15==4


*--2018 lifescore and cluster---

use "D:\真菌分析\data_pre\2018_repeat_datasets_lifescore.dta",clear
merge 1:1 Idind using  "D:\真菌分析\data_pre\repeatDMM_clusterchange.dta"
keep if _merge==3
sort cluster_18
by cluster_18:su lifestyle_score2 // 4.4 4.6 4.65 4.32
gr box lifestyle_score2,over(cluster_18)

*--transport analysis---

use "D:\真菌分析\data_pre\2015_repeat_datasets_lifescore.dta",clear
keep Idind lifestyle_score2
rename lifestyle_score2 life_15
merge 1:1 Idind using "D:\真菌分析\data_pre\2018_repeat_datasets_lifescore.dta"
keep if _merge==3
drop _merge
rename lifestyle_score2 life_18
merge 1:1 Idind using  "D:\真菌分析\data_pre\repeatDMM_clusterchange.dta"
gen life_change=life_18-life_15
sort group
by group:su life_change age
gr box life_change,over(group)

xtile q_score_18=life_18,nq(4)
xtile q_score_15=life_15,nq(4)
gen q_change=q_score_18-q_score_15
sort group
by group:su q_change
order group q_score_15 q_score_18

gen Group=0
replace Group=1 if group==5 | group==9 | group==13 // n to 1
replace Group=2 if group==2  | group==10 | group==14 // n to 2
replace Group=3 if group==3 | group==7 | group==15
replace Group=4 if group==4 | group==8 | group==12 
tab Group
sort Group
by Group:su life_change
mlogit  Group life_change,rr base(4)
mlogit  Group q_change,rr base(4)
mlogit  Group c.q_change,rr
tab Group q_change


*--------------3.7 predict analysis----------

*------age predict analysis----

use "D:\真菌分析\data_pre\2015_all_dataset.dta",clear
drop if age==.
tab sex

preserve
keep if sex==1
keep SampleID-g751 age
save "D:\真菌分析\data_pre\2015_all_age_predict_male.dta",replace
restore

preserve 
keep if sex==2
keep SampleID-g751 age
save "D:\真菌分析\data_pre\2015_all_age_predict_female.dta",replace
restore

keep SampleID-g751 age
save "D:\真菌分析\data_pre\2015_all_age_predict.dta",replace

*-----diversity predict analysis----

use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
gen vegetable=dveg+lveg
keep SampleID age sex BMI  MET city_score city smoke alcohol rice-tuber dveg-othmeat poultry-sauce pastes vegetable faith_pd
rename faith_pd outcome
save "D:\真菌分析\data_pre\2015_all_dataset_diversity_predict.dta",replace

*-----city score predict analysis----

*---micro-cityscore--\

use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
keep SampleID g1-g751 pielou_evenness-shannon_entropy city_score
rename city_score outcome
save "D:\真菌分析\data_pre\2015_all_dataset_city_predict.dta",replace

*---micro+diet-cityscore--\

use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
keep SampleID g1-g751 pielou_evenness-shannon_entropy rice-pastes smoke alc alcohol MET city_score
rename city_score outcome
save "D:\真菌分析\data_pre\2015_all_dataset_city_predict_microdiet.dta",replace

*---diet-cityscore--\

use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
keep SampleID  rice-pastes smoke alc alcohol MET city_score
rename city_score outcome
save "D:\真菌分析\data_pre\2015_all_dataset_city_predict_diet.dta",replace



use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
keep SampleID g1-g751 pielou_evenness-shannon_entropy rice-pastes city_score
rename city_score outcome
save "D:\真菌分析\data_pre\2015_all_dataset_city_predict_.dta",replace



use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
keep SampleID g1-g751 city
rename city outcome
save "D:\真菌分析\data_pre\2015_all_dataset_urban_predict.dta",replace


use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
merge 1:1 SampleID using "D:\真菌分析\data_pre\2015_all_fungi_pc.dta"
keep if _merge==3
rename dim2 outcome
keep SampleID age t2d city-marrige smoke-pet outcome
save "D:\真菌分析\data_pre\2015_all_dataset_pc2_predict.dta",replace

*------------------dietary and lifestyle factors predict fungi--------

*----------------Discovery cohort----------

*---diet data--

use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
drop g1-g751
merge 1:1 SampleID using "D:\真菌分析\data_pre\fungi_clr.dta"
gen vegetable=dveg+lveg
keep SampleID age sex BMI  MET city_score city smoke alcohol rice-tuber dveg-othmeat poultry-sauce pastes vegetable
export excel using "D:\真菌分析\Predict_analysis\diet_fungi.xlsx",  firstrow(variables) sheet("diet") sheetreplace  

*---clr-based---

use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
drop g1-g751
merge 1:1 SampleID using "D:\真菌分析\data_pre\fungi_clr.dta"
keep SampleID pielou_evenness faith_pd observed_features shannon_entropy g1-g751
export excel using "D:\真菌分析\Predict_analysis\diet_fungi.xlsx",  firstrow(variables) sheet("fungi_clr") sheetreplace  

*---raw-micro data--

use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
keep SampleID pielou_evenness faith_pd observed_features shannon_entropy g1-g751
export excel using "D:\真菌分析\Predict_analysis\diet_fungi.xlsx",  firstrow(variables) sheet("fungi_raw") sheetreplace  

*--clr+sd transform---

use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
drop g1-g751
merge 1:1 SampleID using "D:\真菌分析\data_pre\fungi_clr.dta"
foreach v of varlist pielou_evenness faith_pd observed_features shannon_entropy g1-g751{
           egen s_`v' = std(`v')
		   drop `v'
		   rename s_`v' `v'
         }
keep SampleID pielou_evenness faith_pd observed_features shannon_entropy g1-g751
export excel using "D:\真菌分析\Predict_analysis\diet_fungi.xlsx",  firstrow(variables) sheet("fungi_clr_sd") sheetreplace  


import excel "D:\真菌分析\Predict_analysis\diet_fungi.xlsx", sheet("diet") firstrow clear 
drop age city city_score sex BMI
save "D:\真菌分析\data_pre\diet_fungi.dta",replace

*--dataset1---

import excel "D:\真菌分析\Predict_analysis\diet_fungi.xlsx", sheet("fungi_clr") firstrow clear 
merge 1:1 SampleID using "D:\真菌分析\data_pre\diet_fungi.dta"
drop _merge

*--dataset2---

import excel "D:\真菌分析\Predict_analysis\diet_fungi.xlsx", sheet("fungi_raw") firstrow clear 
merge 1:1 SampleID using "D:\真菌分析\data_pre\diet_fungi.dta"
drop _merge

*--dataset3---

import excel "D:\真菌分析\Predict_analysis\diet_fungi.xlsx", sheet("fungi_clr_sd") firstrow clear 
merge 1:1 SampleID using "D:\真菌分析\data_pre\diet_fungi.dta"
drop _merge

*--city/cityscore prediction--

use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
keep SampleID-g751 city_score
save "D:\真菌分析\data_pre\cityscore_predict_raw.dta",replace
 
use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
keep SampleID city_score
save "D:\真菌分析\data_pre\cityscoredata.dta",replace



*----regression explain mapping---

import delimited "D:\真菌分析\Final results\All data\diet_菌群解释度2.csv",  clear 
rename v2 micro
rename v3 number 
rename v4 factor
rename v5 r2
rename value p
merge 1:1 micro using "D:\真菌分析\data_pre\all_g_mapping.dta"
drop if _merge==2
sort r2
replace Genus="diversity" if Genus==""
export excel using "D:\真菌分析\Plot\data\explain_variation.xlsx",  firstrow(variables) sheet("diet") sheetreplace  


*----------------Vlidating results in the repeated cohort----------

*---2015 all and repeated fungi data mapping--

import excel "D:\真菌分析\rawdata\raw_fungi_mapping.xlsx", sheet("repeat") firstrow clear 
save  "D:\真菌分析\data_pre\raw_fungi_mapping.dta",replace

import excel "D:\真菌分析\rawdata\raw_fungi_mapping.xlsx", sheet("all") firstrow clear 
merge 1:1 g using  "D:\真菌分析\data_pre\raw_fungi_mapping.dta"
keep if _merge==3 // 558 matched (202 from all, 74 from repeated without matched)
drop _merge
save  "D:\真菌分析\data_pre\fungi_mapping_all_repeat.dta",replace

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

*---data pre--

use  "D:\真菌分析\data_pre\2015_all_phenotype.dta",clear
rename idind Idind
merge m:m Idind using "D:\真菌分析\data_pre\fugi_15_diversity.dta",force
keep if _merge==3
drop _merge 
merge m:m Idind using "D:\真菌分析\data_pre\fugi_15_g.dta"
keep if _merge==3
save "D:\CNHS\data\Data share\for CDC\2015_repeat_fungi_dataset.dta",replace

drop _merge
drop if observed_features==. | g2==. // 2041 obs
tab district
gen region=substr(SampleID,1,1)
tab region
replace region="Liaoning" if region=="A"
replace region="Jiangsu" if region=="C"
replace region="Hubei" if region=="F"
replace region="Guangxi" if region=="H"
gen vegetable=dveg+lveg

sort region
qui foreach var of varlist age MET income city_score hpc-tc BMI SBP DBP veg rice-pastes  {
     by region: egen `var'_mean= mean(`var')
	 by region: replace `var'=`var'_mean if `var'==.
}

foreach phenotype of varlist t2d dys hyp marrige  smoke alcohol xinjigengse zhongfeng cancer sex hyper_med diabetes_med antibiotic_current antibiotic_6month probiotics kangyan_med kangsuan_med weisuan_med yogurt_drink pet hyper_med city  hyp dys pret2d fuxie gut_disease changdao_shoushu{
	 replace `phenotype'=0 if `phenotype'==9 |`phenotype'==. 
}

misstable sum  age MET income city_score hpc-tc BMI SBP DBP vegetable rice-pastes t2d dys hyp marrige  smoke alcohol xinjigengse zhongfeng cancer sex hyper_med diabetes_med antibiotic_current antibiotic_6month probiotics kangyan_med kangsuan_med weisuan_med yogurt_drink pet hyper_med city  hyp dys pret2d fuxie gut_disease changdao_shoushu

foreach micro of varlist g2-g633 {
	count if `micro'==0 
	if r(N)>1837 {	
	drop  `micro' 
	}
}

save  "D:\真菌分析\data_pre\2015_repeat_datasets.dta",replace

*----2018 repeat data---

use  "D:\真菌分析\data_pre\2018_all_phenotype.dta",clear
rename idind Idind
drop _merge
merge m:m Idind using "D:\真菌分析\data_pre\fugi_18_diversity.dta",force
keep if _merge==3
drop _merge 
merge m:m Idind using "D:\真菌分析\data_pre\fugi_18_g.dta"
keep if _merge==3
save "D:\CNHS\data\Data share\for CDC\2018_repeat_fungi_dataset.dta",replace

drop _merge
drop if observed_features==. | g2==. // 2006 obs
tab district
gen region=substr(SampleID,1,1)
tab region
replace region="Liaoning" if region=="A"
replace region="Jiangsu" if region=="C"
replace region="Hubei" if region=="F"
replace region="Guangxi" if region=="H"
gen vegetable=dveg+lveg

sort region
qui foreach var of varlist age MET income city_score hpc-tc BMI SBP DBP veg rice-pastes  {
     by region: egen `var'_mean= mean(`var')
	 by region: replace `var'=`var'_mean if `var'==.
}

foreach phenotype of varlist t2d dys hyp marrige  smoke alcohol xinjigengse zhongfeng cancer sex hyper_med diabetes_med antibiotic_current antibiotic_6month probiotics kangyan_med kangsuan_med weisuan_med yogurt_drink pet hyper_med city  hyp dys pret2d fuxie gut_disease changdao_shoushu{
	 replace `phenotype'=0 if `phenotype'==9 |`phenotype'==. 
}

misstable sum  age MET income city_score hpc-tc BMI SBP DBP vegetable rice-pastes t2d dys hyp marrige  smoke alcohol xinjigengse zhongfeng cancer sex hyper_med diabetes_med antibiotic_current antibiotic_6month probiotics kangyan_med kangsuan_med weisuan_med yogurt_drink pet hyper_med city  hyp dys pret2d fuxie gut_disease changdao_shoushu

foreach micro of varlist g2-g633 {
	count if `micro'==0 
	if r(N)>1854 {	
	drop  `micro' 
	}
}

save  "D:\真菌分析\data_pre\2018_repeat_datasets.dta",replace


*-----clr-transform------

*---2015---

use  "D:\真菌分析\data_pre\2015_repeat_datasets.dta",clear
keep SampleID g2-g633
save  "D:\真菌分析\data_pre\2015_repeat_forclr.dta",replace

import delimited "D:\真菌分析\Final results\Repeated data\repeat_2015_fungi_clr.csv",  clear 
rename sampleid SampleID 
keep SampleID g2-g633
save "D:\真菌分析\data_pre\2015_repeat_fungi_clr.dta",replace


*---2018---

use  "D:\真菌分析\data_pre\2018_repeat_datasets.dta",clear
keep SampleID g2-g633
save  "D:\真菌分析\data_pre\2018_repeat_forclr.dta",replace

import delimited "D:\真菌分析\Final results\Repeated data\repeat_2018_fungi_clr.csv",  clear 
rename sampleid SampleID 
keep SampleID g2-g633
save "D:\真菌分析\data_pre\2018_repeat_fungi_clr.dta",replace

*---diet data--

use  "D:\真菌分析\data_pre\2015_repeat_datasets.dta",clear
keep SampleID age sex BMI  MET city_score city smoke alcohol rice-tuber dveg-othmeat poultry-sauce pastes vegetable
order SampleID
export excel using "D:\真菌分析\Predict_analysis\repeat_diet_fungi.xlsx",  firstrow(variables) sheet("diet") sheetreplace  

*---clr-based---

use "D:\真菌分析\data_pre\2015_repeat_datasets.dta",clear
drop g2-g633
merge 1:1 SampleID using "D:\真菌分析\data_pre\2015_repeat_fungi_clr.dta"
keep SampleID  pielou_evenness faith_pd observed_features shannon g106	g112	g113	g114	g116	g123	g146	g2	g18	g19	g186	g187	g188	g190	g194	g200	g201	g207	g212	g217	g220	g225	g226	g227	g230	g237	g243	g265	g278	g279	g280	g283	g284	g302	g327	g331	g399	g44	g461	g496	g507	g530	g538	g582	g583	g584	g597	g617	g623	g626	g63	g64	g73
export excel using "D:\真菌分析\Predict_analysis\repeat_diet_fungi.xlsx",  firstrow(variables) sheet("fungi_clr") sheetreplace  

*---validating the discovery associaitons--

use "D:\真菌分析\data_pre\2015_repeat_datasets.dta",clear
drop g2-g633
merge 1:1 SampleID using "D:\真菌分析\data_pre\2015_repeat_fungi_clr.dta"
keep if _merge==3
foreach v of varlist age MET city_score BMI HDL_C LDL_C HbA1c glucose tg tc SBP DBP rice-pastes pielou_evenness-shannon g2-g633 vegetable{
           egen s_`v' = std(`v')
		   drop `v'
		   rename s_`v' `v'
         }


tempname coef
tempfile res
postfile `coef' str200( micro var) float(n rr lul uul p)   using "`res'", replace


foreach outcome of varlist pielou_evenness-shannon g2-g633  {
foreach var of varlist   MET city_score  rice-tuber dveg-othmeat poultry-sauce pastes vegetable{ 
    mixed `outcome' `var' age i.sex BMI ||region:, covariance(unstructured) 
test `var'
 post `coef'   ("`outcome'") ("`var'")  (e(N)) (_b[`var'])  (_b[`var']-1.96*_se[`var'])  (_b[`var']+1.96*_se[`var'])  (chi2tail(1,(_b[`var']/_se[`var'])^2))  
 }
 }

  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\All data\associaiton analysis_repeated.xlsx",  firstrow(variables) sheet("diet_micro_continues") sheetreplace  
restore

*---diet、lifestyle、urbanisation to fungi(catgary vars)--

tempname coef
tempfile res
postfile `coef' str200( micro var) float(n rr lul uul p)   using "`res'", replace

gen rural=0
replace rural=1 if city==1 // 1-city

foreach outcome of varlist pielou_evenness-shannon g2-g633  {
foreach var of varlist  rural smoke alcohol { 
    mixed `outcome' i.`var' age i.sex BMI ||region:, covariance(unstructured) 
test 1.`var'
 post `coef'   ("`outcome'") ("`var'")  (e(N)) (_b[1.`var'])  (_b[1.`var']-1.96*_se[1.`var'])  (_b[1.`var']+1.96*_se[1.`var'])  (chi2tail(1,(_b[1.`var']/_se[1.`var'])^2))  
 }
 }

  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\All data\associaiton analysis_repeated.xlsx",  firstrow(variables) sheet("diet_micro_categray") sheetreplace  
restore


*-----fungi to phenotype------

tempname coef
tempfile res
postfile `coef' str200( var micro ) float(n rr lul uul p) str200(cov)  using "`res'", replace


foreach outcome of varlist  HDL_C LDL_C HbA1c glucose tg tc SBP DBP {
foreach var of varlist pielou_evenness-shannon g2-g633 { 
    mixed `outcome' `var' age i.sex BMI ||region:, covariance(unstructured) 
test `var'
 post `coef'   ("`outcome'") ("`var'")  (e(N)) (_b[`var'])  (_b[`var']-1.96*_se[`var'])  (_b[`var']+1.96*_se[`var'])  (chi2tail(1,(_b[`var']/_se[`var'])^2))  ("${cov`m'}")
 }
 }

  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\All data\associaiton analysis_repeated.xlsx",  firstrow(variables) sheet("fungi_phenotype") sheetreplace  
restore

*----fungi to diseases-------注意探究交互作用

*---model 1: adjust age+sex

xtile =g220,nq(4)
melogit t2d g220 age i.sex  ||region:, covariance(unstructured) or // p=0.067
melogit t2d g220 age i.sex BMI  ||region:, covariance(unstructured) or // 不显著
melogit t2d i.q_g220 age i.sex  ||region:, covariance(unstructured) or // p=0.089

tempname coef
tempfile res
postfile `coef' str200( micro var) float(n rr lul uul p) str200(cov)  using "`res'", replace


foreach outcome of varlist   t2d hyp dys {
foreach var of varlist pielou_evenness-shannon g220 g584 g190 g187 g265 g225 g123{ 
    melogit `outcome' `var' age i.sex  ||region:, covariance(unstructured) 
test `var'
 post `coef'   ("`outcome'") ("`var'")  (e(N)) (exp(_b[`var']))  (exp(_b[`var']-1.96*_se[`var']))  (exp(_b[`var']+1.96*_se[`var']))  (chi2tail(1,(_b[`var']/_se[`var'])^2))  ("${cov`m'}")
 }
 }

  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\All data\associaiton analysis_repeated.xlsx",  firstrow(variables) sheet("fungi_diseases_model1") sheetreplace  
restore

*--------name mapping------


*----continue diet vars--

import excel "D:\真菌分析\Final results\All data\associaiton analysis_repeated.xlsx", sheet("diet_micro_continues") firstrow clear
rename micro code_repeat
merge m:m code_repeat using "D:\真菌分析\data_pre\repeat_g_mapping.dta"
keep if _merge==3
order Phylum Class Order Family Genus micro var p
sort p
drop if p==.
tab var,missing
drop Kingdom _merge
merge m:m var using "D:\真菌分析\data_pre\phenotype_mapping.dta"
drop if p==.
replace Var="Vegetable" if var=="vegetable"
gen micro_var="genus"
replace micro_var="diversity" if Genus==""
replace Genus=code_repeat if Genus==""
keep Phylum Class Order Family Genus Var n rr lul uul p micro_var
order Phylum Class Order Family Genus Var n rr lul uul p micro_var
sort p
save  "D:\真菌分析\data_pre\diet_fungi_continue_repeat.dta",replace

*----categary diet vars--

import excel "D:\真菌分析\Final results\All data\associaiton analysis_repeated.xlsx", sheet("diet_micro_categray") firstrow clear
rename micro code_repeat
merge m:m code_repeat using "D:\真菌分析\data_pre\repeat_g_mapping.dta"
keep if _merge==3
order Phylum Class Order Family Genus micro var p
sort p
drop if p==.
tab var,missing
drop Kingdom _merge
replace var="city" if var=="rural"
merge m:m var using "D:\真菌分析\data_pre\phenotype_mapping.dta"
drop if p==.
gen micro_var="genus"
replace micro_var="diversity" if Genus==""
replace Genus=code_repeat if Genus==""
keep Phylum Class Order Family Genus Var n rr lul uul p micro_var
order Phylum Class Order Family Genus Var n rr lul uul p micro_var
sort p
save  "D:\真菌分析\data_pre\diet_fungi_cat_repeat.dta",replace
 
*---diet vars append-----

use  "D:\真菌分析\data_pre\diet_fungi_continue_repeat.dta",clear
append using   "D:\真菌分析\data_pre\diet_fungi_cat_repeat.dta"
export excel using "D:\真菌分析\Plot\data\mixed_associaitons.xlsx",  firstrow(variables) sheet("diet_city_fungi_repeat") sheetreplace  

*-----fungi-phenotype-

import excel "D:\真菌分析\Final results\All data\associaiton analysis_repeated.xlsx", sheet("fungi_phenotype") firstrow clear
rename micro code_repeat
merge m:m code_repeat using "D:\真菌分析\data_pre\repeat_g_mapping.dta"
keep if _merge==3
order Phylum Class Order Family Genus micro var p
sort p
drop if p==.
tab var,missing
drop Kingdom _merge
merge m:m var using "D:\真菌分析\data_pre\phenotype_mapping.dta"
drop if p==.
gen micro_var="genus"
replace micro_var="diversity" if Genus==""
replace Genus=code_repeat if Genus==""
keep Phylum Class Order Family Genus Var n rr lul uul p micro_var
order Phylum Class Order Family Genus Var n rr lul uul p micro_var
sort p
export excel using "D:\真菌分析\Plot\data\mixed_associaitons.xlsx",  firstrow(variables) sheet("fungi_phenotype_repeat") sheetreplace  

*------fungi-diseases meta analysis------

import excel "D:\真菌分析\Plot\data\mixed_associaitons.xlsx", sheet("fungi_disease_meta") firstrow clear
keep if Disease=="Type 2 diabetes"
sort cohort order
metan rr lul uul, label(namevar=cohort) nooverall  nowt nobox  dp(3) ///
by(Microbiome)fixed effect(rr) null(1)  ///
         boxopt( mcolor(navy8) msymbol(square) ) ///
		  pointopt( msymbol(square) mcolor(navy8) msize(small) ///
            mlabposition(1) ) ///
         ciopt( lcolor(navy8) lwidth(small) ) ///
		 	 xlab(0.7,0.8,0.9,1,1.1,1.2,1.3,1.4) force ///				 
 graphregion(fcolor(white)lcolor(white))  subtitle("")
graph export "D:\真菌分析\Plot\Association analysis\fungi_t2d.pdf", as(pdf) replace

import excel "D:\真菌分析\Plot\data\mixed_associaitons.xlsx", sheet("fungi_disease_meta") firstrow clear
keep if Disease=="Type 2 diabetes"
sort cohort order
metan rr lul uul, label(namevar=cohort) nooverall  nowt nobox  dp(3) ///
by(Microbiome)fixed effect(rr) null(1)  ///
         boxopt( mcolor(navy8) msymbol(square) ) ///
		  pointopt( msymbol(square) mcolor(navy8) msize(small) ///
            mlabposition(1) ) ///
         ciopt( lcolor(navy8) lwidth(small) ) ///
		 	 xlab(0.7,0.8,0.9,1,1.1,1.2,1.3,1.4) force ///				 
 graphregion(fcolor(white)lcolor(white))  subtitle("")
graph export "D:\真菌分析\Plot\Association analysis\fungi_t2d.pdf", as(pdf) replace

import excel "D:\真菌分析\Plot\data\mixed_associaitons.xlsx", sheet("fungi_disease_meta") firstrow clear
keep if Disease=="Dyslipidemia"
sort cohort order
metan rr lul uul, label(namevar=cohort) nooverall  nowt nobox  dp(3) ///
by(Microbiome)fixed effect(rr) null(1)  ///
         boxopt( mcolor(navy8) msymbol(square) ) ///
		  pointopt( msymbol(square) mcolor(navy8) msize(small) ///
            mlabposition(1) ) ///
         ciopt( lcolor(navy8) lwidth(small) ) ///
		 	 xlab(0.8,0.9,1,1.1,1.2,1.3) force ///				 
 graphregion(fcolor(white)lcolor(white))  subtitle("")
graph export "D:\真菌分析\Plot\Association analysis\fungi_dys.pdf", as(pdf) replace

import excel "D:\真菌分析\Plot\data\mixed_associaitons.xlsx", sheet("fungi_disease_meta") firstrow clear
keep if Disease=="Hypertension"
sort cohort order
metan rr lul uul, label(namevar=cohort) nooverall  nowt nobox  dp(3) ///
by(Microbiome)fixed effect(rr) null(1)  ///
         boxopt( mcolor(navy8) msymbol(square) ) ///
		  pointopt( msymbol(square) mcolor(navy8) msize(small) ///
            mlabposition(1) ) ///
         ciopt( lcolor(navy8) lwidth(small) ) ///
		 	 xlab(0.8,0.9,1,1.1,1.2,1.3) force ///				 
 graphregion(fcolor(white)lcolor(white))  subtitle("")
graph export "D:\真菌分析\Plot\Association analysis\fungi_hyp.pdf", as(pdf) replace

*-----3.4-------disease or phenotype predict analysis(匹配vs不匹配策略)----

*----------不匹配策略------不匹配预测更多是年龄与BMI driven的

*---age+sex+BMI

use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
keep SampleID age sex BMI 
export excel using "D:\真菌分析\Predict_analysis\fungi_disease.xlsx",  firstrow(variables) sheet("agesexbmi") sheetreplace  

*---age+sex+BMI+fungi---

use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
drop g1-g751
merge 1:1 SampleID using "D:\真菌分析\data_pre\fungi_clr.dta"
keep SampleID age sex BMI pielou_evenness-shannon_entropy g1-g751
export excel using "D:\真菌分析\Predict_analysis\fungi_disease.xlsx",  firstrow(variables) sheet("agesexbmi_fungi") sheetreplace  

*---disease or phenotype data--

use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
keep SampleID t2d hyp dys pret2d xinjigengse zhongfeng cancer 
tab pret2d
replace pret2d=1 if pret2d==2
inspect
gen disease=0
replace disease=1 if hyp==1 |dys==1 | pret2d==1 | xinjigengse==1 | zhongfeng==1 | cancer==1
tab disease,missing // 6737 cases
export excel using "D:\真菌分析\Predict_analysis\fungi_disease.xlsx",  firstrow(variables) sheet("disease") sheetreplace  

use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
keep SampleID HDL_C LDL_C HbA1c ins glucose tg tc
export excel using "D:\真菌分析\Predict_analysis\fungi_disease.xlsx",  firstrow(variables) sheet("phenotype") sheetreplace  

*--------匹配策略------https://www.ssc.wisc.edu/sscc/pubs/stata_psmatch.htm

*-----匹配策略下的统计分析在一定程度可以延伸出healthy aging的概念（同样的年龄匹配下疾病风险不一样）

*--matched data explore

import delimited "D:\真菌分析\Final results\All data\t2d_case_control_matched.csv",  clear 
gen pair_id=_n
preserve
keep v1 pair_id
gen marker=1
rename v1 SampleID
save "D:\真菌分析\data_pre\t2d_match1.dta",replace // cases
restore

preserve
keep v2 pair_id
gen marker=2
rename v2 SampleID
save "D:\真菌分析\data_pre\t2d_match2.dta",replace // controls
restore

use "D:\真菌分析\data_pre\t2d_match1.dta",clear
append using "D:\真菌分析\data_pre\t2d_match2.dta"
merge m:m SampleID using "D:\真菌分析\data_pre\2015_all_dataset_fill.dta"
keep if _merge==3

tab marker t2d
sort t2d
by t2d: su age BMI
tab t2d sex
egen std_g280=std(g280)
clogit t2d c.std_g280 ,group(pair_id) or //非常显著
keep SampleID pair_id age sex BMI t2d
save "D:\真菌分析\data_pre\t2d_id.dta",replace 
 
*----data for t2d prediciton---不太适合做预测分析，因为无法考虑配对信息-----

use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
drop g1-g751
merge 1:1 SampleID using "D:\真菌分析\data_pre\fungi_clr.dta"
drop _merge
spearman age g280

preserve
keep SampleID age sex   pielou_evenness-shannon_entropy g280 g702 g305 g701 g149 g277 g717 g139 t2d
keep SampleID age sex   pielou_evenness-shannon_entropy g280 g702 t2d
save "D:\真菌分析\data_pre\t2d_all_predict_select.dta",replace 
restore

merge 1:1 SampleID using "D:\真菌分析\data_pre\t2d_id.dta"
keep if _merge==3
egen std_g280=std(g280)
egen std_shannon_entropy=std(shannon_entropy)
clogit t2d c.std_g280 ,group(pair_id) or //非常显著
clogit t2d c.std_shannon_entropy ,group(pair_id) or //非常显著
keep SampleID pair_id age sex BMI  pielou_evenness-shannon_entropy g1-g751 t2d
save "D:\真菌分析\data_pre\pair_fungi_t2d.dta",replace 

keep SampleID age sex BMI  pielou_evenness-shannon_entropy g1-g751 t2d
keep SampleID age sex BMI  pielou_evenness-shannon_entropy g280 g702 g305 g701 g149 g277 g717 g139 t2d
save "D:\真菌分析\data_pre\t2d_match_predict.dta",replace 
keep SampleID t2d age sex BMI
save "D:\真菌分析\data_pre\t2d_match_agesexbmipredict.dta",replace 


*----------------Prediction results pre----------------

import excel "D:\真菌分析\Predict_analysis\results\diet_fungi_abs_shap1.xlsx", sheet("Sheet2") firstrow clear 
foreach micro of varlist pielou_evenness-g751 {
    su `micro'
	replace `micro'=(`micro'-r(min))/(r(max)-r(min))
}
xpose,clear
export excel using "D:\真菌分析\Predict_analysis\results\diet_fungi_abs_shap1.xlsx",  firstrow(variables) sheet("rank") sheetreplace  

import excel "D:\真菌分析\Predict_analysis\results\diet_fungi_abs_shap1.xlsx", sheet("Sheet3") firstrow clear 
rename A micro
merge 1:1 micro using "D:\真菌分析\data_pre\all_g_mapping.dta"
replace Genus=micro if Genus==""
keep if _merge==3
order Phylum Class Order Family Genus micro 
keep Genus age-vegetable
export excel using "D:\真菌分析\Predict_analysis\results\diet_fungi_abs_shap1.xlsx",  firstrow(variables) sheet("rank_") sheetreplace  



*---------------4-----fungi and metabolics------------

*---COIA analysis--

use "D:\真菌分析\data_pre\fugi_15_g_predict1.dta",clear
merge 1:1 SampleID using  "D:\真菌分析\data_pre\metabolism_15.dta"
keep if _merge==3

preserve
keep SampleID g2-g633
export excel using "D:\真菌分析\Final results\microbiota metabolism\raw.xlsx",  firstrow(variables) sheet("fungi") sheetreplace  
restore

preserve
keep SampleID CMPF_N-mwxq04_N
export excel using "D:\真菌分析\Final results\microbiota metabolism\raw.xlsx",  firstrow(variables) sheet("metabolism") sheetreplace  
restore

*-对比细菌--

use "D:\真菌分析\data_pre\bac_15_g_predict.dta",clear
merge 1:1 SampleID using  "D:\真菌分析\data_pre\metabolism_15.dta"
keep if _merge==3

preserve
keep  bac_g32-bac_g921
export excel using "D:\真菌分析\Final results\microbiota metabolism\raw.xlsx",  firstrow(variables) sheet("bac") sheetreplace  
restore

*--------fungi predict metabolism-----

use "D:\真菌分析\data_pre\2015_repeat_fungi_clr.dta",clear
merge 1:1 SampleID using  "D:\真菌分析\data_pre\metabolism_15.dta"
keep if _merge==3

preserve
keep SampleID g2-g633
export excel using "D:\真菌分析\Predict_analysis\fungi_metabolism.xlsx",  firstrow(variables) sheet("fungi") sheetreplace  
restore

preserve
keep SampleID CMPF_N-mwxq04_N
foreach metabolism of varlist CMPF_N-mwxq04_N {
	gen ln_`metabolism'=log(`metabolism') 
	drop `metabolism'
	rename ln_`metabolism' `metabolism'
}
export excel using "D:\真菌分析\Predict_analysis\fungi_metabolism.xlsx",  firstrow(variables) sheet("metabolism") sheetreplace  
restore

*---------------5----Dietary diversity and fungi------

import delimited "D:\真菌分析\Final results\All data\diet_diversity.csv",  clear 
rename sampleid SampleID
merge 1:1 SampleID using  "D:\真菌分析\data_pre\fugi_15all_g_diversity.dta"
keep if _merge==3
drop _merge
merge m:m SampleID using  "D:\真菌分析\data_pre\2015_all_phenotype.dta"
keep if _merge==3
spearman shannon shannon_entropy // -0.0485
logit t2d shannon,or // 1.42
spearman shannon BMI // 0,03

*----------------6-------Species-level analysis--------

*-----s-level and diseases---

*--data pre--

use "D:\真菌分析\data_pre\2015_all_dataset_fill.dta",clear
merge 1:1 SampleID using "D:\真菌分析\data_pre\fugi_15all_s.dta"
keep if _merge==3
count if s335==0
count if s336==0 //显著
count if s337==0
foreach outcome of varlist  HDL_C LDL_C HbA1c glucose tg tc SBP DBP {
foreach var of varlist s372 { 
    mixed `outcome'  age i.sex BMI ||region:, covariance(unstructured) 
 }
 }

drop _merge
foreach micro of varlist s1-s939 {
	count if `micro'==0 
	if r(N)>9626 {	
	drop  `micro' 
	}
}
keep SampleID region s1-s922
save "D:\真菌分析\data_pre\2015_all_s_raw.dta",replace

*--clr-transform--

import delimited "D:\真菌分析\Final results\All data\fungi_clr_s.csv",  clear 
rename sampleid SampleID 
save "D:\真菌分析\data_pre\fungi_clr_s.dta",replace
merge 1:1 SampleID using "D:\真菌分析\data_pre\2015_all_dataset_fill.dta"
keep if _merge==3

foreach v of varlist age MET-tc BMI SBP DBP rice-pastes pielou_evenness-shannon_entropy s1-s922 {
           egen s_`v' = std(`v')
		   drop `v'
		   rename s_`v' `v'
         }
save "D:\真菌分析\data_pre\2015_all_s_phenotype.dta",replace

*-----fungi to phenotype------

use "D:\真菌分析\data_pre\2015_all_s_phenotype.dta",clear
tempname coef
tempfile res
postfile `coef' str200( var micro ) float(n rr lul uul p) str200(cov)  using "`res'", replace


foreach outcome of varlist  HDL_C LDL_C HbA1c glucose tg tc SBP DBP {
foreach var of varlist s1-s922 { 
    mixed `outcome' `var' age i.sex BMI ||region:, covariance(unstructured) 
test `var'
 post `coef'   ("`outcome'") ("`var'")  (e(N)) (_b[`var'])  (_b[`var']-1.96*_se[`var'])  (_b[`var']+1.96*_se[`var'])  (chi2tail(1,(_b[`var']/_se[`var'])^2))  ("${cov`m'}")
 }
 }

  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\All data\associaiton analysis.xlsx",  firstrow(variables) sheet("fungi_s_phenotype") sheetreplace  
restore

*---fungi to diseases--

tempname coef
tempfile res
postfile `coef' str200( micro var) float(n rr lul uul p) str200(cov)  using "`res'", replace


foreach outcome of varlist   t2d hyp dys {
foreach var of varlist pielou_evenness-shannon_entropy g1-g751 { 
    melogit `outcome' `var' age i.sex  ||region:, covariance(unstructured) 
test `var'
 post `coef'   ("`outcome'") ("`var'")  (e(N)) (exp(_b[`var']))  (exp(_b[`var']-1.96*_se[`var']))  (exp(_b[`var']+1.96*_se[`var']))  (chi2tail(1,(_b[`var']/_se[`var'])^2))  ("${cov`m'}")
 }
 }

  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\真菌分析\Final results\All data\associaiton analysis.xlsx",  firstrow(variables) sheet("fungi_diseases_model1") sheetreplace  
restore

*-------------------Validated the results on the N OF 1-------

*--CLR transform---

use "D:\真菌分析\rawdata\N of 1 fungi\WE-MACNUTR1_ITS_rawdata.dta" ,clear
drop p1-f72 s1-s117
save "D:\真菌分析\data_pre\nof1_s_phenotype.dta",replace

foreach micro of varlist g1-g102 {
	count if `micro'==0 
	if r(N)>300 {	
	drop  `micro' 
	}
}

preserve 
keep sampleid g1-g102 
save "D:\真菌分析\data_pre\nof1_s_clr.dta",replace
restore

preserve
drop g1-g102
save "D:\真菌分析\data_pre\nof1_phenotype.dta",replace
restore

import delimited "D:\真菌分析\Final results\nof1_clr.csv",  clear 
merge 1:1 sampleid using "D:\真菌分析\data_pre\nof1_phenotype.dta"
keep if _merge==3
save "D:\真菌分析\data_pre\nof1_phenotype_fungiclr.dta",replace

*--Metabolism data--

import excel "D:\真菌分析\rawdata\N of 1 metabolism\wemac1_metabolites_serum_phase1.xlsx", sheet("sMet") firstrow clear
destring MW0010438-MEDP2333,replace
keep set interv id  MW0010438-MEDN2135
gen SampleID=set+interv+id
save "D:\真菌分析\data_pre\nof1_metabolism.dta",replace

use  "D:\真菌分析\data_pre\nof1_phenotype_fungiclr.dta",clear
gen SampleID=id1+PN
drop _merge
merge 1:1 SampleID using "D:\真菌分析\data_pre\nof1_metabolism.dta"
keep if _merge==3
sort interv
by interv: su g29 g43

sort id
egen code_id=group(id)
sort interv
egen interv_id=group(interv)
xtset code_id interv_id
xtgee  MEDN0378054 g43 age i.sex BMI, family(gaussian ) link(identity) // C18:3-sacc
xtgee  MEDN1269 g43 age i.sex BMI, family(gaussian ) link(identity) // LP-sacc

xtgee  MEDN0378054 g29 age i.sex BMI, family(gaussian ) link(identity) // C18:3-blu
xtgee  MEDP1778110 g29 age i.sex BMI, family(gaussian ) link(identity) // L-Nor-blu
xtgee  MEDP0026110 g29 age i.sex BMI, family(gaussian ) link(identity) // L-va-blu
xtgee  MEDN1840054 g29 age i.sex BMI, family(gaussian ) link(identity) // Pin

mixed MEDN0378054 g43 age i.sex BMI ||interv_id:code_id, covariance(unstructured) 

*----fungal and glucose treat---

import delimited "D:\真菌分析\rawdata\N of 1 fungi\n_of_1_ITS_CGM指标_forwanglong_1214.csv",  clear 
sort pn
egen code_id=group(pn)
sort subgroup
egen interv_id=group(id1)
xtset code_id interv_id
xtgee  pmg_bamean_set treatment_mean_ra_g29_set age i.sex, family(gaussian ) link(identity) // 

fasting_bamean_set hpt_bamean_set auc24_bamean_set mage_bamean_set pmg_bamean_set  treatment_mean_ra_g29_set

*================phase2 dataset========= no sig

import excel "D:\真菌分析\rawdata\N of 1 metabolism\metab_phase2.xlsx", sheet("Sheet2") firstrow clear
save "D:\真菌分析\data_pre\nof1_metabolism2.dta",replace

use  "D:\真菌分析\data_pre\nof1_phenotype_fungiclr.dta",clear
gen SampleID=id1+PN
drop _merge
merge 1:1 SampleID using "D:\真菌分析\data_pre\nof1_metabolism2.dta"
keep if _merge==3

preserve
keep if id1=="S1G1"
regress  MEDN0390 g43 age i.sex BMI  // EPA-sacc // Sig
regress  MEDN0378 g43 age i.sex BMI  // C18:3-sacc // 
regress  MEDN1269*088 g43 age i.sex BMI  // LP-sacc // 
regress  MEDN0378 g29 age i.sex BMI  // C18:3-blu // 
regress  MEDP1778*507 g29 age i.sex BMI  // L-Nor-blu // 
regress  MEDP0026*507 g29 age i.sex BMI  // L-va-blu
regress  MEDN1840 g29 age i.sex BMI   // Pin
restore




sort PN
egen code_id=group(PN)
sort id1
egen interv_id=group(id1)
xtset code_id interv_id
   
xtgee  MEDN0390 g43 age i.sex BMI, family(gaussian ) link(identity) // EPA-sacc // 0.49
xtgee  MEDN0378 g43 age i.sex BMI, family(gaussian ) link(identity) // C18:3-sacc // 0.59
xtgee  MEDN1269*088 g43 age i.sex BMI, family(gaussian ) link(identity) // LP-sacc // 0.27

xtgee  MEDN0378 g29 age i.sex BMI, family(gaussian ) link(identity) // C18:3-blu // 0.035 (负相关)
xtgee  MEDP1778*507 g29 age i.sex BMI, family(gaussian ) link(identity) // L-Nor-blu // 0.25
xtgee  MEDP0026*507 g29 age i.sex BMI, family(gaussian ) link(identity) // L-va-blu
xtgee  MEDN1840 g29 age i.sex BMI, family(gaussian ) link(identity) // Pin

mixed MEDN0378 g43 age i.sex BMI ||interv_id:code_id, covariance(unstructured) 

*=================no-trans----垃圾

use  "D:\真菌分析\data_pre\nof1_s_phenotype.dta",clear
gen SampleID=id1+PN
merge 1:1 SampleID using "D:\真菌分析\data_pre\nof1_metabolism2.dta"
keep if _merge==3

sort PN
egen code_id=group(PN)
sort id1
egen interv_id=group(id1)
xtset code_id interv_id
   
xtgee  MEDN0390 g43 age i.sex BMI, family(gaussian ) link(identity) // EPA-sacc // 0.49
xtgee  MEDN0378 g43 age i.sex BMI, family(gaussian ) link(identity) // C18:3-sacc // 0.59
xtgee  MEDN1269*088 g43 age i.sex BMI, family(gaussian ) link(identity) // LP-sacc // 0.27

xtgee  MEDN0378 g29 age i.sex BMI, family(gaussian ) link(identity) // C18:3-blu // 0.035 (负相关)
xtgee  MEDP1778*507 g29 age i.sex BMI, family(gaussian ) link(identity) // L-Nor-blu // 0.25
xtgee  MEDP0026*507 g29 age i.sex BMI, family(gaussian ) link(identity) // L-va-blu
xtgee  MEDN1840 g29 age i.sex BMI, family(gaussian ) link(identity) // Pin
