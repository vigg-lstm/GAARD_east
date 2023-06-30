# This script identifies the 7 new Cyp6aap CNV alleles identified in this dataset. It cannot be run
# from the GitHub repo because it requires an output from the Ag1000G CNV pipeline that is too large
# to include ('target_regions_analysis.Rdata')

library('Biostrings')
library('stringr')
library('stringi')
library(RColorBrewer)
library(parallel)

load('target_regions_analysis.Rdata')

known.cnvs <- list(ace1   = list(Dup1  = list(FA = matrix(c(3436850,3639550,3437150,3639850), 2, 2),
					                         #BP = data.frame(pos = c(3436926, 3639836), seq = c('GCGAA', 'GGAAT')))
                                              BP = data.frame(pos = c(3436926, 3639836), seq = c('GCGAA', 'TTGTT'))),
                                 Dup2  = list(FA = matrix(c(3447950,3518600,3448250,3518900), 2, 2)),
                                 Del1  = list(FM = matrix(c(3501850,3598750,3502150,3599050), 2, 2)),
                                 Del2  = list(FM = matrix(c(3539300,3573450,3539600,3573750), 2, 2)),
                                 Del3  = list(FM = matrix(c(3535850,3618700,3536150,3619000), 2, 2)),
                                 Del4  = list(FM = matrix(c(3512200,3615990,3512500,3616290), 2, 2))),
                   cyp6   = list(Dup1  = list(FA = matrix(c(28480150, 28483200, 28480450, 28483550), 2, 2)),
                                 Dup1a = list(BP = data.frame(pos = c(28480189, 28483475), 
                                                              seq = c('CGTAG', 'AATTG'))),
                                 Dup1b = list(BP = data.frame(pos = c(28480193, 28483675), 
                                                              seq = c('CTGCT', 'CCTTC'))),
                                 Dup2  = list(FA = matrix(c(28493450, 28497000, 28493750, 28497300), 2, 2),
                                              BP = data.frame(pos = c(28493547, 28497279), 
                                                              seq = c('GCCGC','TTTAA'))),
                                 Dup3  = list(FA = matrix(c(28479350, 28483100, 28479650, 28483400), 2, 2),
                                              BP = data.frame(pos = c(28479407, 28483372),
                                                              seq = c('GCTTA', 'CAAAG'))),
                                 Dup4  = list(FA = matrix(c(28478850, 28482750, 28479150, 28483050), 2, 2),
                                              BP = data.frame(pos = c(28478925, 28483069),
                                                              seq = c('TACTT', 'CATGT'))),
                                 Dup5  = list(FA = matrix(c(28480300, 28484200, 28480600, 28484500), 2, 2), 
                                              BP = data.frame(pos = c(28480372, 28484518),
                                                              seq = c('AAGAG', 'ACAAA'))),
                                 Dup6  = list(FA = matrix(c(28478150, 28483850, 28478450, 28484150), 2, 2),
                                              BP = data.frame(pos = c(28478272, 28484157),
                                                              seq = c('ATCAC', 'CTAGA'))),
                                 Dup7  = list(SS = matrix(c(28478000, 28486000, 28478300, 28486300), 2, 2),
                                              BP = data.frame(pos = c(28478057, 28486036),
                                                              seq = c('AGAGC','TTTTT'))),
                                 Dup8  = list(FA = matrix(c(28475900, 28484700, 28476200, 28485000), 2, 2),
                                              BP = data.frame(pos = c(28475996, 28485005),
                                                              seq = c('AGCGA', 'CAAAT'))),
                                 Dup9  = list(FA = matrix(c(28479100, 28491200, 28479400, 28491500), 2, 2),
                                              BP = data.frame(pos = c(28479181, 28491431),
                                                              seq = c('TGTTC', 'TGTGG'))),
                                 Dup10 = list(FA = matrix(c(28477800, 28490850, 28478100, 28491150), 2, 2),
                                              BP = data.frame(pos = c(28477889, 28491215),
         # It turns out that Dup10 is not quite as simple as made out in the Supp Mat for the Genome Research paper.
         # There is actually some kind of insertion / mutation happening around the breakpoint, and the aligners used
         # for this experiment and in Ag1000G deal with this differently (so we need to check the phase3 data and 
         # onwards to see what happens there since they used a different aligner to phase 2). We therefore need to change
         # the sequence slightly here.
                                                              seq = c('TGTAG','AACTT'))),
                                                             #seq = c('TGTAG','ACTCT'))),
                                 Dup11 = list(FA = matrix(c(28487450, 28517800, 28487750, 28518100), 2, 2),
                                              BP = data.frame(pos = c(28487546, 28518123),
                                                              seq = c('AACAC', 'TTATC'))),
                                 Dup12 = list(FA = matrix(c(28474450, 28519650, 28474750, 28519950), 2, 2),
                                              BP = data.frame(pos = c(28474576, 28520016),
                                                              seq = c('CCGAC', 'ACGGT'))),
                                 Dup13 = list(FA = matrix(c(28472650, 28522350, 28472950, 28522650), 2, 2),
                                              BP = data.frame(pos = c(28472728, 28522671),
                                                              seq = c('ACCGC', 'AGCTG'))),
                                 Dup14 = list(FA = matrix(c(28473800, 28563200, 28474100, 28563500), 2, 2),
                                              BP = data.frame(pos = c(28473874, 28563596),
                                                              seq = c('CCCAC', 'AGTTG'))),
                                 Dup15 = list(FA = matrix(c(28465600, 55958800, 28465900, 55959100), 2, 2),
                                              BP = data.frame(pos = c(28465673, NA),
                                                              seq = c('CAGCC', NA))),
                                 Dup16 = list(FA = matrix(c(28480500, 28483300, 28480800, 28483600), 2, 2),
                                              BP = data.frame(pos = c(28480547, 28483236),
                                                               seq = c('CCATT', 'TTAGT'))),
                                 Dup17 = list(FA = matrix(c(28477500, 28484900, 28477800, 28485200), 2, 2),
                                              BP = data.frame(pos = c(28477540, 28485380),
                                                              seq = c('TGCTG', 'ATCGG'))),
                                 Dup18 = list(FA = matrix(c(28479500, 28494200, 28479800, 28494500), 2, 2),
                                              BP = data.frame(pos = c(28479548, 28494597),
                                                              seq = c('AGTCG', 'TTGTC'))),
                                 Dup19 = list(FA = matrix(c(28475480, 28556250, 28475780, 28556550), 2, 2),
                                              BP = data.frame(pos = c(28475490, 28556726),
                                                              seq = c('AATAG', 'TGTGT'))),
                                 Dup20 = list(FA = matrix(c(28473590, 28794750, 28473890, 28795050), 2, 2),
                                              BP = data.frame(pos = c(28473600, 28795255),
                                                              seq = c('ATACT', 'CAAAA'))),
                                 Dup21 = list(FA = matrix(c(28475100, 28483000, 28475400, 28483300), 2, 2),
                                              BP = data.frame(pos = c(28475128, 28483473),
                                                              seq = c('AGCCG', 'TGCAC'))),
                                 Dup22 = list(FA = matrix(c(28477200, 28484000, 28477500, 28484300), 2, 2),
                                              BP = data.frame(pos = c(28477220, 28484338),
                                                              seq = c('GTGGA', 'CGAGT'))),
                                 Dup23 = list(FA = matrix(c(28497300, 28371800, 28497600, 28372100), 2, 2),
                                              BP = data.frame(pos = c(NA, 28497740),
                                                              seq = c(NA, 'TTGGC'))),
                                 Dup24  = list(SS = matrix(c(28479500, 28483000, 28479800, 28483300), 2, 2),
                                               BP = data.frame(pos = c(28480585, 28483442),
                                                               seq = c('AAACA','TTAAC'))),
				   	             # For 25 and 26, the FA reads would overlap, so we just use the BP reads. 
                                 Dup25 = list(BP = data.frame(pos = c(28480335, 28483384),
                                                              seq = c('GGCGT', 'CATAT'))),
                                 Dup26 = list(BP = data.frame(pos = c(28480166, 28483253),
                                                              seq = c('AACGT', 'TGTGT'))),
                                 Dup27 = list(FA = matrix(c(28496700, 28498900, 28497000, 28499200), 2, 2)),
                                 Dup28 = list(FA = matrix(c(28477700, 28496600, 28478000, 28496900), 2, 2),
                                              BP = data.frame(pos = c(28477710, 28496953),
                                                              seq = c('CTGTA', 'ATTCT'))),
                                 Dup29 = list(FA = matrix(c(28494000, 28496000, 28494300, 28496300), 2, 2),
                                              BP = data.frame(pos = c(28494017, 28496505),
                                                              seq = c('TGGAA', 'TTTGC'))),
                                 Dup30 = list(FA = matrix(c(28478900, 28484700, 28479200, 28485000), 2, 2),
                                              BP = data.frame(pos = c(28478987, 28485033),
                                                              seq = c('AACAG', 'ACGTT'))),
                                 Dup31 = list(FA = matrix(c(28480450, 28492450, 28480750, 28492750), 2, 2),
                                              BP = data.frame(pos = c(28480476, 28492867),  
                                                              seq = c('GCTTC', 'ACGCC'))),
                                 Dup32 = list(FA = matrix(c(28485300, 28494300, 28485600, 28494600), 2, 2),
                                              BP = data.frame(pos = c(28485349, 28494719),  
                                                              seq = c('TATCG', 'CAGAC'))),
                                 Dup33 = list(FA = matrix(c(28479800, 28549500, 28480100, 28549800), 2, 2),
                                              BP = data.frame(pos = c(28479842, 28549886),  
                                                              seq = c('AAAGA', 'ACGGA'))),
                                 Dup34 = list(FA = matrix(c(28480250, 28482800, 28480550, 28483200), 2, 2),
                                              BP = data.frame(pos = c(28480284, 28483348),  
                                                              seq = c('GCCTC', 'ACGCT'))),
                                 Dup35 = list(FA = matrix(c(28472650, 28503550, 28472950, 28503850), 2, 2),
                                              BP = data.frame(pos = c(28472668, 28504018),  
                                                              seq = c('CACGC', 'TGCTT'))),
								 # For Dup35, the end of the CNV is just outside the region where we pulled out
								 # discordant reads, so we can't get the ending breakpoint. 
                                 Dup36 = list(FA = matrix(c(28475400, 28572300, 28475700, 28572600), 2, 2),
                                              BP = data.frame(pos = c(28475403, NA),  
                                                              seq = c('TCTTA', NA))),
                                 Dup37 = list(FA = matrix(c(28477150, 28489000, 28477450, 28489300), 2, 2),
                                              BP = data.frame(pos = c(28477192, 28489472),  
                                                              seq = c('ACTAG', 'GGTCT')))),
                   cyp6mz = list(Dupm1 = list(BP = data.frame(pos = c(6927942, NA),
                                                              seq = c('ATTAT', NA))),
                                 Dupm2 = list(FA = matrix(c(6933100, 6934900, 6933400, 6935200), 2, 2)),
                                 Dupm3 = list(FA = matrix(c(6929600, 6932500, 6929900, 6932800), 2, 2)),
                                 Dupm4 = list(FA = matrix(c(6929900, 6936600, 6930200, 6936900), 2, 2),
                                              BP = data.frame(pos = c(6929933, 6936902),
                                                              seq = c('TTAAA', 'TGTCG'))),
                                 Dupm5 = list(FA = matrix(c(6933800, 6938300, 6934100, 6938600), 2, 2),
                                              BP = data.frame(pos = c(6933972, 6938659),
                                                              seq = c('AAACC', 'GTCGG'))),
                                 Dupz1 = list(FA = matrix(c(6968950, 6979300, 6969250, 6979600), 2, 2),
                                              BP = data.frame(pos = c(6968962, 6979681),
                                                              seq = c('ACGCT', 'AGGTT'))),
                                 Dupz2 = list(FA = matrix(c(6975100, 6977100, 6975400, 6977400), 2, 2),
                                              BP = data.frame(pos = c(NA, 6977514), # clipped reads align at 6975066
                                                              seq = c(NA, 'TAAGA'))),
                                 Dupz3 = list(FA = matrix(c(6971450, 6977800, 6971750, 6978100), 2, 2),
                                              BP = data.frame(pos = c(6971484, NA),  # clipped reads align at 6978084
                                                              seq = c('GCAAA', NA))),
                                 Dupz4 = list(FA = matrix(c(6972700, 6977350, 6973000, 6977650), 2, 2),
                                              BP = data.frame(pos = c(6972775, 6977699),  
                                                              seq = c('GAATG', 'GTCCA'))),
                                 Dupz5 = list(FA = matrix(c(6969800, 6975700, 6970100, 6976000), 2, 2)),
                                 Dupmz1 = list(FA = matrix(c(6982700, 6879900, 6983000, 6880200), 2, 2))),
                   gste   = list(Dup1  = list(FA = matrix(c(6968950, 6979300, 6969250, 6979600), 2, 2),
                                              BP = data.frame(pos = c(28596818, 28598850),
                                                              seq = c('TTTTG', 'CGTTT'))),
                                 # The following definition includes some false positives. 
                                 Dup2  = list(BP.weak = data.frame(pos = c(28596390, 28598923),
                                                                   seq = c('GGGGG', 'TTCCC'))),
                                 Dup3  = list(FA = matrix(c(28590500, 28592950, 28590800, 28593250), 2, 2),
                                              BP = data.frame(pos = c(28590597, 28593254),
                                                              seq = c('TCAAA', 'AGGGC'))),
                                 Dup4  = list(FA = matrix(c(28595050, 28598750, 28595350, 28599050), 2, 2),
                                              BP = data.frame(pos = c(28595162, 28599081),
                                                              seq = c('TTCTA', 'AGAAC'))),
                                 Dup5  = list(BP = data.frame(pos = c(28593122, 28598971),
                                                              seq = c('GTCAT', 'ATTTA'))),
                                 Dup6  = list(FA = matrix(c(28596250, 28601900, 28596550, 28602200), 2, 2),
                                              BP = data.frame(pos = c(28596241, 28602177),
                                                              seq = c('ACAAC', 'GAAGC'))),
                                 # For the XC reads, the first two rows are the CNV start point, and the second two rows are the end point
                                 Dup7  = list(XC = data.frame(c(28597400, 3696450, 28603950, 26597300), c(28597700, 3696750, 28604250, 26597600), c('3R', '2L', '3R', 'UNKN')),
                                              BP = data.frame(pos = c(28597504, 28604250),
                                                              seq = c('GTCCA', 'GCTGT'))),
                                 Dup8  = list(BP = data.frame(pos = c(28594797, 28602349),
                                                              seq = c('GTCCC', 'CAGGG'))),
                                 Dup9  = list(FA = matrix(c(28591050, 28600850, 28591350, 28601150), 2, 2),
                                              BP = data.frame(pos = c(28591140, 28601188),
                                                              seq = c('AGAAG', 'GATGA'))),
                                 Dup10 = list(FA = matrix(c(28593550, 28603350, 28593850, 28603650), 2, 2),
                                              BP = data.frame(pos = c(28593642, 28603786),
                                                              seq = c('TCGCT', 'AAGAC'))),
                                 Dup11 = list(XC = data.frame(c(28581250, 29210650, 28604650, 29210650), c(28581550, 29210950, 28604950, 29210950), c('3R', 'UNKN', '3R', 'UNKN')),
                                              BP = data.frame(pos = c(28581256, 28604994),
                                                              seq = c('CCATT', 'GGTAA'))),
                                 Dup12 = list(FA = matrix(c(28597000, 28599950, 28597300, 28600250), 2, 2),
                                              BP = data.frame(pos = c(28597030, 28600292),
                                                              seq = c('TACTG', 'CATCT'))),
                                 Dup13 = list(FA = matrix(c(28597000, 28598900, 28597300, 28599200), 2, 2),
                                              BP = data.frame(pos = c(28597181, NA), # clipped sequences align at 28599287
                                                              seq = c('TACTC', NA))),
                                 Dup14 = list(FA = matrix(c(28599800, 28607200, 28600100, 28607500), 2, 2),
                                              BP = data.frame(pos = c(28599926, 28607500), 
                                                              seq = c('CGACG', 'ATGCA'))),
                                 Dup15 = list(FA = matrix(c(28596200, 28598500, 28596500, 28598800), 2, 2),
                                              BP = data.frame(pos = c(28596252, 28598948),
                                                              seq = c('TTGGA', 'TTGAC'))),
                                 Dup16 = list(FA = matrix(c(28597300, 28603200, 28597600, 28603500), 2, 2),
                                              BP = data.frame(pos = c(28597383, 28603517), 
                                                              seq = c('ACATT', 'ATTAC')))),
                   cyp9k1 = list(Dup1  = list(FA = matrix(c(15242500, 15244500, 15242800, 15244800), 2, 2),
                                              BP = data.frame(pos = c(15242505,15244812),
                                                            seq = c('GTTTG', 'CATAT'))),
                                 Dup2  = list(FA = matrix(c(15238300, 15240800, 15238600, 15241100), 2, 2),
                                              BP = data.frame(pos = c(15238400, 15241082),
                                                              seq = c('CCGGC',' CGGTA'))),
                                 Dup3  = list(FA = matrix(c(15240300, 15243450, 15240600, 15243750), 2, 2),
                                              BP = data.frame(pos = c(NA, 15243860),
                                                              seq = c(NA, 'TGAAC'))),
                                 Dup4  = list(FA = matrix(c(15240600, 15244200, 15240900, 15244500), 2, 2),
                                              BP = data.frame(pos = c(15240608, 15244503),
                                                              seq = c('ATAAA', 'ACTGG'))),
                                 Dup5  = list(FA = matrix(c(15238800, 15243850, 15239100, 15244150), 2, 2),
                                              BP = data.frame(pos = c(15238911, 15244175),
                                                              seq = c('CACGT', 'AGTAA'))),
                                 Dup6  = list(FA = matrix(c(15236400, 15243250, 15236700, 15243550), 2, 2),
                                              BP = data.frame(pos = c(15236449, 15243646),
                                                              seq = c('TTTTT', 'GTTTT'))),
                                 Dup7  = list(SS = matrix(c(15245400, 15246900, 15245700, 15247200), 2, 2),
                                              BP = data.frame(pos = c(15245768, 15247258),
                                                              seq = c('TTTGT', 'TCTAA'))),
                                 Dup8  = list(FA = matrix(c(15239200, 15247250, 15239500, 15247550), 2, 2),
                                              BP = data.frame(pos = c(15239276, 15247645),
                                                              seq = c('AACAT', 'TTGCT'))),
                                 Dup9  = list(FA = matrix(c(15239100, 15248900, 15239400, 15249200), 2, 2),
                                              BP = data.frame(pos = c(15239184, 15249314),
                                                              seq = c('GCACA', 'AGTAC'))),
                                 Dup10 = list(FA = matrix(c(15234900, 15244750, 15235200, 15245050), 2, 2),
                                              BP = data.frame(pos = c(15234989, 15245128),
                                                              seq = c('GCACC', 'CTGAA'))),
                                 Dup11 = list(FA = matrix(c(15236900, 15246800, 15237200, 15247100), 2, 2),
                                              BP = data.frame(pos = c(15236922, 15247159),
                                                              seq = c('CATTA', 'TATCT'))),
                                 Dup12 = list(FA = matrix(c(15234400, 15244350, 15234700, 15244650), 2, 2),
                                              BP = data.frame(pos = c(15234434, 15244702),
                                                              seq = c('AACAG', 'TACTA'))),
                                 Dup13 = list(FA = matrix(c(15240100, 15250250, 15240400, 15250550), 2, 2),
                                              BP = data.frame(pos = c(15240067, 15250575),
                                                              seq = c('CCTAA', 'GTGTA'))),
                                 # Dup 14 seems to have a different endpoint to Dup15, but the same insertion point way upstream
                                 Dup14 = list(FA = matrix(c(15244200, 9676400, 15244500, 9676700), 2, 2),
                                             # Because Dup14 has the same start pos as Dup15, we don't use the start breakpoint as a diagnostic
                                             #BP = data.frame(pos = c(15233807, 15244936),
                                             #                seq = c('GGGTT', 'CCCAA'))),
                                              BP = data.frame(pos = c(NA, 15244936),
                                                              seq = c(NA, 'CCCAA'))),
                                 Dup15 = list(FA = matrix(c(15246250, 9676400, 15246550, 9676700), 2, 2),
                                             # Because Dup14 has the same start pos as Dup15, we don't use the start breakpoint as a diagnostic
                                             #BP = data.frame(pos = c(15233807, 15246640),
                                             #                seq = c('GGGTT', 'CCCAA'))),
                                              BP = data.frame(pos = c(NA, 15246640),
                                                              seq = c(NA, 'CCCAA'))),
                                 Dup16 = list(FA = matrix(c(15222700, 15244300, 15223000, 15244600), 2, 2),
                                              BP = data.frame(pos = c(NA, 15244755),
                                                              seq = c(NA, 'AAGTA'))),
                                 Dup17 = list(FA = matrix(c(15237150, 15243650, 15237450, 15243950), 2, 2),
                                              BP = data.frame(pos = c(15237138, 15243975),
                                                              seq = c('TTGCT', 'TTTCG'))),
                                 Dup18 = list(FA = matrix(c(15236100, 15243500, 15236400, 15243800), 2, 2),
                                              BP = data.frame(pos = c(NA, 15243915),  # clipped sequence aligns to 15236175
                                                              seq = c(NA, 'CGGCG'))),
                                 Dup19 = list(FA = matrix(c(15238800, 15251100, 15239100, 15251400), 2, 2),
                                              BP = data.frame(pos = c(15238878, 15251503),  
                                                              seq = c('TAAAT', 'GTTAC'))),
                                 Dup20 = list(FA = matrix(c(15237350, 15243100, 15237650, 15243400), 2, 2),
                                              BP = data.frame(pos = c(15237397, 15243514),  
                                                              seq = c('ATGTT', 'TTACG'))),
                                 Dup21 = list(FA = matrix(c(15237450, 15250300, 15237750, 15250600), 2, 2),
                                              BP = data.frame(pos = c(15237482, 15250699),  
                                                              seq = c('CTCTG', 'TTCTC'))),
                                 Dup22 = list(FA = matrix(c(15240650, 15250300, 15240950, 15250600), 2, 2),
                                              BP = data.frame(pos = c(15240680, 15250670),  
                                                              seq = c('TTCCA', 'ATTCT'))),
                                 Dup23 = list(FA = matrix(c(15241800, 15248100, 15242100, 15248400), 2, 2),
                                              BP = data.frame(pos = c(15241929, 15248352),  
                                                              seq = c('AACAA', 'CACGT'))),
                                 Dup24 = list(FA = matrix(c(15238550, 15254800, 15238850, 15255100), 2, 2)),
                                 Dup25 = list(FA = matrix(c(15223500, 15246350, 15223800, 15246650), 2, 2)),
                                 Dup26 = list(FA = matrix(c(15222700, 15247750, 15223000, 15248050), 2, 2)),
                                 Dup27 = list(FA = matrix(c(15237400, 15248350, 15237700, 15248650), 2, 2),
                                              BP = data.frame(pos = c(15237566, NA),  
                                                              seq = c('AATGT', NA))),
                                 Dup28 = list(BP = data.frame(pos = c(NA, 15246640),  
                                                              seq = c(NA, 'TCGAG'))))
)

num.cores = 5
cat('Counting diagnostic reads\n')
diagnostic.read.counts <- mcmapply(count.diagnostic.reads.allsamples, diagnostic.reads, known.cnvs, mc.preschedule = F, mc.cores = num.cores)
cat('Calling read-based CNVs\n')
read.based.cnvs <- mapply(get.read.based.cnvs, cov.cnv.samples, diagnostic.read.counts, 
                          SIMPLIFY = F, MoreArgs = list(threshold.diagnostic.reads = 5))


full.cnv.table <- do.call(cbind, read.based.cnvs)
full.cnv.table <- cbind(rownames(full.cnv.table) %in% high.var.samples, full.cnv.table)
gene.cluster.names <- setNames(c('Ace1', 'Cyp6aap', 'Cyp6mz', 'Gstue', 'Cyp9k1'), 
                               names(diagnostic.read.counts))
# The column names for the full table will be with the gene cluster codes that we used in the phase 2 analysis.
colnames(full.cnv.table) <- c('High.var.sample', unlist(sapply(names(gene.cluster.names), function(x) paste(gene.cluster.names[x], colnames(read.based.cnvs[[x]]), sep = '_'))))
write.table(full.cnv.table, file = 'focal_region_CNV_table_extra_alleles.csv', sep = '\t', col.names = NA, quote = F)


# Plotting functions
n.colours <- 40
duplication.colours <- c(colorRampPalette(brewer.pal(12, 'Paired'))(n.colours), 'grey30', 'grey70')
#duplication.colours <- c(rgb(0.6509804, 0.8078431, 0.8901961), rgb(0.1215686, 0.4705882, 0.7058824), rgb(0.6980392, 0.8745098, 0.5411765), rgb(0.2, 0.627451, 0.172549), rgb(0.9843137, 0.6039216, 0.6), rgb(0.8901961, 0.1019608, 0.1098039), rgb(0.9921569, 0.7490196, 0.4352941), rgb(1, 0.4980392, 0), rgb(0.7921569, 0.6980392, 0.8392157), rgb(0.4156863, 0.2392157, 0.6039216), rgb(1,0,1,0.7), rgb(0.6941176, 0.3490196, 0.1568627), rgb(0.5,0.5,0,0.8), rgb(0,0.5,0.5,0.8), rgb(0.5, 0, 0.5, 0.8), 'yellow', 'lightblue', 'pink', 'orange', 'lightgreen', 'grey30', 'grey70')
names(duplication.colours) <- paste('Dup', c(as.character(1:n.colours), '1a', '1b'), sep = '')

plot.cyp6 <- function(this.sample, list.of.dups = NULL, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$cyp6[1], end.pos = plotting.ranges$cyp6[2]){
	start.index <- which(compact.hmm$cyp6[[this.sample]]$Position >= start.pos)[1]
	end.index <- tail(which(compact.hmm$cyp6[[this.sample]]$Position <= end.pos), 1)
	plot(compact.hmm$cyp6[[this.sample]]$Position[start.index : end.index], compact.hmm$cyp6[[this.sample]]$Normalised_coverage[start.index : end.index], main = this.sample, ylim = c(0,12))
	# We highlight along the x axis where the duplications are as predicted by the FA and SS reads
	if (!is.null(list.of.dups)){
		Dup.order <- paste('Dup', c(as.character(c(37:20, 19, 14, 15, 13:11, 18, 10:7, 17, 6:2)), '1a', '1b'), sep = '')
		list.of.dups <- list.of.dups[Dup.order]
		for (d in names(list.of.dups)[list.of.dups]){
			if (d == 'Dup15'){
				rect(known.cnvs$cyp6[[d]]$BP$pos[1], -0.6, 28555300, 0, col = duplication.colours[d], border = col)
				text(28555300, -0.4, d, cex = 1.2, adj = c(1,0))
			}
			else if (d == 'Dup23'){
				rect(28450000, -0.6, known.cnvs$cyp6[[d]]$BP$pos[2], 0, col = duplication.colours[d], border = col)
				text(known.cnvs$cyp6[[d]]$BP$pos[2], -0.4, d, cex = 1.2, adj = c(1,0))
			}
			else if (d == 'Dup27'){
				rect(28496700, -0.6, 28499200, 0, col = duplication.colours[d], border = col)
				text(28499200, -0.4, d, cex = 1.2, adj = c(1,0))
			}
			else{
				rect(known.cnvs$cyp6[[d]]$BP$pos[1], -0.6, known.cnvs$cyp6[[d]]$BP$pos[2], 0, col = duplication.colours[d], border = col)
				text(known.cnvs$cyp6[[d]]$BP$pos[2], -0.4, d, cex = 1.2, adj = c(1,0))
			}
		}
	}
	discordant.colours <- c(FA = 'blue', SS = 'cyan', FM = 'green', XC = 'red')
	y.ranges <- list(FA = c(0,3), SS = c(3,6), FM = c(6,9), XC = c(9,12))
	for (d in diagnostics){
		if (d == 'BP'){
			add.diagnostics(diagnostic.reads$cyp6[[this.sample]][[d]]$CEP$Position, this.col = 'brown', yrange = c(0,12))
			add.diagnostics(diagnostic.reads$cyp6[[this.sample]][[d]]$CSP$Position, this.col = 'pink', yrange = c(0,12))
		}
		else if (d == 'XC')
			add.diagnostics(diagnostic.reads$cyp6[[this.sample]][[d]]$Position, this.col = discordant.colours[d], yrange = y.ranges[[d]])
		else
			add.diagnostics(diagnostic.reads$cyp6[[this.sample]][[d]], this.col = discordant.colours[d], yrange = y.ranges[[d]])
	}
	these.cnv.states <- compact.hmm$cyp6[[this.sample]]$CNV[start.index : end.index]
	lines(compact.hmm$cyp6[[this.sample]]$Position[start.index : end.index], these.cnv.states , col = 2, lwd = 2)
	# Add the genes to the plot
	these.gene.colours <- gene.colours[1:nrow(gene.coords$cyp6)]
	abline(v = unlist(gene.coords$cyp6[, c('start', 'end')]), col = rep(these.gene.colours, 2))
	text(apply(gene.coords$cyp6[, c('start', 'end')], 1, mean), 9, rownames(gene.coords$cyp6), srt=90, adj = 0, col = these.gene.colours)
}

plot.all.cyp6 <- function(list.of.samples, matrix.of.read.dups = read.based.cnvs$cyp6, diagnostics = c('FA','SS','FM','XC','BP'), start.pos = plotting.ranges$cyp6[1], end.pos = plotting.ranges$cyp6[2]){
	x11()
	x.midpoint <- mean(c(start.pos, end.pos))
	i <- 1
	while(1){
		if (i < 1)
			i <- length(list.of.samples)
		this.sample <- list.of.samples[i]
		if (is.null(matrix.of.read.dups))
			plot.cyp6(this.sample, NULL, diagnostics, start.pos, end.pos)
		else
			plot.cyp6(this.sample, matrix.of.read.dups[this.sample,], diagnostics, start.pos, end.pos)
		x <- locator(1)$x
		if (x <= x.midpoint)
			i <- i-1
		else
			i <- i+1
	}
}

save.image('target_regions_analysis_extra_alleles.Rdata')
