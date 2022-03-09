palette_color = c("chocolate1", "cyan", "gold", "aquamarine", "deepskyblue", 
                  "cyan4", "darkblue", "darkolivegreen1", "darkorchid1",
                  "firebrick1", "firebrick4", "purple", "darksalmon","darkslategray", "darkmagenta")
#define markers
hla_markers = c("HLA-DQB2","HLA-DOA","HLA-DRA","HLA-DMB","HLA-DPB1","HLA-DQA2")
mrp_markers = c("MRPL38","MRPS10","MRPL12","MRPL14","MRPL47","MRPS18A","MRPS24","MRPS36","MRPS5")
nduf_markers = c("NDUFAB1", "NDUFAF3", "NDUFA5", "NDUFS3","ATP5MC1","ATP5MC2","COA3","UQCRC1")
cancer_markers = c("CA9","NDUFA4L2","SLC17A3")
EMT_markers = c("MT2A", "SERPINE1","TM4SF1", "BIRC3", "C3", "CAV1", "TGM2", "MMP7",
                "ANXA2", "LOX", "FSTL3", "LGALS1")
PT_markers = c("NAT8", "SLC3A1", "GPX3", "CYB5A", "SLC25A5", "GSTA1")

palette_color = c("firebrick1", "darkorange", "gold", 
                  "chartreuse", "blue", "mediumorchid1", 
                  "cyan", "lightpink", "sienna4")

#new signatures from Maria
# Complex_V = c("ATP5F1A","ATP5F1C","ATP5F1D","ATP5F1E","ATP5PB",
#               "ATP5MC1","ATP5MC2","ATP5MC3","ATP5PD","ATP5ME","ATP5PF",
#               "ATP5MF","ATP5MG","ATP5O")
Complex_V = c("ATP5A1(ATP5F1A)","ATP5C1(ATP5F1C)","ATP5D(ATP5F1D)","ATP5E(ATP5F1E)",
              "ATP5F1(ATP5PB)","ATP5G1(ATP5MC1)","ATP5G2(ATP5MC2)","ATP5G3(ATP5MC3)",
              "ATP5H(ATP5PD)","ATP5I(ATP5ME)","ATP5J(ATP5PF)","ATP5J2(ATP5MF)",
              "ATP5L(ATP5MG)","ATP5O")
Complex_V = sapply(Complex_V, function(x){
  a = stringr::str_split_fixed(x, "\\(", 2)[1,1]
  return(a)
})
Complex_V = as.character(Complex_V)
Complex_III = c("COQ9","UQCR10","UQCR11","UQCRB","UQCRC1","UQCRC2","UQCRFS1","UQCRQ")
Complex_IV = c("COX10","COX11","COX15","COX17","COX4I1","COX4I2","COX5A","COX5B","COX6A1",
               "COX6A2","COX6B1","COX6B2","COX6C","COX7A1","COX7A2","COX7A2L",
               "COX7B","COX7C","COX8A")
Complex_I = c("NDUFA1","NDUFA10","NDUFA11","NDUFA12","NDUFA13","NDUFA2","NDUFA3",
              "NDUFA4","NDUFA4L2","NDUFA5","NDUFA6","NDUFA7","NDUFA8","NDUFA9",
              "NDUFAB1","NDUFB1","NDUFB10","NDUFB11","NDUFB2","NDUFB3","NDUFB4",
              "NDUFB5","NDUFB6","NDUFB7","NDUFB8","NDUFB9","NDUFC1","NDUFC2",
              "NDUFS1","NDUFS2","NDUFS3","NDUFS4","NDUFS5","NDUFS6","NDUFS7",
              "NDUFS8","NDUFV1","NDUFV2","NDUFV3")
Complex_II = c("SDHA","SDHB","SDHC","SDHD","SUCLA2","SUCLG1","SUCLG2")
glycolyticAndTCA = c("ACO1","ACO2","ALDOA","ALDOB","ALDOC","BPGM","CS","CYB5R1",
                     "CYC1","DLAT","DLD","DLST","ECI1","ENO1","ENO2","ENO3","ETFB",
                     "FH","GAPDH","GCK","GOT1","GOT2","GPD1","GPD2","GPI","HK1",
                     "HK2","HK3","IDH1","IDH2","IDH3A","IDH3B","IDH3G","LDHA",
                     "LDHB","MDH1","MDH2","ME1","ME2","ME3","MICU1","MPC1","MPC2",
                     "MRPL49","OGDH","OGDHL","PC","PCK1","PCK2","PFKL","PFKM","PFKP",
                     "PGAM1","PGAM2","PGK1","PGK2","PKLR","PPA1","PPA2","SCO1",
                     "SLC25A11","SLC25A13","SLC2A1","SLC2A2","SLC2A3","SLC2A4","TPI1")
MRP_positive = c("MRPL16","MRPL49","MRPS35","MRPL46","MRPL10","MRPL39","MRPL19","MRPL33",
                 "MRPS9","MRPL44","MRPL1","MRPS18C","MRPS30","MRPS36","MRPS18A","MRPS18B",
                 "MRPS28","MRPL50","MRPS2")
MRP_negative = c("MRPL12","MRPL14","MRPL15","MRPL17","MRPL18","MRPL23","MRPL24","MRPL34",
                 "MRPL36","MRPL38","MRPL40","MRPL47","MRPL53","MRPL55","MRPS6","MRPS11",
                 "MRPS15","MRPS16","MRPS25")
MT = c("MT-CO1","MT-CO2","MT-CO3","MT-ATP6","MT-ATP8",
       "MT-ND1","MT-ND3","MT-ND4","MT-ND4L")
new_signature = c("OSBPL6","NTM","RFTN1","MAPRE2","CXADR","DYSF","TGFA","EPHA6",
                  "GJB2","FBXL7","DTX4","UNC5B","PKP2","PRSS12","SLC16A2","EPHB1",
                  "MAN1A1","ADM","EDIL3","LGI2","DHRS3","PCDH7","PTPRU","CHSY3",
                  "APLP1","SEMA5B","ATRNL1","UST","GAS6","SLC39A8","SLIT2","EIF4E3",
                  "TESC","PMP22","RBFOX3","WNT5A","FAM171B","AMIGO2","FGF9","PANX2",
                  "TLL2","GPC6","NOG","MAPK11","HLX","RNF182","PIK3AP1","PCDH1","ZIC2",
                  "RYR2","ADAMTS7","MATN3","WT1","HLF","TMEM121","ADORA2B","RTN4RL2",
                  "TMOD2","CXCL16","GJA3","LMO1","SLC25A21","SOBP","RAB20","B4GALNT1",
                  "PTHLH","COL15A1","KIAA1217","KIAA1324","RAB6B","CLDN23","MAL","ITPKA",
                  "ATP8A2","MMP16","RGL1","CACNA1D","TMEM163","HCN4","ZFPM2","KLHL14",
                  "NEGR1","HOXC11","SLFN11","PALM3","ITGB4","ABCA5","SCAMP5","TPPP",
                  "CITED1","CACNB4","SGPP2","NKX3-2","CD83","SPHK1")
CU = c("SOD1","CP","MT-CO1","MT-CO2","DBH","MOXD1","MOXD2P","PAM","SOD3","AFP",
       "ALB","CUTA","S100A5","SPARC","ATP7A","ATP7B","SLC31A1","SLC25A3","SLC31A2",
       "ENOX1","ENOX2","LTF","Cox17","ATOX1","CCS","COMMD1","CUTC","SCO1","SCO2",
       "Cox19","CHCHD7","GLRX","Cox11","MYO5B","COA1","AOC1","AOC3","LOXL3","LOXL4",
       "TYR","TYRP1","F5","MT1A","MT1B","MT1E","MT1F","MT1G","MT1H","MT1M","MT1X",
       "MT2A","MT3")
TCA = c("ACO1","ACO2","DLAT","DLD","DLST","SDHA","SDHB","SDHC","SDHD","SUCLA2",
        "SUCLG1","SUCLG2","FH","IDH2","IDH3A","IDH3B","IDH3G","MDH2","CS","OGDH",
        "OGDHL","PDHX")
MAS = c("MDH1","MDH2","GOT1","GOT2","SLC25A11","SLC25A12","SLC25A13")
HIF1A = c("ADM","ADRA1B","AK3","ALDOA","ALDOC","BNIP3","CA9","CDKN1A","CITED2",
          "CP","EDN1","ENO1","EPO","FLT1","GAPDH","HK1","HK2","HMOX1","IGF2",
          "IGFBP1","IGFBP2","IGFBP3","LDHA","NOS2","P4HA2","PFKL","PGK1","PKM",
          "SERPINE1","SLC2A1","SLC2A3","TF","TFRC","TGFB3","TPI1","VEGFA")
EMT = c("ABI3BP","ACTA2","ADAM12","ANPEP","APLP1","AREG","BASP1","BDNF","BGN",
        "BMP1","CADM1","CALD1","CALU","CAP2","CAPG","CCN1","CCN2","CD44","CD59","CDH11",
        "CDH2","CDH6","COL11A1","COL12A1","COL16A1","COL1A1","COL1A2","COL3A1","COL4A1",
        "COL4A2","COL5A1","COL5A2","COL5A3","COL6A2","COL6A3","COL7A1","COL8A2",
        "COLGALT1","COMP","COPA","CRLF1","CTHRC1","CXCL1","CXCL12","CXCL6","CXCL8",
        "DAB2","DCN","DKK1","DPYSL3","DST","ECM1","ECM2","EDIL3","EFEMP2","ELN",
        "EMP3","ENO2","FAP","FAS","FBLN1","FBLN2","FBLN5","FBN1","FBN2","FERMT2",
        "FGF2","FLNA","FMOD","FN1","FOXC2","FSTL1","FSTL3","FUCA1","FZD8","GADD45A",
        "GADD45B","GAS1","GEM","GJA1","GLIPR1","GPC1","GPX7","GREM1","HTRA1","ID2",
        "IGFBP2","IGFBP3","IGFBP4","IL15","IL32","IL6","INHBA","ITGA2","ITGA5","ITGAV",
        "ITGB1","ITGB3","ITGB5","JUN","LAMA1","LAMA2","LAMA3","LAMC1","LAMC2","LGALS1",
        "LOX","LOXL1","LOXL2","LRP1","LRRC15","LUM","MAGEE1","MATN2","MATN3","MCM7",
        "MEST","MFAP5","MGP","MMP1","MMP14","MMP2","MMP3","MSX1","MXRA5","MYL9","MYLK",
        "NID2","NNMT","NOTCH2","NT5E","NTM","OXTR","P3H1","PCOLCE","PCOLCE2","PDGFRB",
        "PDLIM4","PFN2","PLAUR","PLOD1","PLOD2","PLOD3","PMEPA1","PMP22","POSTN","PPIB",
        "PRRX1","PRSS2","PTHLH","PTX3","PVR","QSOX1","RGS4","RHOB","SAT1","SCG2","SDC1",
        "SDC4","SERPINE1","SERPINE2","SERPINH1","SFRP1","SFRP4","SGCB","SGCD","SGCG",
        "SLC6A8","SLIT2","SLIT3","SNAI2","SNTB1","SPARC","SPOCK1","SPP1","TAGLN","TFPI2",
        "TGFB1","TGFBI","TGFBR3","TGM2","THBS1","THBS2","THY1","TIMP1","TIMP3","TNC",
        "TNFAIP3","TNFRSF11B","TNFRSF12A","TPM1","TPM2","TPM4","VCAM1","VCAN","VEGFA",
        "VEGFC","VIM","WIPF1","WNT5A")
HIF1A_2A = c("ADM","AHNAK2","AK4","AKAP12","ALDOC","ANG","ANGPTL4","ANKZF1","ATXN1",
             "BHLHE40","BNIP3","BNIP3L","CA9","CAV1","CAVIN1","CCN5","CCNG2","CITED2",
             "CSRP2","CXCR4","CYB5A","DPYSL2","DPYSL4","DSC2","DST","EFNA3","EGFR",
             "EGLN1","EGLN3","EGR1","EIF5A","ENO2","ERO1A","FAM13A","FAM162A","FOS",
             "FYN","GBE1","GJA1","GLRX","GPRC5A","GYS1","HEY1","HILPDA","HK2","HLA-DRB1",
             "IGFBP3","IGFBP5","ILVBL","INSIG2","ITPR1","KDM4B","KRT15","LIMCH1","LOX",
             "LOXL1","LOXL2","MAGED4B","MXI1","NDRG1","NOL3","OBSL1","P4HA1","P4HA2",
             "PAM","PDGFB","PDK1","PFKFB3","PGAP1","PGK1","PGM1","PLAC8","PLAUR",
             "PLIN2","PPFIA4","QSOX1","RASA4","RBPJ","RNASE4","RRAGD","S100A4",
             "SAMD4A","SCNN1B","SERPINE1","SFXN3","SH2B2","SH3GL3","SLC2A1",
             "SOX9","SPAG4","SPOCK1","SRD5A3","SRPX","STBD1","STC1","TMEM45A",
             "UPK1A","VEGFA","VEGFC","VLDLR","WSB1","YEATS2","ZNF292","ZNF395","ZNF654")
glycolysis = c("ALDOA","ALDOB","ALDOC","BPGM","ENO1","ENO2","ENO3","GPD1","SLC2A1","SLC2A2",
               "SLC2A3","SLC2A4","PFKL","PFKM","PFKP","PGAM1","PGAM2","PGK1","PGK2","PKLR",
               "PCK1","PCK2","TPI1","LDHA","LDHB","HK1","HK2","HK3","GAPDH")

