# ==================================================================== #
# TITLE                                                                #
# Antimicrobial Resistance (AMR) Analysis                              #
#                                                                      #
# SOURCE                                                               #
# https://gitlab.com/msberends/AMR                                     #
#                                                                      #
# LICENCE                                                              #
# (c) 2019 Berends MS (m.s.berends@umcg.nl), Luz CF (c.f.luz@umcg.nl)  #
#                                                                      #
# This R package is free software; you can freely use and distribute   #
# it for both personal and commercial purposes under the terms of the  #
# GNU General Public License version 2.0 (GNU GPL-2), as published by  #
# the Free Software Foundation.                                        #
#                                                                      #
# This R package was created for academic research and was publicly    #
# released in the hope that it will be useful, but it comes WITHOUT    #
# ANY WARRANTY OR LIABILITY.                                           #
# Visit our website for more info: https://msberends.gitab.io/AMR.     #
# ==================================================================== #

#' @importFrom data.table as.data.table setkey
.onLoad <- function(libname, pkgname) {
  # get new functions not available in older versions of R
  backports::import(pkgname)

  # register data
  if (!all(c("microorganismsDT", "microorganisms.oldDT") %in% ls(envir = asNamespace("AMR")))) {

    microorganisms.oldDT <- as.data.table(AMR::microorganisms.old)
    setkey(microorganisms.oldDT, col_id, fullname)

    assign(x = "microorganismsDT",
           value = make_DT(),
           envir = asNamespace("AMR"))

    assign(x = "microorganisms.oldDT",
           value = microorganisms.oldDT,
           envir = asNamespace("AMR"))

    assign(x = "mo_codes_v0.5.0",
           value = make_trans_tbl(),
           envir = asNamespace("AMR"))
  }
}

#' @importFrom dplyr mutate case_when
#' @importFrom data.table as.data.table setkey
make_DT <- function() {
  microorganismsDT <- AMR::microorganisms %>%
    mutate(prevalence = case_when(
      class == "Gammaproteobacteria"
      | genus %in% c("Enterococcus", "Staphylococcus", "Streptococcus")
      ~ 1,
      phylum %in% c("Proteobacteria",
                    "Firmicutes",
                    "Actinobacteria",
                    "Sarcomastigophora")
      | genus %in% c("Aspergillus",
                     "Bacteroides",
                     "Candida",
                     "Capnocytophaga",
                     "Chryseobacterium",
                     "Cryptococcus",
                     "Elisabethkingia",
                     "Flavobacterium",
                     "Fusobacterium",
                     "Giardia",
                     "Leptotrichia",
                     "Mycoplasma",
                     "Prevotella",
                     "Rhodotorula",
                     "Treponema",
                     "Trichophyton")
      ~ 2,
      TRUE ~ 3
    )) %>%
    as.data.table()
  setkey(microorganismsDT,
         kingdom,
         prevalence,
         fullname)
  microorganismsDT
}

make_trans_tbl <- function() {
  # conversion of old MO codes from v0.5.0 (ITIS) to later versions (Catalogue of Life)
  c(B_ACHRMB = "B_ACHRM", B_ANNMA = "B_ACTNS", B_ACLLS = "B_ALCYC",
    B_AHNGM = "B_ARCHN", B_ARMTM = "B_ARMTMN", B_ARTHRS = "B_ARTHR",
    B_APHLS = "B_AZRHZP", B_BRCHA = "B_BRCHY", B_BCTRM = "B_BRVBCT",
    B_CLRBCT = "B_CLRBC", B_CTRDM = "B_CLSTR", B_CPRMM = "B_CYLND",
    B_DLCLN = "B_DPLCL", B_DMCLM = "B_DSLFT", B_DSLFVB = "B_DSLFV",
    B_FCTRM = "B_FSBCT", B_GNRLA = "B_GRDNR", B_HNRBM = "B_HLNRB",
    B_HPHGA = "B_HNPHGA", B_HCCCS = "B_HYDRC", B_MCRCLS = "B_MCRCL",
    B_MTHYLS = "B_MLSMA", B_MARCLS = "B_MRCLS", B_MGCLS = "B_MSTGC",
    B_MCLLA = "B_MTHYLC", B_MYCPLS = "B_MYCPL", B_NBCTR = "B_NTRBC",
    B_OCLLS = "B_OCNBC", B_PTHRX = "B_PLNKT", B_PCCCS = "B_PRCHL",
    B_PSPHN = "B_PRPHY", B_PDMNS = "B_PSDMN", B_SCCHRP = "B_SCCHR",
    B_SRBCTR = "B_SHRBCTR", B_STRPTC = "B_STRPT", B_SHMNS = "B_SYNTR",
    B_TRBCTR = "B_THRMN", P_ALBMN = "C_ABMNA", F_ACHLY = "C_ACHLY",
    P_ACINT = "C_ACINT", P_ARTCL = "C_ACLNA", P_ACRVL = "C_ACRVL",
    P_ADRCT = "C_ADRCT", P_AMPHS = "C_AHSRS", F_ALBUG = "C_ALBUG",
    P_ALCNT = "C_ALCNT", P_ALFRD = "C_ALFRD", P_ALLGR = "C_ALLGR",
    P_AMPHL = "C_ALPTS", F_ALTHR = "C_ALTHR", P_AMLLA = "C_AMLLA",
    P_ANMLN = "C_AMLNA", P_AMMBC = "C_AMMBC", P_AMMDS = "C_AMMDS",
    P_AMMLG = "C_AMMLG", P_AMMMR = "C_AMMMR", P_AMMMS = "C_AMMMS",
    P_AMMON = "C_AMMON", P_AMMSC = "C_AMMSC", P_AMMSP = "C_AMMSP",
    P_AMMST = "C_AMMST", P_AMMTM = "C_AMMTM", F_AMYCS = "C_AMYCS",
    P_ANARM = "C_ANARM", P_ANGLD = "C_ANGLD", P_ANGLG = "C_ANGLG",
    P_ANNLC = "C_ANNLC", F_ANSLP = "C_ANSLP", F_APDCH = "C_APDCH",
    F_APHND = "C_APHND", F_APLNC = "C_APLNC", F_AQLND = "C_AQLND",
    P_ARCHS = "C_ARCHAS", P_ASTRN = "C_ARNNN", P_ARNPR = "C_ARNPR",
    F_ARSPR = "C_ARSPR", P_ARTST = "C_ARTSTR", P_AMPHC = "C_ARYNA",
    P_ASCHM = "C_ASCHM", P_ASPDS = "C_ASPDS", P_ASTCL = "C_ASTCL",
    P_ASTRG = "C_ASTRGR", P_ASTRM = "C_ASTRMM", P_ASTRR = "C_ASTRR",
    P_ASTRT = "C_ASTRTR", F_ATKNS = "C_ATKNS", F_AYLLA = "C_AYLLA",
    P_BAGGN = "C_BAGGN", P_BCCLL = "C_BCCLL", P_BDLLD = "C_BDLLD",
    P_BGNRN = "C_BGNRN", P_BLCLN = "C_BLCLN", P_BLMND = "C_BLMND",
    P_BLMNL = "C_BLMNL", P_BLPHR = "C_BLPHR", P_BLVNT = "C_BLVNT",
    P_BOLVN = "C_BOLVN", P_BORLS = "C_BORLS", P_BRNNM = "C_BRNNM",
    P_BRSLN = "C_BRSLN", P_BRSRD = "C_BRSRD", F_BRVLG = "C_BRVLG",
    F_BNLLA = "C_BRVLGN", P_BSCCM = "C_BSCCM", F_BSDPH = "C_BSDPH",
    P_BTHYS = "C_BTHYS", P_BTLLN = "C_BTLLN", P_BULMN = "C_BULMN",
    P_CCLDM = "C_CCLDM", P_CDNLL = "C_CDNLL", P_CLPSS = "C_CDNLLP",
    P_CHLDN = "C_CHLDNL", P_CHLST = "C_CHLST", P_CHNLM = "C_CHNLM",
    P_CHRYS = "C_CHRYSL", P_CHTSP = "C_CHTSP", P_CBCDS = "C_CIBCDS",
    P_CLCRN = "C_CLCRN", P_CLMNA = "C_CLMNA", P_CLPDM = "C_CLPDM",
    P_CLPHR = "C_CLPHRY", P_CLVLN = "C_CLVLN", P_CMPNL = "C_CMPNL",
    P_CNCRS = "C_CNCRS", P_CNTCH = "C_CNTCH", F_CNTRM = "C_CNTRMY",
    P_COLPD = "C_COLPD", P_COLPS = "C_COLPS", P_CPRDS = "C_CPRDS",
    P_CRNSP = "C_CPRMA", P_CRBNL = "C_CRBNL", P_CRBRB = "C_CRBRB",
    P_CRBRG = "C_CRBRG", P_CRBRS = "C_CRBRS", P_CRCHS = "C_CRCHS",
    P_CRCLC = "C_CRCLC", P_CRNLC = "C_CRNLC", P_CRNTH = "C_CRNTH",
    P_CRPNT = "C_CRPNT", P_CRSTG = "C_CRSTG", P_CRTHN = "C_CRTHN",
    P_CRTRN = "C_CRTRN", P_CYMBL = "C_CRTTA", P_CRYPT = "C_CRYPT",
    P_CSHMN = "C_CSHMNL", P_CSSDL = "C_CSSDL", P_CLNDS = "C_CSSDLN",
    P_CHRNA = "C_CTHRN", P_CTPSS = "C_CTPSS", P_CUNLN = "C_CUNLN",
    P_CYLND = "C_CVLNA", P_CYCLC = "C_CYCLCB", P_CDNTA = "C_CYCLD",
    P_CYCLG = "C_CYCLG", P_CYCLM = "C_CYCLM", P_CYRTL = "C_CYRTL",
    P_CYSTM = "C_CYSTM", P_DCHLM = "C_DCHLM", P_DCRBS = "C_DCRBS",
    P_DCTYC = "C_DCTYC", P_DIDNM = "C_DIDNM", P_DLPTS = "C_DLPTS",
    P_DNTLN = "C_DNTLN", P_DNTST = "C_DNTST", P_DORTH = "C_DORTH",
    P_DCTYP = "C_DPHMS", F_DPLCY = "C_DPLCY", P_DNDRT = "C_DRTNA",
    P_DSCMM = "C_DSCMM", P_DSCRB = "C_DSCRB", P_DSCRN = "C_DSCRN",
    P_DSCSP = "C_DSCSP", P_DSNBR = "C_DSNBR", P_DYCBC = "C_DYCBC",
    F_DCTYC = "C_DYCHS", F_ECTRG = "C_ECTRG", B_EDWRD = "C_EDWRD",
    P_EGGRL = "C_EGGRL", P_EHLYS = "C_EHLYS", P_EHRNB = "C_EHRNB",
    P_ELPHD = "C_ELPHD", P_ENCHL = "C_ELYDM", P_EPHDM = "C_EPHDM",
    P_EPLTS = "C_EPLTS", P_EPLXL = "C_EPLXL", P_EPNDL = "C_EPNDL",
    P_EPNDS = "C_EPNDS", P_ENLLA = "C_EPSTM", P_EPSTY = "C_EPSTY",
    F_ERYCH = "C_ERYCH", F_ESMDM = "C_ESMDM", P_ESSYR = "C_ESSYR",
    P_FSCHR = "C_FHRNA", P_FLRLS = "C_FLRLS", P_FLNTN = "C_FNTNA",
    P_FRNDC = "C_FRNDC", P_FRNTN = "C_FRNTN", P_FRSNK = "C_FRSNK",
    P_FNLLA = "C_FSCHRN", P_FSSRN = "C_FSSRN", P_FVCSS = "C_FVCSS",
    P_GDRYN = "C_GDRYN", F_GELGN = "C_GELGN", P_GERDA = "C_GERDA",
    P_GLACM = "C_GLACM", P_GLBBL = "C_GLBBL", P_GLBGR = "C_GLBGR",
    P_GLBLN = "C_GLBLN", P_GRTLA = "C_GLBRT", P_GLBTX = "C_GLBTX",
    P_GLLNA = "C_GLLNA", P_GLMSP = "C_GLMSP", P_GLNDL = "C_GLNDL",
    F_GNMCH = "C_GNMCH", P_GOSLL = "C_GOSLL", P_GRNDS = "C_GRNDS",
    P_GRNTA = "C_GRNTA", P_GLBRT = "C_GTLLA", P_GTTLN = "C_GTTLN",
    P_GVLNP = "C_GVLNP", P_GYPSN = "C_GYPSN", P_GYRDN = "C_GYRDN",
    P_HALTR = "C_HALTR", P_HANZW = "C_HANZW", P_HAURN = "C_HAURN",
    P_HELNN = "C_HELNN", P_HLPHR = "C_HHRYA", P_HLNTA = "C_HLNTA",
    F_HLPHT = "C_HLPHT", P_HLSTC = "C_HLSTC", P_HMSPH = "C_HMSPH",
    P_HMTRM = "C_HMTRM", P_HPKNS = "C_HPKNS", P_HPLPH = "C_HPLPH",
    P_HPPCR = "C_HPPCR", P_HNLLA = "C_HPPCRP", P_HRMSN = "C_HRMSN",
    P_HRNLL = "C_HRNLL", F_HRPCH = "C_HRPCH", P_HSTGR = "C_HSTGR",
    P_HSTTL = "C_HSTTL", P_HTRST = "C_HTGNA", P_HTRLL = "C_HTRLL",
    P_HTRPH = "C_HTRPH", F_HYPHC = "C_HYPHC", P_HYPRM = "C_HYPRM",
    P_INTRN = "C_INTRN", P_IRIDI = "C_IRIDI", P_ISLND = "C_ISLND",
    P_JCLLL = "C_JCLLL", P_KHLLL = "C_KHLLL", P_KRNPS = "C_KRNPS",
    P_KRRRL = "C_KRRRL", P_LABOE = "C_LABOE", P_LAGEN = "C_LAGEN",
    P_LBSLL = "C_LBSLL", F_LTHLA = "C_LBYRN", P_LCRYM = "C_LCRYM",
    P_LEMBS = "C_LEMBS", F_LGNDM = "C_LGNDM", P_LGNMM = "C_LGNMM",
    P_LGNPH = "C_LGNPHR", F_LGNSM = "C_LGNSM", P_LGYNP = "C_LGYNP",
    P_LITTB = "C_LITTB", P_LITUL = "C_LITUL", P_LMBDN = "C_LMBDN",
    P_LMRCK = "C_LMRCK", F_LBYRN = "C_LMYXA", P_LNGLN = "C_LNGLN",
    P_LNTCL = "C_LNTCL", P_LOXDS = "C_LOXDS", F_LPTLG = "C_LPTLG",
    F_LNLLA = "C_LPTLGN", F_LPTMT = "C_LPTMT", P_LRYNG = "C_LRYNG",
    P_LTCRN = "C_LTCRN", P_LTHPL = "C_LTHPL", P_LTNTS = "C_LTNTS",
    F_LTRST = "C_LTRST", P_LXPHY = "C_LXPHY", P_MCRTH = "C_MCRTH",
    P_MELNS = "C_MELNS", P_MSDNM = "C_MESDNM", P_METPS = "C_METPS",
    P_MIMSN = "C_MIMSN", P_MINCN = "C_MINCN", P_MLLNL = "C_MLLNL",
    P_MLMMN = "C_MLMMN", F_MNDNL = "C_MNDNL", P_MNLYS = "C_MNLYS",
    P_MNPSS = "C_MNPSS", P_MRGNL = "C_MRGNL", P_MRGNP = "C_MRGNP",
    P_MRSPL = "C_MRSPL", P_MRTNT = "C_MRTNT", P_MSSLN = "C_MSSLN",
    P_MSSSS = "C_MSSSS", P_MTCNT = "C_MTCNT", P_MYCHS = "C_MYCHS",
    P_MYSCH = "C_MYSCH", F_MYZCY = "C_MYZCY", P_NASSL = "C_NASSL",
    P_NBCLN = "C_NBCLN", P_NBCLR = "C_NBCLR", P_NCNRB = "C_NCNRB",
    P_NDBCL = "C_NDBCL", P_NRLLA = "C_NDBCLR", P_NMMLC = "C_NMMLC",
    F_NMTPH = "C_NMTPH", P_NNNLL = "C_NNNLL", P_NODSR = "C_NODSR",
    P_NONIN = "C_NONIN", P_NOURI = "C_NOURI", P_OCLNA = "C_OCLNA",
    P_OGLNA = "C_OGLNA", P_OPHTH = "C_OLMDM", F_OLPDP = "C_OLPDP",
    P_ONYCH = "C_OMPSS", P_OOLIN = "C_OOLIN", P_OPRCL = "C_OPRCL",
    P_ORBLN = "C_ORBLN", F_ORCAD = "C_ORCAD", P_ORDRS = "C_ORDRS",
    P_OPHRY = "C_ORYDM", P_OSNGL = "C_OSNGL", P_OXYTR = "C_OXYTR",
    P_PARRN = "C_PARRN", P_PATRS = "C_PATRS", P_PAVNN = "C_PAVNN",
    P_PTYCH = "C_PCYLS", P_PDPHR = "C_PDPHR", P_PELSN = "C_PELSN",
    F_PHGMY = "C_PHGMY", F_PSDSP = "C_PHRTA", P_PHRYG = "C_PHRYG",
    P_PHYSL = "C_PHYSL", F_PHYTP = "C_PHYTP", P_PLACS = "C_PLACS",
    P_PLCPS = "C_PLCPS", P_PLCPSL = "C_PLCPSL", P_PLCTN = "C_PLCTN",
    P_PLGPH = "C_PLGPH", B_PLGTH = "C_PLGTH", P_PLMRN = "C_PLMRN",
    P_PLNCT = "C_PLNCT", P_PLNDSC = "C_PLNDSC", P_PLNGY = "C_PLNGY",
    P_PLNRBL = "C_PLNLLA", P_PLNLN = "C_PLNLN", P_PLNLR = "C_PLNLR",
    P_PLNRB = "C_PLNRB", P_PLNSP = "C_PLNSPR", P_PLRNM = "C_PLRNM",
    P_PLRST = "C_PLRST", P_PLRTR = "C_PLRTR", F_PLSMD = "C_PLSMD",
    P_PLTYC = "C_PLTYC", P_PSDBL = "C_PLVNA", P_PLYMR = "C_PLYMR",
    P_PLTYN = "C_PNMTM", P_PNRPL = "C_PNRPL", F_PNTSM = "C_PNTSM",
    P_PRCNT = "C_PRCNT", P_PRFSS = "C_PRFSS", P_PRMCM = "C_PRMCUM",
    F_PRNSP = "C_PRNSP", P_PRPND = "C_PRPND", P_PRPYX = "C_PRPYX",
    P_PRRDN = "C_PRRDN", P_PSDDF = "C_PSDDF", P_PSDMC = "C_PSDMC",
    P_PSDND = "C_PSDND", P_PSDNN = "C_PSDNN", P_PSDPL = "C_PSDPLY",
    P_PSMMS = "C_PSMMS", P_PTLLN = "C_PTLLN", P_PTLLND = "C_PTLLND",
    F_PTRSN = "C_PTRSN", P_PULLN = "C_PULLN", P_PUTLN = "C_PUTLN",
    P_PRTTR = "C_PYMNA", P_PYRGL = "C_PYRGL", P_PYRGO = "C_PYRGO",
    P_PYRLN = "C_PYRLN", F_PYTHM = "C_PYTHIM", F_PYTHL = "C_PYTHL",
    P_PYXCL = "C_PYXCL", P_QNQLC = "C_QNQLC", P_RAMLN = "C_RAMLN",
    P_RBRTN = "C_RBRTN", P_RCRVD = "C_RCRVD", P_RCTBL = "C_RCTBL",
    P_RCTCB = "C_RCTCB", P_RCTGL = "C_RCTGL", P_RCTVG = "C_RCTVG",
    P_RDGDR = "C_RDGDR", P_REMNC = "C_REMNC", P_REPHX = "C_REPHX",
    P_RHBDM = "C_RHBDMM", F_RHBDS = "C_RHBDSP", P_RHPDD = "C_RHPDD",
    F_RHPDM = "C_RHPDM", F_RHZDMY = "C_RHZDM", P_RHZMM = "C_RHZMM",
    P_RIVRN = "C_RIVRN", P_ROSLN = "C_ROSLN", P_ROTAL = "C_ROTAL",
    P_RPHDP = "C_RPHDP", P_RPRTN = "C_RPRTN", P_RSSLL = "C_RSSLL",
    P_RTLMM = "C_RTLMM", P_RTYLA = "C_RTYLA", P_RUGID = "C_RUGID",
    F_RZLLP = "C_RZLLP", P_SAGRN = "C_SAGRN", P_SCCMM = "C_SCCMM",
    P_SCCRH = "C_SCCRH", P_SCHLM = "C_SCHLM", F_SCLRS = "C_SCLRS",
    P_SCTLR = "C_SCTLR", P_SEBRK = "C_SEBRK", P_SGMLN = "C_SGMLN",
    P_SGMLP = "C_SGMLP", P_SGMMR = "C_SGMMR", P_SGMVR = "C_SGMVR",
    F_SMMRS = "C_SMMRS", P_SNNDS = "C_SNNDS", P_SORTS = "C_SORTS",
    P_SPHGN = "C_SPHGN", P_SPHNN = "C_SPHNN", P_SNLLA = "C_SPHNNL",
    P_SPHTR = "C_SPHTR", P_SPHTX = "C_SPHTX", P_SPHVG = "C_SPHVG",
    P_SPRDT = "C_SPRDT", P_SPRLC = "C_SPRLC", F_SPRLG = "C_SPRLG",
    P_SPRLL = "C_SPRLL", F_SPRMY = "C_SPRMY", P_SPRPL = "C_SPRPL",
    P_SPRSG = "C_SPRSG", P_SPRST = "C_SPRST", P_SPHNP = "C_SPRTA",
    P_SPRZN = "C_SPRZN", P_SPHRG = "C_SPSNA", P_STHDM = "C_SPTHD",
    P_SRCNR = "C_SRCNR", F_SRLPD = "C_SRLPD", F_SPNGS = "C_SSPRA",
    F_STEIN = "C_STEIN", P_SPTHD = "C_STHDDS", P_STHRP = "C_STHRP",
    P_STNFR = "C_STNFR", P_STNSM = "C_STNSM", P_STNTR = "C_STNTR",
    P_STRBL = "C_STRBL", P_STRMB = "C_STRMB", P_STTSN = "C_STTSN",
    P_STYLN = "C_SYCHA", F_SCHZC = "C_SYTRM", P_TBNLL = "C_TBNLL",
    P_TRCHL = "C_TCHLS", P_TCHNT = "C_TCHNT", P_THRCL = "C_THRCL",
    P_THRMM = "C_THRMM", P_TIARN = "C_TIARN", P_TKPHR = "C_TKPHR",
    P_TLNMA = "C_TLNMA", P_TLYPM = "C_TLYPM", P_TMNDS = "C_TMNDS",
    P_TMNTA = "C_TMNTA", P_TNTNN = "C_TNNDM", P_TTNNS = "C_TNTNN",
    P_TNPSS = "C_TNTNNP", P_TONTN = "C_TONTN", P_TOSAI = "C_TOSAI",
    P_TPHTR = "C_TPHTR", P_TRCHH = "C_TRCHH", P_TRPHS = "C_TRCHLR",
    P_TMMNA = "C_TRCHM", P_TRCHS = "C_TRCHSP", P_TRFRN = "C_TRFRN",
    P_TRLCL = "C_TRLCL", P_TRTXL = "C_TRTXL", P_TRTXS = "C_TRTXS",
    P_TTRHY = "C_TTRHY", F_TTRMY = "C_TTRMY", P_TXTLR = "C_TXTLR",
    F_THRST = "C_TYTRM", P_URLPT = "C_ULPTS", P_UNGLT = "C_UNGLT",
    P_URCNT = "C_URCNT", P_URONM = "C_URONM", P_UROSM = "C_UROSM",
    P_URTRC = "C_URTRC", P_URSTY = "C_UTYLA", P_UVGRN = "C_UVGRN",
    P_VLVLN = "C_VALVLN", P_VGNLN = "C_VGNLN", P_VGNLNP = "C_VGNLNP",
    P_VLNRA = "C_VLVLN", P_VGNCL = "C_VNCLA", P_VRGLN = "C_VRGLN",
    P_VRGLNP = "C_VRGLNP", P_VRTCL = "C_VRTCL", P_WBBNL = "C_WBBNL",
    P_WEBBN = "C_WEBBN", P_WSNRL = "C_WSNRL", P_ZTHMN = "C_ZHMNM",
    B_ZOOGL = "C_ZOOGL", F_DDSCS = "F_DPDSC", F_SCCHR = "F_SMYCS",
    P_AMTRN = "P_ACNTH", F_AMBDM = "P_AMBDM", F_ARCYR = "P_ARCYR",
    F_BADHM = "P_BADHM", F_BDHMP = "P_BDHMP", F_BRBYL = "P_BRBYL",
    F_BRFLD = "P_BRFLD", F_CLMYX = "P_CLMYX", F_CLSTD = "P_CLSTD",
    F_CMTRC = "P_CMTRC", F_CRBRR = "P_CRBRR", F_CRTMY = "P_CRTMY",
    F_CRTRM = "P_CRTRM", F_DCTYD = "P_DCTYD", F_DDYMM = "P_DDYMM",
    F_DIACH = "P_DIACH", F_DIANM = "P_DIANM", F_DIDRM = "P_DIDRM",
    F_ELMYX = "P_ELMYX", F_ESTLM = "P_ESTLM", F_FULIG = "P_FULIG",
    F_HMTRC = "P_HMTRC", F_LCRPS = "P_LCRPS", F_LICEA = "P_LICEA",
    F_LMPRD = "P_LMPRD", F_LPTDR = "P_LPTDR", F_LSTRL = "P_LSTRL",
    F_LYCGL = "P_LYCGL", F_MCBRD = "P_MCBRD", F_MNKTL = "P_MNKTL",
    F_MTTRC = "P_MTTRC", F_MUCLG = "P_MUCLG", F_PHYSR = "P_PHYSR",
    F_PRCHN = "P_PRCHN", F_PRMBD = "P_PRMBD", F_PRTPH = "P_PRTPH",
    F_PSRNA = "P_PSRNA", F_PYSRM = "P_PYSRM", F_RTCLR = "P_RTCLR",
    F_STMNT = "P_STMNT", F_SYMPH = "P_SYMPH", F_TRBRK = "P_TRBRK",
    F_TRICH = "P_TRICH", F_TUBFR = "P_TUBFR")
}
