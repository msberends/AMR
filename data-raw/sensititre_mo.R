codes <- tibble::tribble(
  ~code, ~name,
  "A.DEF", "Abiotrophia defectiva",
  "ABI.SP", "Abiotrophia species",
  "ABS.SP", "Absidia species",
  "A.BRA", "Acholeplasma brassicae",
  "A.PAL", "Acholeplasma palmae",
  "ACO.SP", "Acholeplasma species",
  "A.DEN", "Achrom. xylo ss denitrificans",
  "A.XYL", "Achrom. xylo ss xylosoxidans",
  "A.PIE", "Achromobacter piechaudii",
  "ACH.SP", "Achromobacter species",
  "A.GRPB", "Achromobacter-like group B",
  "A.GRPE", "Achromobacter-like group E",
  "A.GRPF", "Achromobacter-like group F",
  "A.MUL", "Acidiphilium multivorum",
  "ACD.SP", "Acidiphilium species",
  "A.DEL", "Acidovorax delafieldii",
  "A.FAC", "Acidovorax facilis",
  "AVX.SP", "Acidovorax species",
  "A.TEM", "Acidovorax temperans",
  "A.BAU", "Acinetobacter baumannii",
  "A.CALC", "Acinetobacter calcoaceticus",
  "A.HAE", "Acinetobacter haemolyticus",
  "A.JOH", "Acinetobacter johnsonii",
  "A.JUN", "Acinetobacter junii",
  "A.LWO", "Acinetobacter lwoffi",
  "A.RADR", "Acinetobacter radioresistens",
  "A10.SP", "Acinetobacter sp. 10",
  "AC3.SP", "Acinetobacter sp. 3",
  "AC6.SP", "Acinetobacter sp. 6",
  "AC9.SP", "Acinetobacter sp. 9",
  "ACI.SP", "Acinetobacter species",
  "ACR.SP", "Acremonium species",
  "A.EQU", "Actino. equuli ss equuli",
  "A.EQUH", "Actino. equuli ss haemolyticus",
  "A.LIL", "Actinobac. lignieresii-like",
  "A.PLE", "Actinobac. pleuropneumoniae",
  "A.ACM", "Actinobacillus actinomycetemco",
  "A.CAP", "Actinobacillus capsulatus",
  "A.LIG", "Actinobacillus lignieresii",
  "A.SEM", "Actinobacillus seminis",
  "ACT.SP", "Actinobacillus species",
  "AC.SUI", "Actinobacillus suis",
  "A.URE", "Actinobacillus ureae",
  "AC/PA", "Actinobacillus/Pasteurella",
  "ACB.SP", "Actinobaculum species",
  "AM.SUI", "Actinobaculum suis",
  "A.LAT", "Actinomadura latina",
  "A.MAD", "Actinomadura madurae",
  "A.PEL", "Actinomadura pelletieri",
  "ACM.SP", "Actinomadura species",
  "A.BOV", "Actinomyces bovis",
  "A.HOR", "Actinomyces hordeovulneris",
  "A.ISR", "Actinomyces israelii",
  "A.MEY", "Actinomyces meyeri",
  "A.NAE", "Actinomyces naeslundii",
  "A.ODL", "Actinomyces odontolyticus",
  "AMY.SP", "Actinomyces species",
  "A.VIS", "Actinomyces viscosus",
  "A.IRA", "Actinopolyspora iraqiensis",
  "APS.SP", "Actinopolyspora species",
  "A.HYD", "Aero. hydrophila ss hydrophila",
  "A.SAL", "Aero.salmonicid ss salmonicida",
  "AACTIN", "Aerobic Actinomycetes(Inc Noc)",
  "AEC.SP", "Aerococcus species",
  "A.VIR", "Aerococcus viridans",
  "A.FAS", "Aeromicrobium fastidiosum",
  "AMI.SP", "Aeromicrobium species",
  "A.CAV", "Aeromonas caviae",
  "AHYDCM", "Aeromonas hydrophilia complex",
  "A.JAN", "Aeromonas jandaei",
  "A.SCH", "Aeromonas schubertii",
  "A.SOB", "Aeromonas sobria",
  "AEM.SP", "Aeromonas species",
  "A.TRO", "Aeromonas trota",
  "A.VER", "Aeromonas veronii",
  "AEMLSP", "Aeromonas-like species",
  "AGR.SP", "Agrobacterium species",
  "A.FAE", "Alcalig. faecalis ss faecalis",
  "ALC.SP", "Alcaligenes species",
  "ALT.SP", "Alternaria species",
  "ANAERO", "Anaerobic organism",
  "P.PRE", "Anaerococcus prevotii",
  "ANC.SP", "Anaerococcus species",
  "ANY", "ANY ORGANISM",
  "A.PYO", "Arcanobacterium pyogenes",
  "ARC.SP", "Arcanobacterium species",
  "A.CRY", "Arcobacter cryaerophilus",
  "A.NIT", "Arcobacter nitrofigilis",
  "ARB.SP", "Arcobacter species",
  "M.FUL", "Arthroderma fulvum",
  "ART.SP", "Arthroderma species",
  "A.CLA", "Aspergillus calvatus",
  "A.FLA", "Aspergillus flavus",
  "A.FUM", "Aspergillus fumigatus",
  "A.NID", "Aspergillus nidulans",
  "A.NIG", "Aspergillus niger",
  "ASP.SP", "Aspergillus species",
  "A.TER", "Aspergillus terreus",
  "B.ANT", "Bacillus anthracis",
  "B.BAD", "Bacillus badius",
  "B.CER", "Bacillus cereus",
  "B.LAE", "Bacillus laevolacticus",
  "B.LEN", "Bacillus lentus",
  "BAC.SP", "Bacillus species",
  "B.SUB", "Bacillus subtilis ss subtilis",
  "B.THU", "Bacillus thuringiensis",
  "B.CAP", "Bacteroides capillosus",
  "B.DIST", "Bacteroides distasonis",
  "B.EGG", "Bacteroides eggerthii",
  "B.FRA", "Bacteroides fragilis",
  "B.GRPA", "Bacteroides Group 3452A",
  "B.OVA", "Bacteroides ovatus",
  "BCT.SP", "Bacteroides species",
  "B.SPL", "Bacteroides splanchnicus",
  "B.THE", "Bacteroides thetaiotaomicron",
  "B.UNI", "Bacteroides uniformis",
  "B.URE", "Bacteroides ureolyticus",
  "B.VUL", "Bacteroides vulgatus",
  "B.ALP", "Balneatrix alpica",
  "BAL.SP", "Balneatrix species",
  "BER.SP", "Bergeyella species",
  "B.ZOO", "Bergeyella zoohelcum",
  "S.BHG", "Beta Hemolytic Strep Group",
  "B.TREH", "Bibersteinia trehalosi",
  "B.ADO", "Bifidobacterium adolescentis",
  "BIF.SP", "Bifidobacterium species",
  "B.DER", "Blastomyces dermatoides",
  "G.CAP", "Blastoschizomyces capitatus",
  "BLA.SP", "Blastoschizomyces species",
  "B.AVI", "Bordetella avium",
  "B.BRO", "Bordetella bronchiseptica",
  "B.HIN", "Bordetella hinzii",
  "B.HOL", "Bordetella holmesii",
  "B.PAR", "Bordetella parapertussis",
  "B.PER", "Bordetella pertussis",
  "BOR.SP", "Bordetella species",
  "B.TRE", "Bordetella trematum",
  "BOE.SP", "Borrelia species",
  "T.HYO", "Brachyspira hyodysenteriae",
  "BRC.SP", "Brachyspira species",
  "B.DIM", "Brevundimonas diminuta",
  "BRE.SP", "Brevundimonas species",
  "B.VES", "Brevundimonas vesicularis",
  "B.MEL", "Brucella melitensis",
  "BRU.SP", "Brucella species",
  "P.PSEM", "Burkholdera pseudomallei",
  "B.AMBI", "Burkholderia ambifaria",
  "B.ANTH", "Burkholderia anthina",
  "B.CENO", "Burkholderia cencocepacia",
  "B.CEP", "Burkholderia cepacia",
  "B.DOLO", "Burkholderia dolosa",
  "B.GLA", "Burkholderia gladioli",
  "B.MALL", "Burkholderia mallei",
  "B.MULT", "Burkholderia multivorans",
  "B.PSEM", "Burkholderia pseudomallei",
  "B.PYRR", "Burkholderia pyrrocinia",
  "BUR.SP", "Burkholderia species",
  "B.STAB", "Burkholderia stabilis",
  "B.VIET", "Burkholderia vietanamiensis",
  "CEG63", "Buttiauxella ferragutiae",
  "CEG64", "Buttiauxella gaviniae",
  "BUT.SP", "Buttiauxella species",
  "C.JEJ", "C. jejuni ss. jejuni",
  "C.HYO", "C.hyointest ss hyointestinalis",
  "C.SAC", "Caldicellulo. saccharolyticus",
  "CAD.SP", "Caldicellulosiruptor species",
  "C.FER", "Caloramator fervidus",
  "CAO.SP", "Caloramator species",
  "C.FVEN", "Campylo. fetus ss veneralis",
  "C.SBUB", "Campylo. sputorum ss bubulus",
  "C.SSPU", "Campylo. sputorum ss sputorum",
  "C.COL", "Campylobacter coli",
  "C.CON", "Campylobacter concisus",
  "C.NEG", "Campylobacter cult. negative",
  "C.FAE", "Campylobacter faecalis",
  "C.FET", "Campylobacter fetus ss fetus",
  "C.LAR", "Campylobacter lari",
  "C.MUC", "Campylobacter mucosalis",
  "CAM.SP", "Campylobacter species",
  "C.UPS", "Campylobacter upsaliensis",
  "C.ALB", "Candida albicans",
  "C.GLA", "Candida glabrata",
  "T.GLA", "Candida glabrata",
  "C.GUI", "Candida guilliermondii",
  "C.KRU", "Candida krusei",
  "CK6258", "Candida krusei ATCC 6258",
  "C.LIP", "Candida lipolytica",
  "C.LUS", "Candida lusitaniae",
  "C.PARA", "Candida parapsilosis",
  "CP2201", "Candida parapsilosis ATCC22019",
  "T.PIN", "Candida pintolopesii",
  "CA.PSE", "Candida pseudotropicalis",
  "C.RUG", "Candida rugosa",
  "T.CAN", "Candida saitoana",
  "CAN.SP", "Candida species",
  "C.STE", "Candida stellatoidea",
  "C.TRO", "Candida tropicalis",
  "C.CAN", "Capnocytophaga canimorsus",
  "C.CYN", "Capnocytophaga cynodegmi",
  "C.OCH", "Capnocytophaga ochracea",
  "CAP.SP", "Capnocytophaga species",
  "C.SPU", "Capnocytophaga sputigena",
  "C.HOM", "Cardiobacterium hominis",
  "CAR.SP", "Cardiobacterium species",
  "CEG1", "CDC Enteric Gr 1",
  "CEG11", "CDC Enteric Gr 11",
  "CEG41", "CDC Enteric Gr 41",
  "CEG58", "CDC Enteric Gr 58",
  "CEG59", "CDC Enteric Gr 59",
  "CEG60", "CDC Enteric Gr 60",
  "CEG68", "CDC Enteric Gr 68",
  "CEG69", "CDC Enteric Gr 69",
  "CEG76", "CDC Enteric Gr 76",
  "C74/75", "CDC Enteric Group 74 & 75",
  "CDC.EN", "CDC Enteric Groups",
  "C.DAV", "Cedecea davisae",
  "C.LAP", "Cedecea lapagei",
  "C.NET", "Cedecea neteri",
  "CED.SP", "Cedecea species",
  "CHR.SP", "Chromobacterium species",
  "C.VIO", "Chromobacterium violaceum",
  "C.BAL", "Chryseobacterium balustinum",
  "C.GLE", "Chryseobacterium gleum",
  "C.INDG", "Chryseobacterium indologenes",
  "C.INDT", "Chryseobacterium indoltheticum",
  "C.SCO", "Chryseobacterium scophthalmum",
  "CHB.SP", "Chryseobacterium species",
  "C.LUT", "Chryseomonas luteola",
  "CHY.SP", "Chryseomonas species",
  "C.AMA", "Citrobacter amalonaticus",
  "C.AMA1", "Citrobacter amalonaticus bio1",
  "C.BRA", "Citrobacter braakii",
  "C.FAR", "Citrobacter farmeri",
  "C.FRE", "Citrobacter freundii",
  "C.KOS", "Citrobacter koseri",
  "C.SED", "Citrobacter sedlakii",
  "CIT.SP", "Citrobacter species",
  "C.WER", "Citrobacter werkmanii",
  "C.YOU", "Citrobacter youngae",
  "CLA.SP", "Cladosporium species",
  "C.ALG", "Clostridium algidicarnis",
  "CL.PAE", "Clostridium baratii",
  "C.BEI", "Clostridium beijerinckii",
  "C.BIF", "Clostridium bifermentans",
  "C.BOT", "Clostridium botulinum",
  "C.BUT", "Clostridium butyricum",
  "C.CAD", "Clostridium cadaveris",
  "C.CHA", "Clostridium chauvoei",
  "C.DIF", "Clostridium difficile",
  "C.GLY", "Clostridium glycolicum",
  "C.HAS", "Clostridium hastiforme",
  "C.HIS", "Clostridium histolyticum",
  "C.INN", "Clostridium innocuum",
  "C.LIM", "Clostridium limosum",
  "C.NOVA", "Clostridium novyi",
  "CL.PAU", "Clostridium paraputrificum",
  "C.PER", "Clostridium perfringens",
  "C.RAM", "Clostridium ramosum",
  "C.SEP", "Clostridium septicum",
  "C.SOR", "Clostridium sordellii",
  "CLO.SP", "Clostridium species",
  "C.SPO", "Clostridium sporogenes",
  "C.SUB", "Clostridium subterminale",
  "F.SYM", "Clostridium symbiosum",
  "CL.TER", "Clostridium tertium",
  "C.TET", "Clostridium tetani",
  "CNSG", "Coagulase Negative Staph Group",
  "CPSG", "Coagulase Positive Staph Group",
  "C.IMM", "Coccidiodes immitis",
  "EU.AER", "Collinsella aerofaciens",
  "COL.SP", "Collinsella species",
  "COM.SP", "Comamonas species",
  "C.TER", "Comamonas terrigena",
  "C.TES", "Comamonas testosteroni",
  "C.PSED", "Coryne. pseudodiptheriticum",
  "C.PSEG", "Coryne. pseudogenitalium",
  "C.PSET", "Coryne. pseudotuberculosis",
  "C.BOV", "Corynebacterium bovis",
  "C.DIP", "Corynebacterium diptheriae",
  "C.GEN", "Corynebacterium genitalium",
  "C.GLU", "Corynebacterium glutamicum",
  "CGJK", "Corynebacterium Group JK",
  "C.HAE", "Corynebacterium haemolyticum",
  "C.JEI", "Corynebacterium jeikeium",
  "C.KUT", "Corynebacterium kutscheri",
  "C.MAT", "Corynebacterium matruchotii",
  "C.MIN", "Corynebacterium minutissimum",
  "C.MYC", "Corynebacterium mycetoides",
  "C.OVI", "Corynebacterium ovis",
  "C.REN", "Corynebacterium renale",
  "COR.SP", "Corynebacterium species",
  "C.UCL", "Corynebacterium ulcerans",
  "C.URE", "Corynebacterium urealyticum",
  "C.XER", "Corynebacterium xerosis",
  "CT.ALB", "Cryptococcus albidus",
  "C.LAU", "Cryptococcus laurentii",
  "CT.LUT", "Cryptococcus luteolus",
  "C.NEO", "Cryptococcus neoformans",
  "CRY.SP", "Cryptococcus species",
  "CT.TER", "Cryptococcus terrus",
  "CRS.SP", "Cryptosporidium species",
  "CUR.SP", "Curvularia species",
  "D.ACI", "Delftia acidovorans",
  "DEL.SP", "Delftia species",
  "D.CON", "Dermatophilus congolensis",
  "DER.SP", "Dermatophilus species",
  "D.DEH", "Desulfitobact. dehalogenans",
  "DEB.SP", "Desulfitobacterium species",
  "D.BAK", "Desulfuromusa bakii",
  "D.KYS", "Desulfuromusa kysingii",
  "DEM.SP", "Desulfuromusa species",
  "D.SUC", "Desulfuromusa succinoxidans",
  "E.CHM", "E. coli haemolytic mucoid",
  "E.CIN-", "E. coli indole -",
  "E.CNH", "E. coli non-haemolytic",
  "E.CNHM", "E. coli non-haemolytic mucoid",
  "E.NLF", "E. coli non-lactose fermenter",
  "E.SLF", "E. coli slow lactose fermenter",
  "E.COLS", "E. coli sorbitol - *",
  "E.HOS", "Edwardsiella hoshinae",
  "E.ICT", "Edwardsiella ictaluri",
  "EDW.SP", "Edwardsiella species",
  "E.TAR", "Edwardsiella tarda",
  "ET1594", "Edwardsiella tarda ATCC 15947",
  "E.TAR1", "Edwardsiella tarda bio 1",
  "E.LEN", "Eggerthella lenta",
  "EGG.SP", "Eggerthella species",
  "E.COR", "Eikenella corrodens",
  "EIK.SP", "Eikenella species",
  "C.MEN", "Elizabethkingia meningoseptica",
  "E.BRE", "Empedobacter brevis",
  "EMP.SP", "Empedobacter species",
  "EF2921", "Entero. faecalis ATCC 29212",
  "E.ZYM", "Entero. faecalis ss zymogenes",
  "E.AER", "Enterobacter aerogenes",
  "E.AMN", "Enterobacter amnigenus",
  "E.AMN1", "Enterobacter amnigenus bio1",
  "E.AMN2", "Enterobacter amnigenus bio2",
  "E.ASB", "Enterobacter asburiae",
  "E.CAN", "Enterobacter cancerogenus",
  "E.CLO", "Enterobacter cloacae",
  "E.GER", "Enterobacter gergoviae",
  "E.HOR", "Enterobacter hormaechei",
  "E.SAK", "Enterobacter sakazakii",
  "ENT.SP", "Enterobacter species",
  "ENTBAC", "Enterobacteriaceae",
  "E.AVI", "Enterococcus avium",
  "E.AVL", "Enterococcus avium-like",
  "E.CAS", "Enterococcus casseliflavus",
  "E.DUR", "Enterococcus durans",
  "E.FAEL", "Enterococcus faecalis",
  "E.FAEM", "Enterococcus faecium",
  "E.GALM", "Enterococcus gallinarum",
  "E.HIR", "Enterococcus hirae",
  "E.MAL", "Enterococcus malodoratus",
  "E.SOL", "Enterococcus solitarius",
  "ENC.SP", "Enterococcus species",
  "E.FLO", "Epidermophyton floccosum",
  "EPI.SP", "Epidermophyton species",
  "ERW.SP", "Erwinia species",
  "E.RHU", "Erysipelothrix rhusiopathiae",
  "ERY.SP", "Erysipelothrix species",
  "E.BLA", "Escherichia blattae",
  "E.COL", "Escherichia coli",
  "EC2592", "Escherichia coli ATCC 25922",
  "E.CAT", "Escherichia coli atypical",
  "E.BETA", "Escherichia coli Beta",
  "E.H2S", "Escherichia coli H2S+",
  "E.COLH", "Escherichia coli haemolytic",
  "E.K88", "Escherichia coli K88+",
  "E.K987", "Escherichia coli K987P+",
  "E.K99", "Escherichia coli K99+",
  "E.COLM", "Escherichia coli mucoid",
  "E.COLU", "Escherichia coli urea +",
  "E.FER", "Escherichia fergusonii",
  "E.HER", "Escherichia hermannii",
  "ESC.SP", "Escherichia species",
  "E.VUL", "Escherichia vulneris",
  "E.BAR", "Eubacterium barkeri",
  "E.LIM", "Eubacterium limosum",
  "E.MON", "Eubacterium moniliforme",
  "EUB.SP", "Eubacterium species",
  "E.AME", "Ewingella americana",
  "EWI.SP", "Ewingella species",
  "EXO.SP", "Exophialia species",
  "F.NEC", "F. necrophorum ss necrophorum",
  "F.NUC", "F. nucleatum ss nucleatum",
  "F.TUL", "F. tularensis ss tularensis",
  "FIL.SP", "Filifactor species",
  "F.VIL", "Filifactor villosus",
  "P.MAG", "Finegoldia magna",
  "FIN.SP", "Finegoldia species",
  "F.ORY", "Flavimonas oryzihabitans",
  "FLV.SP", "Flavimonas species",
  "S.MIZ", "Flavobacterium mizutaii",
  "FLA.SP", "Flavobacterium species",
  "FRA.SP", "Francisella species",
  "FUNG", "Fungal isolate",
  "FUA.SP", "Fusarium species",
  "F.MOR", "Fusobacterium mortiferum",
  "FUS.SP", "Fusobacterium species",
  "F.VAR", "Fusobacterium varium",
  "PA.ANA", "Gallibacterium anatis",
  "GAL.SP", "Gallibacterium species",
  "GAR.SP", "Gardnerella species",
  "G.VAG", "Gardnerella vaginalis",
  "G.HAE", "Gemella haemolysans",
  "G.MOR", "Gemella morbillorum",
  "GEM.SP", "Gemella species",
  "GEO.SP", "Geobacillus species",
  "B.STE", "Geobacillus stearothermophilus",
  "G.CAN", "Geotrichum candidum",
  "GTC.SP", "Geotrichum species",
  "G.SAN", "Globicatella sanguinis",
  "GLO.SP", "Globicatella species",
  "G.SUL", "Globicatella sulfidifaciens",
  "G.AMA", "Gordona aichiensis",
  "G.AIC", "Gordonia amarae",
  "GOR.SP", "Gordonia species",
  "G.SPU", "Gordonia sputi",
  "GRAM-", "Gram negative organism",
  "GRAM+", "Gram positive organism",
  "A.ADJ", "Granulicatella adiacens",
  "GRA.SP", "Granulicatella species",
  "HACEK", "HACEK Group",
  "A.AMY", "Haemo. actinomycetemcomitans",
  "HP4924", "Haemo. influenzae ATCC 49247",
  "H.AEG", "Haemophilus aegyptius",
  "H.AGN", "Haemophilus agni",
  "H.APH", "Haemophilus aphrophilus",
  "H.INF", "Haemophilus influenzae",
  "H.INFL", "Haemophilus influenzae blact +",
  "H.PARG", "Haemophilus paragallinarum",
  "H.PARI", "Haemophilus parainfluenzae",
  "H.PARP", "Haemophilus paraphrophilus",
  "H.PARS", "Haemophilus parasuis",
  "H.SEG", "Haemophilus segnis",
  "HAE.SP", "Haemophilus species",
  "H.SUT", "Haemophilus suis",
  "HAELSP", "Haemophilus-like species",
  "H.ALV", "Hafnia alvei",
  "HAF.SP", "Hafnia species",
  "H.SSAC", "Hal. sacc. ss saccharolyticum",
  "H.SSEN", "Hal. sacc. ss senegalense",
  "HLN.SP", "Halanaerobium species",
  "H.SAL", "Halococcus salifodinae",
  "HAL.SP", "Halococcus species",
  "HE.CAN", "Helicobacter canis",
  "H.CIN", "Helicobacter cinaedi",
  "H.FEL", "Helicobacter felis",
  "H.FEN", "Helicobacter fennelliae",
  "H.MUS", "Helicobacter mustelae",
  "H.PYL", "Helicobacter pylori",
  "HE.SAL", "Helicobacter salomonis",
  "HEL.SP", "Helicobacter species",
  "H.SOM", "Histophilus somni",
  "HIS.SP", "Histophilus species",
  "H.FOE", "Holophaga foetida",
  "HOL.SP", "Holophaga species",
  "C.LIV", "Janthinobacterium lividium",
  "JAN.SP", "Janthinobacterium species",
  "K.PNER", "K. pneumo. ss rhinoscleromatis",
  "KIN.SP", "Kingella species",
  "K.PNPN", "Kleb. pneumoniae ss pneumoniae",
  "K.GP47", "Klebsiella group 47",
  "K.OXY", "Klebsiella oxytoca",
  "KO8724", "Klebsiella oxytoca ATCC 8724",
  "K.PNEO", "Klebsiella pneumoni ss ozaenae",
  "K.PNEU", "Klebsiella pneumoniae",
  "KLE.SP", "Klebsiella species",
  "K.ASC", "Kluyvera ascorbata",
  "K.CRY", "Kluyvera cryocrescens",
  "E.INT", "Kluyvera intermedia",
  "KLU.SP", "Kluyvera species",
  "K.KRI", "Kocuria kristinae",
  "K.ROS", "Kocuria rosea",
  "MR186", "Kocuria rosea ATCC 186",
  "KOC.SP", "Kocuria species",
  "M.VAR", "Kocuria varians",
  "K.SED", "Kytococcus sedentarius",
  "KYT.SP", "Kytococcus species",
  "L.PNE", "L. pneumophila ss pneumophila",
  "E.SER", "Lacterococcus garvieae",
  "L.ACI", "Lactobacillus acidophilus",
  "L.CAT", "Lactobacillus catenaformis",
  "L.FER", "Lactobacillus fermentum",
  "L.JEN", "Lactobacillus jensenii",
  "LAC.SP", "Lactobacillus species",
  "S.CRE", "Lactococcus lactis ss cremoris",
  "S.LAC", "Lactococcus lactis ss lactis",
  "LCC.SP", "Lactococcus species",
  "L.ADE", "Leclercia adecarboxylata",
  "LEC.SP", "Leclercia species",
  "LEG.SP", "Legionella species",
  "L.GRI", "Leminorella grimontii",
  "L.RIC", "Leminorella richardii",
  "LEM.SP", "Leminorella species",
  "LEP.SP", "Leptospira species",
  "L.CRE", "Leuc. mesenteroide ss cremoris",
  "L.LAC", "Leuconostoc lactis",
  "LEU.SP", "Leuconostoc species",
  "L.GRA", "Listeria grayi",
  "L.INN", "Listeria innocua",
  "L.IVA", "Listeria ivanovii ss ivanovii",
  "L.MON", "Listeria monocytogenes",
  "L.SEE", "Listeria seeligeri",
  "LIS.SP", "Listeria species",
  "L.WEL", "Listeria welshimeri",
  "L.ANG", "Listonella anguillarum",
  "LIT.SP", "Listonella species",
  "S.CASE", "Macrococcus caseolyticus",
  "MAC.SP", "Macrococcus species",
  "M.PAC", "Malassezia pachydermatis",
  "MAL.SP", "Malassezia species",
  "M.HAEA", "Mannheimia haemolytica",
  "M.HAEL", "Mannheimia haemolytica-like",
  "MAN.SP", "Mannheimia species",
  "M.CHA", "Mesoplasma chauliocola",
  "M.COL", "Mesoplasma coleopterae",
  "M.COR", "Mesoplasma corruscae",
  "M.GRA", "Mesoplasma grammopterae",
  "M.PHO", "Mesoplasma photuris",
  "M.PLE", "Mesoplasma pleciae",
  "MES.SP", "Mesoplasma species",
  "M.SYR", "Mesoplasma syrphidae",
  "M.TAB", "Mesoplasma tabanidae",
  "M.BOM", "Methanolobus bombayensis",
  "MTB.SP", "Methanolobus species",
  "M.AMI", "Methylobacterium aminovorans",
  "M.FUJ", "Methylobacterium fujisawaense",
  "M.MES", "Methylobacterium mesophilicum",
  "M.ORG", "Methylobacterium organophilum",
  "M.RAD", "Methylobacterium radiotolerans",
  "M.RHOE", "Methylobacterium rhodesianum",
  "M.RHOI", "Methylobacterium rhodinum",
  "MET.SP", "Methylobacterium species",
  "M.ZAT", "Methylobacterium zatmanii",
  "M.LUT", "Micrococcus luteus",
  "M.LYL", "Micrococcus lylae",
  "MIC.SP", "Micrococcus species",
  "M.DIS", "Microspor canis var distortum",
  "M.AMA", "Microsporum amazonicum",
  "M.AUD", "Microsporum audouinii",
  "M.BOU", "Microsporum boullardi",
  "M.CAN", "Microsporum canis",
  "M.FER", "Microsporum ferrugineum",
  "M.GAL", "Microsporum gallinae",
  "M.GYP", "Microsporum gypseum",
  "M.NAN", "Microsporum nanum",
  "M.PER", "Microsporum persicolor",
  "M.PRA", "Microsporum praecox",
  "M.RAC", "Microsporum racemosum",
  "MIP.SP", "Microsporum species",
  "M.VAN", "Microsporum vanbreuseghemii",
  "MC", "Mixed Culture",
  "MOE.SP", "Moellerella species",
  "M.WIS", "Moellerella wisconsensis",
  "MOO.SP", "Moorella species",
  "M.THAC", "Moorella thermoacetica",
  "M.THAU", "Moorella thermoautotrophica",
  "BRALSP", "Mor (Branhamella)-like species",
  "B.CAT", "Moraxella (Bran) catarrhalis",
  "BH.CAN", "Moraxella (Bran) cuniculi",
  "BH.OVI", "Moraxella (Branhamella) ovis",
  "BRA.SP", "Moraxella (Branhamella) spp",
  "M.ATL", "Moraxella (Moraxella) atlantae",
  "MO.BOV", "Moraxella (Moraxella) bovis",
  "M.LAC", "Moraxella (Moraxella) lacunata",
  "MXM.SP", "Moraxella (Moraxella) spp",
  "M.NLF", "Moraxella nonliquefaciens",
  "M.OSL", "Moraxella osloensis",
  "MOX.SP", "Moraxella species",
  "M.MOR", "Morganella morgani ss morganii",
  "MM2583", "Morganella morganii ATCC 25830",
  "M.MOR1", "Morganella morganii bio 1",
  "MOR.SP", "Morganella species",
  "MUC.SP", "Mucor species",
  "M.AVI", "Myco. avium ss avium",
  "M.BOV", "Myco. bovis ss bovis",
  "M.FOR", "Myco. fortuitum ss fortuitum",
  "M.SCR", "Myco. scrofulaceum",
  "M.PART", "Myco.avium ss paratuberculosis",
  "M.ABS", "Mycobacterium abscessus",
  "M.ASI", "Mycobacterium asiaticum",
  "M.BOK", "Mycobacterium bonickei",
  "M.CHE", "Mycobacterium chelonae",
  "M.GAS", "Mycobacterium gastri",
  "M.GOO", "Mycobacterium goodii",
  "M.GOR", "Mycobacterium gordonae",
  "M.HAE", "Mycobacterium haemophilum",
  "M.HOU", "Mycobacterium houstonense",
  "M.IMM", "Mycobacterium immunogenum",
  "M.INTJ", "Mycobacterium interjectum",
  "M.INT", "Mycobacterium intracellulare",
  "M.KAN", "Mycobacterium kansasii",
  "M.MAG", "Mycobacterium mageritense",
  "M.MAL", "Mycobacterium malmoenese",
  "M.MAR", "Mycobacterium marinum",
  "M.MUC", "Mycobacterium mucogenicum",
  "M.NCH", "Mycobacterium nonchromogenicum",
  "M.PERE", "Mycobacterium peregrinum",
  "M.POR", "Mycobacterium porcinum",
  "M.SEN", "Mycobacterium senegalense",
  "M.SIM", "Mycobacterium simiae",
  "M.SME", "Mycobacterium smegmatis",
  "MYC.SP", "Mycobacterium species",
  "M.SZU", "Mycobacterium szulgai",
  "M.TER", "Mycobacterium terrae",
  "M.TRI", "Mycobacterium triviale",
  "M.TUB", "Mycobacterium tuberculosis",
  "M.WOL", "Mycobacterium wolinskyi",
  "M.XEN", "Mycobacterium xenopi",
  "MYP.SP", "Mycoplasma species",
  "M.ODO", "Myroides odoratus",
  "MYR.SP", "Myroides species",
  "N.DAS", "N. dassonvillei ss albirubida",
  "M.COO", "Namizzia cajetani",
  "NAM.SP", "Namizzia species",
  "N.ENIT", "Nei. elongata ss nitroreducens",
  "N.CAN", "Neisseria canis",
  "N.DEN", "Neisseria denitrificans",
  "N.GON", "Neisseria gonorrhoeae",
  "N.MEN", "Neisseria meningitidis",
  "NEI.SP", "Neisseria species",
  "N.WEA", "Neisseria weaveri",
  "N.SEM", "Nicoletella semolina",
  "NG", "No Growth",
  "NOID", "No Identification Possible",
  "NSG", "No Significant Growth",
  "N.AST", "Nocardia asteroides",
  "N.BRA", "Nocardia brasiliensis",
  "N.FAR", "Nocardia farcinica",
  "N.CAV", "Nocardia otitidiscaviarum",
  "NOC.SP", "Nocardia species",
  "NONENT", "NON-Enterobacteriaceae",
  "NRO", "Non-reactive organism",
  "O.PRO", "Obesumbacterium proteus",
  "OBS.SP", "Obesumbacterium species",
  "O.ANT", "Ochrobactrum anthropi",
  "OCH.SP", "Ochrobactrum species",
  "OLI.SP", "Oligella species",
  "O.UREO", "Oligella ureolytica",
  "O.URET", "Oligella urethralis",
  "O.RHI", "Ornithobact. rhinotracheale",
  "ORN.SP", "Ornithobacterium species",
  "O.OXA", "Oxalophagus oxalicus",
  "OXA.SP", "Oxalophagus species",
  "O.PFE", "Oxobacter pfennigii",
  "OXO.SP", "Oxobacter species",
  "P.PSE", "P.pseudo. ss pseudoalcaligenes",
  "PAE.SP", "Paecilomyces species",
  "P.LAR", "Paen. larvae ss larvae",
  "P.PUL", "Paen. larvae ss pulvifaciens",
  "P.ALV", "Paenibacillus alvei",
  "P.AMY", "Paenibacillus amylolyticus",
  "P.AZO", "Paenibacillus azotofixans",
  "P.MACE", "Paenibacillus macerans",
  "P.MACQ", "Paenibacillus macquariensis",
  "P.PUB", "Paenibacillus pabuli",
  "PAN.SP", "Paenibacillus species",
  "P.GOR", "Paenibacillus validus",
  "P.VAL", "Paenibacillus validus",
  "P.AGG", "Pantoea agglomerans",
  "PAT.SP", "Pantoea species",
  "P.MGAL", "Past. multocida ss gallicida",
  "P.MUL", "Past. multocida ss multocida",
  "P.MSEP", "Past. multocida ss septica",
  "PA.AER", "Pasteurella aerogenes",
  "PA.AEA", "Pasteurella aerogenes atypical",
  "PA.AVI", "Pasteurella avium",
  "P.CAB", "Pasteurella caballi",
  "P.CAN", "Pasteurella canis",
  "P.DAG", "Pasteurella dagmatis",
  "P.GAL", "Pasteurella gallinarum",
  "P.GRPA", "Pasteurella group A",
  "P.GRPB", "Pasteurella group B",
  "P.LAN", "Pasteurella langaaensis",
  "P.PNE", "Pasteurella pneumotropica",
  "PAS.SP", "Pasteurella species",
  "PASTK", "Pasteurella species Taxon K",
  "P.STO", "Pasteurella stomatis",
  "PA.TES", "Pasteurella testudinis",
  "P.TRE", "Pasteurella trehalosi",
  "P.VOL", "Pasteurella volantium",
  "PASLSP", "Pasteurella-like species",
  "PED.SP", "Pediococcus species",
  "PEN.SP", "Penicillium species",
  "PEC.SP", "Peptocococcus species",
  "P.ASA", "Peptoniphilus asaccharolyticus",
  "PEP.SP", "Peptoniphilus species",
  "P.ANA", "Peptostreptococcus anaerobius",
  "PS.BOY", "Peptostreptococcus boydii",
  "P.MIC", "Peptostreptococcus micros",
  "PES.SP", "Peptostreptococcus species",
  "V.DAM", "Photobac. damselae ss damselae",
  "PHO.SP", "Photobacterium species",
  "PLE.SP", "Pleisiomonas speceis",
  "P.SHI", "Plesiomonas shigelloides",
  "B.ASA", "Porphyromonas asaccharolytica",
  "P.CANG", "Porphyromonas cangingivalis",
  "P.CANS", "Porphyromonas cansulci",
  "P.CRE", "Porphyromonas crevioricanis",
  "P.GIN", "Porphyromonas gingivicanis",
  "POR.SP", "Porphyromonas species",
  "B.RRUM", "Prev. ruminicola ss ruminicola",
  "B.BIV", "Prevotella bivia",
  "B.RBRE", "Prevotella brevis",
  "B.DISI", "Prevotella disiens",
  "P.ENO", "Prevotella enoeca",
  "B.MINT", "Prevotella intermedia",
  "B.MMEL", "Prevotella melaninogenica",
  "B.ORA", "Prevotella oralis",
  "PRE.SP", "Prevotella species",
  "P.TAN", "Prevotella tannerae",
  "P.ZOO", "Prevotella zoogleoformans",
  "P.AVI", "Propionibacterium avidum",
  "P.GRA", "Propionibacterium granulosum",
  "A.PRO", "Propionibacterium propionicus",
  "PRP.SP", "Propionibacterium species",
  "P.ACN", "Proprionibacterium acnes",
  "P.MIR", "Proteus mirabilis",
  "P.MYX", "Proteus myxofaciens",
  "P.PEN", "Proteus penneri",
  "PRT.SP", "Proteus species",
  "P.VUL", "Proteus vulgaris",
  "PV6896", "Proteus vulgaris ATCC 6896",
  "P.VUL2", "Proteus vulgaris bio 2",
  "P.VUL3", "Proteus vulgaris bio 3",
  "PROTO", "Protozoa",
  "P.ALCF", "Providencia alcalifaciens",
  "P.RET", "Providencia rettgeri",
  "P.RUS", "Providencia rustigianii",
  "PRV.SP", "Providencia species",
  "P.STU", "Providencia stuartii",
  "P.STU-", "Providencia stuartii ure -",
  "P.STU+", "Providencia stuartii ure +",
  "PSA.SP", "Pseudallescheria species",
  "PA1014", "Pseudo. aeruginosa ATCC 10145",
  "PA2785", "Pseudo. aeruginosa ATCC 27853",
  "P.AER", "Pseudomonas aeruginosa",
  "PA/PP", "Pseudomonas aeruginosa/putida",
  "P.ALCG", "Pseudomonas alcaligenes",
  "P.FLU", "Pseudomonas fluorescens",
  "P.MEN", "Pseudomonas mendocina",
  "P.PUT", "Pseudomonas putida",
  "PSE.SP", "Pseudomonas species",
  "P.STUT", "Pseudomonas stutzeri",
  "P.IMM", "Psychrobacter immobilis",
  "P.PHE", "Psychrobacter phenylpyruvicus",
  "PSY.SP", "Psychrobacter species",
  "RAH.SP", "Rahnella species",
  "P.THO", "Ralstonia mannitolilytica",
  "R.PIC", "Ralstonia pickettii",
  "RAL.SP", "Ralstonia species",
  "R.PLA", "Raoultella planticola",
  "RAO.SP", "Raoultella species",
  "R.TER", "Raoultella terrigena",
  "RMYCO", "Rapid Growing Mycobacteria",
  "RIP", "Reincubate Plate",
  "A.RAD", "Rhizobium radiobacter",
  "RHB.SP", "Rhizobium species",
  "R.ARR", "Rhizopus arrhizus",
  "RHI.SP", "Rhizopus species",
  "R.EQU", "Rhodococcus equi",
  "R.LUT", "Rhodococcus fascians",
  "R.OPA", "Rhodococcus opacus",
  "RHO.SP", "Rhodococcus species",
  "R.ELE", "Rhodoplanes elegans",
  "R.ROS", "Rhodoplanes roseus",
  "RHP.SP", "Rhodoplanes species",
  "RDT.SP", "Rhodotordula species",
  "R.RUB", "Rhodotorula rubra",
  "P.ANT", "Riemerella anatipestifer",
  "RIE.SP", "Riemerella species",
  "R.CER", "Roseomonas cervicalis",
  "R.FAU", "Roseomonas fauriae",
  "R.GIL", "Roseomonas gilardii",
  "RO1.SP", "Roseomonas sp. 1",
  "RO2.SP", "Roseomonas sp. 2",
  "RO3.SP", "Roseomonas sp. 3",
  "ROS.SP", "Roseomonas species",
  "S.MUC", "Rothia mucilaginosa",
  "ROT.SP", "Rothia species",
  "S.HOU", "S. cho. ss houtenae sv houten",
  "S.SALA", "S. choleroesuis ss salamae",
  "S.TCOP", "S. copenhagen",
  "S.GALL", "S.ent.ss enteric sv gallinarum",
  "SE7002", "S.epidermidis ATCC 700296",
  "S.CER", "Saccharomyces cerevisiae",
  "V.COS", "Salin. costicola ss costicola",
  "SLV.SP", "Salinovibrio species",
  "S.CHOC", "Salm. cholerae ss choleraesuis",
  "S.CHOA", "Salm. choleraesuis ss arizonae",
  "S.CSKU", "Salm. cholerasuis b.kunzendorf",
  "S.(A)", "Salm.sp (arizona?) *",
  "S.(SC)", "Salm.sp (cholerae-suis?) *",
  "S.(E)", "Salm.sp (enteritidis?) *",
  "S.(G)", "Salm.sp (gallinarum?) *",
  "S.(PA)", "Salm.sp (paratyphi-A?) *",
  "S.(PB)", "Salm.sp (paratyphi-B?) *",
  "S.(PC)", "Salm.sp (paratyphi-C?) *",
  "S.(P)", "Salm.sp (pullorum?) *",
  "S.(SE)", "Salm.sp (sendai?) *",
  "S.(1)", "Salm.sp (subgenus 1?) *",
  "S.(2)", "Salm.sp (subgenus 2?) *",
  "S(3/A)", "Salm.sp (subgenus 3/ariz?) *",
  "S.(3)", "Salm.sp (subgenus 3?) *",
  "S.(4)", "Salm.sp (subgenus 4?) *",
  "S.(5)", "Salm.sp (subgenus 5?) *",
  "S(1S+)", "Salm.sp (subgrp1 suc+ve?) *",
  "S.(T)", "Salm.sp (typhi?) *",
  "S.(TY)", "Salm.sp (typhimurium?) *",
  "A.HIN", "Salmonella arizonae",
  "S.ENT", "Salmonella enteritidis",
  "S.PTA", "Salmonella paratyphi-A",
  "S.PTB", "Salmonella paratyphi-B",
  "S.PTC", "Salmonella paratyphi-C",
  "S.PUL", "Salmonella pullorum",
  "S.SEN", "Salmonella sendai",
  "SAL.SP", "Salmonella species *",
  "S.SG1", "Salmonella subgenus 1",
  "S.SG2", "Salmonella subgenus 2",
  "S.SG3", "Salmonella subgenus 3",
  "S.SG4", "Salmonella subgenus 4",
  "S.SG5", "Salmonella subgenus 5",
  "S.S1S+", "Salmonella subgrp 1 sucr+ve",
  "S.TYP", "Salmonella typhi",
  "S.TYPM", "Salmonella typhimurium",
  "S.APIO", "Scedosporium apiospermum",
  "SCO.SP", "Scopulariopsis species",
  "SL/SM", "Ser. proteamaculans/marcescens",
  "S.FIC", "Serratia ficaria",
  "S.FON", "Serratia fonticola",
  "S.MAR", "Serratia marcescens",
  "S.MAR1", "Serratia marcescens bio 1",
  "S.ODO", "Serratia odorifera",
  "S.ODO1", "Serratia odorifera biogroup 1",
  "S.ODO2", "Serratia odorifera biogroup 2",
  "S.PLY", "Serratia plymuthica",
  "S.PRO", "Serratia proteamaculans",
  "S.RUB", "Serratia rubidaea",
  "SER.SP", "Serratia species",
  "S.PUTF", "Shewanella putrefaciens",
  "SHE.SP", "Shewanella species",
  "S.(BF)", "Shig.sp (boydii/flex?) *",
  "S.(D1)", "Shig.sp (dysenteriae type1?) *",
  "S.(D)", "Shig.sp (dysenteriae?) *",
  "S.BOY", "Shigella boydii",
  "SB/SF", "Shigella boydii/flexneri",
  "S.DYS", "Shigella dysenteriae",
  "S.D1", "Shigella dysenteriae type 1",
  "S.FLE", "Shigella flexneri",
  "S.SON", "Shigella sonnei",
  "SS2593", "Shigella sonnei ATCC 25931",
  "SHIUT", "Shigella species - untypeable",
  "SHI.SP", "Shigella species *",
  "SHILSP", "Shigella-like species",
  "S.FRE", "Sinorhizobium fredii",
  "S.MELI", "Sinorhizobium meliloti",
  "S.SAH", "Sinorhizobium saheli",
  "SIN.SP", "Sinorhizobium species",
  "S.TER", "Sinorhizobium terangae",
  "SMYCO", "Slow Growing NTB Mycobacteria",
  "S.MUL", "Sphingobacterium multivorum",
  "SPB.SP", "Sphingobacterium species",
  "S.SPI", "Sphingobacterium spiritivorum",
  "S.THA", "Sphingobacterium thalpophilum",
  "S.YAB", "Sphingobacterium yabuuchiae",
  "S.ADH", "Sphingomonas adhaesiva",
  "S.CAPS", "Sphingomonas capsulata",
  "S.PARA", "Sphingomonas parapaucimobilis",
  "S.PAU", "Sphingomonas paucimobilis",
  "SP1.SP", "Sphingomonas sp. 1",
  "SP2.SP", "Sphingomonas sp. 2",
  "SPM.SP", "Sphingomonas species",
  "S.YAN", "Sphingomonas yanoikuyae",
  "SPIRO", "Spirochetes",
  "S.SCHE", "Sporothrix schenckii",
  "SPT.SP", "Sporothrix species",
  "SPO.SP", "Sporotrichum species",
  "S.SAP", "Staph sapro. ss saprophyticus",
  "S.SCOA", "Staph schleiferi ss coagulans",
  "S.SCHL", "Staph schleiferi ss schleiferi",
  "S.SCI", "Staph sciuri ss sciuri",
  "S.CAPI", "Staph. capitis ss capitis",
  "S.CAR", "Staph. carnosus ss carnosus",
  "S.COH", "Staph. cohnii ss cohnii",
  "S.HOM", "Staph. hominis ss hominis",
  "SA2921", "Staph.aur.ss aureus ATCC 29213",
  "SA4330", "Staph.aur.ss.aureus ATCC 43300",
  "S.ARL", "Staphylococcus arlettae",
  "S.AUR", "Staphylococcus aur. ss aureus",
  "S.AURI", "Staphylococcus auricularis",
  "SA.BHA", "Staphylococcus beta haemolytic",
  "S.CAPR", "Staphylococcus caprae",
  "S.CHR", "Staphylococcus chromogenes",
  "S.DEL", "Staphylococcus delphini",
  "S.EPI", "Staphylococcus epidermidis",
  "S.EQUO", "Staphylococcus equorum",
  "S.GAL", "Staphylococcus gallinarum",
  "S.HAE", "Staphylococcus haemolyticus",
  "S.HYI", "Staphylococcus hyicus",
  "S.INT", "Staphylococcus intermedius",
  "S.KLO", "Staphylococcus kloosii",
  "S.LEN", "Staphylococcus lentus",
  "S.LUG", "Staphylococcus lugdunensis",
  "S.PSE", "Staphylococcus pseudintermediu",
  "S.SAC", "Staphylococcus saccharolyticus",
  "S.SCH", "Staphylococcus schleiferi",
  "S.SIM", "Staphylococcus simulans",
  "STA.SP", "Staphylococcus species",
  "S.WAR", "Staphylococcus warneri",
  "S.XYL", "Staphylococcus xylosus",
  "S.AFR", "Stenotrophomonas africana",
  "S.MALT", "Stenotrophomonas maltophilia",
  "STE.SP", "Stenotrophomonas species",
  "S.DEQU", "Str.dysgalact. ss equisimilis",
  "S.MIL", "Strep milleri",
  "S.PAR", "Strep parasanguinis",
  "S.SAN", "Strep sanguinis",
  "V.STR", "Strep Viridans Group",
  "S.CON", "Strep. con. ss constellatus",
  "S.DYSG", "Strep. dys ss dysgalactiae",
  "S.EQUS", "Strep. dys. ss equisimilis",
  "S.EQUI", "Strep. equi ss equi",
  "S.ZOO", "Strep. equi ss zooepidemicus",
  "SP4961", "Strep. pneumoniae ATCC 49619",
  "S.ACI", "Streptococcus acidominimus",
  "S.AGA", "Streptococcus agalactiae",
  "S.ANG", "Streptococcus anginosus",
  "S.GRPA", "Streptococcus beta group - A",
  "S.GRPB", "Streptococcus beta group - B",
  "S.GRPC", "Streptococcus beta group - C",
  "S.GRPD", "Streptococcus beta group - D",
  "S.GRPE", "Streptococcus beta group - E",
  "S.GRPF", "Streptococcus beta group - F",
  "S.GRPG", "Streptococcus beta group - G",
  "S.GRPL", "Streptococcus beta group - L",
  "S.GRPN", "Streptococcus beta group - N",
  "S.BHA", "Streptococcus beta haemolytic",
  "S.BOV", "Streptococcus bovis",
  "S.BOVV", "Streptococcus bovis variens",
  "S.BOV1", "Streptococcus bovis variens 1",
  "S.BOV2", "Streptococcus bovis variens 2",
  "S.CAN", "Streptococcus canis",
  "S.CRI", "Streptococcus cristatus",
  "S.DIF", "Streptococcus difficilis",
  "S.EQUN", "Streptococcus equinus",
  "S.EQUC", "Streptococcus equisimilis GP C",
  "S.GOR", "Streptococcus gordonii",
  "S.INFA", "Streptococcus infantis",
  "S.INI", "Streptococcus iniae",
  "S.INTM", "Streptococcus intermedius",
  "S.MIL2", "Streptococcus milleri II",
  "S.MITS", "Streptococcus mitis",
  "S.MUT", "Streptococcus mutans",
  "S.ORA", "Streptococcus oralis",
  "S.PER", "Streptococcus peroris",
  "S.PHO", "Streptococcus phocae",
  "S.PNE", "Streptococcus pneumoniae",
  "S.POR", "Streptococcus porcinus",
  "S.PYO", "Streptococcus pyogenes",
  "S.SALI", "Streptococcus salivarius",
  "S.SAN1", "Streptococcus sanguinis I",
  "S.SOB", "Streptococcus sobrinus",
  "STR.SP", "Streptococcus species",
  "S.SUI", "Streptococcus suis",
  "S.THE", "Streptococcus thermophilus",
  "S.UBE", "Streptococcus uberis",
  "S.ANU", "Streptomyces anulatus",
  "S.SOM", "Streptomyces somaliensis",
  "STM.SP", "Streptomyces species",
  "T.PTY", "Tatumella ptyseos",
  "TA.EQU", "Taylorella equigenitalis",
  "TAY.SP", "Taylorella species",
  "T.THEC", "Th. thermocopriae",
  "T.THES", "Th. thermosaccharolyticum",
  "T.KIV", "Thermoanaerobacter kivui",
  "THA.SP", "Thermoanaerobacter species",
  "T.PRO", "Thermococcus profundus",
  "TMC.SP", "Thermococcus species",
  "T.COM", "Thermodesulfobacterium commune",
  "THD.SP", "Thermodesulfobacterium species",
  "TOR.SP", "Torulopsis species",
  "T.FOE", "Trichomonas foetus",
  "TRM.SP", "Trichomonas species",
  "T.AJE", "Trichophyton ajelloi",
  "T.CON", "Trichophyton concentricum",
  "TP.EQU", "Trichophyton equinum",
  "T.MEN", "Trichophyton mentagrophytes",
  "T.RUB", "Trichophyton rubrum",
  "T.SCH", "Trichophyton schoenleinii",
  "T.SOU", "Trichophyton soudanense",
  "TRI.SP", "Trichophyton species",
  "T.TER", "Trichophyton terrestre",
  "T.TON", "Trichophyton tonsurans",
  "T.VER", "Trichophyton verrucosum",
  "T.BEI", "Trichosperon beigelii",
  "TRC.SP", "Trichosperon species",
  "TRO.SP", "Tropheryma species",
  "T.WHI", "Tropheryma whipplei",
  "T.INC", "Tsukamurella inchonensis",
  "T.PAU", "Tsukamurella paurometabola",
  "T.PULM", "Tsukamurella pulmonis",
  "TSU.SP", "Tsukamurella species",
  "T.STR", "Tsukamurella strandjordae",
  "T.TYR", "Tsukamurella tyrosinosolvens",
  "V.PRV", "Veillonella parvula",
  "VEI.SP", "Veillonella species",
  "V.ALG", "Vibrio alginolyticus",
  "V.CHO1", "Vibrio cholerae 0.1",
  "VC-S-", "Vibrio cholerae 0.1-sucrose-*",
  "VC-S+", "Vibrio cholerae 0.1-sucrose+*",
  "VC+ET", "Vibrio cholerae 01+ El Tor*",
  "VC-ET", "Vibrio cholerae 01-El Tor*",
  "V.CHO", "Vibrio cholerae*",
  "V.DIA", "Vibrio diazotrophicus",
  "V.FLU", "Vibrio fluvialis",
  "V.FUR", "Vibrio furnissii",
  "V.HAR", "Vibrio harveyi",
  "V.HOL", "Vibrio hollisae",
  "V.MET", "Vibrio metschnikovii",
  "V.MIM", "Vibrio mimicus",
  "V.NER", "Vibrio nereis",
  "V.PAR", "Vibrio parhaemolyticus",
  "V.PRO", "Vibrio proteolyticus",
  "VIB.SP", "Vibrio species",
  "V.VUL", "Vibrio vulnificus",
  "VIBLSP", "Vibrio-like species",
  "WEE.SP", "Weeksella species",
  "W.VIR", "Weeksella virosa",
  "WEI.SP", "Weisella species",
  "L.VIR", "Weissella viridescens",
  "XAN.SP", "Xanthomomas species",
  "YEAST", "Yeast",
  "Y.ENT", "Yer. ent. ss enterocolitica",
  "Y.FRE", "Yersinia frederiksenii",
  "Y.INT", "Yersinia intermedia",
  "Y.KRI", "Yersinia kristensenii",
  "Y.PEST", "Yersinia pestis",
  "Y.PSE", "Yersinia pseudotuberculosis",
  "Y.RUC", "Yersinia ruckeri",
  "YER.SP", "Yersinia species",
  "Y.PES", "Yersinia spp. (pestis?)",
  "Y.REG", "Yokenella regensburgei",
  "YOK.SP", "Yokenella species"
)

codes$mo_name <- mo_name(codes$name, keep_synonyms = FALSE, minimum_matching_score = 0.650)

import <- codes |> filter(mo_name != "(unknown name)")

import$mo <- as.mo(import$mo_name)

microorganisms.codes <- microorganisms.codes |>
  bind_rows(
    tibble(code = toupper(import$code),
           mo = import$mo) |>
      distinct()
  ) |>
  arrange(code)
class(microorganisms.codes$mo) <- c("mo", "character")
usethis::use_data(microorganisms.codes, overwrite = TRUE, compress = "xz", version = 2)
rm(microorganisms.codes)
devtools::load_all()

