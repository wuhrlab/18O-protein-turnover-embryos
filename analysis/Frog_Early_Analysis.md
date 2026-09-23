Frog Early Analysis
================
2026-01-07

``` r
rm(list=ls(all=T))
library(lsa)
```

    ## Loading required package: SnowballC

``` r
library(ggplot2)
library(fgsea)
library(stats)
library(GO.db)
```

    ## Loading required package: AnnotationDbi

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## Loading required package: generics

    ## 
    ## Attaching package: 'generics'

    ## The following objects are masked from 'package:base':
    ## 
    ##     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
    ##     setequal, union

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
    ##     unsplit, which.max, which.min

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: IRanges

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:grDevices':
    ## 
    ##     windows

    ## 

``` r
library(tidyr)
```

    ## 
    ## Attaching package: 'tidyr'

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     expand

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following object is masked from 'package:AnnotationDbi':
    ## 
    ##     select

    ## The following objects are masked from 'package:IRanges':
    ## 
    ##     collapse, desc, intersect, setdiff, slice, union

    ## The following objects are masked from 'package:S4Vectors':
    ## 
    ##     first, intersect, rename, setdiff, setequal, union

    ## The following object is masked from 'package:Biobase':
    ## 
    ##     combine

    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     combine, intersect, setdiff, setequal, union

    ## The following object is masked from 'package:generics':
    ## 
    ##     explain

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(stringr)
library(KEGGREST)
library(rrvgo)
library(org.Hs.eg.db)
```

    ## 

``` r
library(clusterProfiler)
```

    ## clusterProfiler v4.18.4 Learn more at https://yulab-smu.top/contribution-knowledge-mining/
    ## 
    ## Please cite:
    ## 
    ## G Yu. Thirteen years of clusterProfiler. The Innovation. 2024,
    ## 5(6):100722

    ## 
    ## Attaching package: 'clusterProfiler'

    ## The following object is masked from 'package:AnnotationDbi':
    ## 
    ##     select

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     slice

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     rename

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

``` r
library(minpack.lm)
library(Biostrings)
```

    ## Loading required package: XVector

    ## Loading required package: Seqinfo

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

``` r
library(Peptides)
library(ggpubr)
library(forcats)
library(effsize)
library(patchwork)

knitr::opts_chunk$set(fig.path = "figures/Frog_Early_Analysis/")
```

``` r
Frog_T5_all.min <- c(161, 221, 281, 401, 521, 641,
                     881, 1605, 1961, 2918, 4405) - 161
Frog_T5_light.min <- Frog_T5_all.min[c(1:2, 4:5, 7:8, 10:11)]

Frog_T6_all.min <- c(150, 210, 270, 390, 510, 635,
                     930, 1591, 1950, 2910, 4394) - 150
Frog_T6_light.min <- Frog_T6_all.min[c(1:2, 4:5, 7:8, 10:11)]

avg_light_time <- (Frog_T5_light.min + Frog_T6_light.min)/2
avg_heavy_time <- (Frog_T5_all.min + Frog_T6_all.min)/2

light.channels <- c("T0_L", "N1", "N3", "N4",
                    "N7", "N9",  "N11", "N12")
heavy.channels <- c("T0_H", "O1", "O2", "O3", "O4", "O5",
                    "O7", "O9", "O10", "O11", "O12")

XLA_fits_df <- read.csv("Files/Fits/Dyn-Model/XLA_O18_T5-6_FinalFits_ALL.csv")
head(XLA_fits_df)
```

    ##                   Protein_ID   XLA_Gene Human_Gene         kD    HL_Hrs
    ## 1  XBmRNA10569|XBXL10_1g5638    ccnb1.S      CCNB1 0.05817981 0.1985646
    ## 2 XBmRNA35507|XBXL10_1g19214  ccnb1.2.L      CCNB1 0.05236107 0.2206306
    ## 3 XBmRNA34126|XBXL10_1g18541  LOC398156      PTTG1 0.05055962 0.2284917
    ## 4 XBmRNA40199|XBXL10_1g21682  ccnb1.2.S      CCNB1 0.04860139 0.2376980
    ## 5   XBmRNA1720|XBXL10_1g1100    ccna2.L      CCNA2 0.03471313 0.3327977
    ## 6 XBmRNA20588|XBXL10_1g10943    ccna1.S      CCNA1 0.03404127 0.3393661
    ##   Estimate      T5_kD      T6_kD deltaBIC      T0_L        N1        N3
    ## 1 interval 0.05196905 0.06439057 59.78006 1.7147274 1.0893752 1.1090663
    ## 2    point         NA 0.05236107 19.02400 0.5715232 0.3167322 0.4122321
    ## 3    point         NA 0.05055962 53.52402 0.7593631 0.2741591 0.3558655
    ## 4    point 0.04860139         NA 16.69609 0.8122257 0.5816428 0.5083152
    ## 5    point         NA 0.03471313 50.29648 0.3719585 0.2490067 0.3498724
    ## 6 interval 0.02838204 0.03970050 51.62111 2.1104195 1.5523841 1.7135697
    ##          N4        N7         N9        N11          N12     T0_H        O1
    ## 1 1.2411665 1.8280703 0.92944429 0.06468336 2.346668e-02 9.479349 0.6087783
    ## 2 0.4270137 1.8856945 1.67675998 1.13389467 1.576150e+00 6.415712 0.3157785
    ## 3 0.4436930 1.9537817 1.95445907 0.69835989 1.560319e+00 9.303156 0.5306161
    ## 4 0.3789010 1.0024334 1.50784145 1.53953923 1.669101e+00 6.061506 0.5321921
    ## 5 0.3315833 1.4536714 1.20658237 1.33815842 2.699167e+00 8.538092 1.2239487
    ## 6 1.7205897 0.8069156 0.07777375 0.01829269 5.503508e-05 8.399209 1.3292753
    ##           O2        O3        O4         O5        O7         O9       O10
    ## 1 0.00000000 0.2046160 0.1690556 0.15700874 0.1380140 0.06531466 0.1457002
    ## 2 0.05760827 0.4069183 0.6042256 0.44871175 0.5156761 0.49982006 0.2410755
    ## 3 0.00000000 0.1234161 0.2365997 0.18184335 0.1569108 0.00000000 0.0000000
    ## 4 0.20740221 0.4029011 0.3367262 0.26181593 0.7716691 0.50081417 0.4388957
    ## 5 0.00000000 0.2604006 0.2245385 0.15593400 0.2258017 0.17048824 0.1304815
    ## 6 0.24473731 0.3145507 0.2879517 0.08938614 0.1435368 0.00000000 0.0000000
    ##          O11       O12 mClass Number_Pep_Fits
    ## 1 0.03216316 0.0000000    Deg               2
    ## 2 0.74512284 0.7493515    Deg               1
    ## 3 0.18012809 0.2873295    Deg               1
    ## 4 0.48312743 1.0029505    Deg               1
    ## 5 0.07031518 0.0000000    Deg               1
    ## 6 0.11168100 0.0796721    Deg               3

``` r
max_fit_kd <- log(2)/(tail(avg_heavy_time, 1)*3)
max_fit_hl <- log(2)/max_fit_kd/60

cat("Max Half-life (hrs):\t", max_fit_hl) 
```

    ## Max Half-life (hrs):  212.2

``` r
human_map <- read.csv("Files/Reference/Xen10_to_Human_Name_Assign.csv")
human_map["Protein_ID"] <- sapply(human_map$Protein_ID, function(x){
  split_id <- strsplit(x, "\\|")[[1]]
  reformat <- paste0(split_id[3], "|", split_id[2])

  return(reformat) })

human_map <- human_map[human_map$Protein_ID %in% XLA_fits_df$Protein_ID,]
head(human_map)
```

    ##                    Protein_ID    XLA_Gene           Human_FASTA Human_Gene
    ## 10 XBmRNA18680|XBXL10_1g10005   akirin1.S sp|Q9H9L7|AKIR1_HUMAN    AKIRIN1
    ## 13 XBmRNA18689|XBXL10_1g10009    srsf10.S sp|O75494|SRS10_HUMAN     SRSF10
    ## 21 XBmRNA18700|XBXL10_1g10017    zbtb8b.S sp|Q8NAP8|ZBT8B_HUMAN     ZBTB8B
    ## 23   XBmRNA1559|XBXL10_1g1002  arhgap10.L sp|A1A4S6|RHG10_HUMAN   ARHGAP10
    ## 29 XBmRNA18740|XBXL10_1g10027    thrap3.S sp|Q9Y2W1|TR150_HUMAN     THRAP3
    ## 30 XBmRNA18747|XBXL10_1g10028    map7d1.S sp|Q3KQU3|MA7D1_HUMAN     MAP7D1
    ##                                          Description
    ## 10                                         Akirin-1 
    ## 13          Serine/arginine-rich splicing factor 10 
    ## 21 Zinc finger and BTB domain-containing protein 8B 
    ## 23                 Rho GTPase-activating protein 10 
    ## 29    Thyroid hormone receptor-associated protein 3 
    ## 30                 MAP7 domain-containing protein 1 
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                AA
    ## 10                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    MACGATLKRSMEFEALMSPQSPKRRRCAPLPGSPATPSPQRCAIRPEMQQGQQQPLSQLGGDRRLTPEQILQNIKQEYTRYQRRRQLEGAFNQCEAGALNEVQASCSSLTAPSSPGSLVKKDQPTFSLRQVGILCERLLKDHEDKIREEYEQILNIKLAEQYESFVKFTHDQIMRRYGARPASYVS*
    ## 13                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            MSRYLRPPNSSLFVRNIADDIRSEDLRREFGRYGPIVDVYVPLDYYNRRPRGFAYVQFEDVRDAEDALHNLDKKWICGRQIEIQFAQGDRKTPHQMKAKEGSSTYGSSRYDDDRYNRRSRSRSYERRRSRSRSIEQNYGRSNSPRGGRGAERLRRSRSRSDHGRFNRRNRSRSRSGSNSRSRSKSEPKKTVREQGSGSRSNSRGHSKADSKSRCRESSKYNRESRREEQASKSPSRSVSRSRSKSRSRSWNSHKSSGH*
    ## 21                                                                                                                                                                                                                                                                                                                                                                                                                                             MHSWSNLPVCNNCTRTEMDVNCQSAGRRHEEDGCDELRCCDSSPRARRYFCVGAMEMSSYYTKLLGELNEQRKRDFFCDCSIIVEGRIFKAHRNVLFANSGYFRAMLVHYIQDSGRHSTASLDIVTSEAFSTILDFLYSGKLNVCGENVIEVMSAASYLQMTDVVNFCKAYIHSSLDICRKIEKESSFGQADSGSDGINSGREAELGASPEKENDSDCQKDPPGGDCSSCNSIELMVKHHPTDGSNESKSSNKAVEPKEEFDSDVVEVSEGVQTYHIPTGLEQGEEGLHSASGVDIACNNYHMKQFLEALLRNSSSQRKDDAVHHFPRDFESRQEDTGVPMSSMMDIQGDWFGEDTGDVLVVPIKLHKCPFCPYTAKQKGILKRHIRSHTGERPYPCDVCGKRFTRQEHLRSHALTVHRSNRPIICKGCRRTFTSNLSPGLRRFGLCDSCTCVTTTHDDSNHPGGTEAASESMDKEEEVDGDWPIYIESGDENEGPEEEEIDDKDQIHREVLM*
    ## 23                                                                                                                                                         MGLQPLEFSDCYLDSPWLRERIRAHEAELDRTNRFIKDLLKDGKNLIAATKSLSAAQRKFAHSLRDFKFEFIGDAETDDERCIDASLQEFSTFLQNLEEQREIMALNVNETLVKPLERFRKEQLGGVKELKKKFDKETEKNYSLLDKHLSLSAKKKEPHLQEADVQVEQNRQHFYQLSLDYVCKLQEIQERKKFEFVEPMLSFFQGIFTFYHQGYELANDFNHYKLDLQINIQNTRNLFDGTRSEVEDLLRKIRRNPQEHKRASPFTMEGYLYVQEKRPAPFGSSWVKHYCMYKKDSKSFTMFPYEHRSGGKIGDGDSFSLKSCIKRHTDTIDRRFCFDIEASERPGVLTTMQSFSEDNRCLWMEVLDGKESRFLSLSRSVSRPEGGARLDKYGFAIVKNCIRAIETRGINDQGLYRVVGVSSKVQRLLSLLIDVKTCCDVDLDSSEEWEVKTVTSALKLYLRSLPEPLMTHELHDQFVNPAKSGSPESRVTSIHHLIHQLPEKNREMLDILITHLANVARHAKQNLMTVANLGVVFGPTLMRPQEETVAAIMDLKFQNIVVEILIENHEKIFKNPPSGLEAESDSLSISPPNAPPRQSRRQHRKRPVPVYNLSLELDNADMHLTPREDTPTGSTDSLSPQSLTPTPPHTNPTADERNHITANIGSPADWSASERNLIHSPLVLWINDESPAASAVTTSQAAEATPVKDPEEPLTLPAHLLTTSDTSAESSAKRKAKAVYPCEAEHSSELSFEVGAVFEDVHHSREPGWLEGTLNGRRGLIPENYVHFL*
    ## 29                        MSKNKSKSRSPSSRSMSRSRSRSFSKSRSRSRSLSHSRKRRHSSRSRSRSYSPSYNRERNHPRVYQNRDFRGHNRGYRRPYYFRGRGRGFYPRGQHNRGGYGNYRPNWQNYNRQPYSPRRGRSRSRSPKRRSGSPRSRSRSRNSDKSTSDRSRSSSSRSSSNHSRVEVSKRKSKKEKKSHSKDQRGSAQPVDDDSKEQSGSGGADTGSKTEGSKAWQDVEAYDTSPVPQDSPAERAPTVKSSVQSVVVRRCSPRPSPVQKTSPPSPSVLQSSGAYRASSRQSPFDQSLSPPRKSPLAKSPPLGSLYGNQKEEGRSAEMAAPIGTGYKRFMEEQRNKAAELEKENKDKDKMSSLERMKDRTSPTDFAMNELEKAYRKSQSPKRFKMRDGFEKLKLAELRFAKEEAEQEKKGRSRKDSDSDSKNQDPYDPAKWEELSFIPSAKEKRRKSEDMDDDPYIERPKKEEKSSKRESAHKGFLPEKNFKVTGFKSVKERSESPPVIKSAEVRDKEKVSAKEELAFTKSTFCISREAGPSVRLDSFDEDLARPSGVLAQERKLSRDLVHSNKKDQEFKSIFQHIQSAQPQRSPSELFAQHIVTIVHHVKEHHFGPSSMTLSERFAKYLKRAKDQESSKSKKSPEIHRRIDISPSAFRKHGFLHEEAKHSKETGLKGDGKYRDEPSDLRQDIERRKKHKDKDSKRDHSRDSGDSRDSSRSRDRSSEKTEKSRKSSKKHKKHRKVRERSRSSSSSSHSSHSLRAPGNEEYPAESEEKEDGAGGFDKSRMANKDFQGNSERGRARGTLQFRARGTRPWGRGAFPGNNNNNNNDFQKRNRDEEWDPEYTPKSKKYYLHDDREGENEEKWMNRGRGRANFARGRGRFLFRKPGISPKWAHDKFSGEEGEIEEDESGAENKDEKGNIMSIAE*
    ## 30 MENQAVGNYKEQQQKIPDSLLNGEAEDVPQNFTETPTCNTEPRLQADQQKIRTETIAPTSVAHSHSSPPQTDGLSSAMADNGATFPLHATGQKSVPAEQKPVSPSQTIAPCFTIVEGPRLSPEGQSSSPQTEGITPFSEALGSPKLDQNIKYDGSSPPAVSSPISAVSPRQKPDTQKAEQRQKQAKERREERAKYLAAKKSVWLEKEEKAKQLREKQLQDRRKKLEEQRLKAEKKRVLLEERQRVKLEKNKERYESAVNRNSKKTWAEIRQQRWSWAGAFNPNSPGHKGGANRCSVSAVNLPRHVDSIINKRLSKSSATLWNTPSRNRSLQLSPWESSIVDRLMTPTLSFLARSRSVATLPGNGRDTNSHVCPRSASASPLTPCNGHKHRCPERRRTIGSSLDVTPKKRSDSSPKKKEKKEKDRENEREKNALCRERVLKKRQTLPSTKTRLSTTPSEIKSLKSKNRPSSPSTPTRRPASPSPAASISLTPKLTSPKSTQGAVKIKPKTEKSKEEKTPSKIKEKKEELKIEKEQQSSSTPKSEEPTEQKVTECTEQSTDDDVTSAAPTAASPKPALLPVPVVLVTPSPEPVFSKPMAAAPAAVIAPSTTPSTSKPVVSTPAAAIASSPAPSPSKSIAPTPVAVIPPSPAPSPGKPMAGTTDKEEAARLLSEKRRQAREQREREEQERREREEQKRREREEKARKEAEERLQREEEARRLEEQQKAEELERRQNEEEQQRKEEERKAKSEQEEMERLQKQREEAEAKAREEAEKQRLEREQHFQKEERERLERKKRLEEIMKRTRKSDATDKKKEDKQALNGKETNLESIPSKDEAKSSEDIIKIEHPEKKTMWLQDRLGTSEPSQGVTGGILNGVQPPKQENGFSSKDSSSVFEEVITLAERSGSNGEKGIPPADPIIAFSGKDPFSKKSSVQPHQVTEVL*
    ##       E_LtoH    E_HtoL Human_KB
    ## 10  3.89e-81  1.67e-82   Q9H9L7
    ## 13  3.47e-89        NA   O75494
    ## 21  0.00e+00        NA   Q8NAP8
    ## 23  0.00e+00  0.00e+00   A1A4S6
    ## 29  0.00e+00  0.00e+00   Q9Y2W1
    ## 30 3.40e-115 9.87e-132   Q3KQU3
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          GO_BP
    ## 10                                 myoblast migration involved in skeletal muscle regeneration [GO:0014839]; negative regulation of satellite cell differentiation [GO:1902725]; negative regulation of skeletal muscle satellite cell proliferation [GO:1902723]; positive regulation of lamellipodium assembly [GO:0010592]; positive regulation of macrophage chemotaxis [GO:0010759]; positive regulation of myoblast differentiation [GO:0045663]; positive regulation of transcription by RNA polymerase II [GO:0045944]
    ## 13                                                                                          cytosolic transport [GO:0016482]; mRNA splice site selection [GO:0006376]; mRNA splicing, via spliceosome [GO:0000398]; negative regulation of mRNA splicing, via spliceosome [GO:0048025]; regulation of mRNA splicing, via spliceosome [GO:0048024]; regulation of transcription, DNA-templated [GO:0006355]; RNA splicing, via transesterification reactions [GO:0000375]; spliceosomal tri-snRNP complex assembly [GO:0000244]
    ## 21                                                                                                                                                                                                                                                                                                                                                                                                                                                               regulation of transcription by RNA polymerase II [GO:0006357]
    ## 23                                                                                                                                                                                                                                                                                                                       cytoskeleton organization [GO:0007010]; negative regulation of apoptotic process [GO:0043066]; regulation of small GTPase mediated signal transduction [GO:0051056]; signal transduction [GO:0007165]
    ## 29 circadian rhythm [GO:0007623]; mRNA processing [GO:0006397]; mRNA stabilization [GO:0048255]; nuclear-transcribed mRNA catabolic process [GO:0000956]; positive regulation of circadian rhythm [GO:0042753]; positive regulation of mRNA splicing, via spliceosome [GO:0048026]; positive regulation of transcription by RNA polymerase II [GO:0045944]; positive regulation of transcription, DNA-templated [GO:0045893]; regulation of alternative mRNA splicing, via spliceosome [GO:0000381]; RNA splicing [GO:0008380]
    ## 30                                                                                                                                                                                                                                                                                                                                                                                                                                                                          microtubule cytoskeleton organization [GO:0000226]
    ##                                                                                                                                                                                                           GO_CC
    ## 10                                                                                                        chromatin [GO:0000785]; nuclear membrane [GO:0031965]; nucleoplasm [GO:0005654]; nucleus [GO:0005634]
    ## 13 axon terminus [GO:0043679]; cytoplasm [GO:0005737]; cytosol [GO:0005829]; dendrite [GO:0030425]; neuronal cell body [GO:0043025]; nuclear speck [GO:0016607]; nucleoplasm [GO:0005654]; nucleus [GO:0005634]
    ## 21                                                                                                                                                                 chromatin [GO:0000785]; nucleus [GO:0005634]
    ## 23                                                                                                             cytosol [GO:0005829]; perinuclear region of cytoplasm [GO:0048471]; plasma membrane [GO:0005886]
    ## 29                                                                extracellular exosome [GO:0070062]; mediator complex [GO:0016592]; nuclear speck [GO:0016607]; nucleoplasm [GO:0005654]; nucleus [GO:0005634]
    ## 30                                                                                                                          cytoplasm [GO:0005737]; microtubule cytoskeleton [GO:0015630]; spindle [GO:0005819]
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                  GO_MF
    ## 10                                                                                                                                                                                                                                                                                                                                                                                                     transcription coregulator activity [GO:0003712]
    ## 13                                                                                                                                                                                                                                                                                                                                                     RNA binding [GO:0003723]; RS domain binding [GO:0050733]; unfolded protein binding [GO:0051082]
    ## 21                                                                                                                                                                                                                                    DNA-binding transcription factor activity, RNA polymerase II-specific [GO:0000981]; metal ion binding [GO:0046872]; RNA polymerase II transcription regulatory region sequence-specific DNA binding [GO:0000977]
    ## 23                                                                                                                                                                                                                                                                                                                                                                                                              GTPase activator activity [GO:0005096]
    ## 29 ATP binding [GO:0005524]; DNA binding [GO:0003677]; nuclear receptor coactivator activity [GO:0030374]; phosphoprotein binding [GO:0051219]; RNA binding [GO:0003723]; RNA polymerase II cis-regulatory region sequence-specific DNA binding [GO:0000978]; thyroid hormone receptor binding [GO:0046966]; transcription coactivator activity [GO:0003713]; transcription coregulator activity [GO:0003712]; vitamin D receptor binding [GO:0042809]
    ## 30                                                                                                                                                                                                                                                                                                                                                                                                                                                 N/A
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             GO_All
    ## 10                                                                                                                                                                                                                                                                                                                                                                                                                                                                             chromatin [GO:0000785]; nuclear membrane [GO:0031965]; nucleoplasm [GO:0005654]; nucleus [GO:0005634]; transcription coregulator activity [GO:0003712]; myoblast migration involved in skeletal muscle regeneration [GO:0014839]; negative regulation of satellite cell differentiation [GO:1902725]; negative regulation of skeletal muscle satellite cell proliferation [GO:1902723]; positive regulation of lamellipodium assembly [GO:0010592]; positive regulation of macrophage chemotaxis [GO:0010759]; positive regulation of myoblast differentiation [GO:0045663]; positive regulation of transcription by RNA polymerase II [GO:0045944]
    ## 13                                                                                                                                                                                                                                                                                                                                                                               axon terminus [GO:0043679]; cytoplasm [GO:0005737]; cytosol [GO:0005829]; dendrite [GO:0030425]; neuronal cell body [GO:0043025]; nuclear speck [GO:0016607]; nucleoplasm [GO:0005654]; nucleus [GO:0005634]; RNA binding [GO:0003723]; RS domain binding [GO:0050733]; unfolded protein binding [GO:0051082]; cytosolic transport [GO:0016482]; mRNA splice site selection [GO:0006376]; mRNA splicing, via spliceosome [GO:0000398]; negative regulation of mRNA splicing, via spliceosome [GO:0048025]; regulation of mRNA splicing, via spliceosome [GO:0048024]; regulation of transcription, DNA-templated [GO:0006355]; RNA splicing, via transesterification reactions [GO:0000375]; spliceosomal tri-snRNP complex assembly [GO:0000244]
    ## 21                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   chromatin [GO:0000785]; nucleus [GO:0005634]; DNA-binding transcription factor activity, RNA polymerase II-specific [GO:0000981]; metal ion binding [GO:0046872]; RNA polymerase II transcription regulatory region sequence-specific DNA binding [GO:0000977]; regulation of transcription by RNA polymerase II [GO:0006357]
    ## 23                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 cytosol [GO:0005829]; perinuclear region of cytoplasm [GO:0048471]; plasma membrane [GO:0005886]; GTPase activator activity [GO:0005096]; cytoskeleton organization [GO:0007010]; negative regulation of apoptotic process [GO:0043066]; regulation of small GTPase mediated signal transduction [GO:0051056]; signal transduction [GO:0007165]
    ## 29 extracellular exosome [GO:0070062]; mediator complex [GO:0016592]; nuclear speck [GO:0016607]; nucleoplasm [GO:0005654]; nucleus [GO:0005634]; ATP binding [GO:0005524]; DNA binding [GO:0003677]; nuclear receptor coactivator activity [GO:0030374]; phosphoprotein binding [GO:0051219]; RNA binding [GO:0003723]; RNA polymerase II cis-regulatory region sequence-specific DNA binding [GO:0000978]; thyroid hormone receptor binding [GO:0046966]; transcription coactivator activity [GO:0003713]; transcription coregulator activity [GO:0003712]; vitamin D receptor binding [GO:0042809]; circadian rhythm [GO:0007623]; mRNA processing [GO:0006397]; mRNA stabilization [GO:0048255]; nuclear-transcribed mRNA catabolic process [GO:0000956]; positive regulation of circadian rhythm [GO:0042753]; positive regulation of mRNA splicing, via spliceosome [GO:0048026]; positive regulation of transcription by RNA polymerase II [GO:0045944]; positive regulation of transcription, DNA-templated [GO:0045893]; regulation of alternative mRNA splicing, via spliceosome [GO:0000381]; RNA splicing [GO:0008380]
    ## 30                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         cytoplasm [GO:0005737]; microtubule cytoskeleton [GO:0015630]; spindle [GO:0005819]; microtubule cytoskeleton organization [GO:0000226]
    ##                                                                                                                                                                                                                                                                                                        GO_IDs
    ## 10                                                                                                                                                             GO:0000785; GO:0003712; GO:0005634; GO:0005654; GO:0010592; GO:0010759; GO:0014839; GO:0031965; GO:0045663; GO:0045944; GO:1902723; GO:1902725
    ## 13                                                                         GO:0000244; GO:0000375; GO:0000398; GO:0003723; GO:0005634; GO:0005654; GO:0005737; GO:0005829; GO:0006355; GO:0006376; GO:0016482; GO:0016607; GO:0030425; GO:0043025; GO:0043679; GO:0048024; GO:0048025; GO:0050733; GO:0051082
    ## 21                                                                                                                                                                                                                                     GO:0000785; GO:0000977; GO:0000981; GO:0005634; GO:0006357; GO:0046872
    ## 23                                                                                                                                                                                                             GO:0005096; GO:0005829; GO:0005886; GO:0007010; GO:0007165; GO:0043066; GO:0048471; GO:0051056
    ## 29 GO:0000381; GO:0000956; GO:0000978; GO:0003677; GO:0003712; GO:0003713; GO:0003723; GO:0005524; GO:0005634; GO:0005654; GO:0006397; GO:0007623; GO:0008380; GO:0016592; GO:0016607; GO:0030374; GO:0042753; GO:0042809; GO:0045893; GO:0045944; GO:0046966; GO:0048026; GO:0048255; GO:0051219; GO:0070062
    ## 30                                                                                                                                                                                                                                                             GO:0000226; GO:0005737; GO:0005819; GO:0015630

``` r
enrich_kd <- XLA_fits_df[XLA_fits_df$Protein_ID %in% human_map$Protein_ID,
                         c("Protein_ID", "kD", "O12")]
enrich_kd <- enrich_kd %>% arrange(kD, desc(O12))

enrich_kd["Rank"] <- seq(1, nrow(enrich_kd))
enrich_kd["Signed_Rank"] <- enrich_kd$Rank - mean(enrich_kd$Rank)
enrich_kd["Signed_Rank"] <- enrich_kd$Signed_Rank / max(enrich_kd$Signed_Rank)
enrich_kd <- enrich_kd %>%
  select(Signed_Rank, Protein_ID, kD) %>%
  arrange(desc(Signed_Rank))

protein_stats <- setNames(enrich_kd$Signed_Rank, enrich_kd$Protein_ID)

head(protein_stats)
```

    ##  XBmRNA10569|XBXL10_1g5638 XBmRNA35507|XBXL10_1g19214 
    ##                  1.0000000                  0.9997659 
    ## XBmRNA34126|XBXL10_1g18541 XBmRNA40199|XBXL10_1g21682 
    ##                  0.9995318                  0.9992978 
    ##   XBmRNA1720|XBXL10_1g1100 XBmRNA20588|XBXL10_1g10943 
    ##                  0.9990637                  0.9988296

``` r
cat("\n\n")
```

``` r
tail(protein_stats)
```

    ## XBmRNA75660|XBXL10_1g40210 XBmRNA32250|XBXL10_1g17494 
    ##                 -0.9988296                 -0.9990637 
    ## XBmRNA36209|XBXL10_1g19569 XBmRNA34226|XBXL10_1g18592 
    ##                 -0.9992978                 -0.9995318 
    ## XBmRNA81187|XBXL10_1g43183 XBmRNA61399|XBXL10_1g32545 
    ##                 -0.9997659                 -1.0000000

``` r
human_GO <- human_map
human_GO <- apply(human_GO, 1, function(x){

  BP_terms <- trimws(strsplit(x["GO_BP"], ";")[[1]])
  BP_id <- sub(".*\\[(GO:\\d+)\\].*", "\\1", BP_terms)
  BP_desc <- sub("\\s*\\[GO:.*\\]", "", BP_terms)

  MF_terms <- trimws(strsplit(x["GO_MF"], ";")[[1]])
  MF_id <- sub(".*\\[(GO:\\d+)\\].*", "\\1", MF_terms)
  MF_desc <- sub("\\s*\\[GO:.*\\]", "", MF_terms)

  CC_terms <- trimws(strsplit(x["GO_CC"], ";")[[1]])
  CC_id <- sub(".*\\[(GO:\\d+)\\].*", "\\1", CC_terms)
  CC_desc <- sub("\\s*\\[GO:.*\\]", "", CC_terms)


  out_df <- data.frame(Protein_ID=rep(x["Protein_ID"],
                                      length(c(BP_id, MF_id, CC_id))),
                       ID=c(BP_id, MF_id, CC_id),
                       Desc=c(BP_desc, MF_desc, CC_desc),
                       Type=c(rep("BP", length(BP_id)),
                              rep("MF", length(MF_id)),
                              rep("CC", length(CC_id))))

  return(out_df) }) %>% bind_rows

#This dataframe is solely to count terms and remove <10 annotations
human_GO_counts <- data.frame(table(human_GO$ID))
human_GO_counts <- human_GO_counts[human_GO_counts$Freq > 20,]

#df of GO terms for identified and mapped proteins
human_GO_annot <- unique(human_GO[c("ID", "Type", "Desc")])
human_GO_annot <- human_GO_annot[human_GO_annot$ID %in% human_GO_counts$Var1,]
row.names(human_GO_annot) <- NULL

head(human_GO_annot)
```

    ##           ID Type                                                      Desc
    ## 1 GO:0010592   BP             positive regulation of lamellipodium assembly
    ## 2 GO:0045663   BP           positive regulation of myoblast differentiation
    ## 3 GO:0045944   BP positive regulation of transcription by RNA polymerase II
    ## 4 GO:0003712   MF                        transcription coregulator activity
    ## 5 GO:0000785   CC                                                 chromatin
    ## 6 GO:0031965   CC                                          nuclear membrane

``` r
filtered_human_GO <- human_GO
filtered_human_GO <- filtered_human_GO[filtered_human_GO$ID %in% human_GO_annot$ID,]

human_GO_ref <- lapply(unique(filtered_human_GO$ID),
       function(x){

         term <- human_GO[human_GO$ID == x,]$Protein_ID

         return(term) })

names(human_GO_ref) <- unique(filtered_human_GO$ID)

human_GO_ref[1]
```

    ## $`GO:0010592`
    ##  [1] "XBmRNA18680|XBXL10_1g10005" "XBmRNA19567|XBXL10_1g10465"
    ##  [3] "XBmRNA26742|XBXL10_1g14180" "XBmRNA31760|XBXL10_1g17133"
    ##  [5] "XBmRNA36775|XBXL10_1g19878" "XBmRNA41236|XBXL10_1g22177"
    ##  [7] "XBmRNA44838|XBXL10_1g23972" "XBmRNA46308|XBXL10_1g24809"
    ##  [9] "XBmRNA5537|XBXL10_1g3061"   "XBmRNA58364|XBXL10_1g30904"
    ## [11] "XBmRNA59208|XBXL10_1g31302" "XBmRNA5677|XBXL10_1g3138"  
    ## [13] "XBmRNA75906|XBXL10_1g40372" "XBmRNA76142|XBXL10_1g40474"
    ## [15] "XBmRNA76762|XBXL10_1g40761" "XBmRNA77331|XBXL10_1g40994"
    ## [17] "XBmRNA81312|XBXL10_1g43254" "XBmRNA81462|XBXL10_1g43328"
    ## [19] "XBmRNA81967|XBXL10_1g43561" "XBmRNA82324|XBXL10_1g43733"
    ## [21] "XBmRNA82626|XBXL10_1g43906" "XBmRNA10457|XBXL10_1g5579" 
    ## [23] "XBmRNA10576|XBXL10_1g5642"  "XBmRNA12953|XBXL10_1g6920"

``` r
# Run enrichment
set.seed(123)   # reproducible fgsea permutations
res <- fgseaMultilevel(
  pathways = human_GO_ref, #Human GO
  stats    = protein_stats,
  minSize  = 10,
  maxSize  = 500)

#Order your results by significance
res <- res[order(res$padj), ]

# fgsea's built-in collapsing function
# Compares the leading edge genes of nested pathways and drops redundant ones
collapsed_pathways <- collapsePathways(res, human_GO_ref, protein_stats)

#Filter main results to keep only the 'parent' representative terms
res <- res[res$pathway %in% collapsed_pathways$mainPathways, ]
res <- res[res$padj < 0.05,] %>% as.data.frame
res <- merge(res, human_GO_annot, by.x="pathway",  by.y="ID")
res <- res %>% arrange(desc(NES))
res <- res[c("pathway", "NES", "padj", "size", "Type", "Desc")]

head(res)
```

    ##      pathway      NES         padj size Type
    ## 1 GO:0051301 3.021002 4.029105e-26  293   BP
    ## 2 GO:0000981 2.964091 1.785756e-14  104   MF
    ## 3 GO:0001228 2.915806 2.729617e-10   54   MF
    ## 4 GO:0007052 2.829938 3.551169e-10   62   BP
    ## 5 GO:0003700 2.805185 2.085899e-09   56   MF
    ## 6 GO:0000978 2.713094 9.838673e-12  135   MF
    ##                                                                       Desc
    ## 1                                                            cell division
    ## 2    DNA-binding transcription factor activity, RNA polymerase II-specific
    ## 3 DNA-binding transcription activator activity, RNA polymerase II-specific
    ## 4                                             mitotic spindle organization
    ## 5                                DNA-binding transcription factor activity
    ## 6    RNA polymerase II cis-regulatory region sequence-specific DNA binding

``` r
custom_GO_terms <- human_GO_annot[human_GO_annot$ID %in% 
  c("GO:0051301", "GO:0000981", "GO:0004712", "GO:0061630", "GO:0098609",
    "GO:0005759", "GO:0006457", "GO:0005840", "GO:0006635", "GO:0000502",
    "GO:0008017"),]
custom_GO_terms <- custom_GO_terms[order(custom_GO_terms$ID),]

custom_GO_terms["Custom_Name"] <-
  c("Proteasome complex", "Transcription factor activity", "Protein kinase activity",
    "Mitochondrial matrix", "Ribosome", "Protein folding", "Fatty acid β-oxidation",
    "Microtubule binding",
    "Cell division", "Ubiquitin protein ligase activity", "Cell-cell adhesion")
custom_GO_terms <- custom_GO_terms[-c(2,3)]
custom_GO_terms <- merge(res, custom_GO_terms,
                         by.x="pathway", by.y="ID") %>% arrange(desc(NES))

custom_GO_terms
```

    ##       pathway       NES         padj size Type
    ## 1  GO:0051301  3.021002 4.029105e-26  293   BP
    ## 2  GO:0000981  2.964091 1.785756e-14  104   MF
    ## 3  GO:0098609  2.676564 1.939590e-08   55   BP
    ## 4  GO:0008017  2.308142 5.512485e-08  159   MF
    ## 5  GO:0061630  2.171335 1.997782e-06  141   MF
    ## 6  GO:0004712  2.034145 8.925974e-07  221   MF
    ## 7  GO:0005840 -1.892944 3.233633e-03   54   CC
    ## 8  GO:0006457 -2.219182 9.956499e-07  154   BP
    ## 9  GO:0000502 -2.556083 5.198599e-08   79   CC
    ## 10 GO:0006635 -2.603559 5.930066e-07   48   BP
    ## 11 GO:0005759 -3.069076 3.984653e-34  397   CC
    ##                                                                     Desc
    ## 1                                                          cell division
    ## 2  DNA-binding transcription factor activity, RNA polymerase II-specific
    ## 3                                                     cell-cell adhesion
    ## 4                                                    microtubule binding
    ## 5                                      ubiquitin protein ligase activity
    ## 6                      protein serine/threonine/tyrosine kinase activity
    ## 7                                                               ribosome
    ## 8                                                        protein folding
    ## 9                                                     proteasome complex
    ## 10                                             fatty acid beta-oxidation
    ## 11                                                  mitochondrial matrix
    ##                          Custom_Name
    ## 1                      Cell division
    ## 2      Transcription factor activity
    ## 3                 Cell-cell adhesion
    ## 4                Microtubule binding
    ## 5  Ubiquitin protein ligase activity
    ## 6            Protein kinase activity
    ## 7                           Ribosome
    ## 8                    Protein folding
    ## 9                 Proteasome complex
    ## 10            Fatty acid β-oxidation
    ## 11              Mitochondrial matrix

``` r
p1 <- ggplot(custom_GO_terms, aes(x = NES, size = size,
                                  y = reorder(Custom_Name, NES))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(alpha = 1) +
  geom_point(shape = 1, color = "black", stroke=0.8) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 50)) +  # wrap for display only
  scale_color_continuous(low = "grey80", high = "grey30", trans = "reverse") +
  labs(
    x = "Normalized enrichment score",
    y = "GO term",
    size = "Set size"
  ) +
  theme_bw() +
  theme(
    legend.text = element_text(size = 12, colour = "black"),
    legend.title = element_text(size = 12, colour = "black"),    
    axis.text = element_text(size = 13, colour = "black"),    
    plot.title = element_text(size = 18, hjust = 0.5),
    axis.title.x = element_text(size=16),
    axis.title.y = element_blank(),    
    panel.border = element_rect(linewidth = 1),
    plot.margin = margin(t = 5, r = 10, b = 15, l = 7, unit = "pt"),
    legend.position = c(0.25,0.75),
    legend.background = element_rect(colour = "black", linewidth = 0.5,
                                     fill = "white")
  ) +
  coord_cartesian(xlim=c(-3.25,3.25))

# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=6, height=4, res=300)

p1
```

![](figures/Frog_Early_Analysis/GO%20term%20plot-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
XLA_deg_proteins <- XLA_fits_df$Protein_ID

interproscan_XLA <- 
  read.csv("Systems/Frog_Interproscan/XLA_10p1_longest-transcript_interproscan.tsv",
           sep="\t")
interproscan_XLA <- interproscan_XLA[interproscan_XLA$Protein_ID %in%
                                       XLA_deg_proteins,]
interproscan_XLA <- unique(interproscan_XLA[c("Protein_ID", "Signature_Accession",
                                              "Signature_Description")])

sub_IP_XLA <- data.frame(table(interproscan_XLA$Signature_Accession)) %>%
  arrange(Freq)
sub_IP_XLA <- sub_IP_XLA[sub_IP_XLA$Freq>=5,]$Var1
sub_IP_XLA <- sapply(sub_IP_XLA, as.character)

length(sub_IP_XLA)
```

    ## [1] 722

``` r
sub_IP_XLA[1:10]
```

    ##  [1] "PF00026" "PF00053" "PF00055" "PF00058" "PF00084" "PF00094" "PF00117"
    ##  [8] "PF00156" "PF00200" "PF00210"

``` r
domain_sub <- XLA_fits_df[c("Protein_ID", "XLA_Gene", "Human_Gene",
                            "kD", "HL_Hrs", "mClass")]
domain_sub <- merge(domain_sub,
                    interproscan_XLA[c("Protein_ID", "Signature_Accession",
                                       "Signature_Description")],
                    by="Protein_ID") %>% arrange(HL_Hrs)
# domain_sub <-  domain_sub[domain_sub$Signature_Accession %in%
#                             sub_IP_XLA,] %>% arrange(desc(Avg_kD))

head(domain_sub)
```

    ##                   Protein_ID   XLA_Gene Human_Gene         kD    HL_Hrs mClass
    ## 1  XBmRNA10569|XBXL10_1g5638    ccnb1.S      CCNB1 0.05817981 0.1985646    Deg
    ## 2  XBmRNA10569|XBXL10_1g5638    ccnb1.S      CCNB1 0.05817981 0.1985646    Deg
    ## 3 XBmRNA35507|XBXL10_1g19214  ccnb1.2.L      CCNB1 0.05236107 0.2206306    Deg
    ## 4 XBmRNA35507|XBXL10_1g19214  ccnb1.2.L      CCNB1 0.05236107 0.2206306    Deg
    ## 5 XBmRNA34126|XBXL10_1g18541  LOC398156      PTTG1 0.05055962 0.2284917    Deg
    ## 6 XBmRNA40199|XBXL10_1g21682  ccnb1.2.S      CCNB1 0.04860139 0.2376980    Deg
    ##   Signature_Accession                         Signature_Description
    ## 1             PF00134                     Cyclin, N-terminal domain
    ## 2             PF02984                     Cyclin, C-terminal domain
    ## 3             PF02984                     Cyclin, C-terminal domain
    ## 4             PF00134                     Cyclin, N-terminal domain
    ## 5             PF04856 Securin sister-chromatid separation inhibitor
    ## 6             PF02984                     Cyclin, C-terminal domain

``` r
pfam_domains <- c("PF00225", "PF01437", "PF00096", "PF00595", "PF00567",
                  "PF02984", "PF00028", "PF00176", "PF02207", "PF00439",
                  "PF13445", "PF00632", "PF00069", "PF07525", "PF00046", 
                  "PF00022", "PF00472", "PF00191", "PF01399", "PF00012",
                  "PF00227", "PF00118", "PF00106")

names(pfam_domains) <-
  c("Kinesin motor", "Plexin repeat", "Znf-C2H2", "PDZ", "Tudor",
    "Cyclin-Cterm", "Cadherin", "SNF2-related", "UBR box", "Bromo",
    "RING-Znf", "HECT", "Protein kinase", "SOCS box", "Homeo",
    "Actin", "RF-1", "Annexin", "PCI", "Hsp70",
    "Proteasome subunit", "TCP1-cnp60 chaperonin", "Short chain dehydrogenase")
```

``` r
domain_sub[domain_sub$Signature_Accession == "PF00028",]
```

    ##                      Protein_ID      XLA_Gene Human_Gene           kD   HL_Hrs
    ## 427  XBmRNA22763|XBXL10_1g12166       pcdh1.L      PCDH1 0.0006215640 18.58610
    ## 542  XBmRNA31101|XBXL10_1g16561       pcdh1.S      PCDH1 0.0005495646 21.02110
    ## 655  XBmRNA34590|XBXL10_1g18758        cdh3.L       CDH1 0.0004871643 23.71367
    ## 793  XBmRNA52360|XBXL10_1g27887        dsc3.L       DSC1 0.0004386308 26.33753
    ## 1078 XBmRNA39444|XBXL10_1g21313        cdh3.S       CDH1 0.0003320322 34.79317
    ## 1484 XBmRNA55520|XBXL10_1g29508        dsg2.S       DSG2 0.0002603570 44.37159
    ## 1527 XBmRNA55517|XBXL10_1g29505        dsc3.S       DSC2 0.0002529830 45.66494
    ## 3079    XBmRNA1264|XBXL10_1g860  LOC108719599       FAT1 0.0001328455 86.96156
    ##      mClass Signature_Accession Signature_Description
    ## 427     Deg             PF00028       Cadherin domain
    ## 542     Deg             PF00028       Cadherin domain
    ## 655     Deg             PF00028       Cadherin domain
    ## 793     Deg             PF00028       Cadherin domain
    ## 1078    Deg             PF00028       Cadherin domain
    ## 1484    Deg             PF00028       Cadherin domain
    ## 1527   Flat             PF00028       Cadherin domain
    ## 3079   Flat             PF00028       Cadherin domain

``` r
domain_all_data <- lapply(pfam_domains, function(pfam){

  sub_df <- domain_sub[domain_sub$Signature_Accession == pfam,]
  sub_df$Signature_Description <- rep(names(pfam_domains[pfam_domains==pfam]),
                                      nrow(sub_df))

  return(sub_df)}) %>% bind_rows
domain_all_data <- domain_all_data[domain_all_data$Signature_Accession %in% sub_IP_XLA,]

domain_all_data <- domain_all_data[!domain_all_data$Signature_Description %in%
    c("Actin", "RF-1", "Annexin", "Hsp70", "Proteasome subunit",
      "TCP1-cnp60 chaperonin", "Short chain dehydrogenase"),]

head(domain_all_data)
```

    ##                   Protein_ID  XLA_Gene Human_Gene           kD     HL_Hrs
    ## 1 XBmRNA79044|XBXL10_1g42043   kif22.L      KIF22 0.0142095807  0.8130045
    ## 2 XBmRNA83366|XBXL10_1g44357   kif22.S      KIF22 0.0113377716  1.0189351
    ## 3 XBmRNA65420|XBXL10_1g34623  kif14l.L      KIF14 0.0015913238  7.2596494
    ## 4 XBmRNA23190|XBXL10_1g12407  kif20a.L     KIF20A 0.0014123316  8.1797031
    ## 5 XBmRNA25526|XBXL10_1g13446   kif23.L      KIF23 0.0011160482 10.3512136
    ## 6   XBmRNA7535|XBXL10_1g4091       N/A        N/A 0.0008141075 14.1903293
    ##   mClass Signature_Accession Signature_Description
    ## 1    Deg             PF00225         Kinesin motor
    ## 2    Deg             PF00225         Kinesin motor
    ## 3    Deg             PF00225         Kinesin motor
    ## 4    Deg             PF00225         Kinesin motor
    ## 5    Deg             PF00225         Kinesin motor
    ## 6    Deg             PF00225         Kinesin motor

``` r
domain_all_data["Signature_Description"] <-
  factor(domain_all_data$Signature_Description,
         levels=c("PCI", "Protein kinase", "HECT", "Bromo",
                  "PDZ", "Homeo", "SNF2-related", "RING-Znf",
                  "Znf-C2H2", "Tudor", "UBR box", "Kinesin motor",
                  "SOCS box", "Cadherin", "Plexin repeat", "Cyclin-Cterm"))

deg_domains <- ggplot(domain_all_data, aes(x = Signature_Description,
                                           y = HL_Hrs)) +

  # geom_hline(yintercept = max_fit_hl, linetype = "dashed",
  #            color = "black", linewidth=1) +
  geom_boxplot(outliers = FALSE, fill = "grey90", color = "black", alpha = 0.7) +
  stat_summary(fun = median, geom = "point", 
               shape = 21, size = 3, fill = "red3", color = "black") +
  labs(title = "Protein Half-Lives by InterPro Domain",
       x = "Protein domain", y = "Half-life (hours)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15,colour="black",
                                   angle=35, hjust=1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=16,colour="black"),
        axis.title.y = element_text(size=16),
        plot.title = element_blank(),
        panel.border = element_rect(linewidth=1.1),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
        legend.position = "none",
        plot.margin = margin(t = 10, r = 25, b = 10,
                             l = 10, unit = "pt")) 

# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=9, height=4.25, res=300)

deg_domains
```

![](figures/Frog_Early_Analysis/Domain%20plot-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
XLA_iupred <- read.csv("Systems/iupred2a/XLA_10p1_IUPred2A-results.csv")

XLA_idr_res <- XLA_fits_df[c("Protein_ID", "XLA_Gene", "Human_Gene",
                             "kD", "HL_Hrs", "mClass")]
XLA_idr_res <- merge(XLA_idr_res, XLA_iupred, by="Protein_ID") %>%
  arrange(desc(kD))

XLA_idr_res["IDR_Class"] <- "N/A"
XLA_idr_res[XLA_idr_res$Mean_IUPred>=0.5,"IDR_Class"] <- "Disordered"
XLA_idr_res[XLA_idr_res$Mean_IUPred<0.5,"IDR_Class"] <- "Ordered"

# --- Are disordered proteins more likely to be degrading AT ALL? ------------
XLA_idr_res$IDR_Class <- factor(XLA_idr_res$IDR_Class,
                                levels = c("Ordered", "Disordered"))

prop_tab <- table(IDR_Class = XLA_idr_res$IDR_Class,
                  Degrading = XLA_idr_res$mClass)

# chi-square (Fisher is needlessly slow at this n; result is identical in spirit)
prop_test <- chisq.test(prop_tab)

print(prop_test)
```

    ## 
    ##  Pearson's Chi-squared test with Yates' continuity correction
    ## 
    ## data:  prop_tab
    ## X-squared = 76.823, df = 1, p-value < 2.2e-16

``` r
cat("\n\n\n")
```

``` r
# --- AMONG degraders, do disordered ones degrade faster? ------------

deg <- XLA_idr_res %>% filter(mClass=="Deg") 

rate_test <- t.test(log10(kD) ~ IDR_Class, data = deg)
print(rate_test)
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  log10(kD) by IDR_Class
    ## t = 0.10853, df = 739, p-value = 0.9136
    ## alternative hypothesis: true difference in means between group Ordered and group Disordered is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.04102863  0.04583034
    ## sample estimates:
    ##    mean in group Ordered mean in group Disordered 
    ##                -3.807253                -3.809654

``` r
IDR_boxplot_df <- cbind(XLA_idr_res[c("IDR_Class", "HL_Hrs")],
                        data.frame(Set=rep("All", nrow(XLA_idr_res))))
IDR_boxplot_df <- rbind(IDR_boxplot_df,
                        cbind(deg[c("IDR_Class", "HL_Hrs")],
                              data.frame(Set=rep("Degrading", nrow(deg)))))

p1 <- ggplot(IDR_boxplot_df,
             aes(x = IDR_Class, y = HL_Hrs, fill = Set,
                 group=interaction(IDR_Class, Set))) +
  geom_boxplot(width = 0.5, color = "black",
               outliers = FALSE, alpha=0.7,
               position=position_dodge(width = 0.6)) +
  stat_summary(fun = median, geom = "point",
               shape = 21, size = 3, fill = "red3", color = "black",
               position = position_dodge(width = 0.6)) +
  scale_fill_manual(values = c(All = "#D55E00", Degrading = "#5E4FA2")) +
  theme_bw() +
  theme(axis.text.y = element_text(size=16, colour="black"),
        axis.text.x = element_text(size = 15, colour="black"),
        axis.title.y = element_text(size=16),
        axis.title.x = element_blank(),
        plot.title = element_blank(),
        panel.border = element_rect(linewidth=1.2),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
        legend.position = "none",
        plot.margin = margin(t = 0, r = 25, b = 0, l = 10, unit = "pt")) +
  labs(y=expression("Half-life (hours)")) +
  coord_cartesian(ylim=c(0,210)) +
  scale_y_continuous(breaks=c(0, 50, 100, 150, 200))
```

``` r
dist_idr_df <- 
  rbind(data.frame(table(XLA_idr_res[XLA_idr_res$IDR_Class=="Ordered",]$mClass)),
        data.frame(table(XLA_idr_res[XLA_idr_res$IDR_Class=="Disordered",]$mClass)))
dist_idr_df["IDR"] <- c("Ordered", "Ordered", "Disordered", "Disordered")
dist_idr_df[1:2, "Freq"] <- dist_idr_df$Freq[1:2] / 7617
dist_idr_df[3:4, "Freq"] <- dist_idr_df$Freq[3:4] / 1187
dist_idr_df$IDR <- factor(dist_idr_df$IDR, levels=c("Ordered", "Disordered"))
dist_idr_df$Var1 <- factor(dist_idr_df$Var1, levels=c("Flat", "Deg"))

dist_idr_df
```

    ##   Var1      Freq        IDR
    ## 1  Deg 0.2753052    Ordered
    ## 2 Flat 0.7246948    Ordered
    ## 3  Deg 0.4001685 Disordered
    ## 4 Flat 0.5998315 Disordered

``` r
p2 <- ggplot(dist_idr_df, aes(x = IDR, y = Freq*100, fill = Var1)) +
  geom_col(width = 0.6, colour = "black", alpha=0.7) +
  scale_fill_manual(values = c(Deg = "#5E4FA2",
                               Flat = "#767676"),
                    labels = c("No measurable degradation", "Degrading")) +
  labs(y = "% of prot.") +
  theme_bw() +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y  = element_text(size = 16, colour = "black"),
        axis.title.y = element_text(size = 16),
        panel.border = element_rect(linewidth = 1.2),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        plot.margin = margin(t = 0, r = 25,
                             b = 0, l = 0, unit = "pt")) +
  scale_y_continuous(breaks = c(0,50,100), expand = expansion(mult = c(0, 0)))


# --- stack, aligned x-axis ---
p_combined <- p2 / p1 + plot_layout(heights = c(1, 2.5))

# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=4, height=3.5, res=300)

p_combined
```

![](figures/Frog_Early_Analysis/Final%20disorder%20plot-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
target_bins <- data.frame(
  Macro_Compartment = c("Nucleus", "Cytoplasm", "Mitochondrion", "Endoplasmic reticulum", 
                        "Plasma membrane", "Golgi apparatus", "Extracellular region", 
                        "Chromatin", "Vesicle"),
  Macro_GO_ID = c("GO:0005634", "GO:0005737", "GO:0005739", "GO:0005783", 
                  "GO:0005886", "GO:0005794", "GO:0005576", 
                  "GO:0000785", "GO:0031982"))

#---------------------------

offspring_map <- as.list(GOCCOFFSPRING[target_bins$Macro_GO_ID])

mapping_df <- lapply(names(offspring_map),
                               function(parent_id) {
  data.frame(Macro_GO_ID = parent_id,
             Specific_GO_ID = offspring_map[[parent_id]],
             stringsAsFactors = FALSE) }) %>% bind_rows

mapping_df <- rbind(data.frame(Macro_GO_ID = target_bins$Macro_GO_ID,
                               Specific_GO_ID = target_bins$Macro_GO_ID),
                    mapping_df)

mapping_df <- merge(mapping_df, target_bins, by="Macro_GO_ID")
head(mapping_df)
```

    ##   Macro_GO_ID Specific_GO_ID Macro_Compartment
    ## 1  GO:0000785     GO:0033503         Chromatin
    ## 2  GO:0000785     GO:0000785         Chromatin
    ## 3  GO:0000785     GO:0000123         Chromatin
    ## 4  GO:0000785     GO:0000124         Chromatin
    ## 5  GO:0000785     GO:0000786         Chromatin
    ## 6  GO:0000785     GO:0000791         Chromatin

``` r
XLA_macro_CC <- XLA_fits_df[c("Protein_ID", "XLA_Gene", "Human_Gene",
                              "kD", "HL_Hrs", "mClass")]
XLA_macro_CC <- merge(XLA_macro_CC,
                      unique(human_GO[c("Protein_ID", "ID")]),
                      by="Protein_ID")
XLA_macro_CC <- XLA_macro_CC[XLA_macro_CC$ID %in%
                               mapping_df$Specific_GO_ID,]
XLA_macro_CC <- merge(XLA_macro_CC, mapping_df[2:3],
                      by.x="ID", by.y="Specific_GO_ID")
XLA_macro_CC["ID"] <- NULL
XLA_macro_CC <- unique(XLA_macro_CC) %>% arrange(Protein_ID)

CC_prot_counts <- data.frame(table(XLA_macro_CC$Protein_ID))

XLA_macro_CC <- XLA_macro_CC[XLA_macro_CC$Protein_ID %in%
                               CC_prot_counts[CC_prot_counts$Freq<=2,]$Var1,]

all_prot_CC <- unique(XLA_macro_CC[-7])
all_prot_CC["Macro_Compartment"] <- "All proteins"

XLA_macro_CC <- rbind(XLA_macro_CC, all_prot_CC)

head(XLA_macro_CC)
```

    ##                   Protein_ID    XLA_Gene Human_Gene           kD    HL_Hrs
    ## 1  XBmRNA10003|XBXL10_1g5362     sbno1.S      SBNO1 2.718905e-04  42.48935
    ## 9  XBmRNA10035|XBXL10_1g5375    vps33a.S     VPS33A 1.000000e-06 212.20000
    ## 10 XBmRNA10035|XBXL10_1g5375    vps33a.S     VPS33A 1.000000e-06 212.20000
    ## 18 XBmRNA10040|XBXL10_1g5381     psmd9.S      PSMD9 1.999634e-05 212.20000
    ## 19 XBmRNA10040|XBXL10_1g5381     psmd9.S      PSMD9 1.999634e-05 212.20000
    ## 20 XBmRNA10054|XBXL10_1g5386  tmem120b.S   TMEM120B 1.192172e-05 212.20000
    ##    mClass Macro_Compartment
    ## 1     Deg           Nucleus
    ## 9    Flat         Cytoplasm
    ## 10   Flat           Vesicle
    ## 18   Flat           Nucleus
    ## 19   Flat         Cytoplasm
    ## 20   Flat           Nucleus

``` r
CC_ordering <- lapply(unique(XLA_macro_CC$Macro_Compartment),
       function(x){
         
         sub_df <- XLA_macro_CC[XLA_macro_CC$Macro_Compartment == x,]
         
         max_hl_rows <- nrow(sub_df[sub_df$HL_Hrs == max_fit_hl,])

         percent_collapsed <- max_hl_rows/nrow(sub_df)
         median_hl <- median(sub_df$HL_Hrs)
          
         out_line <- data.frame(CC=x, HL=median_hl,
                                Collapse_Percent = round(percent_collapsed*100,0),
                                Deg_Percent = 100-round(percent_collapsed*100,0))

         return(out_line) }) %>% bind_rows %>%
  arrange(CC != "All proteins", -HL, -Collapse_Percent)

XLA_macro_CC$Macro_Compartment <- factor(XLA_macro_CC$Macro_Compartment,
                                         levels=CC_ordering$CC)

CC_ordering
```

    ##                       CC       HL Collapse_Percent Deg_Percent
    ## 1           All proteins 212.2000               63          37
    ## 2          Mitochondrion 212.2000               91           9
    ## 3  Endoplasmic reticulum 212.2000               81          19
    ## 4              Cytoplasm 212.2000               69          31
    ## 5                Vesicle 212.2000               61          39
    ## 6        Golgi apparatus 212.2000               59          41
    ## 7                Nucleus 212.2000               55          45
    ## 8        Plasma membrane 202.0717               49          51
    ## 9              Chromatin 160.4056               40          60
    ## 10  Extracellular region 146.2835               41          59

``` r
p_kd_CC <- ggplot(data=XLA_macro_CC,
                  aes(x = Macro_Compartment, y = HL_Hrs)) +

  geom_boxplot(outliers=FALSE, width = 0.7, color = "black", fill = "#5E4FA2", alpha=0.7) +
  
  geom_vline(xintercept = 1.5, linetype="dashed",
             color = "black", linewidth = 1) +

  stat_summary(fun = median, geom = "point", 
               shape = 21, size = 3, fill = "red3", color = "black") +

  # geom_text(data = CC_ordering,
  #           aes(x = CC, y = max_fit_hl, label = paste0(Deg_Percent, "%")),
  #           vjust = -0.6, size = 5, inherit.aes = FALSE) +
    
  theme_bw() +
  labs(x = "Cellular Compartment",
       y = "Half-life (hours)") +
  
  theme(axis.text.x = element_text(angle = 35, size=15, hjust = 1, colour="black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=16,colour="black"),
        axis.title.y = element_text(size=16, hjust=0.2),        
        plot.title = element_blank(),
        panel.border = element_rect(linewidth=1.2),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
        plot.margin = margin(t = 10, r = 25, b = 10,
                             l = 10, unit = "pt")) +
  coord_cartesian(ylim = c(0,230)) +
  scale_y_continuous(breaks=seq(0,200,50))
  
p_kd_CC
```

![](figures/Frog_Early_Analysis/Bottom%20of%20Compartment%20plot-1.png)<!-- -->

``` r
prop_df <- CC_ordering %>%
  pivot_longer(cols = -c(CC, HL),
               names_to = "Type",
               values_to = "Percent")
prop_df$CC <- factor(prop_df$CC, levels=CC_ordering$CC)

# --- top track: 100% stacked proportion bar ---
p_bar <- ggplot(prop_df, aes(x = CC, y = Percent, fill = Type)) +
  geom_col(width = 0.7, colour = "black", alpha=0.7) +
  geom_vline(xintercept = 1.5, linetype = "dashed", colour = "black", linewidth = 1) +
  scale_fill_manual(values = c(Deg_Percent = "#5E4FA2",
                               Collapse_Percent = "#767676"),
                    labels = c("No measurable degradation", "Degrading")) +
  labs(y = "% of prot.") +
  theme_bw() +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y  = element_text(size = 16, colour = "black"),
        axis.title.y = element_text(size = 16),
        panel.border = element_rect(linewidth = 1.2),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        plot.margin = margin(t = 15, r = 25,
                             b = 0, l = 0, unit = "pt")) +
  scale_y_continuous(breaks = c(0,50,100), expand = expansion(mult = c(0, 0)))

# --- bottom: your existing boxplot, x-title/labels live here ---
p_box <- p_kd_CC + theme(plot.margin = margin(t = 0, r = 25, b = 0, l = 10, unit = "pt"))

# --- stack, aligned x-axis ---
p_combined <- p_bar / p_box + plot_layout(heights = c(1, 4.3))

# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=8, height=5, res=300)

p_combined
```

![](figures/Frog_Early_Analysis/Complete%20compartment%20plot-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
norm_set <- read.csv("Data/XLA_Norm/XLA-O18_YolkNormSet_T8-NYS-Decay.csv")$Protein_ID
```

``` r
all_tracks_df <- XLA_fits_df
all_tracks_df["Light_FC"] <- apply(all_tracks_df, 1, function(row){
  
  light_fc <- as.numeric(row["N12"]) / as.numeric(row["T0_L"])
  
  return(light_fc) })

inc_min_fc <- quantile(all_tracks_df[all_tracks_df$Protein_ID %in% norm_set,]$Light_FC, 0.99)
dec_min_fc <- quantile(all_tracks_df[all_tracks_df$Protein_ID %in% norm_set,]$Light_FC, 0.01)

inc_min_fc
```

    ##      99% 
    ## 1.439371

``` r
dec_min_fc
```

    ##        1% 
    ## 0.6218052

``` r
all_tracks_df["Light_Type"] <- rep("Unchanging", nrow(all_tracks_df))
all_tracks_df[all_tracks_df$Light_FC > inc_min_fc, "Light_Type"] <- "Increasing"
all_tracks_df[all_tracks_df$Light_FC < dec_min_fc, "Light_Type"] <- "Decreasing"

all_tracks_df["Merged_Class"] <- paste0(all_tracks_df$Light_Type, "_", all_tracks_df$mClass)

table(all_tracks_df$Merged_Class)
```

    ## 
    ##  Decreasing_Deg Decreasing_Flat  Increasing_Deg Increasing_Flat  Unchanging_Deg 
    ##             258              30            1605            3692             709 
    ## Unchanging_Flat 
    ##            2512

``` r
summary.data <- data.frame(
  category = c("No measurable degradation", "Balanced Protein Turnover",
               "Protein Degradation Only", "Net Protein Increase with Degradation",
               "Protein Synthesis Only"),
  count = c(2512+30, 709, 258, 1605, 3692))

sum(summary.data$count)
```

    ## [1] 8806

``` r
summary.data$category
```

    ## [1] "No measurable degradation"            
    ## [2] "Balanced Protein Turnover"            
    ## [3] "Protein Degradation Only"             
    ## [4] "Net Protein Increase with Degradation"
    ## [5] "Protein Synthesis Only"

``` r
cat("\n")
```

``` r
round(summary.data$count / sum(summary.data$count) * 100,1)
```

    ## [1] 28.9  8.1  2.9 18.2 41.9

``` r
summary.data["category"] <- factor(summary.data$category,
                                   levels=c("No measurable degradation",
                                            "Protein Synthesis Only",
                                            "Net Protein Increase with Degradation",
                                            "Balanced Protein Turnover",
                                            "Protein Degradation Only"))

# Create a pie chart
deg.pie <- ggplot(summary.data, aes(x = 1, y = count, fill = category)) +
  geom_bar(width = 1, stat = "identity", color = "black", alpha=0.8,
           linewidth=0.25) +
  coord_polar(theta = "y") +
  xlim(0.5, 1.5) +
  scale_fill_manual(values = c("grey100", "grey80", "grey55", "grey20", "black")) +
  theme_void()

deg.pie
```

![](figures/Frog_Early_Analysis/Pie%20chart%20of%20kinetic%20classes-1.png)<!-- -->

``` r
ggsave('test.png', bg="transparent", plot=deg.pie,width=4,height=4,units="in")
```

``` r
plot_protein <- function(protein_id){

  protein_sub <- XLA_fits_df[XLA_fits_df$Protein_ID==protein_id,]  

  print(protein_sub$HL_Hrs)
  
  light_data <- as.numeric(protein_sub[10:17])
  light_data <- light_data / light_data[1]
  
  heavy_data <- as.numeric(protein_sub[18:28])
  heavy_data <- heavy_data / heavy_data[1]
      
  plot_data <- data.frame(Time=c(avg_light_time, avg_heavy_time),
                          Values=c(light_data, heavy_data),
                          Group=c(rep("Control", length(avg_light_time)),
                                  rep("O18", length(avg_heavy_time))))  
  plot_data["Time"] <- plot_data$Time/60
  

  p1 <- ggplot() +
    geom_line(data=plot_data[plot_data$Group=="O18",],
              aes(x=Time, y=Values, group=Group),
              color="#36A893", size=3) +
    geom_line(data=plot_data[plot_data$Group=="Control",],
              aes(x=Time, y=Values, group=Group),
              color="#767676", size=3) +    
    theme_bw() +
    theme(axis.text=element_text(size=18,colour="black"),
          plot.title = element_text(size=22, hjust=0.5),
          axis.title = element_blank(),
          panel.border = element_rect(linewidth=2),
          panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
          panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
          legend.position = "none",
          aspect.ratio = 1,
          plot.margin = margin(t = 5, r = 5, b = 5,
                               l = 5, unit = "pt")) +    
    geom_hline(yintercept = 1, linetype="dashed", linewidth=1) +
    coord_cartesian(ylim=c(0,2.25)) +
    scale_y_continuous(breaks=c(0,1,2))
  
  return(p1) }

# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=10, height=5, res=300)

plot_protein("XBmRNA25690|XBXL10_1g13517") |
plot_protein("XBmRNA30988|XBXL10_1g16504") |
plot_protein("XBmRNA74431|XBXL10_1g39492") |
plot_protein("XBmRNA32706|XBXL10_1g17749") |
plot_protein("XBmRNA79420|XBXL10_1g42248")
```

    ## [1] 212.2

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once per session.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## [1] 212.2
    ## [1] 7.526297
    ## [1] 55.37235
    ## [1] 33.23004

![](figures/Frog_Early_Analysis/Example%20proteins-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
sessionInfo()
```

    ## R version 4.5.3 (2026-03-11 ucrt)
    ## Platform: x86_64-w64-mingw32/x64
    ## Running under: Windows 11 x64 (build 26200)
    ## 
    ## Matrix products: default
    ##   LAPACK version 3.12.1
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.utf8 
    ## [2] LC_CTYPE=English_United States.utf8   
    ## [3] LC_MONETARY=English_United States.utf8
    ## [4] LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.utf8    
    ## 
    ## time zone: America/Los_Angeles
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] patchwork_1.3.2        effsize_0.8.1          forcats_1.0.1         
    ##  [4] ggpubr_0.6.3           Peptides_2.4.6         Biostrings_2.78.0     
    ##  [7] Seqinfo_1.0.0          XVector_0.50.0         minpack.lm_1.2-4      
    ## [10] clusterProfiler_4.18.4 org.Hs.eg.db_3.22.0    rrvgo_1.22.0          
    ## [13] KEGGREST_1.50.0        stringr_1.6.0          dplyr_1.2.0           
    ## [16] tidyr_1.3.2            GO.db_3.22.0           AnnotationDbi_1.72.0  
    ## [19] IRanges_2.44.0         S4Vectors_0.48.1       Biobase_2.70.0        
    ## [22] BiocGenerics_0.56.0    generics_0.1.4         fgsea_1.36.2          
    ## [25] ggplot2_4.0.3          lsa_0.73.4             SnowballC_0.7.1       
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RColorBrewer_1.1-3      rstudioapi_0.19.0       jsonlite_2.0.0         
    ##   [4] tidydr_0.0.6            umap_0.2.10.0           magrittr_2.0.4         
    ##   [7] ggtangle_0.1.2          farver_2.1.2            rmarkdown_2.31         
    ##  [10] ragg_1.5.2              fs_2.1.0                vctrs_0.7.1            
    ##  [13] memoise_2.0.1           ggtree_4.0.5            askpass_1.2.1          
    ##  [16] rstatix_0.7.3           htmltools_0.5.9         broom_1.0.13           
    ##  [19] Formula_1.2-5           gridGraphics_0.5-1      htmlwidgets_1.6.4      
    ##  [22] plyr_1.8.9              cachem_1.1.0            igraph_2.3.2           
    ##  [25] mime_0.13               lifecycle_1.0.5         pkgconfig_2.0.3        
    ##  [28] gson_0.1.0              Matrix_1.7-5            R6_2.6.1               
    ##  [31] fastmap_1.2.0           shiny_1.13.0            digest_0.6.39          
    ##  [34] aplot_0.2.9             enrichplot_1.30.5       colorspace_2.1-2       
    ##  [37] ggnewscale_0.5.2        RSpectra_0.16-2         textshaping_1.0.5      
    ##  [40] RSQLite_3.53.2          labeling_0.4.3          abind_1.4-8            
    ##  [43] polyclip_1.10-7         httr_1.4.8              compiler_4.5.3         
    ##  [46] bit64_4.8.2             fontquiver_0.2.1        withr_3.0.3            
    ##  [49] backports_1.5.1         S7_0.2.1                BiocParallel_1.44.0    
    ##  [52] carData_3.0-6           DBI_1.3.0               ggforce_0.5.0          
    ##  [55] R.utils_2.13.0          ggsignif_0.6.4          MASS_7.3-65            
    ##  [58] openssl_2.4.2           rappdirs_0.3.4          tools_4.5.3            
    ##  [61] otel_0.2.0              ape_5.8-1               scatterpie_0.2.6       
    ##  [64] httpuv_1.6.17           R.oo_1.27.1             glue_1.8.0             
    ##  [67] nlme_3.1-169            GOSemSim_2.36.0         promises_1.5.0         
    ##  [70] grid_4.5.3              gridBase_0.4-7          cluster_2.1.8.2        
    ##  [73] reshape2_1.4.5          snow_0.4-4              gtable_0.3.6           
    ##  [76] R.methodsS3_1.8.2       data.table_1.18.4       car_3.1-5              
    ##  [79] xml2_1.5.2              ggrepel_0.9.8           pillar_1.11.1          
    ##  [82] yulab.utils_0.2.4       later_1.4.8             splines_4.5.3          
    ##  [85] tweenr_2.0.3            treeio_1.34.0           lattice_0.22-9         
    ##  [88] bit_4.6.0               tidyselect_1.2.1        fontLiberation_0.1.0   
    ##  [91] tm_0.7-18               knitr_1.51              fontBitstreamVera_0.1.1
    ##  [94] NLP_0.3-2               xfun_0.58               pheatmap_1.0.13        
    ##  [97] stringi_1.8.7           lazyeval_0.2.3          ggfun_0.2.0            
    ## [100] yaml_2.3.12             evaluate_1.0.5          codetools_0.2-20       
    ## [103] wordcloud_2.6           gdtools_0.5.1           tibble_3.3.1           
    ## [106] qvalue_2.42.0           ggplotify_0.1.3         cli_3.6.5              
    ## [109] xtable_1.8-8            reticulate_1.46.0       systemfonts_1.3.2      
    ## [112] treemap_2.4-4           Rcpp_1.1.1-1.1          png_0.1-9              
    ## [115] parallel_4.5.3          blob_1.3.0              DOSE_4.4.0             
    ## [118] slam_0.1-55             tidytree_0.4.7          ggiraph_0.9.6          
    ## [121] scales_1.4.0            purrr_1.2.2             crayon_1.5.3           
    ## [124] rlang_1.1.7             cowplot_1.2.0           fastmatch_1.1-8
