Calculating absolute abundances for frog and fly
================
Edward Cruz
2025-09-24

``` r
rm(list=ls(all=T))

library(Biostrings)
```

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

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:grDevices':
    ## 
    ##     windows

    ## Loading required package: XVector

    ## Loading required package: Seqinfo

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

``` r
library(stringi)
library(arrow)
```

    ## 
    ## Attaching package: 'arrow'

    ## The following object is masked from 'package:utils':
    ## 
    ##     timestamp

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

    ## The following objects are masked from 'package:Biostrings':
    ## 
    ##     collapse, intersect, setdiff, setequal, union

    ## The following object is masked from 'package:Seqinfo':
    ## 
    ##     intersect

    ## The following object is masked from 'package:XVector':
    ## 
    ##     slice

    ## The following objects are masked from 'package:IRanges':
    ## 
    ##     collapse, desc, intersect, setdiff, slice, union

    ## The following objects are masked from 'package:S4Vectors':
    ## 
    ##     first, intersect, rename, setdiff, setequal, union

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
library(Peptides)
library(ggplot2)

knitr::opts_chunk$set(fig.path = "figures/turnover_abs_quant/")
```

``` r
#FASTA file used to search data
frog_fasta <- readAAStringSet("Files/FASTA/XENLA_longest-transcript_10p1_GFY-Fix.fasta")

#Converting to a dataframe for merging
frog_AA.df <- data.frame(
  Protein_ID = names(frog_fasta),
  Sequence = as.character(frog_fasta),
  stringsAsFactors = FALSE )
row.names(frog_AA.df) <- NULL
frog_AA.df["Protein_ID"] <- sapply(strsplit(frog_AA.df$`Protein_ID`, split = " "),
                                   "[", 1)

#FASTA file used to search data
fly_fasta <- readAAStringSet("Files/Reference/Dmel_Uniprot_Proteome_092722.fasta")

#Converting to a dataframe for merging
fly_AA.df <- data.frame(
  Protein_ID = names(fly_fasta),
  Sequence = as.character(fly_fasta),
  stringsAsFactors = FALSE )
row.names(fly_AA.df) <- NULL
fly_AA.df["Protein_ID"] <- sapply(strsplit(fly_AA.df$`Protein_ID`, split = " "),
                                  "[", 1)
fly_AA.df["Protein_ID"] <- sapply(strsplit(fly_AA.df$Protein_ID, split = "\\|"),
                                  "[", 2)
#---- Calculating MW ----------

frog_AA.df["Sequence"] <- str_remove_all(frog_AA.df$Sequence,
                                         "(?i)[^ARNDCEQGHILKMFPSTWYV]")
fly_AA.df["Sequence"] <- str_remove_all(fly_AA.df$Sequence,
                                        "(?i)[^ARNDCEQGHILKMFPSTWYV]")

frog_AA.df["kDa"] <- sapply(frog_AA.df$Sequence, function(seq){
    return(mw(seq)/1000) })
fly_AA.df["kDa"] <- sapply(fly_AA.df$Sequence, function(seq){
    return(mw(seq)/1000) })

frog_AA.df <- frog_AA.df %>% select(Protein_ID, kDa, Sequence)
fly_AA.df <- fly_AA.df %>% select(Protein_ID, kDa, Sequence)
```

``` r
tryptic_digest <- function(fasta_file) {

  #FASTA file used to search data
  protein_fasta <- 
    readAAStringSet(fasta_file)
  
  #Converting to a dataframe for merging
  protein_seq.df <- data.frame(
    Protein_ID = names(protein_fasta),
    Sequence = as.character(protein_fasta),
    stringsAsFactors = FALSE )
  row.names(protein_seq.df) <- NULL
  
  #Setting names and sequences for contaminants
  contams_seq.df <- protein_seq.df[grepl("contaminant",
                                         protein_seq.df$Protein_ID),]
  contams_seq.df["Protein_ID"] <- sapply(strsplit(contams_seq.df$Protein_ID, split = "\\|"),
                                         "[", 2)
  contams_seq.df["Description"] <- rep("contaminant", nrow(contams_seq.df))
  
  #Setting names for fwd sequences
  protein_seq.df <- protein_seq.df[!grepl("contaminant",
                                          protein_seq.df$Protein_ID),]
  protein_seq.df["Description"] <- rep("Temp", nrow(protein_seq.df))
  
  #Merging contaminants and organizing data.frame
  protein_seq.df <- rbind(contams_seq.df, protein_seq.df)
  
  #Determining number of tryptic peptides in sample
  protein_seq.df <- cbind(protein_seq.df,
                          lapply(protein_seq.df$Sequence, function(seq){
    
    peps <- stri_split_regex(seq, "(?<=[RK])")[[1]]  
    peps <- peps[nchar(peps) >= 7 & nchar(peps) <= 30]
    
    theo_peps <- length(peps)
    peps <- paste(peps, collapse = ",")
    
    outline <- data.frame(Seq_Length = nchar(seq),
                          Theo_Peps = theo_peps,
                          Fragments = peps)
    
    return(outline) }) %>% bind_rows )
  
  protein_seq.df <- protein_seq.df %>% 
    select(Protein_ID, Description, Seq_Length, Theo_Peps,
           Sequence, Fragments)
  
  return(protein_seq.df) }
```

``` r
XLA_Digest <- 
  tryptic_digest("Files/FASTA/DIA-NN/XENLA_longest-transcript_10p1_GFY-Fix_w_contams.fasta")

XLA_Digest[!XLA_Digest$Description=="contaminant", "Description"] <-
  sapply(strsplit(XLA_Digest[!XLA_Digest$Description=="contaminant",]$`Protein_ID`,
                  split = " "), "[", 2)
XLA_Digest[is.na(XLA_Digest$Description), "Description"] <- "N/A"

XLA_Digest[!XLA_Digest$Description=="contaminant", "Protein_ID"] <-
  sapply(strsplit(XLA_Digest[!XLA_Digest$Description=="contaminant",]$`Protein_ID`,
                  split = " "), "[", 1)
```

``` r
Dmel_Digest <- 
  tryptic_digest("Files/FASTA/DIA-NN/Dmel_Uniprot_Proteome_092722_w_contams.fasta")

Dmel_Digest[!Dmel_Digest$Description=="contaminant", "Description"] <-
  sapply(Dmel_Digest[!Dmel_Digest$Description=="contaminant",]$Protein_ID,
         function(x){
           description <- str_match(x, " ([^\r\n]*?)(?=\\s*OS=)")[,2]
         return(description) })
Dmel_Digest[!Dmel_Digest$Description=="contaminant", "Protein_ID"] <- 
  sapply(strsplit(Dmel_Digest[!Dmel_Digest$Description=="contaminant",]$Protein_ID,
                  split = "\\|"), "[", 2)

head(Dmel_Digest)
```

    ##   Protein_ID Description Seq_Length Theo_Peps
    ## 1     P20930 contaminant       4061       210
    ## 2     P15924 contaminant       2871       163
    ## 3     Q86YZ3 contaminant       2850        85
    ## 4     Q5D862 contaminant       2391        74
    ## 5     P23172 contaminant       1512        74
    ## 6     Q87022 contaminant       1505        69
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        Sequence
    ## 1 MSTLLENIFAIINLFKQYSKKDKNTDTLSKKELKELLEKEFRQILKNPDDPDMVDVFMDHLDIDHNKKIDFTEFLLMVFKLAQAYYESTRKENLPISGHKHRKHSHHDKHEDNKQEENKENRKRPSSLERRNNRKGNKGRSKSPRETGGKRHESSSEKKERKGYSPTHREEEYGKNHHNSSKKEKNKTENTRLGDNRKRLSERLEEKEDNEEGVYDYENTGRMTQKWIQSGHIATYYTIQDEAYDTTDSLLEENKIYERSRSSDGKSSSQVNRSRHENTSQVPLQESRTRKRRGSRVSQDRDSEGHSEDSERHSGSASRNHHGSAWEQSRDGSRHPRSHDEDRASHGHSADSSRQSGTRHAETSSRGQTASSHEQARSSPGERHGSGHQQSADSSRHSATGRGQASSAVSDRGHRGSSGSQASDSEGHSENSDTQSVSGHGKAGLRQQSHQESTRGRSGERSGRSGSSLYQVSTHEQPDSAHGRTGTSTGGRQGSHHEQARDSSRHSASQEGQDTIRGHPGSSRGGRQGSHHEQSVNRSGHSGSHHSHTTSQGRSDASHGQSGSRSASRQTRNEEQSGDGTRHSGSRHHEASSQADSSRHSQVGQGQSSGPRTSRNQGSSVSQDSDSQGHSEDSERWSGSASRNHHGSAQEQSRDGSRHPRSHHEDRAGHGHSADSSRKSGTRHTQNSSSGQAASSHEQARSSAGERHGSRHQLQSADSSRHSGTGHGQASSAVRDSGHRGSSGSQATDSEGHSEDSDTQSVSGHGQAGHHQQSHQESARDRSGERSRRSGSFLYQVSTHKQSESSHGWTGPSTGVRQGSHHEQARDNSRHSASQDGQDTIRGHPGSSRRGRQGSHHEQSVDRSGHSGSHHSHTTSQGRSDASRGQSGSRSASRTTRNEEQSRDGSRHSGSRHHEASSHADISRHSQAGQGQSEGSRTSRRQGSSVSQDSDSEGHSEDSERWSGSASRNHRGSAQEQSRHGSRHPRSHHEDRAGHGHSADSSRQSGTPHAETSSGGQAASSHEQARSSPGERHGSRHQQSADSSRHSGIPRRQASSAVRDSGHWGSSGSQASDSEGHSEESDTQSVSGHGQDGPHQQSHQESARDWSGGRSGRSGSFIYQVSTHEQSESAHGRTRTSTGRRQGSHHEQARDSSRHSASQEGQDTIRAHPGSRRGGRQGSHHEQSVDRSGHSGSHHSHTTSQGRSDASHGQSGSRSASRQTRKDKQSGDGSRHSGSRHHEAASWADSSRHSQVGQEQSSGSRTSRHQGSSVSQDSDSERHSDDSERLSGSASRNHHGSSREQSRDGSRHPGFHQEDRASHGHSADSSRQSGTHHTESSSHGQAVSSHEQARSSPGERHGSRHQQSADSSRHSGIGHRQASSAVRDSGHRGSSGSQVTNSEGHSEDSDTQSVSAHGQAGPHQQSHKESARGQSGESSGRSRSFLYQVSSHEQSESTHGQTAPSTGGRQGSRHEQARNSSRHSASQDGQDTIRGHPGSSRGGRQGSYHEQSVDRSGHSGYHHSHTTPQGRSDASHGQSGPRSASRQTRNEEQSGDGSRHSGSRHHEPSTRAGSSRHSQVGQGESAGSKTSRRQGSSVSQDRDSEGHSEDSERRSESASRNHYGSAREQSRHGSRNPRSHQEDRASHGHSAESSRQSGTRHAETSSGGQAASSQEQARSSPGERHGSRHQQSADSSTDSGTGRRQDSSVVGDSGNRGSSGSQASDSEGHSEESDTQSVSAHGQAGPHQQSHQESTRGQSGERSGRSGSFLYQVSTHEQSESAHGRTGPSTGGRQRSRHEQARDSSRHSASQEGQDTIRGHPGSSRGGRQGSHYEQSVDSSGHSGSHHSHTTSQERSDVSRGQSGSRSVSRQTRNEKQSGDGSRHSGSRHHEASSRADSSRHSQVGQGQSSGPRTSRNQGSSVSQDSDSQGHSEDSERWSGSASRNHLGSAWEQSRDGSRHPGSHHEDRAGHGHSADSSRQSGTRHTESSSRGQAASSHEQARSSAGERHGSHHQLQSADSSRHSGIGHGQASSAVRDSGHRGYSGSQASDSEGHSEDSDTQSVSAQGKAGPHQQSHKESARGQSGESSGRSGSFLYQVSTHEQSESTHGQSAPSTGGRQGSHYDQAQDSSRHSASQEGQDTIRGHPGPSRGGRQGSHQEQSVDRSGHSGSHHSHTTSQGRSDASRGQSGSRSASRKTYDKEQSGDGSRHSGSHHHEASSWADSSRHSLVGQGQSSGPRTSRPRGSSVSQDSDSEGHSEDSERRSGSASRNHHGSAQEQSRDGSRHPRSHHEDRAGHGHSAESSRQSGTHHAENSSGGQAASSHEQARSSAGERHGSHHQQSADSSRHSGIGHGQASSAVRDSGHRGSSGSQASDSEGHSEDSDTQSVSAHGQAGPHQQSHQESTRGRSAGRSGRSGSFLYQVSTHEQSESAHGRTGTSTGGRQGSHHKQARDSSRHSTSQEGQDTIHGHPGSSSGGRQGSHYEQLVDRSGHSGSHHSHTTSQGRSDASHGHSGSRSASRQTRNDEQSGDGSRHSGSRHHEASSRADSSGHSQVGQGQSEGPRTSRNWGSSFSQDSDSQGHSEDSERWSGSASRNHHGSAQEQLRDGSRHPRSHQEDRAGHGHSADSSRQSGTRHTQTSSGGQAASSHEQARSSAGERHGSHHQQSADSSRHSGIGHGQASSAVRDSGHRGYSGSQASDNEGHSEDSDTQSVSAHGQAGSHQQSHQESARGRSGETSGHSGSFLYQVSTHEQSESSHGWTGPSTRGRQGSRHEQAQDSSRHSASQDGQDTIRGHPGSSRGGRQGYHHEHSVDSSGHSGSHHSHTTSQGRSDASRGQSGSRSASRTTRNEEQSGDGSRHSGSRHHEASTHADISRHSQAVQGQSEGSRRSRRQGSSVSQDSDSEGHSEDSERWSGSASRNHHGSAQEQLRDGSRHPRSHQEDRAGHGHSADSSRQSGTRHTQTSSGGQAASSHEQARSSAGERHGSHHQQSADSSRHSGIGHGQASSAVRDSGHRGYSGSQASDNEGHSEDSDTQSVSAHGQAGSHQQSHQESARGRSGETSGHSGSFLYQVSTHEQSESSHGWTGPSTRGRQGSRHEQAQDSSRHSASQYGQDTIRGHPGSSRGGRQGYHHEHSVDSSGHSGSHHSHTTSQGRSDASRGQSGSRSASRTTRNEEQSGDSSRHSVSRHHEASTHADISRHSQAVQGQSEGSRRSRRQGSSVSQDSDSEGHSEDSERWSGSASRNHRGSVQEQSRHGSRHPRSHHEDRAGHGHSADRSRQSGTRHAETSSGGQAASSHEQARSSPGERHGSRHQQSADSSRHSGIPRGQASSAVRDSRHWGSSGSQASDSEGHSEESDTQSVSGHGQAGPHQQSHQESARDRSGGRSGRSGSFLYQVSTHEQSESAHGRTRTSTGRRQGSHHEQARDSSRHSASQEGQDTIRGHPGSSRRGRQGSHYEQSVDRSGHSGSHHSHTTSQGRSDASRGQSGSRSASRQTRNDEQSGDGSRHSWSHHHEASTQADSSRHSQSGQGQSAGPRTSRNQGSSVSQDSDSQGHSEDSERWSGSASRNHRGSAQEQSRDGSRHPTSHHEDRAGHGHSAESSRQSGTHHAENSSGGQAASSHEQARSSAGERHGSHHQQSADSSRHSGIGHGQASSAVRDSGHRGSSGSQASDSEGHSEDSDTQSVSAHGQAGPHQQSHQESTRGRSAGRSGRSGSFLYQVSTHEQSESAHGRAGPSTGGRQGSRHEQARDSSRHSASQEGQDTIRGHPGSRRGGRQGSYHEQSVDRSGHSGSHHSHTTSQGRSDASHGQSGSRSASRETRNEEQSGDGSRHSGSRHHEASTQADSSRHSQSGQGESAGSRRSRRQGSSVSQDSDSEAYPEDSERRSESASRNHHGSSREQSRDGSRHPGSSHRDTASHVQSSPVQSDSSTAKEHGHFSSLSQDSAYHSGIQSRGSPHSSSSYHYQSEGTERQKGQSGLVWRHGSYGSADYDYGESGFRHSQHGSVSYNSNPVVFKERSDICKASAFGKDHPRYYATYINKDPGLCGHSSDISKQLGFSQSQRYYYYE
    ## 2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       MSCNGGSHPRINTLGRMIRAESGPDLRYEVTSGGGGTSRMYYSRRGVITDQNSDGYCQTGTMSRHQNQNTIQELLQNCSDCLMRAELIVQPELKYGDGIQLTRSRELDECFAQANDQMEILDSLIREMRQMGQPCDAYQKRLLQLQEQMRALYKAISVPRVRRASSKGGGGYTCQSGSGWDEFTKHVTSECLGWMRQQRAEMDMVAWGVDLASVEQHINSHRGIHNSIGDYRWQLDKIKADLREKSAIYQLEEEYENLLKASFERMDHLRQLQNIIQATSREIMWINDCEEEELLYDWSDKNTNIAQKQEAFSIRMSQLEVKEKELNKLKQESDQLVLNQHPASDKIEAYMDTLQTQWSWILQITKCIDVHLKENAAYFQFFEEAQSTEAYLKGLQDSIRKKYPCDKNMPLQHLLEQIKELEKEREKILEYKRQVQNLVNKSKKIVQLKPRNPDYRSNKPIILRALCDYKQDQKIVHKGDECILKDNNERSKWYVTGPGGVDMLVPSVGLIIPPPNPLAVDLSCKIEQYYEAILALWNQLYINMKSLVSWHYCMIDIEKIRAMTIAKLKTMRQEDYMKTIADLELHYQEFIRNSQGSEMFGDDDKRKIQSQFTDAQKHYQTLVIQLPGYPQHQTVTTTEITHHGTCQDVNHNKVIETNRENDKQETWMLMELQKIRRQIEHCEGRMTLKNLPLADQGSSHHITVKINELKSVQNDSQAIAEVLNQLKDMLANFRGSEKYCYLQNEVFGLFQKLENINGVTDGYLNSLCTVRALLQAILQTEDMLKVYEARLTEEETVCLDLDKVEAYRCGLKKIKNDLNLKKSLLATMKTELQKAQQIHSQTSQQYPLYDLDLGKFGEKVTQLTDRWQRIDKQIDFRLWDLEKQIKQLRNYRDNYQAFCKWLYDAKRRQDSLESMKFGDSNTVMRFLNEQKNLHSEISGKRDKSEEVQKIAELCANSIKDYELQLASYTSGLETLLNIPIKRTMIQSPSGVILQEAADVHARYIELLTRSGDYYRFLSEMLKSLEDLKLKNTKIEVLEEELRLARDANSENCNKNKFLDQNLQKYQAECSQFKAKLASLEELKRQAELDGKSAKQNLDKCYGQIKELNEKITRLTYEIEDEKRRRKSVEDRFDQQKNDYDQLQKARQCEKENLGWQKLESEKAIKEKEYEIERLRVLLQEEGTRKREYENELAKVRNHYNEEMSNLRNKYETEINITKTTIKEISMQKEDDSKNLRNQLDRLSRENRDLKDEIVRLNDSILQATEQRRRAEENALQQKACGSEIMQKKQHLEIELKQVMQQRSEDNARHKQSLEEAAKTIQDKNKEIERLKAEFQEEAKRRWEYENELSKVRNNYDEEIISLKNQFETEINITKTTIHQLTMQKEEDTSGYRAQIDNLTRENRSLSEEIKRLKNTLTQTTENLRRVEEDIQQQKATGSEVSQRKQQLEVELRQVTQMRTEESVRYKQSLDDAAKTIQDKNKEIERLKQLIDKETNDRKCLEDENARLQRVQYDLQKANSSATETINKLKVQEQELTRLRIDYERVSQERTVKDQDITRFQNSLKELQLQKQKVEEELNRLKRTASEDSCKRKKLEEELEGMRRSLKEQAIKITNLTQQLEQASIVKKRSEDDLRQQRDVLDGHLREKQRTQEELRRLSSEVEALRRQLLQEQESVKQAHLRNEHFQKAIEDKSRSLNESKIEIERLQSLTENLTKEHLMLEEELRNLRLEYDDLRRGRSEADSDKNATILELRSQLQISNNRTLELQGLINDLQRERENLRQEIEKFQKQALEASNRIQESKNQCTQVVQERESLLVKIKVLEQDKARLQRLEDELNRAKSTLEAETRVKQRLECEKQQIQNDLNQWKTQYSRKEEAIRKIESEREKSEREKNSLRSEIERLQAEIKRIEERCRRKLEDSTRETQSQLETERSRYQREIDKLRQRPYGSHRETQTECEWTVDTSKLVFDGLRKKVTAMQLYECQLIDKTTLDKLLKGKKSVEEVASEIQPFLRGAGSIAGASASPKEKYSLVEAKRKKLISPESTVMLLEAQAATGGIIDPHRNEKLTVDSAIARDLIDFDDRQQIYAAEKAITGFDDPFSGKTVSVSEAIKKNLIDRETGMRLLEAQIASGGVVDPVNSVFLPKDVALARGLIDRDLYRSLNDPRDSQKNFVDPVTKKKVSYVQLKERCRIEPHTGLLLLSVQKRSMSFQGIRQPVTVTELVDSGILRPSTVNELESGQISYDEVGERIKDFLQGSSCIAGIYNETTKQKLGIYEAMKIGLVRPGTALELLEAQAATGFIVDPVSNLRLPVEEAYKRGLVGIEFKEKLLSAERAVTGYNDPETGNIISLFQAMNKELIEKGHGIRLLEAQIATGGIIDPKESHRLPVDIAYKRGYFNEELSEILSDPSDDTKGFFDPNTEENLTYLQLKERCIKDEETGLCLLPLKEKKKQVQTSQKNTLRKRRVVIVDPETNKEMSVQEAYKKGLIDYETFKELCEQECEWEEITITGSDGSTRVVLVDRKTGSQYDIQDAIDKGLVDRKFFDQYRSGSLSLTQFADMISLKNGVGTSSSMGSGVSDDVFSSSRHESVSKISTISSVRNLTIRSSSFSDTLEESSPIAAIFDTENLEKISITEGIERGIVDSITGQRLLEAQACTGGIIHPTTGQKLSLQDAVSQGVIDQDMATRLKPAQKAFIGFEGVKGKKKMSAAEAVKEKWLPYEAGQRFLEFQYLTGGLVDPEVHGRISTEEAIRKGFIDGRAAQRLQDTSSYAKILTCPKTKLKISYKDAINRSMVEDITGLRLLEAASVSSKGLPSPYNMSSAPGSRSGSRSGSRSGSRSGSRSGSRRGSFDATGNSSYSYSYSFSSSSIGH
    ## 3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            MPKLLQGVITVIDVFYQYATQHGEYDTLNKAELKELLENEFHQILKNPNDPDTVDIILQSLDRDHNKKVDFTEYLLMIFKLVQARNKIIGKDYCQVSGSKLRDDTHQHQEEQEETEKEENKRQESSFSHSSWSAGENDSYSRNVRGSLKPGTESISRRLSFQRDFSGQHNSYSGQSSSYGEQNSDSHQSSGRGQCGSGSGQSPNYGQHGSGSGQSSSNDTHGSGSGQSSGFSQHKSSSGQSSGYSQHGSGSGHSSGYGQHGSRSGQSSRGERHRSSSGSSSSYGQHGSGSRQSLGHGRQGSGSRQSPSHVRHGSGSGHSSSHGQHGSGSSYSYSRGHYESGSGQTSGFGQHESGSGQSSGYSKHGSGSGHSSSQGQHGSTSGQASSSGQHGSSSRQSSSYGQHESASRHSSGRGQHSSGSGQSPGHGQRGSGSGQSPSSGQHGTGFGRSSSSGPYVSGSGYSSGFGHHESSSEHSSGYTQHGSGSGHSSGHGQHGSRSGQSSRGERQGSSAGSSSSYGQHGSGSRQSLGHSRHGSGSGQSPSPSRGRHESGSRQSSSYGPHGYGSGRSSSRGPYESGSGHSSGLGHQESRSGQSSGYGQHGSSSGHSSTHGQHGSTSGQSSSCGQHGATSGQSSSHGQHGSGSSQSSRYGQQGSGSGQSPSRGRHGSDFGHSSSYGQHGSGSGWSSSNGPHGSVSGQSSGFGHKSGSGQSSGYSQHGSGSSHSSGYRKHGSRSGQSSRSEQHGSSSGLSSSYGQHGSGSHQSSGHGRQGSGSGHSPSRVRHGSSSGHSSSHGQHGSGTSCSSSCGHYESGSGQASGFGQHESGSGQGYSQHGSASGHFSSQGRHGSTSGQSSSSGQHDSSSGQSSSYGQHESASHHASGRGRHGSGSGQSPGHGQRGSGSGQSPSYGRHGSGSGRSSSSGRHGSGSGQSSGFGHKSSSGQSSGYTQHGSGSGHSSSYEQHGSRSGQSSRSEQHGSSSGSSSSYGQHGSGSRQSLGHGQHGSGSGQSPSPSRGRHGSGSGQSSSYGPYRSGSGWSSSRGPYESGSGHSSGLGHRESRSGQSSGYGQHGSSSGHSSTHGQHGSTSGQSSSCGQHGASSGQSSSHGQHGSGSSQSSGYGRQGSGSGQSPGHGQRGSGSRQSPSYGRHGSGSGRSSSSGQHGSGLGESSGFGHHESSSGQSSSYSQHGSGSGHSSGYGQHGSRSGQSSRGERHGSSSGSSSHYGQHGSGSRQSSGHGRQGSGSGHSPSRGRHGSGLGHSSSHGQHGSGSGRSSSRGPYESRSGHSSVFGQHESGSGHSSAYSQHGSGSGHFCSQGQHGSTSGQSSTFDQEGSSTGQSSSYGHRGSGSSQSSGYGRHGAGSGQSPSRGRHGSGSGHSSSYGQHGSGSGWSSSSGRHGSGSGQSSGFGHHESSSWQSSGCTQHGSGSGHSSSYEQHGSRSGQSSRGERHGSSSGSSSSYGQHGSGSRQSLGHGQHGSGSGQSPSPSRGRHGSGSGQSSSYSPYGSGSGWSSSRGPYESGSSHSSGLGHRESRSGQSSGYGQHGSSSGHSSTHGQHGSTSGQSSSCGQHGASSGQSSSHGQHGSGSSQSSGYGRQGSGSGQSPGHGQRGSGSRQSPSYGRHGSGSGRSSSSGQHGSGLGESSGFGHHESSSGQSSSYSQHGSGSGHSSGYGQHGSRSGQSSRGERHGSSSRSSSRYGQHGSGSRQSSGHGRQGSGSGQSPSRGRHGSGLGHSSSHGQHGSGSGRSSSRGPYESRSGHSSVFGQHESGSGHSSAYSQHGSGSGHFCSQGQHGSTSGQSSTFDQEGSSTGQSSSHGQHGSGSSQSSSYGQQGSGSGQSPSRGRHGSGSGHSSSYGQHGSGSGWSSSSGRHGSGSGQSSGFGHHESSSWQSSGYTQHGSGSGHSSSYEQHGSRSGQSSRGEQHGSSSGSSSSYGQHGSGSRQSLGHGQHGSGSGQSPSPSRGRHGSGSGQSSSYGPYGSGSGWSSSRGPYESGSGHSSGLGHRESRSGQSSGYGQHGSSSGHSSTHGQHGSASGQSSSCGQHGASSGQSSSHGQHGSGSSQSSGYGRQGSGSGQSPGHGQRGSGSRQSPSYGRHGSGSGRSSSSGQHGPGLGESSGFGHHESSSGQSSSYSQHGSGSGHSSGYGQHGSRSGQSSRGERHGSSSGSSSRYGQHGSGSRQSSGHGRQGSGSGHSPSRGRHGSGSGHSSSHGQHGSGSGRSSSRGPYESRSGHSSVFGQHESGSGHSSAYSQHGSGSGHFCSQGQHGSTSGQSSTFDQEGSSTGQSSSHGQHGSGSSQSSSYGQQGSGSGQSPSRGRHGSGSGHSSSYGQHGSGSGWSSSSGRHGSGSGQSSGFGHHESSSWQSSGYTQHGSGSGHSSSYEQHGSRSGQSSRGERHGSSSGSSSSYGQHGSGSRQSLGHGQHGSGSGQSPSPSRGRHGSGSGQSSSYSPYGSGSGWSSSRGPYESGSGHSSGLGHRESRSGQSSGYGQHGSSSGHSSTHGQHGSTSGQSSSCGQHGASSGQSSSHGQHGSGSSQSSGYGRQGSGSGQSPGHGQRGSGSRQSPSYGRHGSGSGRSSSSGQHGSGLGESSGFGHHESSSGQSSSYSQHGSGSGHSSGYGQHGSRSGQSSRGERHGSSSGSSSHYGQHGSGSRQSSGHGRQGSGSGQSPSRGRHGSGLGHSSSHGQHGSGSGRSSSRGPYESRLGHSSVFGQHESGSGHSSAYSQHGSGSGHFCSQGQHGSTSGQSSTFDQEGSSTGQSSSYGHRGSGSSQSSGYGRHGAGSGQSLSHGRHGSGSGQSSSYGQHGSGSGQSSGYSQHGSGSGQDGYSYCKGGSNHDGGSSGSYFLSFPSSTSPYEYVQEQRCYFYQ
    ## 4                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       MTDLLRSVVTVIDVFYKYTKQDGECGTLSKGELKELLEKELHPVLKNPDDPDTVDVIMHMLDRDHDRRLDFTEFLLMIFKLTMACNKVLSKEYCKASGSKKHRRGHRHQEEESETEEDEEDTPGHKSGYRHSSWSEGEEHGYSSGHSRGTVKCRHGSNSRRLGRQGNLSSSGNQEGSQKRYHRSSCGHSWSGGKDRHGSSSVELRERINKSHISPSRESGEEYESGSGSNSWERKGHGGLSCGLETSGHESNSTQSRIREQKLGSSCSGSGDSGRRSHACGYSNSSGCGRPQNASSSCQSHRFGGQGNQFSYIQSGCQSGIKGGQGHGCVSGGQPSGCGQPESNPCSQSYSQRGYGARENGQPQNCGGQWRTGSSQSSCCGQYGSGGSQSCSNGQHEYGSCGRFSNSSSSNEFSKCDQYGSGSSQSTSFEQHGTGLSQSSGFEQHVCGSGQTCGQHESTSSQSLGYDQHGSSSGKTSGFGQHGSGSGQSSGFGQCGSGSGQSSGFGQHGSVSGQSSGFGQHGSVSGQSSGFGQHESRSRQSSYGQHGSGSSQSSGYGQYGSRETSGFGQHGLGSGQSTGFGQYGSGSGQSSGFGQHGSGSGQSSGFGQHESRSGQSSYGQHSSGSSQSSGYGQHGSRQTSGFGQHGSGSSQSTGFGQYGSGSGQSSGFGQHVSGSGQSSGFGQHESRSGHSSYGQHGFGSSQSSGYGQHGSSSGQTSGFGQHELSSGQSSSFGQHGSGSGQSSGFGQHGSGSGQSSGFGQHESRSGQSSYGQHSSGSSQSSGYGQHGSRQTSGFGQHGSGSSQSTGFGQYGSGSGQSAGFGQHGSGSGQSSGFGQHESRSHQSSYGQHGSGSSQSSGYGQHGSSSGQTSGFGQHRSSSGQYSGFGQHGSGSGQSSGFGQHGTGSGQYSGFGQHESRSHQSSYGQHGSGSSQSSGYGQHGSSSGQTFGFGQHRSGSGQSSGFGQHGSGSGQSSGFGQHESGSGKSSGFGQHESRSSQSNYGQHGSGSSQSSGYGQHGSSSGQTTGFGQHRSSSGQYSGFGQHGSGSDQSSGFGQHGTGSGQSSGFGQYESRSRQSSYGQHGSGSSQSSGYGQHGSNSGQTSGFGQHRPGSGQSSGFGQYGSGSGQSSGFGQHGSGTGKSSGFAQHEYRSGQSSYGQHGTGSSQSSGCGQHESGSGPTTSFGQHVSGSDNFSSSGQHISDSGQSTGFGQYGSGSGQSTGLGQGESQQVESGSTVHGRQETTHGQTINTTRHSQSGQGQSTQTGSRVTRRRRSSQSENSDSEVHSKVSHRHSEHIHTQAGSHYPKSGSTVRRRQGTTHGQRGDTTRHGHSGHGQSTQTGSRTSGRQRFSHSDATDSEVHSGVSHRPHSQEQTHSQAGSQHGESESTVHERHETTYGQTGEATGHGHSGHGQSTQRGSRTTGRRGSGHSESSDSEVHSGGSHRPQSQEQTHGQAGSQHGESGSTVHGRHGTTHGQTGDTTRHAHYHHGKSTQRGSSTTGRRGSGHSESSDSEVHSGGSHTHSGHTHGQSGSQHGESESIIHDRHRITHGQTGDTTRHSYSGHEQTTQTGSRTTGRQRTSHSESTDSEVHSGGSHRPHSREHTYGQAGSQHEEPEFTVHERHGTTHGQIGDTTGHSHSGHGQSTQRGSRTTGRQRSSHSESSDSEVHSGVSHTHTGHTHGQAGSQHGQSESIVPERHGTTHGQTGDTTRHAHYHHGLTTQTGSRTTGRRGSGHSEYSDSEGYSGVSHTHSGHTHGQARSQHGESESIVHERHGTIHGQTGDTTRHAHSGHGQSTQTGSRTTGRRSSGHSEYSDSEGHSGFSQRPHSRGHTHGQAGSQHGESESIVDERHGTTHGQTGDTSGHSQSGHGQSTQSGSSTTGRRRSGHSESSDSEVHSGGSHTHSGHTHSQARSQHGESESTVHKRHQTTHGQTGDTTEHGHPSHGQTIQTGSRTTGRRGSGHSEYSDSEGPSGVSHTHSGHTHGQAGSHYPESGSSVHERHGTTHGQTADTTRHGHSGHGQSTQRGSRTTGRRASGHSEYSDSEGHSGVSHTHSGHAHGQAGSQHGESGSSVHERHGTTHGQTGDTTRHAHSGHGQSTQRGSRTAGRRGSGHSESSDSEVHSGVSHTHSGHTYGQARSQHGESGSAIHGRQGTIHGQTGDTTRHGQSGHGQSTQTGSRTTGRQRSSHSESSDSEVHSEASPTHSGHTHSQAGSRHGQSGSSGHGRQGTTHGQTGDTTRHAHYGYGQSTQRGSRTTGRRGSGHSESSDSEVHSWGSHTHSGHIQGQAGSQQRQPGSTVHGRLETTHGQTGDTTRHGHSGYGQSTQTGSRSSRASHFQSHSSERQRHGSSQVWKHGSYGPAEYDYGHTGYGPSGGSRKSISNSHLSWSTDSTANKQLSRH
    ## 5                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      MSSLLNSLLPEYFKPKTNLNINSSRVQYGFNARIDMQYEDDSGTRKGSRPNAFMSNTVAFIGNYEGIIVDDIPILDGLRADIFDTHGDLDMGLVEDALSKSTMIRRNVPTYTAYASELLYKRNLTSLFYNMLRLYYIKKWGSIKYEKDAIFYDNGHACLLNRQLFPKSRDASLESSLSLPEAEIAMLDPGLEFPEEDVPAILWHGRVSSRATCILGQACSEFAPLAPFSIAHYSPQLTRKLFVNAPAGIEPSSGRYTHEDVKDAITILVSANQAYTDFEAAYLMLAQTLVSPVPRTAEASAWFINAGMVNMPTLSCANGYYPALTNVNPYHRLDTWKDTLNHWVAYPDMLFYHSVAMIESCYVELGNVARVSDSDAINKYTFTELSVQGRPVMNRGIIVDLTLVAMRTGREISLPYPVSCGLTRTDALLQGTEIHVPVVVKDIDMPQYYNAIDKDVIEGQETVIKVKQLPPAMYPIYTYGINTTEFYSDHFEDQVQVEMAPIDNGKAVFNDARKFSKFMSIMRMMGNDVTATDLVTGRKVSNWADNSSGRFLYTDVKYEGQTAFLVDMDTVKARDHCWVSIVDPNGTMNLSYKMTNFRAAMFSRNKPLYMTGGSVRTIATGNYRDAAERLRAMDETLRLKPFKITEKLDFSCSSLRDTKFVGQQYAILTPSGTTTDIRSGRGTNQSYRRGRTSTGYRIGVEDDEDLDIGTVKYIVPLYLNGDNVAQNCLEATHVLIKACSIANRIVDDGEGHCFTQQGLAQQWIFHRGEMIFVKAVRIGQLNAYYVDYKNVTNYSLKTAAQVGATISNNLRHGFVDNQQDAYTRLVANYSDTRKWIRDNFTYNYNMEKEKYRITQYHHTHVRLKDLFPSRKIVKLEGYEALLAMMLDRFNNIESTHVTFFTYLRALPDREKEVFISLVLNYNGLGREWLKSEGVRAKQAQGTVKYDMSKLFELNVLENGVDEEVDWEKEKRNRSDIKTVNISYAKVLEHCRELFIMARAEGKRPMRMKWQEYWRQRAVIMPGGSVHSQHPVEQDVIRVLPREIRSKKGVASVMPYKEQKYFTSRRPEIHAYTSTKYEWGKVRALYGCDFSSHTMADFGLLQCEDTFPGFVPTGSYANEDYVRTRIAGTHSLIPFCYDFDDFNSQHSKEAMQAVIDAWISVYHDKLTDDQIEAAKWTRNSVDRMVAHQPNTGETYDVKGTLFSGWRLTTFFNTALNYCYLANAGINSLVPTSLHNGDDVFAGIRTIADGISLIKNAAATGVRANTTKMNIGTIAEFLRVDMRAKNSTGSQYLTRGIATFTHSRVESDAPLTLRNLVSAYKTRYDEILARGASIDNMKPLYRKQLFFARKLFNVEKDIVDNLITMDISCGGLQEKGRVSEMVLQEVDIENIDSYRKTRMIAKLIDKGVGDYTAFLKTNFSEIADAITRETRVESVTKAYNVKKKTVVRAFRDLSAAYHERAVRHAWKGMSGLHIVNRIRMGVSNLVMVVSKINPAKANVLAKSGDPTKWLAVLT
    ## 6                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             MLRFVTKNSQDKSSDLFSICSDRGTFVAHNRVRTDFKFDNLVFNRVYGVSQKFTLVGNPTVCFNEGSSYLEGIAKKYLTLDGGLAIDNVLNELRSTCGIPGNAVASHAYNITSWRWYDNHVALLMNMLRAYHLQVLTEQGQYSAGDIPMYHDGHVKIKLPVTIDDTAGPTQFAWPSDRSTDSYPDWAQFSESFPSIDVPYLDVRPLTVTEVNFVLMMMSKWHRRTNLAIDYEAPQLADKFAYRHALTVQDADEWIEGDRTDDQFRPPSSKVMLSALRKYVNHNRLYNQFYTAAQLLAQIMMKPVPNCAEGYAWLMHDALVNIPKFGSIRGRYPFLLSGDAALIQATALEDWSAIMAKPELVFTYAMQVSVALNTGLYLRRVKKTGFGTTIDDSYEDGAFLQPETFVQAALACCTGQDAPLNGMSDVYVTYPDLLEFDAVTQVPITVIEPAGYNIVDDHLVVVGVPVACSPYMIFPVAAFDTANPYCGNFVIKAANKYLRKGAVYDKLEAWKLAWALRVAGYDTHFKVYGDTHGLTKFYADNGDTWTHIPEFVTDGDVMEVFVTAIERRARHFVELPRLNSPAFFRSVEVSTTIYDTHVQAGAHAVYHASRINLDYVKPVSTGIQVINAGELKNYWGSVRRTQQGLRSGRSYDASCNAYRRTYSWRCPRRVDRTGGQCFSRVNVIEPSHGPRPTRYILQEPGTYPAWIRFRNRVQAVSRQKATHFLFDIVPAAVISDFTTSDTSSFAYKSHTYAVNVTALRFSDTYALYVQTDTNMTILSPAARRQASATYSQVAGFCYNTPTVMDSLANILDVDRNIRPKHFKGLRLYTRSKVTAQHHTHLRPDELVEAAAKVSPRRKYYLMCVVELLANLQVDLEAAVATILAYVLTLSEKFVPIFLDSRAIWVGEPGPDALTARLKASSGQIKSIHTADYEPLTELFELAVLMNRGVGHVSWQAEKDHRLNPDVAVVDQARLYSCVRDMFEGSKQTYKYPFMTWDDYTANRWEWVPGGSVHSQYEEDNDYIYPGQYTRNKFITVNKMPKHKISRMIASPPEVRAWTSTKYEWGKQRAIYGTDLRSTLITNFAMFRCEDVLTHKFPVGDQAEAAKVHKRVNMMLDGASSFCFDYDDFNSQHSIASMYTVLCAFRDTFSRNMSDEQAEAMNWVCESVRHMWVLDPDTKEWYRLQGTLLSGWRLTTFMNTVLNWAYMKLAGVFDLDDVQDSVHNGDDVMISLNRVSTAVRIMDAMHRINARAQPAKCNLFSISEFLRVEHGMSGGDGLGAQYLSRSCATLVHSRIESNEPLSVVRVMEADQARLRDLANRTRVQSAVTAIKEQLDKRVTKIFGVGDDVVRDIHTAHRVCGGISTDTWAPVETKIITDNEAYEIPYEIDDPSFWPGVNDYAYKVWKNFGERLEFNKIKDAVARGSRSTIALKRKARITSKKNEFANKSEWERTMYKAYKGLAVSYYANLSKFMSIPPMANIEFGQARYAMQAALDSSDPLRALQVIL
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               Fragments
    ## 1 MSTLLENIFAIINLFK,NTDTLSK,NPDDPDMVDVFMDHLDIDHNK,IDFTEFLLMVFK,LAQAYYESTR,ENLPISGHK,HESSSEK,GYSPTHR,NHHNSSK,EDNEEGVYDYENTGR,WIQSGHIATYYTIQDEAYDTTDSLLEENK,SSSQVNR,HENTSQVPLQESR,DSEGHSEDSER,HSGSASR,NHHGSAWEQSR,ASHGHSADSSR,HAETSSR,GQTASSHEQAR,HGSGHQQSADSSR,GQASSAVSDR,GSSGSQASDSEGHSENSDTQSVSGHGK,QQSHQESTR,SGSSLYQVSTHEQPDSAHGR,TGTSTGGR,QGSHHEQAR,HSASQEGQDTIR,GHPGSSR,QGSHHEQSVNR,SGHSGSHHSHTTSQGR,SDASHGQSGSR,NEEQSGDGTR,HHEASSQADSSR,HSQVGQGQSSGPR,NQGSSVSQDSDSQGHSEDSER,WSGSASR,NHHGSAQEQSR,AGHGHSADSSR,HTQNSSSGQAASSHEQAR,HQLQSADSSR,HSGTGHGQASSAVR,SGSFLYQVSTHK,QSESSHGWTGPSTGVR,QGSHHEQAR,HSASQDGQDTIR,GHPGSSR,QGSHHEQSVDR,SGHSGSHHSHTTSQGR,HHEASSHADISR,HSQAGQGQSEGSR,QGSSVSQDSDSEGHSEDSER,WSGSASR,GSAQEQSR,AGHGHSADSSR,QSGTPHAETSSGGQAASSHEQAR,HQQSADSSR,QASSAVR,SGSFIYQVSTHEQSESAHGR,QGSHHEQAR,HSASQEGQDTIR,QGSHHEQSVDR,SGHSGSHHSHTTSQGR,SDASHGQSGSR,QSGDGSR,HHEAASWADSSR,HSQVGQEQSSGSR,HQGSSVSQDSDSER,HSDDSER,LSGSASR,NHHGSSR,HPGFHQEDR,ASHGHSADSSR,QSGTHHTESSSHGQAVSSHEQAR,HQQSADSSR,HSGIGHR,QASSAVR,GQSGESSGR,SFLYQVSSHEQSESTHGQTAPSTGGR,HSASQDGQDTIR,GHPGSSR,QGSYHEQSVDR,SGHSGYHHSHTTPQGR,SDASHGQSGPR,NEEQSGDGSR,HHEPSTR,HSQVGQGESAGSK,QGSSVSQDR,DSEGHSEDSER,NHYGSAR,ASHGHSAESSR,HAETSSGGQAASSQEQAR,HQQSADSSTDSGTGR,QDSSVVGDSGNR,SGSFLYQVSTHEQSESAHGR,TGPSTGGR,HSASQEGQDTIR,GHPGSSR,QGSHYEQSVDSSGHSGSHHSHTTSQER,QSGDGSR,HHEASSR,HSQVGQGQSSGPR,NQGSSVSQDSDSQGHSEDSER,WSGSASR,NHLGSAWEQSR,HPGSHHEDR,AGHGHSADSSR,HTESSSR,GQAASSHEQAR,HGSHHQLQSADSSR,HSGIGHGQASSAVR,GYSGSQASDSEGHSEDSDTQSVSAQGK,AGPHQQSHK,GQSGESSGR,SGSFLYQVSTHEQSESTHGQSAPSTGGR,QGSHYDQAQDSSR,HSASQEGQDTIR,GHPGPSR,QGSHQEQSVDR,SGHSGSHHSHTTSQGR,EQSGDGSR,HSGSHHHEASSWADSSR,HSLVGQGQSSGPR,GSSVSQDSDSEGHSEDSER,NHHGSAQEQSR,AGHGHSAESSR,QSGTHHAENSSGGQAASSHEQAR,HGSHHQQSADSSR,HSGIGHGQASSAVR,SGSFLYQVSTHEQSESAHGR,TGTSTGGR,HSTSQEGQDTIHGHPGSSSGGR,QGSHYEQLVDR,SGHSGSHHSHTTSQGR,SDASHGHSGSR,NDEQSGDGSR,HHEASSR,ADSSGHSQVGQGQSEGPR,NWGSSFSQDSDSQGHSEDSER,WSGSASR,NHHGSAQEQLR,AGHGHSADSSR,HTQTSSGGQAASSHEQAR,HGSHHQQSADSSR,HSGIGHGQASSAVR,HEQAQDSSR,HSASQDGQDTIR,GHPGSSR,QGYHHEHSVDSSGHSGSHHSHTTSQGR,NEEQSGDGSR,HHEASTHADISR,HSQAVQGQSEGSR,QGSSVSQDSDSEGHSEDSER,WSGSASR,NHHGSAQEQLR,AGHGHSADSSR,HTQTSSGGQAASSHEQAR,HGSHHQQSADSSR,HSGIGHGQASSAVR,HEQAQDSSR,HSASQYGQDTIR,GHPGSSR,QGYHHEHSVDSSGHSGSHHSHTTSQGR,NEEQSGDSSR,HHEASTHADISR,HSQAVQGQSEGSR,QGSSVSQDSDSEGHSEDSER,WSGSASR,GSVQEQSR,AGHGHSADR,HAETSSGGQAASSHEQAR,HQQSADSSR,GQASSAVR,SGSFLYQVSTHEQSESAHGR,QGSHHEQAR,HSASQEGQDTIR,GHPGSSR,QGSHYEQSVDR,SGHSGSHHSHTTSQGR,NDEQSGDGSR,HSWSHHHEASTQADSSR,HSQSGQGQSAGPR,NQGSSVSQDSDSQGHSEDSER,WSGSASR,GSAQEQSR,HPTSHHEDR,AGHGHSAESSR,QSGTHHAENSSGGQAASSHEQAR,HGSHHQQSADSSR,HSGIGHGQASSAVR,SGSFLYQVSTHEQSESAHGR,AGPSTGGR,HSASQEGQDTIR,QGSYHEQSVDR,SGHSGSHHSHTTSQGR,SDASHGQSGSR,NEEQSGDGSR,HHEASTQADSSR,HSQSGQGESAGSR,QGSSVSQDSDSEAYPEDSER,NHHGSSR,HPGSSHR,DTASHVQSSPVQSDSSTAK,EHGHFSSLSQDSAYHSGIQSR,GSPHSSSSYHYQSEGTER,GQSGLVWR,HGSYGSADYDYGESGFR,HSQHGSVSYNSNPVVFK,YYATYINK,DPGLCGHSSDISK,QLGFSQSQR
    ## 2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       MSCNGGSHPR,AESGPDLR,YEVTSGGGGTSR,GVITDQNSDGYCQTGTMSR,HQNQNTIQELLQNCSDCLMR,AELIVQPELK,YGDGIQLTR,ELDECFAQANDQMEILDSLIR,QMGQPCDAYQK,LLQLQEQMR,GGGGYTCQSGSGWDEFTK,HVTSECLGWMR,AEMDMVAWGVDLASVEQHINSHR,GIHNSIGDYR,SAIYQLEEEYENLLK,QLQNIIQATSR,EIMWINDCEEEELLYDWSDK,NTNIAQK,QEAFSIR,MSQLEVK,QESDQLVLNQHPASDK,IEAYMDTLQTQWSWILQITK,CIDVHLK,ENAAYFQFFEEAQSTEAYLK,GLQDSIR,NMPLQHLLEQIK,QVQNLVNK,GDECILK,IEQYYEAILALWNQLYINMK,SLVSWHYCMIDIEK,TIADLELHYQEFIR,NSQGSEMFGDDDK,IQSQFTDAQK,QETWMLMELQK,QIEHCEGR,NLPLADQGSSHHITVK,SVQNDSQAIAEVLNQLK,DMLANFR,YCYLQNEVFGLFQK,LENINGVTDGYLNSLCTVR,ALLQAILQTEDMLK,LTEEETVCLDLDK,SLLATMK,AQQIHSQTSQQYPLYDLDLGK,VTQLTDR,DNYQAFCK,QDSLESMK,FGDSNTVMR,NLHSEISGK,IAELCANSIK,DYELQLASYTSGLETLLNIPIK,TMIQSPSGVILQEAADVHAR,YIELLTR,FLSEMLK,IEVLEEELR,DANSENCNK,FLDQNLQK,YQAECSQFK,LASLEELK,QAELDGK,LTYEIEDEK,NDYDQLQK,ENLGWQK,VLLQEEGTR,EYENELAK,NHYNEEMSNLR,YETEINITK,LNDSILQATEQR,AEENALQQK,ACGSEIMQK,QHLEIELK,QSLEEAAK,AEFQEEAK,WEYENELSK,NNYDEEIISLK,NQFETEINITK,TTIHQLTMQK,EEDTSGYR,AQIDNLTR,SLSEEIK,NTLTQTTENLR,VEEDIQQQK,ATGSEVSQR,QQLEVELR,QSLDDAAK,CLEDENAR,VQYDLQK,ANSSATETINK,VQEQELTR,VEEELNR,TASEDSCK,LEEELEGMR,ITNLTQQLEQASIVK,DVLDGHLR,LSSEVEALR,QLLQEQESVK,LQSLTENLTK,EHLMLEEELR,LEYDDLR,SEADSDK,NATILELR,SQLQISNNR,TLELQGLINDLQR,QALEASNR,NQCTQVVQER,LEDELNR,STLEAETR,QQIQNDLNQWK,ETQSQLETER,ETQTECEWTVDTSK,LVFDGLR,VTAMQLYECQLIDK,SVEEVASEIQPFLR,GAGSIAGASASPK,YSLVEAK,LISPESTVMLLEAQAATGGIIDPHR,LTVDSAIAR,DLIDFDDR,QQIYAAEK,AITGFDDPFSGK,TVSVSEAIK,LLEAQIASGGVVDPVNSVFLPK,NFVDPVTK,VSYVQLK,IEPHTGLLLLSVQK,SMSFQGIR,QPVTVTELVDSGILR,PSTVNELESGQISYDEVGER,DFLQGSSCIAGIYNETTK,LGIYEAMK,PGTALELLEAQAATGFIVDPVSNLR,LPVEEAYK,GLVGIEFK,AVTGYNDPETGNIISLFQAMNK,LLEAQIATGGIIDPK,LPVDIAYK,GYFNEELSEILSDPSDDTK,GFFDPNTEENLTYLQLK,DEETGLCLLPLK,QVQTSQK,VVIVDPETNK,EMSVQEAYK,GLIDYETFK,ELCEQECEWEEITITGSDGSTR,TGSQYDIQDAIDK,SGSLSLTQFADMISLK,NGVGTSSSMGSGVSDDVFSSSR,ISTISSVR,SSSFSDTLEESSPIAAIFDTENLEK,ISITEGIER,GIVDSITGQR,LLEAQACTGGIIHPTTGQK,LSLQDAVSQGVIDQDMATR,AFIGFEGVK,MSAAEAVK,WLPYEAGQR,FLEFQYLTGGLVDPEVHGR,ISTEEAIR,LQDTSSYAK,SMVEDITGLR,LLEAASVSSK,GLPSPYNMSSAPGSR,GSFDATGNSSYSYSYSFSSSSIGH
    ## 3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    LLQGVITVIDVFYQYATQHGEYDTLNK,ELLENEFHQILK,NPNDPDTVDIILQSLDR,VDFTEYLLMIFK,DYCQVSGSK,DDTHQHQEEQEETEK,QESSFSHSSWSAGENDSYSR,PGTESISR,DFSGQHNSYSGQSSSYGEQNSDSHQSSGR,SSSGQSSGYSQHGSGSGHSSGYGQHGSR,SSSGSSSSYGQHGSGSR,QSLGHGR,QSPSHVR,HGSGSGHSSSHGQHGSGSSYSYSR,GHYESGSGQTSGFGQHESGSGQSSGYSK,QSSSYGQHESASR,GQHSSGSGQSPGHGQR,GSGSGQSPSSGQHGTGFGR,QGSSAGSSSSYGQHGSGSR,QSLGHSR,HGSGSGQSPSPSR,QSSSYGPHGYGSGR,GPYESGSGHSSGLGHQESR,YGQQGSGSGQSPSR,SGSGQSSGYSQHGSGSSHSSGYR,SEQHGSSSGLSSSYGQHGSGSHQSSGHGR,QGSGSGHSPSR,HGSGSGQSPGHGQR,GSGSGQSPSYGR,HGSGSGR,HGSGSGQSSGFGHK,SSSGQSSGYTQHGSGSGHSSSYEQHGSR,SEQHGSSSGSSSSYGQHGSGSR,QSLGHGQHGSGSGQSPSPSR,HGSGSGQSSSYGPYR,SGSGWSSSR,GPYESGSGHSSGLGHR,QGSGSGQSPGHGQR,QSPSYGR,HGSGSGR,HGSSSGSSSHYGQHGSGSR,QSSGHGR,QGSGSGHSPSR,HGSGLGHSSSHGQHGSGSGR,GSGSSQSSGYGR,HGAGSGQSPSR,HGSGSGHSSSYGQHGSGSGWSSSSGR,HGSSSGSSSSYGQHGSGSR,QSLGHGQHGSGSGQSPSPSR,HGSGSGQSSSYSPYGSGSGWSSSR,GPYESGSSHSSGLGHR,QGSGSGQSPGHGQR,QSPSYGR,HGSGSGR,YGQHGSGSR,QSSGHGR,QGSGSGQSPSR,HGSGLGHSSSHGQHGSGSGR,HGSGSGHSSSYGQHGSGSGWSSSSGR,GEQHGSSSGSSSSYGQHGSGSR,QSLGHGQHGSGSGQSPSPSR,HGSGSGQSSSYGPYGSGSGWSSSR,GPYESGSGHSSGLGHR,QGSGSGQSPGHGQR,QSPSYGR,HGSGSGR,HGSSSGSSSR,YGQHGSGSR,QSSGHGR,QGSGSGHSPSR,HGSGSGHSSSHGQHGSGSGR,HGSGSGHSSSYGQHGSGSGWSSSSGR,HGSSSGSSSSYGQHGSGSR,QSLGHGQHGSGSGQSPSPSR,HGSGSGQSSSYSPYGSGSGWSSSR,GPYESGSGHSSGLGHR,QGSGSGQSPGHGQR,QSPSYGR,HGSGSGR,HGSSSGSSSHYGQHGSGSR,QSSGHGR,QGSGSGQSPSR,HGSGLGHSSSHGQHGSGSGR,GSGSSQSSGYGR,HGAGSGQSLSHGR
    ## 4                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        SVVTVIDVFYK,QDGECGTLSK,ELHPVLK,NPDDPDTVDVIMHMLDR,LDFTEFLLMIFK,LTMACNK,HQEEESETEEDEEDTPGHK,HSSWSEGEEHGYSSGHSR,QGNLSSSGNQEGSQK,SSCGHSWSGGK,HGSSSVELR,SHISPSR,ESGEEYESGSGSNSWER,GHGGLSCGLETSGHESNSTQSR,LGSSCSGSGDSGR,SHACGYSNSSGCGR,PQNASSSCQSHR,FGGQGNQFSYIQSGCQSGIK,ENGQPQNCGGQWR,FSNSSSSNEFSK,QSSYGQHGSGSSQSSGYGQYGSR,SGQSSYGQHSSGSSQSSGYGQHGSR,SGQSSYGQHSSGSSQSSGYGQHGSR,SSGFGQHESR,SSGFAQHEYR,QETTHGQTINTTR,HSQSGQGQSTQTGSR,SSQSENSDSEVHSK,HSEHIHTQAGSHYPK,QGTTHGQR,HGHSGHGQSTQTGSR,FSHSDATDSEVHSGVSHR,PHSQEQTHSQAGSQHGESESTVHER,HETTYGQTGEATGHGHSGHGQSTQR,GSGHSESSDSEVHSGGSHR,PQSQEQTHGQAGSQHGESGSTVHGR,HGTTHGQTGDTTR,HAHYHHGK,GSSTTGR,ITHGQTGDTTR,HSYSGHEQTTQTGSR,TSHSESTDSEVHSGGSHR,EHTYGQAGSQHEEPEFTVHER,HGTTHGQIGDTTGHSHSGHGQSTQR,HGTTHGQTGDTTR,HAHYHHGLTTQTGSR,GSGHSEYSDSEGYSGVSHTHSGHTHGQAR,SQHGESESIVHER,HGTIHGQTGDTTR,HAHSGHGQSTQTGSR,SSGHSEYSDSEGHSGFSQR,GHTHGQAGSQHGESESIVDER,SGHSESSDSEVHSGGSHTHSGHTHSQAR,SQHGESESTVHK,HQTTHGQTGDTTEHGHPSHGQTIQTGSR,HGTTHGQTADTTR,HGHSGHGQSTQR,HGTTHGQTGDTTR,HAHSGHGQSTQR,GSGHSESSDSEVHSGVSHTHSGHTYGQAR,SQHGESGSAIHGR,QGTIHGQTGDTTR,HGQSGHGQSTQTGSR,SSHSESSDSEVHSEASPTHSGHTHSQAGSR,HGQSGSSGHGR,QGTTHGQTGDTTR,HAHYGYGQSTQR,QPGSTVHGR,LETTHGQTGDTTR,HGHSGYGQSTQTGSR,ASHFQSHSSER,HGSSQVWK,HGSYGPAEYDYGHTGYGPSGGSR,SISNSHLSWSTDSTANK
    ## 5                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           MSSLLNSLLPEYFK,TNLNINSSR,VQYGFNAR,IDMQYEDDSGTR,PNAFMSNTVAFIGNYEGIIVDDIPILDGLR,ADIFDTHGDLDMGLVEDALSK,NVPTYTAYASELLYK,NLTSLFYNMLR,DAIFYDNGHACLLNR,ATCILGQACSEFAPLAPFSIAHYSPQLTR,LFVNAPAGIEPSSGR,YTHEDVK,VSDSDAINK,YTFTELSVQGR,GIIVDLTLVAMR,EISLPYPVSCGLTR,TDALLQGTEIHVPVVVK,DIDMPQYYNAIDK,DVIEGQETVIK,AVFNDAR,MMGNDVTATDLVTGR,VSNWADNSSGR,FLYTDVK,YEGQTAFLVDMDTVK,DHCWVSIVDPNGTMNLSYK,PLYMTGGSVR,TIATGNYR,AMDETLR,LDFSCSSLR,FVGQQYAILTPSGTTTDIR,GTNQSYR,IGVEDDEDLDIGTVK,YIVPLYLNGDNVAQNCLEATHVLIK,ACSIANR,IVDDGEGHCFTQQGLAQQWIFHR,GEMIFVK,IGQLNAYYVDYK,NVTNYSLK,TAAQVGATISNNLR,HGFVDNQQDAYTR,LVANYSDTR,DNFTYNYNMEK,ITQYHHTHVR,LEGYEALLAMMLDR,FNNIESTHVTFFTYLR,EVFISLVLNYNGLGR,QAQGTVK,LFELNVLENGVDEEVDWEK,TVNISYAK,ELFIMAR,AVIMPGGSVHSQHPVEQDVIR,GVASVMPYK,PEIHAYTSTK,IAGTHSLIPFCYDFDDFNSQHSK,EAMQAVIDAWISVYHDK,LTDDQIEAAK,MVAHQPNTGETYDVK,GTLFSGWR,TIADGISLIK,NAAATGVR,MNIGTIAEFLR,NSTGSQYLTR,GIATFTHSR,VESDAPLTLR,NLVSAYK,YDEILAR,GASIDNMK,DIVDNLITMDISCGGLQEK,VSEMVLQEVDIENIDSYR,GVGDYTAFLK,TNFSEIADAITR,DLSAAYHER,GMSGLHIVNR,MGVSNLVMVVSK
    ## 6                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             SSDLFSICSDR,GTFVAHNR,FDNLVFNR,VYGVSQK,FTLVGNPTVCFNEGSSYLEGIAK,YLTLDGGLAIDNVLNELR,STCGIPGNAVASHAYNITSWR,WYDNHVALLMNMLR,AYHLQVLTEQGQYSAGDIPMYHDGHVK,LPVTIDDTAGPTQFAWPSDR,STDSYPDWAQFSESFPSIDVPYLDVR,PLTVTEVNFVLMMMSK,TNLAIDYEAPQLADK,HALTVQDADEWIEGDR,VMLSALR,LYNQFYTAAQLLAQIMMK,PVPNCAEGYAWLMHDALVNIPK,YPFLLSGDAALIQATALEDWSAIMAK,PELVFTYAMQVSVALNTGLYLR,VAGYDTHFK,VYGDTHGLTK,HFVELPR,LNSPAFFR,SVEVSTTIYDTHVQAGAHAVYHASR,INLDYVK,PVSTGIQVINAGELK,NYWGSVR,SYDASCNAYR,TGGQCFSR,VNVIEPSHGPR,YILQEPGTYPAWIR,ATHFLFDIVPAAVISDFTTSDTSSFAYK,SHTYAVNVTALR,FSDTYALYVQTDTNMTILSPAAR,VTAQHHTHLR,PDELVEAAAK,FVPIFLDSR,AIWVGEPGPDALTAR,ASSGQIK,SIHTADYEPLTELFELAVLMNR,GVGHVSWQAEK,LNPDVAVVDQAR,DMFEGSK,YPFMTWDDYTANR,WEWVPGGSVHSQYEEDNDYIYPGQYTR,MIASPPEVR,AIYGTDLR,STLITNFAMFR,CEDVLTHK,FPVGDQAEAAK,NMSDEQAEAMNWVCESVR,HMWVLDPDTK,LQGTLLSGWR,LTTFMNTVLNWAYMK,LAGVFDLDDVQDSVHNGDDVMISLNR,IMDAMHR,CNLFSISEFLR,VEHGMSGGDGLGAQYLSR,SCATLVHSR,IESNEPLSVVR,VMEADQAR,VQSAVTAIK,IFGVGDDVVR,DIHTAHR,VCGGISTDTWAPVETK,IITDNEAYEIPYEIDDPSFWPGVNDYAYK,GLAVSYYANLSK,FMSIPPMANIEFGQAR,YAMQAALDSSDPLR

``` r
read_DIANN <- function(file_name) {

  #Selecting proteins that meet q-value threshold
  df <- read_parquet(file_name)
  df <- df[df$Lib.PG.Q.Value <= 0.01,]

  # For a Protein.Group to be listed by itself, it therefore requires at
  #     least one peptide that can UNIQUELY map to it.  
  df <- df[!grepl("contam", df$Protein.Names),]
  df <- df[!is.na(df$Precursor.Quantity),]
  df <- df[!df$Precursor.Quantity==0,]

  df <- df[c("Run", "Protein.Group", "Genes", "Stripped.Sequence",
             "Precursor.Quantity")]
  df <- df %>% 
    group_by(Run, `Protein.Group`, Genes, `Stripped.Sequence`) %>%
    dplyr::summarise(`Precursor.Quantity` = sum(`Precursor.Quantity`),
                     .groups = "drop") %>%
    as.data.frame

  return(df) }


XLA_Abs_DIANN <- read_DIANN("Data/DIA/Frog_AbsQ_report.parquet")
Dmel_Abs_DIANN <- read_DIANN("Data/DIA/Fly_AbsQ_report.parquet")

head(XLA_Abs_DIANN)
```

    ##                     Run Protein.Group         Genes Stripped.Sequence
    ## 1 ORC_05940_Frog_T5-RT0    XBgroup104 XBproteins104      AIIIFVPVPQLK
    ## 2 ORC_05940_Frog_T5-RT0    XBgroup104 XBproteins104      AQLRELNITAAK
    ## 3 ORC_05940_Frog_T5-RT0    XBgroup104 XBproteins104       DVVFEFPEFQL
    ## 4 ORC_05940_Frog_T5-RT0    XBgroup104 XBproteins104          EIEVGAGR
    ## 5 ORC_05940_Frog_T5-RT0    XBgroup104 XBproteins104         EIEVGAGRK
    ## 6 ORC_05940_Frog_T5-RT0    XBgroup104 XBproteins104          ELNITAAK
    ##   Precursor.Quantity
    ## 1          633052416
    ## 2            2469481
    ## 3            7774395
    ## 4           35786228
    ## 5           38811700
    ## 6          120480145

``` r
head(Dmel_Abs_DIANN)
```

    ##                    Run Protein.Group        Genes   Stripped.Sequence
    ## 1 ORC_05948_Fly_T6-T0A    A0A021WW64  FBgn0058460          EALYQGIFHR
    ## 2 ORC_05948_Fly_T6-T0A    A0A023GRW3 Dmel\\CG4896       AISDESDYVDFQK
    ## 3 ORC_05948_Fly_T6-T0A    A0A023GRW3 Dmel\\CG4896 DSFGATAAMPISSTNVGSR
    ## 4 ORC_05948_Fly_T6-T0A    A0A023GRW3 Dmel\\CG4896     LNDYVPEAGPPAISK
    ## 5 ORC_05948_Fly_T6-T0A    A0A023GRW3 Dmel\\CG4896         MGWSEGQGLGK
    ## 6 ORC_05948_Fly_T6-T0A    A0A0B4JCU3        ReepA       ERGYSAVLQLGSK
    ##   Precursor.Quantity
    ## 1           309537.9
    ## 2           195800.7
    ## 3           301966.9
    ## 4           928276.4
    ## 5           515928.5
    ## 6          1553789.6

``` r
merge_DIA_FASTA <- function(data_df, fasta_df){
  
  all_data <- data_df %>%
    #Exploding dataframe by treating all ";" as individual entries
    mutate(ProteinToken = strsplit(`Protein.Group`, ";", fixed = TRUE)) %>%
    tidyr::unnest(ProteinToken) %>% #Unnesting exploded df
    #Merging individual tokens with FASTA file
    inner_join(fasta_df[c("Protein_ID", "Theo_Peps", "Seq_Length")],
               by = c("ProteinToken" = "Protein_ID")) %>% 
    
    #Grouping by what we started with and determining med,min,max
    #     for all peptides in Protein.Group
    group_by(across(all_of(names(data_df)))) %>%
    dplyr::summarise(
      Theo_Peps_Median = median(Theo_Peps, na.rm = TRUE),
      Theo_Peps_Min = min(Theo_Peps, na.rm = TRUE),
      Theo_Peps_Max = max(Theo_Peps, na.rm = TRUE),
      Length_Median = median(Seq_Length, na.rm = TRUE),
      Length_Min = min(Seq_Length, na.rm = TRUE),
      Length_Max = max(Seq_Length, na.rm = TRUE),
      .groups = "drop") %>%
    
    #Getting protein groups that are within 20% spread in theoretical peptides
    filter(
      Theo_Peps_Max <= 1.2 * Theo_Peps_Min,
      Length_Max <= 1.2 * Length_Min) %>%

    #Returning back original columns and modified theo_peps and sequence length    
    transmute(
      !!!syms(names(data_df)),
      Theo_Peps = Theo_Peps_Median,
      Length = Length_Median)    
  
  return(all_data) }

XLA_Abs_FASTA_merge <- merge_DIA_FASTA(XLA_Abs_DIANN, XLA_Digest)
Dmel_Abs_FASTA_merge <- merge_DIA_FASTA(Dmel_Abs_DIANN, Dmel_Digest)

head(XLA_Abs_FASTA_merge)
```

    ## # A tibble: 6 × 7
    ##   Run         Protein.Group Genes Stripped.Sequence Precursor.Quantity Theo_Peps
    ##   <chr>       <chr>         <chr> <chr>                          <dbl>     <dbl>
    ## 1 ORC_05940_… XBgroup104    XBpr… AIIIFVPVPQLK              633052416         10
    ## 2 ORC_05940_… XBgroup104    XBpr… AQLRELNITAAK                2469481.        10
    ## 3 ORC_05940_… XBgroup104    XBpr… DVVFEFPEFQL                 7774395         10
    ## 4 ORC_05940_… XBgroup104    XBpr… EIEVGAGR                   35786228         10
    ## 5 ORC_05940_… XBgroup104    XBpr… EIEVGAGRK                  38811700         10
    ## 6 ORC_05940_… XBgroup104    XBpr… ELNITAAK                  120480145.        10
    ## # ℹ 1 more variable: Length <dbl>

``` r
head(Dmel_Abs_FASTA_merge)
```

    ## # A tibble: 6 × 7
    ##   Run         Protein.Group Genes Stripped.Sequence Precursor.Quantity Theo_Peps
    ##   <chr>       <chr>         <chr> <chr>                          <dbl>     <dbl>
    ## 1 ORC_05948_… A0A021WW64    "FBg… EALYQGIFHR                   309538.        12
    ## 2 ORC_05948_… A0A023GRW3    "Dme… AISDESDYVDFQK                195801.        44
    ## 3 ORC_05948_… A0A023GRW3    "Dme… DSFGATAAMPISSTNV…            301967.        44
    ## 4 ORC_05948_… A0A023GRW3    "Dme… LNDYVPEAGPPAISK              928276.        44
    ## 5 ORC_05948_… A0A023GRW3    "Dme… MGWSEGQGLGK                  515928.        44
    ## 6 ORC_05948_… A0A0B4JCU3    "Ree… ERGYSAVLQLGSK               1553790.        40
    ## # ℹ 1 more variable: Length <dbl>

``` r
fly_yolk_set <- c("P02843", "P02844", "P06607") #Yp1, Yp2, Yp3
names(fly_yolk_set) <- c("Yp1", "Yp2", "Yp3")

frog_yolk_set <- c("XBmRNA35014|XBXL10_1g18990", "XBmRNA39801|XBXL10_1g21498",
                   "XBmRNA35013|XBXL10_1g18989", "XBmRNA35012|XBXL10_1g18988",
                   "XBmRNA67639|XBXL10_1g35803")
names(frog_yolk_set) <- c("vtgb1.L", "vtgb1.S", "vtga1.L", "vtga2.L", "serpina")
```

``` r
Abs_Conc_per_run <- function(data_df, run_no, cell_conc) {
  sub_df <- data_df[data_df$Run == run_no,]
  sub_df <- sub_df %>% group_by(`Protein.Group`, Genes, Theo_Peps) %>%
    dplyr::summarise(Sum_Pep_Int = sum(`Precursor.Quantity`),
                     .groups = "drop")

  sub_df <- sub_df[!sub_df$Theo_Peps == 0,] #Removing weird proteins
  
  #Estimating molar units based on peptides identified    
  sub_df["iBAQ_Value"] <- sub_df$Sum_Pep_Int / sub_df$Theo_Peps

  return(sub_df) }

unnest_data <- function(data_df) {
  
  data_df <- data_df %>%
      #Exploding dataframe by treating all ";" as individual entries
      mutate(Ind_Protein_ID = strsplit(`Protein.Group`, ";", fixed = TRUE)) %>%
      unnest(Ind_Protein_ID) %>% #Unnesting exploded df
      select(Ind_Protein_ID, everything())
  
  return(data_df) }
```

``` r
XLA_T5_RT0_Abs <- Abs_Conc_per_run(XLA_Abs_FASTA_merge, 'ORC_05940_Frog_T5-RT0', 2)
XLA_T5_N12_Abs <- Abs_Conc_per_run(XLA_Abs_FASTA_merge, 'ORC_05942_Frog_T5-N12', 2)
XLA_T6_RT0_Abs <- Abs_Conc_per_run(XLA_Abs_FASTA_merge, 'ORC_05944_Frog_T6-RT0', 2)
XLA_T6_N12_Abs <- Abs_Conc_per_run(XLA_Abs_FASTA_merge, 'ORC_05946_Frog_T6-N12', 2)

XLA_T5_RT0_Abs <- unnest_data(XLA_T5_RT0_Abs)
XLA_T5_N12_Abs <- unnest_data(XLA_T5_N12_Abs)
XLA_T6_RT0_Abs <- unnest_data(XLA_T6_RT0_Abs)
XLA_T6_N12_Abs <- unnest_data(XLA_T6_N12_Abs)

head(XLA_T5_RT0_Abs)
```

    ## # A tibble: 6 × 6
    ##   Ind_Protein_ID Protein.Group Genes         Theo_Peps Sum_Pep_Int iBAQ_Value
    ##   <chr>          <chr>         <chr>             <dbl>       <dbl>      <dbl>
    ## 1 XBgroup104     XBgroup104    XBproteins104        10 1687260692. 168726069.
    ## 2 XBgroup106     XBgroup106    XBproteins106        10    7522030.    752203.
    ## 3 XBgroup109     XBgroup109    XBproteins109        15  202097062.  13473137.
    ## 4 XBgroup111     XBgroup111    XBproteins111        15     908997.     60600.
    ## 5 XBgroup112     XBgroup112    XBproteins112        18    2612846.    145158.
    ## 6 XBgroup113     XBgroup113    XBproteins113         4   63020564.  15755141.

``` r
#-------------------

Dmel_T6_T0A_Abs <- Abs_Conc_per_run(Dmel_Abs_FASTA_merge, 'ORC_05948_Fly_T6-T0A', 1.6)
Dmel_T6_T0B_Abs <- Abs_Conc_per_run(Dmel_Abs_FASTA_merge, 'ORC_05950_Fly_T6-T0B', 1.6)
Dmel_T6_N16_Abs <- Abs_Conc_per_run(Dmel_Abs_FASTA_merge, 'ORC_05952_Fly_T6-N16', 1.6)
Dmel_T7_T0A_Abs <- Abs_Conc_per_run(Dmel_Abs_FASTA_merge, 'ORC_05954_Fly_T7-T0A', 1.6)
Dmel_T7_T0B_Abs <- Abs_Conc_per_run(Dmel_Abs_FASTA_merge, 'ORC_05956_Fly_T7-T0B', 1.6)
Dmel_T7_N16_Abs <- Abs_Conc_per_run(Dmel_Abs_FASTA_merge, 'ORC_05958_Fly_T7-N16', 1.6)

Dmel_T6_T0A_Abs <- unnest_data(Dmel_T6_T0A_Abs)
Dmel_T6_T0B_Abs <- unnest_data(Dmel_T6_T0B_Abs)
Dmel_T6_N16_Abs <- unnest_data(Dmel_T6_N16_Abs)
Dmel_T7_T0A_Abs <- unnest_data(Dmel_T7_T0A_Abs)
Dmel_T7_T0B_Abs <- unnest_data(Dmel_T7_T0B_Abs)
Dmel_T7_N16_Abs <- unnest_data(Dmel_T7_N16_Abs)

head(Dmel_T6_T0A_Abs)
```

    ## # A tibble: 6 × 6
    ##   Ind_Protein_ID Protein.Group Genes            Theo_Peps Sum_Pep_Int iBAQ_Value
    ##   <chr>          <chr>         <chr>                <dbl>       <dbl>      <dbl>
    ## 1 A0A021WW64     A0A021WW64    "FBgn0058460"           12     309538.     25795.
    ## 2 A0A023GRW3     A0A023GRW3    "Dmel\\CG4896"          44    1941972.     44136.
    ## 3 A0A0B4JCU3     A0A0B4JCU3    "ReepA"                 40   18476376.    461909.
    ## 4 A0A0B4JCZ0     A0A0B4JCZ0    "CG15684"               44    9425127.    214207.
    ## 5 A0A0B4JCZ3     A0A0B4JCZ3    "pre-mod(mdg4)-…         4   10140653.   2535163.
    ## 6 A0A0B4JCZ8     A0A0B4JCZ8    "mim"                   64   17841774.    278778.

Citation for frog absolute proteomics:

Gurdon and Wickens // Evi’s Paper

280 ug yolk protein 45 ug non-yolk protein

``` r
split_frog_yolk <- function(df, yolk_amount, nonyolk_amount) {

  df <- merge(df, frog_AA.df, by.x="Ind_Protein_ID", by.y="Protein_ID")

  sub_yolk <- df[df$Ind_Protein_ID %in% frog_yolk_set,]
  sub_nonyolk <- df[!df$Ind_Protein_ID %in% frog_yolk_set,]
  
  #Calculate amount of protein in milligrams
  sub_yolk["iBAQ_WeightFraction"] <- sub_yolk$iBAQ_Value * sub_yolk$kDa #mult by kda for weight
  #Calculate mass fraction
  sub_yolk["iBAQ_WeightFraction"] <- sub_yolk$iBAQ_WeightFraction / sum(sub_yolk$iBAQ_WeightFraction)  
  sub_yolk["iBAQ_Weight_mg"] <- sub_yolk$iBAQ_WeightFraction * (yolk_amount/1000) #mg calc
  #Calculate mg/mL then divide by kDa -> mM
  sub_yolk["iBAQ_Conc"] <- sub_yolk$iBAQ_Weight_mg / 0.001  
  sub_yolk["iBAQ_Conc"] <- sub_yolk$iBAQ_Conc / sub_yolk$kDa  
      
  #--------------

  #Calculate amount of protein in milligrams
  sub_nonyolk["iBAQ_WeightFraction"] <- sub_nonyolk$iBAQ_Value * sub_nonyolk$kDa #mult by kda for weight
  #Calculate mass fraction
  sub_nonyolk["iBAQ_WeightFraction"] <- sub_nonyolk$iBAQ_WeightFraction /
    sum(sub_nonyolk$iBAQ_WeightFraction)  
  sub_nonyolk["iBAQ_Weight_mg"] <- sub_nonyolk$iBAQ_WeightFraction * (nonyolk_amount/1000) #mg calc
  #Calculate mg/mL then divide by kDa -> mM
  sub_nonyolk["iBAQ_Conc"] <- sub_nonyolk$iBAQ_Weight_mg / 0.001  
  sub_nonyolk["iBAQ_Conc"] <- sub_nonyolk$iBAQ_Conc / sub_nonyolk$kDa  

  #--------------
  
  recalc_df <- rbind(sub_yolk, sub_nonyolk)
  recalc_df <- recalc_df %>% arrange(desc(iBAQ_Conc))
  recalc_df[c("iBAQ_WeightFraction", "iBAQ_Weight_mg")] <- NULL

  recalc_df <- recalc_df %>% 
    relocate(Sequence, .after = last_col())  
    
  return(recalc_df) }

#These are from Evi's paper assuming:
#  ~ non-yolk protein amount doesn't increase until after NF30
XLA_T5_RT0_Abs <- split_frog_yolk(XLA_T5_RT0_Abs, 280, 45)
XLA_T5_N12_Abs <- split_frog_yolk(XLA_T5_N12_Abs, 280, 50)
XLA_T6_RT0_Abs <- split_frog_yolk(XLA_T6_RT0_Abs, 280, 45)
XLA_T6_N12_Abs <- split_frog_yolk(XLA_T6_N12_Abs, 280, 50)

head(XLA_T5_RT0_Abs)
```

    ##               Ind_Protein_ID              Protein.Group        Genes Theo_Peps
    ## 1 XBmRNA67639|XBXL10_1g35803 XBmRNA67639|XBXL10_1g35803   serpina6.L        27
    ## 2 XBmRNA35014|XBXL10_1g18990 XBmRNA35014|XBXL10_1g18990      vtgb1.L        87
    ## 3 XBmRNA35013|XBXL10_1g18989 XBmRNA35013|XBXL10_1g18989      vtga1.L       105
    ## 4 XBmRNA35012|XBXL10_1g18988 XBmRNA35012|XBXL10_1g18988      vtga2.L       101
    ## 5 XBmRNA39801|XBXL10_1g21498 XBmRNA39801|XBXL10_1g21498      vtgb1.S        91
    ## 6                  XBgroup52                  XBgroup52 XBproteins52        20
    ##   Sum_Pep_Int iBAQ_Value       kDa  iBAQ_Conc
    ## 1 16575471163  613906339  49.52268 1.18752305
    ## 2 20835884583  239492926 203.11782 0.46326834
    ## 3 18288864612  174179663 202.34768 0.33692821
    ## 4 11020786564  109116699 201.54506 0.21107225
    ## 5  3834001832   42131888 200.92837 0.08149873
    ## 6 15707480683  785374034  41.98385 0.01719761
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    Sequence
    ## 1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      MHLLVYLSLFFALALASVTEISLDNKHRHRHEQQGHHDSAKHGHQKDKQQQEQIKNDEGKLTKEEKILSEENSDFSVNLFNQLSTESKRSPRKNIFFSPISISAAFYMLALGAKSETHQQILKGLSFNKKKLSESQVHEAFKRLIEDSNNPMKAHQFTIGNALFVEQTVNILKGFEENVKHYYQAGVFPMNFKDPDNAKKQLNNYVKDKTHGVIQEMIRELDSNTEMVLVNYVLYKGEWANNFNPTLTQKSLFSVDKNTNVTVQMMNRLGLYRTYQDDDCKIIELPYKNDTAMLLVVPQLGKIQELVLTSKLINHWYESLATSIVDLYMPTFSISGKVVLKDTLRKMGISDIFTDKADLTGISEQIKLKVSMASHNAVLNVNEFGTEAVGATSAQASPTKLFPPFLIDSPFLVMIYSRTLGSQLFMGKVMDPTNAQ
    ## 2 MRGIILALLLALAGSEKSQYEPFFSESKTYVYNYEGIILNGIPENGLARSGIKLNCKAEISGYAQRSYMLKIRDPEIKEFNGLWPKDPFTRATKLSQILAEHLTRPVKFEYRNGQVGNIFAAEDVSETALNIQRGILNLLQLTIKKSQNVYDLQEDSISGVCHTKYTIQEDKKAERLTVMKSRDFDRCESRVEKKIGVAYLETCTQCKKKVQHLRGSATFTYKLQNKDQGALIVEASAKQVYRYAPYYELYGAAVMEARQILTLAETKSNRVSESQVQLTKRGGLNYHFWQEEHQMPFHLLKTKNAESQVVETLQHLVQNNQELVHGESTSKFFELVQLLRTTNQENIQAIWKQFSDRPQHRHWLVNAIPMAGSVDALKFIKQMIHNQDMTHQEAATAILMAMHYTRASQRAVEIAADLVTDSRVKNAPSLYKSALLAYGSLVYKYCVDSQSCPEAALQPLHEFAAEAANNAHQEEITLALKAIGNAGQPASIKRIQKFLPGFSSGASRLPVRLQTDAVMALRNIAIKDPRKVQEILLQLFMDRDVHPEVRMITCVALFETNPGIAIVTAVANAVLRESKDNMQLASFTYSHIKSLSKSMMPEFYNLASACRVALNIMNPVFDQLSFRYSKAFHVDTFNYPMMGGLALNVFFINSPNTAFPVFALIKVRECSFGASTDVFEVALRAEGLQEVIRKQNIQFSEFSMHKKIVKILKALMAYKEAPANVPLVSASIKLLGQEVNFNQATRDTIQNAMKLLNEPAERHTMVKKILNKLLNGVTGQWSHPMLVAEARHIVPTCVGLPLEIATHSAAVANAALNVDMKVSPSPSNDFSLSQLLESNIQLHTDMSPSVLIHLVNTMGVNTPLIQSGLEFHGKVRIQTPAKFTAKLDMKDKNFKIETPPCQREVELVTVRHQMFAVTRNVGQQDSEKRMLVMPKQHGPSTMEQHFEKSGRTSAESASMMEVSSEMGIRNKNAQHQHILNNVESYRYCFRLSKLNINACLQSKLNTVQGDRRSFLYENIGEHEIKFLLKPTQNEPAIEKLQVEITAGPRASSKIIALAEVEEKEGENLDESAVQKRLKIILGINEEEYSARNSSLVRKASKKRKVKKTKDAGETDVLEAERKWQQSSSSSSSSSSPSSSSSSSSSSSSSSSSSSSSSSKSSSSSSSSSPQRNSRNKGIKVDVLSRVSKRQQEKNKDTHGQTGRNKHRASSSESSSSESSSSESSRSSESSSSSSSESSSSSESSRRHGKSGKQQQGDREQQREKQTHQNKQKQQKHRSRHPQDPSQMPPQLYKLKFQRPYKQQYKGMSGQSESSESSSSSSSSSSSSESARSPHHKKFLGDKKSPTVVITAQAVRNDNVKQGYQVILYTDDTTRPKMQAYLVDITGKSRWRACAEAFVRNPHEAQAELKWGQNCKQYKMEFIAQTGIFGNQPALRLRVNWPRIPSKLKSVGEEVGEFLPGAAYMMGLNTKFQRNPSKQAKFIFSLSSPRTMNAVIQVPWATCYYQAINLPFAVTAAYHLHSSDIQAPTWNIFANSYSSMVDRLSGECTMSSDQIKTFNEVKFNYTLPGSCFHVLAQDCSTELKFLVLAKKSPESSELRDINIKLADYDIDMYTSSEEIRLKVNGLEVPKEKLPYIASTNPIVEIKMEEGALVLTAPVFGIEQLYFNGKLTKLRTAVWMRGGTCGICGRHDDETEREYQMPDGYAARDENSFTHSWILPEEPCSGTCKLQHTSVKLEKTIHIDDQKSKCQSVRPVLRCVKGCSPTSVTPVTVGFHCLPFDASTDMLESQNIMGKTEDLVETVDAHTACSCEALKCAA
    ## 3  MRGIILALLLAIAGSERTQIEPVFSESKTSVYNYEAVILNGFPESGLSRAGIKINCKVEISAYAQRSYFLKIQSPELKEYNGVWPKDPFTRSSKLTQALAEQLTKPARFQYSNGRVGDIFVSDDVSDTVLNIYRGILNLLQVTIKKSQNVYDLQEPSIGGICHTRYVLQEDSRGDRISIIKATDFNNCHEKVFKGIGFELTESCESCKQFHRNIRGTATYTYKLKGRDQGSVIMEVTARQVLQFTPFAERNGAATMESRQILVWTGSKSGPIHPPQIQLKNRGNLHYQFASELHQMPFHLMKTKSPEAQAVEVLQHLMKDTQQQIREDTPVKFLQLIQLLRSSDFEDLQALWKQFAQRTQYRRWLLDAIPMAGTVDCLKFVKQLIHNEELNPYEAAVTITLALRSARPSQRAFQISTDFVQDSKVQKYSTVHKAAILAYGTMVKKYCDQLSSCPEQALEPLHDLAAEAADKGHAEIIALALKALGNAGQRESIKRIQKFLPGFSSSAYQLPVRIQTDAVMALRNIAKQDPQKVQEILLQIFMDRDVRTEVRMMACLALFETKPHLATVTTIANVAARESKTNLQLASFTYSQLKALAKSSVPHLEPLAAACSVALKILNPSLDNLGYRYSKVMRVDTFKYNLMAGGAAKVYVMNSANTMFPVFILAKFREYLAGVESDILEIGIRGEGIEEILRKQNIQFADFPMRKKISQILKTLLGFKGLPSQMPLISGYVKVLGQEIAFTEVNKEVIQNIIRALNEPAERHTMIRNILNKFLGGVAGHYSQALMTGEYRYLIPTTVGLPAEFSLYHSAIVNAAVNSDIKVKPSPSKDFSVAQLMESQIQLNADVNPSFMFYKVGTMGINSRLIQAGFEFHGKVHARLPAKFTAYLDMKDRNFKIETPPCQQENHLVELRAQTFAVTRNIEDLDAARKNLVLPRNSEQNILKKQFQSTGRSSAEGASMMEDSSEMGPRKYSTEPGHPQFAPNINSYDACTKLSKAGVHFCIQCKTHNAASHRNTFLYHIIGEHEFKLTMKPAHTEGAIEKLQLEIAAGPKAASKLMGLIEVEETEGEKMGDSAVQKRLKLILGIDDSVMNTNETAVLRSKKSKKERKEHKRPHDAEVVEAKRQQSSSSSSSSPSSSSSSSSSSSSSSSSSSSSSSSSSSSYRDRPNRRRQKSGQQGESSSSSSQKQDKHKPQENRKHGQKAQSSSSSSSSSSSSSSSSSSSSSSRSSSSSSSESSPSRSHNKQHENEQRKKQPKQNQQNQQNKNKYGSSSSSSGSSSSSSSSEMWTKKKHHRHFYDLHFRRAARTKDKKQRGSQWSSSSESSSSSSESANTNRRKTNFLGDKESPVFVATFRAVLNDNTKQGYQMVVYQDQHSSKQQIQAYVMDISKSRWATCLNAVVNNQYEAQASLKWGQNCQDYKINFKAQTGNFGNQPALKVTANWPKIPSKLKSTGKYVAEYVPGAMYMMGFQGEYKRNSQRQVKLVFALSSPRTCDVVIQTPRVTVYYRALRLPVPIPVGHHTKENVLQTPTWNIFAEAPQIIMDSLQGECKVAQDQITTFNGVDLVSALPENCYHVLAQDCSPEMKFMVLMRNSKESPNHKAINVKLGTYDIDMHYSSEMLKMKINGIELSEERLPYKSFEDPTVELKKKGNGVSLSAPEYGINSLDYDGLTFKLKAAIWMKGKTCGICGHNDDENEKELQMPDGSVAKGPMSFIHSWILPAESCSEGCNLKRSLVQLEKEIDGAMAKCYSVQPVLRCAKGCSPLKTEKVSTGFHCLPTDASLDLPTSQTRLEKSEDFSESVDAHTACSCETSQCAA
    ## 4           MKGIVLALLLALAGSERTHIEPVFSESKISVYNYEAVILNGFPESGLSRAGIKINCKVEISAYAQRSYFLKIQSPEIKEYNGVWPKDPFTRSSKLTQALAEQLTKPARFEYSNGRVGDIFVADDVSDTVANIYRGILNLLQVTIKKSQDVYDLQESSVGGICHTRYVIQEDKRGDQIRIIKSTDFNNCQDKVSKTIGLELAEFCHSCKQLNRVIQGAATYTYKLKGRDQGTVIMEVTARQVLQVTPFAERHGAATMESRQVLAWVGSKSGQLTPPQIQLKNRGNLHYQFASELHQMPIHLMKTKSPEAQAVEVLQHLVQDTQQHIREDAPAKFLQLVQLLRASNFENLQALWKQFAQRTQYRRCLLDALPMAGTVDCLKFIKQLIHNEELTTQEAAVLITFAMRSARPGQRNFQISADLVQDSKVQKYSTVHKAAILAYGTMVRRYCDQLSSCPEHALEPLHELAAEAANKGHYEDIALALKALGNAGQPESIKRIQKFLPGFSSSADQLPVRIQTDAVMALRNIAKEDPRKVQEILLQIFMDRDVRTEVRMMACLALFETRPGLATVTAIANVAARESKTNLQLASFTFSQMKALSKSSVPHLEPLAAACSVALKILNPSLDNLGYRYSKVMRVDTFKYNLMAGAAAKVFIMNSANTMFPVFILAKFREYTSLVENDDIEIGIRGEGIEEFLRKQNIQFANFPMRKKISQIVKSLLGFKGLPSQVPLISGYIKLFGQEIAFTELNKEVIQNTIQALNQPAERHTMIRNVLNKLLNGVVGQYARRWMTWEYRHIIPTTVGLPAELSLYQSAIVHAAVNSDVKVKPTPSGDFSAAQLLESQIQLNGEVKPSVLVHTVATMGINSPLFQAGIEFHGKVHAHLPAKFTAFLDMKDRNFKIETPPFQQENHLVEIRAQTFAFTRNIADLDSARKTLVVPRNNEQNILKKHFETTGRTSAEGASMMEDSSEMGPKKYSAEPGHHQYAPNINSYDACTKFSKAGVHLCIQCKTHNAASRRNTIFYQAVGEHDFKLTMKPAHTEGAIEKLQLEITAGPKAASKIMGLVEVEGTEGEPMDETAVTKRLKMILGIDESRKDTNETALYRSKQKKKNKIHNRRLDAEVVEARKQQSSLSSSSSSSSSSSSSSSSSSSSSSSSSPSSSSSSSYSKRSKRREHNPHHQRESSSSSSQEQNKKRNLQENRKHGQKGMSSSSSSSSSSSSSSSSSSSSSSSSSSSSEENRPHKNRQHDNKQAKMQSNQHQQKKNKFSESSSSSSSSSSSEMWNKKKHHRNFYDLNFRRTARTKGTEHRGSRLSSSSESSSSSSESAYRHKAKFLGDKEPPVLVVTFKAVRNDNTKQGYQMVVYQEYHSSKQQIQAYVMDISKTRWAACFDAVVVNPHEAQASLKWGQNCQDYKINMKAETGNFGNQPALRVTANWPKIPSKWKSTGKVVGEYVPGAMYMMGFQGEYKRNSQRQVKLVFALSSPRTCDVVIRIPRLTVYYRALRLPVPIPVGHHAKENVLQTPTWNIFAEAPKLIMDSIQGECKVAQDQITTFNGVDLASALPENCYNVLAQDCSPEMKFMVLMRNSKESPNHKDINVKLGEYDIDMYYSADAFKMKINNLEVSEEHLPYKSFNYPTVEIKKKGNGVSLSASEYGIDSLDYDGLTFKFRPTIWMKGKTCGICGHNDDESEKELQMPDGSVAKDQMRFIHSWILPAESCSEGCNLKHTLVKLEKAIATDGAKAKCYSVQPVLRCAKGCSPVKTVEVSTGFHCLPSDVSLDLPEGQIRLEKSEDFSEKVEAHTACSCETSPCAA
    ## 5                  MRGIILALLLALAGCEKSQYEPFFSESKTYVYNYEGIILNGIPENGLARSGIKLNCKVELSGYAQRSYMLKIRNPEIKEFNGLWPKDPFTRATKLSQILAEQLTRPVKFEYRNGQVGNIFAPEDVSETVLNIQRGILNMLQVTIKKSQNVYDLQEDSISGVCHTKYTVQEDKKAERLIVIKSRDFENCESRVEKKIGVAYMETCTQCKKRAKHLQGSATFTYKLQNKDQGALIVEASAKQVYRYTPYYELYGAAVMEARQILTLAETKSSRNSESQVQLTKRGGLNYHFWQEVHQMPFHLLKTKNAESQVVETLQHLVQNNQELVHGESTSKFFQLVQLLRATDQENIQAIWKQFSDRPQHRHWLVNAITMAGSVDALKFIKQMIHNQDMTHQEAATAILMAMHYTRASQRVVEIAADLLTDSKVQNSPTLYKSAMLAYGSLVYKYCVDSQSCPEAALQPLHEFAAEAANKAHQEDITLALKAIGNAGQPASIKRIQKFLPGFSSGASRFPVRVQTDAVMALRNIAKKDPRKVQEVLLQLFMDRDVHPEVRMITCVALFETNPGIAIVTAVANAVLRESKDNLQLASFTYSHIKSLTKSTMPEFHNLASACRVALNILNPVFDQLSFRYSKAFHVDTFSYPMMGGLALNVFFINSPRTAFPVFALLKVRESFAGASTDVLEVALRAEGLQEVIRKQNIQFSEYSMRKKIVKILKALLAYKEAPANVPLVSASIKMLGQEVNFNQVTRDTIQNAMKLLNEPAERHTVVKRVLNKLLSGATGQWSHPMLLTEARHIIPTCVGLPLEIATHSAAVATTALNVDMKVTPSPSNDFSLSQLLESNIQLHTDMSPSVLIHVVETMGVNTPLIQSGLEFHGKLRIQTPVKFTAKLDMKDKNFKMEFPPCQREVEILTARHQMFAVTRNVGQLDSEKRMLVVPKQHGPSTLEQHFEKSGRTSAESASMMEVSSEMGIRKNAYSGGGAQHQRIPNHVQSYHFCLRLSKLNINACLQSKLNTAHGDMQSFLYEQIGEHETKVILKPVQNDAAIEKLQLEITTGPRASSKIIALTEVEVKEGEQLDESAVQKRLKIILGINEEEYRARNRSLVGKESKKAKIKKNKGAGETDVLEAGRNWKQSSSSSSSSPSSSSSSSSSSSSSSSSSSSSSSSSESSSSSSQRNRRNHGRKVDILSRVSQQAMNKETHGQTGRQKQRDSSSSSSSSQSSNKESSSSSSSSSSESSSSSESSRGKHQQNSKQQTRQNKQKQQKHHGSHPQDPSQMPPQMYKLRFQRPYKQQYKGMSSQSSESSSSSSSSSSSESARSQHHNKFLGDKKPPTVVITAQAVRNDNVKQGYQVILYTDDTTRPKMQAYVVDITGKSRWRACAEAYIRNPHEAQAELKWGQNCKQYKMEFIAQTGIFGNQPALRFKMNWPRIPSKLKSAGEEVAEFLPGAAYMMGFNTKNQRNPSKQAKIILSLSSPRTMNAVIQVPWATCYYQAVNLPFAVTSAYHLHSSDIQAPTWNIFADTYSSLVDRFSGECTMSRDQIKTFNEVKFNYTLPGSCFHVLAQDCSSELKFLVMAKKSPESSELRDINIKLADYDIDMYTSSEEIRLKVNGLEVPKEKLPYIASADPAVEIKMEEKAVILSAPDFGIEQLYYNGKLTKVRTAVWMRGSTCGICGQHDDETEREYQMPGGHAARDENRFTHSWILPEEPCSGTCKLQHTSVKLERTVHIDDQESKCYSVRPVLRCVKGCSPTSVTPVTVGFHCLPADASTDLVESQNIMGKTEDLVDTVDAHTACSCEALKCAA
    ## 6                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 MCDDDETTALVCDNGSGLVKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQSKRGILTLKYPIEHGIITNWDDMEKIWHHTFYNELRVAPEEHPTLLTEAPLNPKANREKMTQIMFETFNVPAMYVAIQAVLSLYASGRTTGIVLDSGDGVTHNVPIYEGYALPHAIQRLDLAGRDLTDYLMKILTERGYSFVTTAEREIVRDIKEKLAYVALDFENEMATAASSSSLEKSYELPDGQVITIGNERFRCPETLFQPSFIGMESAGIHETTYNSIMKCDIDIRKDLYANNVLSGGTTMYPGIADRMQKEITALAPSTMKIKIIAPPERKYSVWIGGSILASLSTFQQMWITKQEYDEAGPSIVHRKCF

Citation for fly absolute concentrations:

110 mg of protein in mature oocyte 80 mg at Stage 5 60 mg at hatching

<https://ijdb.ehu.eus/article/pdf/2518452>
<https://link.springer.com/content/pdf/10.1007/BF00848689.pdf>

-\> Total protein decreases by ~25% at hatching from window start

1.8 ug total protein at start 1.35 ug total protein at end

Trying: 0.8 ug total protein at start 0.6 ug total protein at end

``` r
fly_yolk_set <- c("P02843", "P02844", "P06607") #Yp1, Yp2, Yp3
names(fly_yolk_set) <- c("Yp1", "Yp2", "Yp3")

split_fly_yolk <- function(df, protein_amount) {

  df <- merge(df, fly_AA.df, by.x="Ind_Protein_ID", by.y="Protein_ID")

  #Calculate amount of protein in milligrams
  df["iBAQ_WeightFraction"] <- df$iBAQ_Value * df$kDa #mult by kda for weight
  #Calculate mass fraction
  df["iBAQ_WeightFraction"] <- df$iBAQ_WeightFraction / sum(df$iBAQ_WeightFraction)
  df["iBAQ_Weight_mg"] <- df$iBAQ_WeightFraction * (protein_amount/1000) #mg calc
  #Calculate mg/mL then divide by kDa -> mM
  df["iBAQ_Conc"] <- df$iBAQ_Weight_mg / 0.00001
  df["iBAQ_Conc"] <- df$iBAQ_Conc / df$kDa
  df <- df %>% arrange(desc(iBAQ_Conc))

  df[c("iBAQ_WeightFraction", "iBAQ_Weight_mg")] <- NULL
  df <- df %>% relocate(Sequence, .after = last_col())

  return(df) }

Dmel_T6_T0A_Abs <- split_fly_yolk(Dmel_T6_T0A_Abs, 0.8)
Dmel_T6_T0B_Abs <- split_fly_yolk(Dmel_T6_T0B_Abs, 0.8)
Dmel_T6_N16_Abs <- split_fly_yolk(Dmel_T6_N16_Abs, 0.6) #total decrease by 25%
Dmel_T7_T0A_Abs <- split_fly_yolk(Dmel_T7_T0A_Abs, 0.8)
Dmel_T7_T0B_Abs <- split_fly_yolk(Dmel_T7_T0B_Abs, 0.8)
Dmel_T7_N16_Abs <- split_fly_yolk(Dmel_T7_N16_Abs, 0.6) #total decrease by 25%
```

``` r
merge_datasets <- function(data1, data2, l1, l2) {

  protein_counts <- c(data1$Ind_Protein_ID, data2$Ind_Protein_ID)  
  protein_counts <- data.frame(table(protein_counts))
  
  single_ids <- protein_counts[protein_counts$Freq==1,]
  multiple_ids <- protein_counts[protein_counts$Freq>1,]
    
  no_reps <- rbind(data1[data1$Ind_Protein_ID %in% single_ids[[1]],],
                   data2[data2$Ind_Protein_ID %in% single_ids[[1]],])
  no_reps <- no_reps[c("Ind_Protein_ID", "iBAQ_Conc", "Sequence")]
  no_reps["Type"] <- "Point"
  
  reps <- merge(data1[data1$Ind_Protein_ID %in% multiple_ids[[1]],
                      c("Ind_Protein_ID", "iBAQ_Conc")],
                data2[data2$Ind_Protein_ID %in% multiple_ids[[1]],
                      c("Ind_Protein_ID", "iBAQ_Conc", "Sequence")],
                by="Ind_Protein_ID", suffixes=c(paste("_", c(l1, l2), sep="")),
                all=TRUE)

  reps["iBAQ_Conc"] <- (reps[[2]] + reps[[3]])/2
  reps["Type"] <- "Interval"

  out_df <- rbind(reps[c("Ind_Protein_ID", "iBAQ_Conc",
                         "Type", "Sequence")],
                  no_reps[c("Ind_Protein_ID", "iBAQ_Conc",
                            "Type", "Sequence")])

  return(out_df) }


XLA_RT0_Abs <- merge_datasets(XLA_T5_RT0_Abs, XLA_T6_RT0_Abs, "T5", "T6")
XLA_N12_Abs <- merge_datasets(XLA_T5_N12_Abs, XLA_T6_N12_Abs, "T5", "T6")

Dmel_T0_Abs <-
  merge_datasets(merge_datasets(Dmel_T6_T0A_Abs, Dmel_T6_T0B_Abs, "A", "B"),
                 merge_datasets(Dmel_T7_T0A_Abs, Dmel_T7_T0B_Abs, "A", "B"),
                 "T6", "T7")
Dmel_N16_Abs <- merge_datasets(Dmel_T6_N16_Abs, Dmel_T7_N16_Abs, "T6", "T7")
```

``` r
write.csv(XLA_RT0_Abs, "Files/AbsC/XLA_RT0_Abs-Concentrations.csv",
          row.names = FALSE)
write.csv(XLA_N12_Abs, "Files/AbsC/XLA_N12_Abs-Concentrations.csv",
          row.names = FALSE)

write.csv(Dmel_T0_Abs, "Files/AbsC/Dmel_T0_Abs-Concentrations.csv",
          row.names = FALSE)
write.csv(Dmel_N16_Abs, "Files/AbsC/Dmel_N16_Abs-Concentrations.csv",
          row.names = FALSE)
```

``` r
human_map <- read.csv("Files/Reference/Xen10_to_Human_Name_Assign.csv")
human_map["Protein_ID"] <- sapply(human_map$Protein_ID, function(x){
  split_id <- strsplit(x, "\\|")[[1]]
  reformat <- paste0(split_id[3], "|", split_id[2])

  return(reformat) })

human_map <- human_map[c("Protein_ID", "XLA_Gene", "Human_Gene", "Description")]
```

``` r
Dmel.genes <- read.csv("Files/Reference/Dmel_GeneTable_Descriptions.csv")
colnames(Dmel.genes)[1] <- "Protein_ID"
Dmel.genes <- Dmel.genes[c("Protein_ID", "Gene_Symbol", "Name", "Description")]
```

``` r
supp_df <- XLA_RT0_Abs
supp_df["iBAQ_Conc_uM"] <- supp_df$iBAQ_Conc * 1000

supp_df <- merge(human_map, supp_df[c("Ind_Protein_ID", "iBAQ_Conc_uM")],
                 by.x="Protein_ID", by.y="Ind_Protein_ID") %>%
  arrange(desc(iBAQ_Conc_uM))

dir.create("Files/Supp_Tables", showWarnings = FALSE, recursive = TRUE)
write.csv(supp_df, "Files/Supp_Tables/Frog_Start-AbsQ_Supplementary_Table.csv",
          row.names = FALSE)
```

``` r
supp_df <- Dmel_T0_Abs
supp_df["iBAQ_Conc_uM"] <- supp_df$iBAQ_Conc * 1000
supp_df <- merge(Dmel.genes, supp_df[c("Ind_Protein_ID", "iBAQ_Conc_uM")],
                 by.x="Protein_ID", by.y="Ind_Protein_ID") %>%
  arrange(desc(iBAQ_Conc_uM))

dir.create("Files/Supp_Tables", showWarnings = FALSE, recursive = TRUE)
write.csv(supp_df, "Files/Supp_Tables/Fly_Start-AbsQ_Supplementary_Table.csv",
          row.names = FALSE)
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
    ##  [1] ggplot2_4.0.3       Peptides_2.4.6      stringr_1.6.0      
    ##  [4] dplyr_1.2.0         tidyr_1.3.2         arrow_25.0.1       
    ##  [7] stringi_1.8.7       Biostrings_2.78.0   Seqinfo_1.0.0      
    ## [10] XVector_0.50.0      IRanges_2.44.0      S4Vectors_0.48.1   
    ## [13] BiocGenerics_0.56.0 generics_0.1.4     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] bit_4.6.0          gtable_0.3.6       compiler_4.5.3     crayon_1.5.3      
    ##  [5] Rcpp_1.1.1-1.1     tidyselect_1.2.1   assertthat_0.2.1   scales_1.4.0      
    ##  [9] yaml_2.3.12        fastmap_1.2.0      R6_2.6.1           knitr_1.51        
    ## [13] tibble_3.3.1       RColorBrewer_1.1-3 pillar_1.11.1      tzdb_0.5.0        
    ## [17] rlang_1.1.7        utf8_1.2.6         xfun_0.58          S7_0.2.1          
    ## [21] bit64_4.8.2        otel_0.2.0         cli_3.6.5          withr_3.0.3       
    ## [25] magrittr_2.0.4     grid_4.5.3         digest_0.6.39      rstudioapi_0.19.0 
    ## [29] lifecycle_1.0.5    vctrs_0.7.1        evaluate_1.0.5     glue_1.8.0        
    ## [33] farver_2.1.2       rmarkdown_2.31     purrr_1.2.2        tools_4.5.3       
    ## [37] pkgconfig_2.0.3    htmltools_0.5.9
