``

# Mito_genes_data_analysis
code for data analysis completed on mitochondrial associated genes found in possum. 


# To pull out relevant human gene stuff through BASH
cut -f1 possum_mito_gene_matches.txt | grep -v "NA" | grep -v "human_gene" > to_grep_human.txt
cut -f2 mito_genes_to_human_name.txt >> to_grep_human.txt

lines=`wc -l to_grep_human.txt | awk '{ print $1 }'`

for i in `seq 1 1 $lines`;
  do gene=`head -n $i to_grep_human.txt | tail -n 1`;
  grep "ID=gene-"$gene";" GCF_000001405.39_GRCh38.p13_genomic.gff >> human_mito_gene.gff;
done

# Download those to local computer

# loading libraries
library(tidyverse)
library(GGally)

# Read in our summary file
input <- read_tsv("possum_mito_gene_matches.txt")
human_input <- read_tsv("human_mito_gene.gff",col_names = FALSE)

# Derived from NCBI's genome summaries
possum_genome_info <- read_tsv("possum_genome_info.txt")
human_genome_info <- read_tsv("human_genome_info.txt")

# BASIC RESULTS
# input has 1739 lines
input

# How many total genes in mito_genes
input %>% select(gene_symbol) %>% distinct %>% dim()
# 1705 genes in total in mito_gene

# How many of these genes did we NOT find
input %>% filter(is.na(human_gene)) %>% select(gene_symbol) %>% distinct %>% dim()
# We found all but 66

# 1705-66 1639 of the genes have been found in possum

# Number of genes found more than once
input %>% filter(!(is.na(human_gene))) %>% count(human_gene) %>% count(n)
#n    nn
#<int> <int>
#  1     1  1617
#2     2    13
#3     3     6
#4     4     3

# SUMMARISING WHICH LOCI ARE PUTATIVELY MULTI-COPY
input %>% filter(!(is.na(human_gene))) %>% count(human_gene) %>% filter(n>1) %>% arrange(human_gene) %>% print(n=22)
# If pseudogenes present for these in ncbi gene, that may help with interpetation

# COPIES FOUND ON EACH CHROMOSOME IN POSSUMS
# Gives counts per Chromosome
input %>% count(seqid) %>% filter(!is.na(seqid)) %>% 
  mutate(RefSeq=ifelse(seqid %in% possum_genome_info$RefSeq,seqid,"Unknown"))
## A tibble: 17 × 3
#seqid              n RefSeq     
#<chr>          <int> <chr>      
#  1 NC_003039.1       13 NC_003039.1
#2 NC_050573.1      267 NC_050573.1
#3 NC_050574.1      217 NC_050574.1
#4 NC_050575.1      238 NC_050575.1
#5 NC_050576.1      215 NC_050576.1
#6 NC_050577.1      158 NC_050577.1
#7 NC_050578.1      130 NC_050578.1
#8 NC_050579.1      156 NC_050579.1
#9 NC_050580.1      150 NC_050580.1
#10 NC_050581.1       95 NC_050581.1
#11 NC_050582.1       17 NC_050582.1
#12 NW_023494377.1     7 Unknown    
#13 NW_023494378.1     4 Unknown    
#14 NW_023494416.1     2 Unknown    
#15 NW_023494454.1     1 Unknown    
#16 NW_023494531.1     1 Unknown    
#17 NW_023494553.1     2 Unknown  

# Total number of "unknowns" are:
7+4+2+1+1+2

chrom_rep <- input %>% count(seqid) %>% filter(!is.na(seqid)) %>% 
  mutate(RefSeq=ifelse(seqid %in% possum_genome_info$RefSeq,seqid,"Unknown")) %>% 
  filter(RefSeq!="Unknown")

chrom_rep <- rbind(chrom_rep,(c("Various",17,NA)))
chrom_rep <- full_join(chrom_rep,possum_genome_info,by="RefSeq")
chrom_rep <- chrom_rep %>% mutate(n=as.numeric(n))

# Correlation between scaffold size and number of mitochondrially associated genes
ggplot(chrom_rep) + geom_point(mapping=aes(x=`Size (Mb)`,y=n)) +
  geom_smooth(method='lm',mapping=aes(x=`Size (Mb)`,y=n))

# What other stuff is number of mitochondrially genes associated with
# Looks like with everything that is also correlated with length
ggpairs(chrom_rep %>% select(n,`Size (Mb)`,`GC%`,Protein,rRNA,tRNA,`Other RNA`,Gene,Pseudogene))

# I'm not really sure how to fit relationships using a known (0,0) intercept
# so cheated by "printing out" n and Size (Mb) and copying that into excel
chrom_rep %>% select(n,`Size (Mb)`)

# From excel
# Allowing intercept to vary
# y = 0.4809x + 4.8759R² = 0.9782

# Fixing intercept to 0
# y = 0.4934xR² = 0.9939

# Can then look at difference between expected and actual numbers of genes
chrom_rep %>% mutate(expected_genes=0.4934*`Size (Mb)`,obs_exp=n/expected_genes) %>% 
  select(RefSeq, Name, obs_exp)
# Not surprisingly, mtDNA has way more mitochondrial genes than expected for its length
# Interestingly, the X chromosome has about 40% fewer genes than expected for its length
# This might be due to the X chromosome not being great, so conservatively we could add 
# the X chromosome and unknowns together and see if that makes a difference
(chrom_rep$n[which(chrom_rep$Name=="X")]+chrom_rep$n[which(chrom_rep$Name=="Un")])/((chrom_rep$`Size (Mb)`[which(chrom_rep$Name=="X")]+chrom_rep$`Size (Mb)`[which(chrom_rep$Name=="Un")])*0.4934)
# Still has about 26% fewer genes than expected even if including unassigned contigs as X contgs

# creating a "new input" that has the chromosome names
input <- input %>% rowwise() %>% mutate(chrom=ifelse((seqid %in% chrom_rep$seqid),chrom_rep$Name[which(chrom_rep$seqid==seqid)],"Unplaced"))

# COPIES FOUND ON EACH CHROMOSOME IN HUMANS
# Gives counts per Chromosome
human_chrom_rep <- human_input %>%  mutate(X1=gsub("GCF_000001405.39_GRCh38.p13_genomic.gff:","",X1)) %>% 
  count(X1) %>% filter(!is.na(X1)) %>%
  mutate(RefSeq=ifelse(X1 %in% human_genome_info$RefSeq,X1,"Unknown"))

human_chrom_rep <- full_join(human_chrom_rep,human_genome_info,by="RefSeq")
human_chrom_rep <- human_chrom_rep %>% mutate(n=as.numeric(n))

# Correlation between scaffold size and number of mitochondrially associated genes
ggplot(human_chrom_rep) + geom_point(mapping=aes(x=`Size (Mb)`,y=n)) +
  geom_smooth(method='lm',mapping=aes(x=`Size (Mb)`,y=n))
# To me it looks like a much less tight relationship in humans than in possums!
# Although we haven't chased down matches in humans to the same extent
# so that is a limitation that should be stated

# What other stuff is number of mitochondrially genes associated with
# Looks like with everything that is also correlated with length
ggpairs(human_chrom_rep %>% select(n,`Size (Mb)`,`GC%`,Protein,rRNA,tRNA,`Other RNA`,Gene,Pseudogene))
# Looks like gene content is actually more important in humans
# correlation is much higher than size
# 0.903 (gene) vs 0.674 (length)

# Perhaps this is because there is more variation
# in length in possum chromosomes (so more explanatory power)
# this seems to be supported by comparing the lengths
# of the chromosomes
ggplot() + geom_point(mapping=aes(x="human",y=human_chrom_rep$`Size (Mb)`)) +
  geom_point(mapping=aes(x="possum",y=chrom_rep$`Size (Mb)`)) 

# However, let's check out patterns in humans based on # of protein per chromosome
human_chrom_rep %>% select(n,Protein) %>% print(n=27)
# y = 0.012x R² = 0.9731
# y = 0.0114x + 3.7232R² = 0.87
# Can then look at difference between expected and actual numbers of genes
human_chrom_rep %>% mutate(expected_genes=0.012*Protein,obs_exp=n/expected_genes) %>% 
  select(RefSeq, Name, obs_exp) %>% print(n=27)

# However, let's check out patterns in humans based on length per chromosome
human_chrom_rep %>% select(n,`Size (Mb)`) %>% print(n=27)
# y = 0.3671x + 19.032
# y = 0.4897x
# Can then look at difference between expected and actual numbers of genes
human_chrom_rep %>% mutate(expected_genes=0.4897*`Size (Mb)`,obs_exp=n/expected_genes) %>% 
  select(RefSeq, Name, obs_exp) %>% print(n=27)
# Lots more variation in humans

# creating a new inputs that has the chromosome names for downstream plotting
input <- input %>% rowwise() %>% mutate(chrom=ifelse((seqid %in% chrom_rep$seqid),chrom_rep$Name[which(chrom_rep$seqid==seqid)],"Unplaced"))
human_input <- human_input %>% rowwise() %>% mutate(chrom=ifelse((X1 %in% human_chrom_rep$RefSeq),human_chrom_rep$Name[which(human_chrom_rep$RefSeq==X1)],"Unplaced"))

# SLIDING WINDOW STUFF POSSUMS
# scaffold names
chrom_names <- unique(input$chrom)[!is.na(unique(input$chrom))]
# excluding unplaced
chrom_names <- chrom_names[chrom_names!="Unplaced"]
chrom_names <- c(sort(as.numeric(chrom_names)),"X","MT")

#initialising output
output <- NULL
# could tweak these - window width is the distance across which the number of mitochondrial genes is counted
window_width <- 10000000
# The windows "slide" by advance_ratio * window_width
advance_ratio <- 0.1
# intiialising the output index
output_index <- 1
# Stepping through each of the seqid_names
for (i in chrom_names) {
  # filtering just to this scaffold
  temp <- input %>% filter(chrom==i) %>% arrange(start)
  # starting at 0
  slide_start <- 0
  slide_end <- slide_start + window_width
  # calculating things for the first window
  temp_output <- c((slide_start+slide_end)/2,(dim(temp %>% filter(start <= slide_end & start >= slide_start))[1]))
  # moving the window along
  slide_start <- slide_start + round((advance_ratio*window_width),digits=0)
  slide_end <- slide_start + window_width
  # Keeping on marching along til we get within the window_width of teh end
  while (slide_start <= (temp$start[dim(temp)[1]]-window_width)) {
    # Add each of the windows and their counts to temp_output
    temp_output <- rbind(temp_output,c((slide_start+slide_end)/2,(dim(temp %>% filter(start <= slide_end & start >= slide_start))[1])))
    # Slide things along
    slide_start <- slide_start + round((advance_ratio*window_width),digits=0)
    slide_end <- slide_start + window_width
  }
  # Record the outputs from that scaffold 
  output[[output_index]] <- temp_output
  # count up the output index so the next scaffold can go in the next slot
  output_index <- output_index+1
}

# Converting output list into one big output
output_combined <- NULL

for (i in 1:length(output)) {
  output[[i]] <- as_tibble(output[[i]])
  if ((ncol(output[[i]])==1)) {
    output[[i]] <- c(output[[i]]$value[1],output[[i]]$value[2],i,chrom_names[i])
    names(output[[i]]) <- c("V1", "V2","yvalue","chrom")
  } else {
    output[[i]] <- output[[i]] %>% mutate(yvalue=i,chrom=chrom_names[i])
  }  
  output_combined <- rbind(output_combined,output[[i]])
}
  
output_combined  <- output_combined %>% mutate(V1=as.numeric(V1), V2=as.numeric(V2))

output_combined$chrom <- factor(output_combined$chrom, levels = chrom_names)

ggplot() + geom_point(data=output_combined,mapping=aes(x=V1,y=as.factor(chrom),colour=V2), shape="|", size=5) +
  scale_colour_viridis_c(name="Number of genes") +
  theme_bw() +
  ylab("Scaffold") +
  xlab("Location (bp)")
# save as a pdf 12 inches by 3 inches
  
# SLIDING WINDOW STUFF HUMANS
# scaffold names
chrom_names <- unique(human_input$chrom)[!is.na(unique(human_input$chrom))]
# excluding unplaced
chrom_names <- chrom_names[chrom_names!="Unplaced"]
chrom_names <- c(sort(as.numeric(chrom_names)),"X","MT")

#initialising output
output <- NULL
# could tweak these - window width is the distance across which the number of mitochondrial genes is counted
window_width <- 10000000
# The windows "slide" by advance_ratio * window_width
advance_ratio <- 0.1
# intiialising the output index
output_index <- 1
# Stepping through each of the seqid_names
for (i in chrom_names) {
  # filtering just to this scaffold
  temp <- human_input %>% filter(chrom==i) %>% arrange(X4)
  # starting at 0
  slide_start <- 0
  slide_end <- slide_start + window_width
  # calculating things for the first window
  temp_output <- c((slide_start+slide_end)/2,(dim(temp %>% filter(X4 <= slide_end & X4 >= slide_start))[1]))
  # moving the window along
  slide_start <- slide_start + round((advance_ratio*window_width),digits=0)
  slide_end <- slide_start + window_width
  # Keeping on marching along til we get within the window_width of teh end
  while (slide_start <= (temp$X4[dim(temp)[1]]-window_width)) {
    # Add each of the windows and their counts to temp_output
    temp_output <- rbind(temp_output,c((slide_start+slide_end)/2,(dim(temp %>% filter(X4 <= slide_end & X4 >= slide_start))[1])))
    # Slide things along
    slide_start <- slide_start + round((advance_ratio*window_width),digits=0)
    slide_end <- slide_start + window_width
  }
  # Record the outputs from that scaffold 
  output[[output_index]] <- temp_output
  # count up the output index so the next scaffold can go in the next slot
  output_index <- output_index+1
}

# Converting output list into one big output
output_combined <- NULL

for (i in 1:length(output)) {
  output[[i]] <- as_tibble(output[[i]])
  if ((ncol(output[[i]])==1)) {
    output[[i]] <- c(output[[i]]$value[1],output[[i]]$value[2],i,chrom_names[i])
    names(output[[i]]) <- c("V1", "V2","yvalue","chrom")
  } else {
    output[[i]] <- output[[i]] %>% mutate(yvalue=i,chrom=chrom_names[i])
  }  
  output_combined <- rbind(output_combined,output[[i]])
}

output_combined  <- output_combined %>% mutate(V1=as.numeric(V1), V2=as.numeric(V2))

output_combined$chrom <- factor(output_combined$chrom, levels = chrom_names)

ggplot() + geom_point(data=output_combined,mapping=aes(x=V1,y=as.factor(chrom),colour=V2), shape="|", size=5) +
  scale_colour_viridis_c(name="Number of genes") +
  theme_bw() +
  ylab("Scaffold") +
  xlab("Location (bp)")
# save as a pdf 6 inches by 6 inches

# IMPROVING THE POSSUM GENOME ANNOTATION
# How many of these are genes where the gene name in possum differs from mito_genes
input %>% mutate(possum_name=gsub(";.*","",gsub("ID=gene-","",ID))) %>% filter(gene_symbol!=possum_name)
# 301 in total where the name differs between possum and mito_genes

# Of those, 205 instances where it appears a gene has not yet been annotated in possums with the correct gene name
input %>% mutate(possum_name=gsub(";.*","",gsub("ID=gene-","",ID))) %>% filter(gene_symbol!=possum_name) %>% filter(grepl("^LOC",possum_name)) %>% print(n=205)

# These 205 instances are found across 171 distinct genes (i.e. some genes in mito_genes map to more than one possum locus)
input %>% mutate(possum_name=gsub(";.*","",gsub("ID=gene-","",ID))) %>% filter(gene_symbol!=possum_name) %>% filter(grepl("^LOC",possum_name)) %>% distinct(human_gene)

Loc_possum_genes <- input %>% mutate(possum_name=gsub(";.*","",gsub("ID=gene-","",ID))) %>% filter(gene_symbol!=possum_name) %>% filter(grepl("^LOC",possum_name)) %>% distinct(human_gene) %>% arrange(human_gene)

# These are the 22 genes without specified annotation in possum, where there are multiple possum matches (therefore that might be the reason they are not annotated)
dup_Loc_possum_genes <- Loc_possum_genes %>% rowwise() %>% mutate(copies=sum(input$human_gene==human_gene,na.rm=TRUE)) %>% filter(copies>1)

# Presumably the other 171-22 "LOC" genes represent ones where the annotation could be updated
input %>% mutate(possum_name=gsub(";.*","",gsub("ID=gene-","",ID))) %>% filter(gene_symbol!=possum_name) %>% filter(grepl("^LOC",possum_name)) %>% filter(!(human_gene %in% dup_Loc_possum_genes$human_gene)) %>% arrange(human_gene)

# LENGTH DISTRIBUTION
ggplot((input %>% mutate(length=abs(end-start)) %>% summarise(length)),aes(length)) +
  geom_histogram()

ggplot((input %>% mutate(length=abs(end-start)) %>% filter(length<100000)),aes(length)) +
  geom_histogram()

# Comparing location in humans to location in possums
# Just keeping locations where gene name is same
match_human <- human_input %>% mutate(matchup=gsub(";.*","",gsub("ID=gene-","",X9)))
match_human_poss <- inner_join(x=match_human,y=input,by=c("matchup" = "human_gene"))

match_human_poss$chrom.y <- factor(match_human_poss$chrom.y,levels = c(sort(as.numeric(unique(match_human_poss$chrom.y))),"X","MT","Unplaced"))
match_human_poss$chrom.x <- factor(match_human_poss$chrom.x,levels = c(sort(as.numeric(unique(match_human_poss$chrom.x))),"X","MT","Unplaced"))

# Plotting human location coloured by possums
ggplot(match_human_poss) + geom_jitter(mapping=aes(x=X4,y=chrom.x,fill=chrom.y),color="black",shape=21,width=0,height=0.2,size=2) +
  scale_fill_brewer(type = "qual",palette = "Paired") +
  theme_bw()
# Export 6 x 6
# Human X made up of multiple possum chromosomes

to_humans <- match_human_poss %>% group_by(chrom.x,chrom.y) %>% count() %>% arrange(chrom.y,desc(n))

# Scatter of how many genes on each chromosome assign to different human
# chromosomes
ggplot(to_humans) + geom_point(mapping=aes(x=chrom.y,y=n,colour=chrom.x))
# Possum X and possum mtDNA only have matches to human X and human Y
# suggests low number of mtDNA genes on possum X might be because
# in other species they are located not on sex chromosome
# potentially due to stronger selective pressure through inclusive fitness?
# e.g. marsupial sperm being way faster
