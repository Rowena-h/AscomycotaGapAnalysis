######Ascomycota genome gap analysis###########
##This script uses some webpage scraping, so results will differ as new data becomes available##

library(ape)
library(DescTools)
library(rvest)
library(ggplot2)
library(ggpubr)
library(ggstance)
library(ggtree)
library(grid)
library(gridExtra)
library(gtable)
library(multcompView)
library(RCurl)
library(scales)
library(stringr)


##MAIN FIGURE##

##Generate Ascomycota order-level phylogeny for side by side plot##

#Read in Ascomycota taxonomy (from Wijayawardene et al. 2018)
tax.df <- read.csv("Ascomycota outline 2017.csv", stringsAsFactors=TRUE)
#Add phylum column
tax.df$Phylum <- as.factor("Ascomycota")
#Remove incertae sedis orders
tax.df <- tax.df[!grepl("incertae sedis", tax.df$Order, ignore.case=TRUE),]
#Remove duplicate orders
tax.df2 <- tax.df[!duplicated(tax.df$Order),c(1:2,5)]
#Make order level tree with dataframe
asc.tree <- as.phylo(~Phylum/Class/Order, data=tax.df2)

##Cytometric genome size data (from http://www.zbi.ee/fungal-genomesize/, not from assemblies)##

#Read in genome size data
df <- read.csv("fungi_genome_sizes.csv")
#Subset genome size dataframe for just Ascomycota
asc.df <- subset(df, PHYLUM == "Ascomycota")
#Add estimates from Le Cam et al 2019
lecam <- read.csv("lecam_etal_2019.csv")
asc.df <- rbind(asc.df, lecam)
#Fix genus and species names for unknown (sp.) taxa
for (i in 1:length(asc.df$GENUS)) {
  if(length(grep("\\bsp\\.", asc.df$GENUS[i])) > 0) {
    asc.df$GENUS[i] <- sub(" sp\\.", "", asc.df$GENUS[i])
    asc.df$SPECIES[i] <- "sp."
  }
}
#Find genera that don't match current classification to check for species synonyms
sort(unique(paste(asc.df$GENUS, asc.df$SPECIES)[is.na(match(asc.df$GENUS, tax.df$Genus))]))
#Read in classification corrections
corrections <- read.csv("name_check.csv")
#For each species...
for (i in 1:length(asc.df$SPECIES)) {
  #If the species name matches the list of name to correct..
  if (length(grep(paste(asc.df$GENUS[i], asc.df$SPECIES[i]), corrections$Old)) > 0) {
    idx <- grep(paste(asc.df$GENUS[i], asc.df$SPECIES[i]), corrections$Old)
    #Replace the genus and species names with the current names
    asc.df$GENUS[i] <- unlist(str_split(corrections$Current[idx], " "))[1]
    asc.df$SPECIES[i] <- unlist(str_split(corrections$Current[idx], " "))[2]
  }
}
#Create vector of method to exclude (genome assembly, unreliable or unknown methods)
exclude <- c("CS", "genomic reconstruction", "", "CHEF gel electrophoresis", "DAPI-PC", "DAPI-IC", "PFGE", "CS and PFGE", "quantitative real-time PCR", "Re-association kinetics", "Integrated Physical/Genetic Map")
#Remove genome size data for these methods
no.CS.df <- subset(asc.df, !METHOD %in% exclude)
#Update taxonomy data
no.CS.df$ORDER <- tax.df$Order[match(no.CS.df$GENUS, tax.df$Genus)]
no.CS.df$CLASS <- tax.df$Class[match(no.CS.df$GENUS, tax.df$Genus)]
#Remove rows with no order classification
no.CS.df <- no.CS.df[!is.na(no.CS.df$ORDER),]
#Make vector with unique orders
no.CS.orders <- as.vector(unique(no.CS.df$ORDER))
#Extract genome size measurements
no.CS.res <- list()
for (i in 1:length(no.CS.orders)) {
  no.CS.res[[i]] <- no.CS.df$X1C.in.Mbp[no.CS.df$ORDER==no.CS.orders[i]]
}
#Add order names
no.CS.res <- setNames(no.CS.res,no.CS.orders)
#Remove genome sizes of 0 or NA
no.CS.res <- lapply(no.CS.res, function(x) x[x!=0])
no.CS.res <- lapply(no.CS.res, function(x) x[!is.na(x)])
#Filter for orders with more than 5 genome size measurements
#no.CS.res <- no.CS.res[lapply(no.CS.res, length) >= 5]
#Convert results list to vector
no.CS.taxa <- setNames(unlist(no.CS.res, recursive=TRUE, use.names=FALSE), rep(names(no.CS.res), lengths(no.CS.res)))
#Make dataframe of genome sizes for orders
no.CS.df2 <- data.frame(order=names(no.CS.taxa), size=no.CS.taxa)
#Identify outliers
#outliers <- boxplot(no.CS.df2$size, plot=FALSE)$out
#no.CS.df2.out <- no.CS.df2
#Make dataframe excluding outliers
#no.CS.df2 <- no.CS.df2[-which(no.CS.df2$size %in% outliers),]

#Mean for cytometric methods
mean(no.CS.df2$size)
#Number of species included
length(unique(paste0(no.CS.df$GENUS, no.CS.df$SPECIES)))
#Number/proportion of orders without cytometric genome size data
paste0(length(tax.df3$order[tax.df3$noCS == "N"]), "/",length(tax.df3$order), ", ", (length(tax.df3$order[tax.df3$noCS == "N"])) / length(tax.df3$order) * 100,"%")
#Number/proportion of classes without cytometric genome size data
paste0(length(class.df[-20,]$class[class.df[-20,]$noCS == "N"]), "/",length(class.df[-20,]$class), ", ", (length(class.df[-20,]$class[class.df[-20,]$noCS == "N"])) / length(class.df[-20,]$class) * 100,"%")


##Assembly-based genome size data##

#Download and read in ncbi genome data (< 3 MB file)
download.file("ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/eukaryotes.txt", destfile=paste0(Sys.Date(), "_eukaryotes.txt"))
ncbi <- read.csv(paste0(Sys.Date(), "_eukaryotes.txt"), header=TRUE, sep="\t")
#Filter for Ascomycota
ncbi <- ncbi[ncbi$SubGroup == "Ascomycetes",]
#Filter for only 'whole' genomes
ncbi <- ncbi[ncbi$Status != "Chromosome",]
#Remove duplicate biosamples
ncbi <- ncbi[!duplicated(ncbi$BioSample.Accession[ncbi$BioSample.Accession != "-"]),]
#Remove genomes which are too small to be credible
ncbi <- ncbi[ncbi$Size..Mb. > 1,]

#Extract assembly method information
#Download and read in file with ftp links to assemblies (< 300 MB file)
download.file("ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt", destfile=paste0(Sys.Date(), "_assembly_summary_genbank.txt"))
assembly.sum <- read.csv(paste0(Sys.Date(), "_assembly_summary_genbank.txt"), skip=1, header=TRUE, sep="\t", quote="")
#Match to ncbi data
assembly.sum <- assembly.sum[match(ncbi$Assembly.Accession, assembly.sum$X..assembly_accession),]

#Scrape Mycocosm website for data
mycocosm <- read_html("https://mycocosm.jgi.doe.gov/ascomycota/ascomycota.info.html")
#Make dataframe
myc <- as.data.frame(mycocosm %>%
                       html_nodes(xpath="/html/body/div[4]/div/table") %>%
                       html_table())
#Add portal code
myc$Portal <- sub("/", "", mycocosm %>%
  html_nodes("td:nth-child(2) a") %>%
  html_attr("href"))
#Remove ncbi assemblies duplicated in mycocosm
assembly.sum <- assembly.sum[is.na(match(assembly.sum$asm_name, myc$Portal)),]
ncbi <- ncbi[match(assembly.sum$X..assembly_accession, ncbi$Assembly.Accession),]

#Remove square brackets from ncbi names
ncbi$X.Organism.Name <- sub('\\[', "", ncbi$X.Organism.Name)
ncbi$X.Organism.Name <- sub('\\]', "", ncbi$X.Organism.Name)

#Add column for genus
ncbi$Genus <- gsub(" .*", "", ncbi$X.Organism.Name)
myc$Genus <- gsub(" .*", "", myc$Name)
#Correct Mycocosm genome size units
myc$Assembly.Length <- as.numeric(gsub(",","", myc$Assembly.Length)) / 1000000
#Add taxonomy data to ncbi and mycocosm dataframes
for (i in c("Family", "Order", "Class")) {
  myc[[i]] <- tax.df[match(myc$Genus, tax.df$Genus), i]
  ncbi[[i]] <- tax.df[match(ncbi$Genus, tax.df$Genus), i]
}

#Combine Mycocosm and NCBI to make dataframe of genome sizes
CS.df <- data.frame(order=c(as.character(ncbi$Order), as.character(myc$Order)), size=c(as.numeric(ncbi$Size..Mb.), as.numeric(myc$Assembly.Length)))
#Remove incertae sedis and NA orders
CS.df2 <- CS.df[!is.na(CS.df$order),]
CS.df2 <- CS.df2[!grepl("incertae sedis", CS.df2$order),]
#Make vector of orders
CS.orders <- as.vector(unique(CS.df2$order))
#Identify outliers
#outliers <- boxplot(CS.df2$size, plot=FALSE)$out
#CS.df2.out <- CS.df2
#Make dataframe excluding outliers
#CS.df2 <- CS.df2[-which(CS.df2$size %in% outliers),]

#Mean for assembly methods
mean(CS.df2$size)
#Number of strains included
length(unique(c(ncbi$X.Organism.Name, myc$Name)))
#Number/proportion of orders without assembly-based genome size data
paste0(length(tax.df3$order[tax.df3$CS == "N"]), "/",length(tax.df3$order), ", ", (length(tax.df3$order[tax.df3$CS == "N"])) / length(tax.df3$order) * 100,"%")
#Number/proportion of classes without assembly-based genome size data
paste0(length(class.df[-20,]$class[class.df[-20,]$CS == "N"]), "/",length(class.df[-20,]$class), ", ", (length(class.df[-20,]$class[class.df[-20,]$CS == "N"])) / length(class.df[-20,]$class) * 100,"%")


##Test significant difference of mean for each order##

#Filter for orders with 3 or more genome size measurements (for normality testing)
no.CS.filt <- no.CS.df2
CS.filt <- CS.df2

for (i in 1:length(CS.orders)) {
  if (length(CS.filt$order[grepl(CS.orders[i], CS.filt$order)]) < 3) {
    CS.filt <- CS.filt[CS.filt$order != CS.orders[i],]
  }
}

for (i in 1:length(no.CS.orders)) {
  if (length(no.CS.filt$order[grepl(no.CS.orders[i], no.CS.filt$order)]) < 3) {
    no.CS.filt <- no.CS.filt[no.CS.filt$order != no.CS.orders[i],]
  }
}

#Create vector of orders with both CS and no CS data
both.filt <- intersect(no.CS.filt$order, CS.filt$order)
#Create dataframe for results
sig <- data.frame(order=both.filt, pvalue=NA, noCS=NA, CS=NA)

#Significance testing
for (i in 1:length(both.filt)) {
  print(both.filt[i])
  x <- no.CS.df2$size[no.CS.df2$order == both.filt[i]]
  y <- CS.df2$size[CS.df2$order == both.filt[i]]
  
  #Add columns with means
  sig$noCS[i] <- mean(x)
  sig$CS[i] <- mean(y)
  
  #Test for normality of data
  shapiro.x <- shapiro.test(x)
  shapiro.y <- shapiro.test(y)
  
  #Reject normality if p-value < 0.05 and do Wilcoxon, otherwise do t-test
  if (shapiro.x$p.value < 0.05 || shapiro.y$p.value < 0.05) {
    print("Reject normality, doing Wilcoxon test")
    wilcox <- wilcox.test(x=x, y=y, mu=0)
    sig$pvalue[sig$order == both.filt[i]] <- wilcox$p.value
  } else {
    print("Normal data, doing t-test")
    ttest <- t.test(x=x, y=y, mu=0)
    sig$pvalue[sig$order == both.filt[i]] <- ttest$p.value
  }
}

#Which orders have significantly bigger mean genome estimates from assembly
sig[sig$CS > sig$noCS & sig$pvalue < 0.05,]


##Create plotting dataframes

#Create dataframe for branch labels
branches <- data.frame(node=asc.tree$edge[,2], 
                       edge_num=1:nrow(asc.tree$edge), 
                       branch_label=rep(NA,nrow(asc.tree$edge)))
#Read in tree node labels for classes
classnodes <- read.csv("classnodes.csv", header=TRUE)
#Add class
for (i in 1:length(classnodes$class)) {
  branches$branch_label[branches$edge_num == classnodes$node[i]] <- as.character(classnodes$class[i])
}
#Order by node
branches <- branches[order(branches$node),]
#Add colour
branches$colour <- classnodes$colour[match(branches$branch_label, classnodes$class)]

#Make dataframe matching taxonomy data with the Ascomycota tree tip labels for colours
tax.df3 <- data.frame(order=asc.tree$tip.label, 
                      class=tax.df2$Class[match(asc.tree$tip.label, unique(tax.df2$Order))], 
                      CS="N", 
                      noCS="N")
#Add to CS column for whether the order has CS genome size data
tax.df3$CS[!is.na(match(tax.df3$order, CS.df$order))] <- "Y"
#Add to noCS column for whether the order has noCS genome size data
tax.df3$noCS[!is.na(match(tax.df3$order, no.CS.df2$order))] <- "Y"

#Make dataframe of class data to colour class branches
class.df <- data.frame(class=unique(tax.df3$class),
                       CS="N",
                       noCS="N")
#Add genome data to class dataframe
class.df$CS[match(unique(tax.df3$class[which(tax.df3$CS == "Y")]), class.df$class)] <- "Y"
class.df$noCS[match(unique(tax.df3$class[which(tax.df3$noCS == "Y")]), class.df$class)] <- "Y"
#Add genome data to branches dataframe
branches$CS.fontcol <- class.df$CS[match(branches$branch_label, class.df$class)]
branches$noCS.fontface <- class.df$noCS[match(branches$branch_label, class.df$class)]

#Dataframes for vertical grid lines
box.lines <- data.frame(x=seq(0, 230, 10), .panel="Genome size (Mbp/1C)", stringsAsFactors=TRUE)
noCSbox.lines.out <- data.frame(x=seq(0, RoundTo(max(c(no.CS.df2$size, CS.df2$size)), 100, FUN=ceiling), 100), .panel="Genome size (Mbp/1C)", stringsAsFactors=TRUE)
#Dataframe for mean lines for cytometric and assembly-based methods
no.CS.mean <- data.frame(x=mean(no.CS.df2$size), .panel="Genome size (Mbp/1C)", stringsAsFactors=TRUE)
CS.mean <- data.frame(x=mean(CS.df2$size), .panel="Genome size (Mbp/1C)", stringsAsFactors=TRUE)

#Add class data to genome size dataframe for plot colours
no.CS.df3 <- no.CS.df2
no.CS.df3$class <- tax.df2$Class[match(no.CS.df3$order,tax.df2$Order)]

#Add class data to genome size dataframe (INCLUDING OUTLIERS)
#no.CS.df3.out <- no.CS.df2.out
#no.CS.df3.out["class"] <- tax.df2$Class[match(no.CS.df3.out$order,tax.df2$Order)]
#Create a list of dataframes with means for each class
#means.out <- list()
#for (i in sort(truenoCS$class)) {
#  means.out[[i]] <- data.frame(x=(mean(no.CS.df3.out[no.CS.df3.out$class==i,]$size)), .panel='Genome size (Mbp/1C)')
#}

#Create dataframe for adding sample size to plot
counts.df <- tax.df3
#Remove orders without genome size data
counts.df <- counts.df[which(counts.df$noCS == "Y" | counts.df$CS == "Y"),]
#Add column with count of sample size for no CS/CS
counts.df$noCScount <- table(unlist(no.CS.df2$order))[match(counts.df$order, names(table(unlist(no.CS.df2$order))))]
counts.df$CScount <- table(unlist(CS.df2$order))[match(counts.df$order,names(table(unlist(CS.df2$order))))]
#Put brackets around count numbers
counts.df$noCScount <- paste0("(", counts.df$noCScount,")")
counts.df$CScount <- paste0("(", counts.df$CScount,")")
#Replace NA with 0
counts.df$noCScount[counts.df$noCScount == "(NA)"] <- "(0)"
counts.df$CScount[counts.df$CScount == "(NA)"] <- "(0)"

#Add column with asterisks for orders with significant difference of means
#0.05 > x > 0.01
for (i in 1:length(sig$order[sig$pvalue < 0.05 & sig$pvalue > 0.01])) {
  counts.df$sig[counts.df$order == sig$order[sig$pvalue < 0.05 & sig$pvalue > 0.01][i]] <- "*"
}
#0.01 > x > 0.001
for (i in 1:length(sig$order[sig$pvalue < 0.01 & sig$pvalue > 0.001])) {
  counts.df$sig[counts.df$order == sig$order[sig$pvalue < 0.01 & sig$pvalue > 0.001][i]] <- "**"
}
#0.001 > x
for (i in 1:length(sig$order[sig$pvalue < 0.001])) {
  counts.df$sig[counts.df$order == sig$order[sig$pvalue < 0.001][i]] <- "***"
}

#Add column with the upper limit per order (to position labels)
counts.df$noCSmax <- NA
for (i in unique(no.CS.df2$order)) {
  counts.df[counts.df$order==i,]$noCSmax <- max(no.CS.df2[no.CS.df2$order==i & no.CS.df2$size < 220,]$size)
}
counts.df$noCSmax.out <- NA
for (i in unique(no.CS.df2$order)) {
  counts.df[counts.df$order==i,]$noCSmax.out <- max(no.CS.df2[no.CS.df2$order==i,]$size)
}
counts.df$CSmax.out <- NA
for (i in unique(CS.df2$order)) {
  counts.df[counts.df$order==i,]$CSmax.out <- max(CS.df2[CS.df2$order==i,]$size)
}

##Number of genome assemblies##

#Create data frame for number of genome assemblies per order
num.df <- data.frame(order=tax.df3$order)

#Add column with number of genomes in mycocosm and ncbi
for (i in 1:length(num.df$order)) {
  num.df$num[i] <- length(grepl(num.df$order[i], myc$Order)[grepl(num.df$order[i], myc$Order) == TRUE]) + length(grepl(num.df$order[i], ncbi$Order)[grepl(num.df$order[i], ncbi$Order) == TRUE])
}

#Delete rows with no genomes
num.df <- num.df[num.df$num > 0,]

#Dataframes for vertical grid lines
bar.lines <- data.frame(x=seq(0, RoundTo(max(num.df$num), 200, FUN=ceiling), 50), .panel="Number of genome assemblies", stringsAsFactors=TRUE)

#Make dataframe of tree tip order for row shading
tips.df <- subset(fortify(asc.tree), isTip)
tips.df <- with(tips.df, label[order(y, decreasing=T)])
tips.df <- data.frame(tip=tips.df, min=seq(length(tips.df))-0.5, max=seq(length(tips.df))+0.5, col=NA)
tips.df$col <- rep_len(c(0,1),length(tips.df$tip))

#Function to remove unwanted elements from plots (https://stackoverflow.com/questions/36779537/ggplot2-facet-wrap-y-axis-scale-on-the-first-row-only)
gtable_filter_remove <- function (x, name, trim = TRUE){
  matches <- !(x$layout$name %in% name)
  x$layout <- x$layout[matches, , drop = FALSE]
  x$grobs <- x$grobs[matches]
  if (trim) 
    x <- gtable_trim(x)
  x
}


#Plot tree against number of genomes and CS/noCS genome sizes 

#Dummy plot for correct scale of boxplot (not showing extreme outliers)
gg.dummy <- ggtree(asc.tree) +
  geom_rect(data=tips.df,
            aes(x=NULL,
                y=NULL,
                ymin=min,
                ymax=max,
                fill=as.factor(col),
                xmin=-Inf, 
                xmax=+Inf),
            alpha=0.15) +
  geom_tree()
gg.dummy1 <- gg.dummy %<+% tax.df3 +
  geom_tiplab(size=2.8,
              aes(subset=noCS != "Y", colour=CS)) +
  geom_tiplab(size=2.8,
              aes(subset=noCS == "Y", colour=CS), fontface="bold.italic")
gg.dummy2 <- gg.dummy1 %<+% branches + 
  geom_label2(aes(x=branch, label=branch_label, colour=CS.fontcol, subset=noCS.fontface == "N"),
             fill="white",
             size=2.8,
             label.padding=unit(0.15, "lines"),
             label.size=0) +
  geom_label2(aes(x=branch, label=branch_label, colour=CS.fontcol, subset=noCS.fontface == "Y"),
             fontface="bold.italic",
             fill="white",
             size=2.8,
             label.padding=unit(0.15, "lines"),
             label.size=0) +
  scale_colour_manual(values=c("darkgrey", "black"))
gg.dummy3 <- gg.dummy2 +
  geom_vline(data=box.lines,
             color="grey94",
             aes(xintercept=x)) +
  geom_vline(data=no.CS.mean,
             color="black",
             linetype="dashed",
             aes(xintercept=x)) +
  geom_vline(data=CS.mean,
             color="darkgrey",
             linetype="dashed",
             aes(xintercept=x))
gg.dummy4 <- facet_plot(gg.dummy3,
                  panel="Number of genome assemblies",
                  data=num.df,
                  geom=geom_barh,
                  aes(x=num, fill=class),
                  stat='identity')
gg.dummy5 <- facet_plot(gg.dummy4,
                  panel="Genome size (Mbp/1C)",
                  data=no.CS.df2[no.CS.df2$size < 220,],
                  geom=geom_boxploth,
                  outlier.size=1,
                  aes(x=size, group=label, fill=class))
gg.dummy5 <- facet_plot(gg.dummy5,
                  panel="Genome size (Mbp/1C)",
                  data=CS.df2,
                  geom=geom_boxploth,
                  outlier.size=1,
                  colour="darkgrey",
                  linetype="dotted",
                  alpha=0.3,
                  aes(x=size, group=label, fill=class))
gg.dummy5 <- facet_plot(gg.dummy5,
                  panel="Genome size (Mbp/1C)",
                  data=counts.df,
                  geom=geom_text,
                  size=2,
                  nudge_x=10,
                  aes(x=noCSmax, label=noCScount))
gg.dummy5 <- facet_plot(gg.dummy5,
                  panel="Genome size (Mbp/1C)",
                  data=counts.df,
                  geom=geom_text,
                  size=2,
                  nudge_x=10,
                  nudge_y=0.5,
                  aes(x=noCSmax, label=sig)) +
  scale_x_continuous(breaks=seq(0,230,by=20), 
                     limits=c(0,230),
                     expand=c(0,0),
                     position="top",
                     sec.axis = dup_axis()) +
  scale_y_continuous(expand=expansion(mult=0.005)) +
  scale_fill_manual(values=c("snow3","white",as.vector(classnodes$colour[classnodes$colour != ""]))) +
  theme_classic() +
  theme(panel.border=element_blank(),
        legend.position="none",
        strip.placement="outside",
        axis.text.x.top=element_text(angle=45, vjust=0.3),
        axis.text.x.bottom=element_text(angle=45, hjust=1),
        panel.spacing=unit(1, "lines"),
        strip.background=element_blank())

#Make plot into table
dummy <- ggplotGrob(gg.dummy5)
#Change widths of bar and box panels
dummy$widths[7] <- 0.70*dummy$widths[7]
dummy$widths[9] <- 0.70*dummy$widths[9]


#Main plot

gg.main <- ggtree(asc.tree) +
  geom_rect(data=tips.df,
            aes(x=NULL,
                y=NULL,
                ymin=min,
                ymax=max,
                fill=as.factor(col),
                xmin=-Inf, 
                xmax=+Inf),
            alpha=0.15)

for (i in 1:length(na.omit(branches[branches$colour != "",])$node)) {
  gg.main <- gg.main +
    geom_hilight(node=na.omit(branches[branches$colour != "",])$node[i],
                 extend=1.2,
                 alpha=0.1, fill=na.omit(branches[branches$colour != "",])$colour[i])
}

gg.main <- gg.main + geom_tree()
gg.main1 <- gg.main %<+% tax.df3 +
  geom_tiplab(size=2.8,
              aes(subset=noCS != "Y", colour=CS)) +
  geom_tiplab(size=2.8,
              aes(subset=noCS == "Y", colour=CS), fontface="bold.italic")
gg.main2 <- gg.main1 %<+% branches +
  geom_label2(aes(x=branch, label=branch_label, colour=CS.fontcol, subset=noCS.fontface == "N"),
              fill="white",
              size=2.8,
              label.padding=unit(0.15, "lines"),
              label.size=0) +
  geom_label2(aes(x=branch, label=branch_label, colour=CS.fontcol, subset=noCS.fontface == "Y"),
              fontface="bold.italic",
              fill="white",
              size=2.8,
              label.padding=unit(0.15, "lines"),
              label.size=0) +
  scale_colour_manual(values=c("darkgrey", "black"))
gg.main3 <- gg.main2 +
  geom_vline(data=box.lines,
             color="grey94",
             aes(xintercept=x)) +
  geom_vline(data=no.CS.mean,
             color="black",
             linetype="dashed",
             aes(xintercept=x)) +
  geom_vline(data=CS.mean,
             color="darkgrey",
             linetype="dashed",
             aes(xintercept=x)) +
  geom_vline(data=bar.lines,
             color="grey94",
             aes(xintercept=x))
gg.main4 <- facet_plot(gg.main3 + xlim_tree(3),
                  panel="Number of genome assemblies",
                  data=num.df,
                  geom=geom_barh,
                  aes(x=num, fill=class),
                  stat='identity')
gg.main4 <- facet_plot(gg.main4,
                  panel="Number of genome assemblies",
                  data=num.df,
                  geom=geom_text,
                  size=2,
                  nudge_x=50,
                  aes(x=num, label=paste0("(",num,")")))
gg.main5 <- facet_plot(gg.main4,
                  panel="Genome size (Mbp/1C)",
                  data=no.CS.df2,
                  geom=geom_boxploth,
                  #outlier.shape=NA,
                  outlier.size=1,
                  aes(x=size, group=label, fill=class))
gg.main5 <- facet_plot(gg.main5,
                  panel="Genome size (Mbp/1C)",
                  data=CS.df2,
                  geom=geom_boxploth,
                  #outlier.shape=NA,
                  outlier.size=1,
                  colour="darkgrey",
                  linetype="dotted",
                  alpha=0.3,
                  aes(x=size, group=label, fill=class))
gg.main5 <- facet_plot(gg.main5,
                  panel="Genome size (Mbp/1C)",
                  data=counts.df,
                  geom=geom_text,
                  size=2,
                  nudge_x=5,
                  aes(x=noCSmax, label=noCScount)) +
  scale_x_continuous(breaks=pretty_breaks(10),
                     expand=c(0,0),
                     position="top",
                     sec.axis=dup_axis()) +
  scale_y_continuous(expand=expansion(mult=0.005)) +
  scale_fill_manual(values=c("snow3","white",as.vector(classnodes$colour[classnodes$colour != ""]))) +
  theme_classic() +
  theme(panel.border=element_blank(),
        legend.position="none",
        strip.placement="outside",
        axis.text.x.top=element_text(angle=45, vjust=0.3),
        axis.text.x.bottom=element_text(angle=45, hjust=1),
        panel.spacing=unit(1, "lines"),
        strip.background=element_blank())

#Make plot into table
gg.main.tab <- ggplotGrob(gg.main5)
#Change widths of bar and box panels
gg.main.tab$widths[7] <- 0.70*gg.main.tab$widths[7]
gg.main.tab$widths[9] <- 0.70*gg.main.tab$widths[9]
#Make vector of elements in table
elements <- gg.main.tab$layout$name
#Remove unwanted axes on tree panel and entire third panel with bad scale
gg.main6 <- gtable_filter_remove(gg.main.tab, name=elements[c(5,8,11,13,4,7,10)], trim=FALSE)

#Replace third panel with dummy boxplot and axes
element.replace <- c("panel-1-3", "axis-t-3", "axis-b-3")
for (i in element.replace) {
  pos <- c(subset(dummy$layout, name == i, se=t:r))
  gg.main6 <- gtable_add_grob(gg.main6, dummy$grobs[[which(dummy$layout$name == i)]], pos$t, pos$l, pos$b, pos$r, name=i)
}

tiff(file=paste0("allgenomesplot_", Sys.Date(), ".tiff"), height=15, width=10, units="in", res=300)

grid.draw(gg.main6)

dev.off()


#Supplementary figure including extreme outliers

gg.supp <- ggtree(asc.tree) +
  geom_rect(data=tips.df,
            aes(x=NULL,
                y=NULL,
                ymin=min,
                ymax=max,
                fill=as.factor(col),
                xmin=-Inf, 
                xmax=+Inf),
            alpha=0.15)

for (i in 1:length(na.omit(branches[branches$colour != "",])$node)) {
  gg.supp <- gg.supp +
    geom_hilight(node=na.omit(branches[branches$colour != "",])$node[i],
                 extend=1.2,
                 alpha=0.1, fill=na.omit(branches[branches$colour != "",])$colour[i])
}

gg.supp <- gg.supp + 
  geom_tree()

gg.supp1 <- gg.supp %<+% tax.df3 +
  geom_tiplab(size=2.8,
              aes(subset=noCS != "Y", colour=CS)) +
  geom_tiplab(size=2.8,
              aes(subset=noCS == "Y", colour=CS), fontface="bold.italic") +
  scale_colour_manual(values=c("darkgrey", "black"))

gg.supp2 <- gg.supp1 %<+% branches + 
  geom_label2(aes(x=branch, label=branch_label, colour=CS.fontcol, subset=noCS.fontface == "N"),
              fill="white",
              size=2.8,
              label.padding=unit(0.15, "lines"),
              label.size=0) +
  geom_label2(aes(x=branch, label=branch_label, colour=CS.fontcol, subset=noCS.fontface == "Y"),
              fontface="bold.italic",
              fill="white",
              size=2.8,
              label.padding=unit(0.15, "lines"),
              label.size=0)

gg.supp3 <- gg.supp2 +
  geom_vline(data=noCSbox.lines.out,
             color="grey94",
             aes(xintercept=x))
gg.supp3 <- facet_plot(gg.supp3 + xlim_tree(3),
                       panel="Genome size (Mbp/1C)",
                       data=no.CS.df2,
                       geom=geom_boxploth,
                       outlier.size=1,
                       aes(x=size, group=label, fill=class))
gg.supp3 <- facet_plot(gg.supp3,
                       panel="Genome size (Mbp/1C)",
                       data=CS.df2,
                       geom=geom_boxploth,
                       outlier.size=1,
                       colour="darkgrey",
                       linetype="dotted",
                       alpha=0.3,
                       aes(x=size, group=label, fill=class))

gg.supp4 <- facet_plot(gg.supp3,
                       panel="Sample size",
                       data=counts.df[counts.df$CS == "Y",],
                       geom=geom_text,
                       size=2,
                       hjust=1,
                       aes(x=500, label=CScount))
gg.supp4 <- facet_plot(gg.supp4,
                       panel="Sample size",
                       data=counts.df[counts.df$noCS == "Y",],
                       geom=geom_text,
                       fontface="bold.italic",
                       size=2,
                       hjust=0,
                       aes(x=600, label=noCScount))  +
  xlim_expand(c(0, 1000), panel="Sample size") +
  scale_x_continuous(breaks=pretty_breaks(10),
                     expand=c(0, 0),
                     position="top",
                     sec.axis=dup_axis()) +
  scale_y_continuous(expand=expansion(mult=0.005)) +
  scale_fill_manual(values=c("snow3","white",as.vector(classnodes$colour[classnodes$colour != ""]))) +
  coord_cartesian(clip="off") +
  theme_classic() +
  theme(panel.border=element_blank(),
        legend.position="none",
        panel.spacing=unit(0, "lines"),
        plot.margin=unit(c(0,0,5,0), "mm"),
        strip.placement="outside",
        strip.background=element_blank())

#Remove unwanted axes on tree panel
gg.supp.tab <- ggplotGrob(gg.supp4)
elements <- gg.supp.tab$layout$name
gg.supp5 <- gtable_filter_remove(gg.supp.tab, name=elements[c(5, 7, 8, 10, 11, 13, 15)], trim=FALSE)
#Change widths of bar and box panels
gg.supp5$widths[9] <- 0.1*gg.supp5$widths[9]

tiff(file=paste0("genomesizeplot_outliers_", Sys.Date(), ".tiff"), height=15, width=10, units="in", res=300)

grid.draw(gg.supp5)

dev.off()



## SPECIES-LEVEL COMPARISON ##
## FIGURE 2a ##

##Identify case study species with both cytometric and assembly-based measurements
#Add name field to genome size dataframe 
no.CS.df$name <- paste0(no.CS.df$GENUS," ",no.CS.df$SPECIES)

#Create an empty list for species with both cytometric and assembly-based 
species.comp <- list()

#For each species with cytometric data...
for (i in 1:length(unique(no.CS.df$name))) {
  #If a genome assembly is also in NCBI...
  if (length(grep(unique(no.CS.df$name)[i], ncbi$X.Organism.Name)) > 0) {
    #Add the number of genome assemblies to the list
    species.comp[[unique(no.CS.df$name)[i]]] <- length(grep(unique(no.CS.df$name)[i], ncbi$X.Organism.Name))
  }
}

#Remove unknown species (i.e. sp.)
species.comp <- species.comp[-grep("\\bsp\\b", names(species.comp))]
#Make vector of abbreviations
abb <- abbreviate(names(species.comp), minlength=3)

#For each species-level case study...
for (i in 1:length(species.comp)) {
  
  #Get NCBI assembly links
  assemblies <- assembly.sum[agrep(names(species.comp)[i], assembly.sum$organism_name),]
  ftp.links <- paste0(assemblies$ftp_path, "/", assemblies$X..assembly_accession,"_", assemblies$asm_name, "_assembly_report.txt")
  
  #Create a results dataframe of genome sizes and assembly methods
  ncbi.methods.df <- data.frame(species=rep(names(species.comp)[i], length(ftp.links)), size=ncbi$Size..Mb.[agrep(names(species.comp)[i], ncbi$X.Organism.Name)], method=NA, type="CS") 
  
    #For each assembly...
    for (j in 1:length(ftp.links)) {
      
      #Try to download the assembly report
      report <- NULL
      report <- tryCatch(unlist(strsplit(getURL(ftp.links[j]), "\\r*\\n")), error = function(e) {e$message})
      
      #If the assembly method is recorded...
      if (length(report[grep("Assembly method", report)]) > 0) {
      #Extract the assembly method and add to dataframe
      ncbi.methods.df$method[j] <- sub("# Assembly method: ", "", report[grep("Assembly method", report)])
      }
      
    }
  
  #Create a dataframe of genome sizes and cytometric methods from the fungal genome size database
  no.CS.methods.df <- data.frame(species=no.CS.df$name[grep(names(species.comp)[i], no.CS.df$name)], size=no.CS.df$X1C.in.Mbp[grep(names(species.comp)[i], no.CS.df$name)], method=no.CS.df$METHOD[grep(names(species.comp)[i], no.CS.df$name)], type="noCS")
  
  #Combine the assembly-based and cytometric dataframes
  methods.df <- rbind(ncbi.methods.df, no.CS.methods.df)
  #Remove any rows without methods
  methods.df <- methods.df[!is.na(methods.df$method),]
  #Make species name uniform
  methods.df$species <- names(species.comp)[i]
  
  #Rename with abbreviation
  assign(paste0(abb[i],".methods.df"), methods.df)
}

#Print number of species-level comparisons possible
length(species.comp)

rm(no.CS.methods.df)
rm(ncbi.methods.df)
#Combine all methods dataframes
all.methods.df <- do.call("rbind", mget(objects(pattern="*.methods.df")))
rownames(all.methods.df) <- NULL

#Correct method names
all.methods.df$method[grep("allpaths lg", all.methods.df$method, ignore.case=TRUE)] <- "ALLPATHS-LG"
all.methods.df$method[grep("smrt", all.methods.df$method, ignore.case=TRUE)] <- "SMRT Analysis"
all.methods.df$method[grep("//bgs", all.methods.df$method, ignore.case=TRUE)] <- "Newbler"
all.methods.df$method <- gsub("_", " ", all.methods.df$method)

#Make method names uniform
assemblers <- read.csv("assemblers.csv", header=FALSE)$V1
for (i in assemblers) {
  all.methods.df$method[grep(paste0("\\b", i), all.methods.df$method, ignore.case=TRUE)] <- i
}

#Add asterisk if reported method isn't known
if (length(all.methods.df$method[is.na(match(all.methods.df$method, assemblers))]) > 0) {
  all.methods.df$method[is.na(match(all.methods.df$method, assemblers)) & all.methods.df$type == "CS"] <- paste0(all.methods.df$method[is.na(match(all.methods.df$method, assemblers)) & all.methods.df$type == "CS"], "*")
}

#Abbreviate methods
#methods.df$method <-word(methods.df$method, 1)

#For each species-level case study...
for (i in 1:length(species.comp)) {
  
  #Get the dataframe for the species
  methods.df <- all.methods.df[grep(names(species.comp)[i], all.methods.df$species),]
  
  #If there are at least 3 unique methods with more than 1 measurement...
  if (length(unique(methods.df$method)) > 2 & length(methods.df$method) != length(unique(methods.df$method))) {
    
    #Remove hyphens for Tukey testing
    methods.df$method <- gsub("-", " ", methods.df$method)
    #Tukey significance testing
    tukey <- TukeyHSD(aov(lm(size ~ method, data=methods.df)))
    #Make dataframe for ggplot with tukey groups
    sig.df <- data.frame(multcompLetters(tukey[["method"]][,4])["Letters"])
    sig.df <- data.frame(Treatment=rownames(sig.df), Letters=sig.df$Letters)
    
    if (length(sig.df$Treatment) > 0) {
      sig.df$species <- names(species.comp)[i]
    }
    
    #Rename with abbreviation
    assign(paste0(abb[names(abb) == names(species.comp)[i]],".methods.df"), methods.df)
    assign(paste0(abb[names(abb) == names(species.comp)[i]],".tukey.df"), sig.df)
  }
}

#Find species for which there were enough different methods to perform Tukey testing
species.mainfig <- sub(".tukey.df", "", objects(pattern="*.tukey.df"))
#Combine species dataframes
spec.df <- do.call("rbind", mget(paste0(species.mainfig, ".methods.df")))
rownames(spec.df) <- NULL
spec.df$size <- as.numeric(spec.df$size)
#Combine Tukey dataframe
tukey.df <- do.call("rbind", mget(paste0(species.mainfig, ".tukey.df")))

#Create dataframe for plot labels
labels.df <- unique(spec.df[c(1,3,4)])

for (i in 1:length(labels.df$species)) {
  #Add field with max y position for label
  labels.df$max[i] <- max(spec.df$size[spec.df$species == labels.df$species[i] & spec.df$method == labels.df$method[i]])
  #Add field with sample size
  labels.df$count[i] <- paste0("n=",length(spec.df$size[spec.df$species == labels.df$species[i] & spec.df$method == labels.df$method[i]]))
}

#Create dataframe for mean lines
mean.df <- unique(labels.df[1])
mean.df$type <- "CS"

#Add mean of assembly methods for each species
for (i in 1:length(mean.df$species)) {
  mean.df$mean[i] <- mean(spec.df$size[spec.df$species == mean.df$species[i] & spec.df$type == "CS"])
}

#Add class to dataframe
spec.df$class <- NA

for (i in 1:length(unique(spec.df$species))) {
  spec.df$class[which(spec.df$species == unique(spec.df$species)[i])] <- as.character(tax.df$Class[grep(word(unique(spec.df$species), 1)[i], tax.df$Genus)])
}

#Facet boxplot of genome sizes for each species
gg.spec <- ggplot(spec.df, aes(method, size)) +
  geom_hline(data=mean.df, 
             linetype="dashed", 
             aes(yintercept=mean)) +
  geom_violin(aes(color=type),
              show.legend=FALSE) +
  geom_boxplot(aes(color=type),
               width=0.1,
               outlier.size=1) +
  geom_point(aes(fill=class),
             shape=22,
             colour="white",
             size=5,
             alpha=0) +
  facet_grid(. ~ species,
             scales="free",
             space="free",
             labeller=label_wrap_gen(width=8)) +
  geom_text(data=tukey.df,
            aes(x=Treatment, y=Inf, label=Letters),
            family="mono",
            hjust=1,
            angle=90,
            size=3) +
  geom_text(data=labels.df,
            vjust=-1,
            size=2,
            aes(x=method, y=max, label=count)) +
  scale_y_continuous(breaks=seq(20, max(spec.df$size) + 10, 20),
                     expand=expansion(mult=c(0.05,0.2))) +
  scale_color_manual(labels=c("Genome assembly", "Cytometric"), values=c("black", "red")) +
  scale_fill_manual(values=classnodes$colour[match(sort(unique(spec.df$class)), classnodes$class)]) +
  guides(fill=guide_legend(override.aes=list(alpha=0.5))) +
  labs(subtitle=expression(bold("a")), x="", y="Genome size (Mbp)", col="", fill="", size=2) +
  theme(strip.text=element_text(face="bold.italic", size=8),
        axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)),
        axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
        axis.text.x=element_text(size=rel(0.7), angle=45, hjust=1),
        legend.position="top",
        legend.box="vertical",
        legend.margin=margin(0,0,0,0),
        legend.key=element_blank(),
        plot.title.position="plot",
        plot.margin=unit(c(0,0,-5,0), "mm"))
  
#Create list for colours of species depending on class
spec.cols <- list()

for (i in 1:length(unique(spec.df$species))) {
  spec.cols[[i]] <- classnodes$colour[match(as.character(tax.df$Class[grep(word(unique(spec.df$species), 1)[i], tax.df$Genus)]), classnodes$class)]
}

#Convert list to vector
spec.cols <- unlist(spec.cols)
#Adjust transparency
spec.cols <- adjustcolor(spec.cols, alpha.f = 0.5)

#Convert plot to gtable
gg.spec1 <- ggplot_gtable(ggplot_build(gg.spec))
#Find facet strips
stripr <- which(grepl('strip-t', gg.spec1$layout$name))
#Replace strip colours according to class vector
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', gg.spec1$grobs[[i]]$grobs[[1]]$childrenOrder))
  gg.spec1$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- spec.cols[k]
  k <- k+1
}


## FIGURE 2b ##

#Identify case study species with a large range of sizes from assemblies
species.comp.x <- list()

for (i in 1:length(unique(ncbi$X.Organism.Name))) {
  if (length(ncbi$X.Organism.Name[which(ncbi$X.Organism.Name == unique(ncbi$X.Organism.Name)[i])]) > 2) {
    species.comp.x[[unique(ncbi$X.Organism.Name)[i]]] <- (max(as.numeric(ncbi$Size..Mb.[which(ncbi$X.Organism.Name == unique(ncbi$X.Organism.Name)[i])])) - min(as.numeric(ncbi$Size..Mb.[which(ncbi$X.Organism.Name == unique(ncbi$X.Organism.Name)[i])])))
    }
}

#Print number of species with extreme genome size disparity between different assemblies (> 20Mbp)
paste0(length(species.comp.x[species.comp.x > 20]), "/", length(species.comp.x), ", ", round(length(species.comp.x[species.comp.x > 20]) / length(species.comp.x) * 100), "%")

#Filter for species with range > 20Mbp
species.comp.x <- species.comp.x[species.comp.x > 20]
#Remove duplicate species from previous figure
species.comp.x <- species.comp.x[is.na(match(names(species.comp.x), unique(spec.df$species)))]
#Create abbrevation
abb.x <- abbreviate(names(species.comp.x), minlength=3)

#For each species-level case study...
for (i in 1:length(species.comp.x)) {
  
  #Search for the species name in NCBI assemblies
  assemblies <- assembly.sum[agrep(names(species.comp.x)[i], assembly.sum$organism_name),]
  
  #If assembly links were found for the name...
  if (length(assemblies$X..assembly_accession) > 0) {
    
    #Get NCBI assembly links
    ftp.links <- paste0(assemblies$ftp_path, "/", assemblies$X..assembly_accession,"_", assemblies$asm_name, "_assembly_report.txt")
    
    #Create a results dataframe of genome sizes and assembly methods
    methods.df <- data.frame(species=rep(names(species.comp.x)[i], length(ftp.links)), size=ncbi$Size..Mb.[agrep(names(species.comp.x)[i], ncbi$X.Organism.Name)], method=NA) 
    
    #For each assembly...
    for (j in 1:length(ftp.links)) {
      
      #Try to download the assembly report
      report <- NULL
      report <- tryCatch(unlist(strsplit(getURL(ftp.links[j]), "\\r*\\n")), error = function(e) {e$message})
      
      #If the assembly method is recorded...
      if (length(report[grep("Assembly method", report)]) > 0) {
        #Extract the assembly method and add to dataframe
        methods.df$method[j] <- sub("# Assembly method: ", "", report[grep("Assembly method", report)])
      }
      
    }
    
    #Remove any rows without methods
    methods.df <- methods.df[!is.na(methods.df$method),]
    #Make species name uniform
    methods.df$species <- names(species.comp.x)[i]
    
    #Rename with abbreviation
    assign(paste0(abb.x[i],".x.methods.df"), methods.df)
  }
  
}

#Combine all methods dataframes
all.x.methods.df <- do.call("rbind", mget(objects(pattern="*.x.methods.df")))
rownames(all.x.methods.df) <- NULL

#Correct method names
all.x.methods.df$method[grep("allpaths lg", all.x.methods.df$method, ignore.case=TRUE)] <- "ALLPATHS-LG"
all.x.methods.df$method[grep("\\bgs", all.x.methods.df$method, ignore.case=TRUE)] <- "Newbler"
all.x.methods.df$method[grep("smrt", all.x.methods.df$method, ignore.case=TRUE)] <- "SMRT Analysis"
all.x.methods.df$method <- gsub("_", " ", all.x.methods.df$method)

#Make method names uniform
for (i in assemblers) {
  all.x.methods.df$method[grep(paste0("\\b", i), all.x.methods.df$method, ignore.case=TRUE)] <- i
}

#Add asterisk if reported method isn't known
if (length(all.x.methods.df$method[is.na(match(all.x.methods.df$method, assemblers))]) > 0) {
  all.x.methods.df$method[is.na(match(all.x.methods.df$method, assemblers))] <- paste0(all.x.methods.df$method[is.na(match(all.x.methods.df$method, assemblers))], "*")
}

#For each species-level case study...
for (i in 1:length(species.comp.x)) {
  
  #Get the dataframe for the species
  methods.df <- all.x.methods.df[grep(names(species.comp.x)[i], all.x.methods.df$species),]
  
  #If there are at least 3 unique methods with more than 1 measurement
  if (length(unique(methods.df$method)) > 2 & length(methods.df$method) != length(unique(methods.df$method))) {
    
    #Remove hyphens for Tukey testing
    methods.df$method <- gsub("-", " ", methods.df$method)
    #Tukey significance testing
    tukey <- TukeyHSD(aov(lm(size ~ method, data=methods.df)))
    #Make dataframe for ggplot with tukey groups
    sig.df <- data.frame(multcompLetters(tukey[["method"]][,4])["Letters"])
    sig.df <- data.frame(Treatment=rownames(sig.df), Letters=sig.df$Letters)
    
    if (length(sig.df$Treatment) > 0) {
      sig.df$species <- names(species.comp.x)[i]
    }
    
    #Rename with abbreviation
    assign(paste0(abb.x[names(abb.x) == names(species.comp.x)[i]],".x.methods.df"), methods.df)
    assign(paste0(abb.x[names(abb.x) == names(species.comp.x)[i]],".x.tukey.df"), sig.df)
  }
}

#Find species for which there were enough different methods to perform Tukey testing
species.ncbi.mainfig <- sub(".x.tukey.df", "", objects(pattern="*.x.tukey.df"))
#Combine species dataframes
spec.x.df <- do.call("rbind", mget(paste0(species.ncbi.mainfig, ".x.methods.df")))
rownames(spec.x.df) <- NULL
spec.x.df$size <- as.numeric(spec.x.df$size)
#Combine Tukey dataframe
tukey.x.df <- do.call("rbind", mget(paste0(species.ncbi.mainfig, ".x.tukey.df")))

#Create dataframe for plot labels
labels.x.df <- unique(spec.x.df[c(1,3)])

for (i in 1:length(labels.x.df$species)) {
  #Add field with max y position for label
  labels.x.df$max[i] <- max(spec.x.df$size[spec.x.df$species == labels.x.df$species[i] & spec.x.df$method == labels.x.df$method[i]])
  #Add field with sample size
  labels.x.df$count[i] <- paste0("n=",length(spec.x.df$size[spec.x.df$species == labels.x.df$species[i] & spec.x.df$method == labels.x.df$method[i]]))
}

#Add class to dataframe
spec.x.df$class <- NA

for (i in 1:length(unique(spec.x.df$species))) {
  spec.x.df$class[which(spec.x.df$species == unique(spec.x.df$species)[i])] <- as.character(tax.df$Class[grep(word(unique(spec.x.df$species), 1)[i], tax.df$Genus)])
}

#Plot boxplots
gg.spec.x <- ggplot(spec.x.df, aes(method, size)) +
  geom_violin(position="dodge") +
  geom_boxplot(width=0.1,
               outlier.size=1) +
  facet_grid(. ~ species,
             scales="free",
             space="free",
             labeller=label_wrap_gen(width=7)) +
  geom_point(aes(fill=class),
             shape=22,
             colour="white",
             size=5,
             alpha=0) +
  geom_text(data=tukey.x.df,
            aes(x=Treatment, y=Inf, label=Letters),
            family="mono",
            hjust=1,
            angle=90,
            size=3) +
  geom_text(data=labels.x.df,
            vjust=-1,
            size=2,
            aes(x=method, y=max, label=count)) +
  scale_y_continuous(breaks=seq(20, max(spec.x.df$size) + 10, 20),
                     expand=expansion(mult=c(0.05,0.15))) +
  scale_fill_manual(values=classnodes$colour[match(sort(unique(spec.x.df$class)), classnodes$class)]) +
  guides(fill=guide_legend(override.aes=list(alpha=0.5))) +
  labs(x="Method of genome size inference", y="Genome size (Mbp)", col="", fill="", size=2) +
  theme(strip.text=element_text(face="bold.italic", size=8),
        axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)),
        axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
        axis.text.x=element_text(size=rel(0.7), angle=45, hjust=1),
        legend.position="top",
        legend.box="vertical",
        legend.margin=margin(0,0,0,0),
        legend.key=element_blank(),
        plot.title.position="plot",
        plot.margin=unit(c(0,0,0,0), "mm")) +
  labs(subtitle=expression(bold("b")))

#Create list for colours of species depending on class
spec.x.cols <- list()

for (i in 1:length(unique(spec.x.df$species))) {
  spec.x.cols[[i]] <- classnodes$colour[match(as.character(tax.df$Class[grep(word(unique(spec.x.df$species), 1)[i], tax.df$Genus)]), classnodes$class)]
}

#Convert list to vector
spec.x.cols <- unlist(spec.x.cols)
#Adjust transparency
spec.x.cols <- adjustcolor(spec.x.cols, alpha.f = 0.5)

#Convert plot to gtable
gg.spec.x1 <- ggplot_gtable(ggplot_build(gg.spec.x))
#Find facet strips
stripr <- which(grepl('strip-t', gg.spec.x1$layout$name))
#Replace strip colours according to class vector
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', gg.spec.x1$grobs[[i]]$grobs[[1]]$childrenOrder))
  gg.spec.x1$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- spec.x.cols[k]
  k <- k+1
}


#Plot species-level comparison together
tiff(file=paste0("speccomparisonfig_", Sys.Date(), ".tiff"), height=8, width=10, units="in", res=300)
grid.arrange(gg.spec1, gg.spec.x1, heights=c(1.1, 1))
dev.off()


## SUPPLEMENTARY FIGURE 2 ##

#Combine species dataframes not shown in Fig2a with both cytometric and assembly-based methods
spec.supp.df <- do.call("rbind", mget(paste0(abb[is.na(match(abb, species.mainfig))], ".methods.df")))
rownames(spec.supp.df) <- NULL
species.suppfig <- unique(spec.supp.df$species)

#Correct method names
spec.supp.df$method[grep("allpaths lg", spec.supp.df$method, ignore.case=TRUE)] <- "ALLPATHS-LG"
spec.supp.df$method[grep("\\bgs", spec.supp.df$method, ignore.case=TRUE)] <- "Newbler"
spec.supp.df$method[grep("smrt", spec.supp.df$method, ignore.case=TRUE)] <- "SMRT Analysis"
spec.supp.df$method[grep("light microscopy", spec.supp.df$method, ignore.case=TRUE)] <- "LM"
spec.supp.df$method <- gsub("_", " ", spec.supp.df$method)

#Make method names uniform
for (i in assemblers) {
  spec.supp.df$method[grep(paste0("\\b", i), spec.supp.df$method, ignore.case=TRUE)] <- i
}

#Remove extralong pipeline
spec.supp.df <- spec.supp.df[-which(sapply(strsplit(spec.supp.df$method, " "), length) > 2),]

#Add asterisk if reported method isn't known
if (length(spec.supp.df$method[is.na(match(spec.supp.df$method, assemblers))]) > 0) {
  spec.supp.df$method[is.na(match(spec.supp.df$method, assemblers)) & spec.supp.df$type == "CS"] <- paste0(spec.supp.df$method[is.na(match(spec.supp.df$method, assemblers)) & spec.supp.df$type == "CS"], "*")
}

#Remove species missing known assembly-based methods
for (i in 1:length(species.suppfig)) {
  if (length(unique(spec.supp.df$type[which(spec.supp.df$species == species.suppfig[i])])) < 2) {
    spec.supp.df <- spec.supp.df[-which(spec.supp.df$species == species.suppfig[i]),]
  }
}

#Make sizes numeric
spec.supp.df$size <- as.numeric(spec.supp.df$size)
#Make dataframe of species and their sizes to sort by size
spec.supp.order <- data.frame(species=unique(spec.supp.df$species, max=NA))

for (i in 1:length(spec.supp.order$species)) {
  spec.supp.order$max[i] <- max(spec.supp.df$size[spec.supp.df$species == spec.supp.order$species[i]])
}

#Sort species in order of size (low to high)
spec.supp.order <- spec.supp.order[order(spec.supp.order$max),]
#Sort main dataframe by size
spec.supp.df <- spec.supp.df[order(factor(spec.supp.df$species, levels=unique(spec.supp.order$species))),]

#Assign species to rows of 5 each for plot
rows <- rep(c(1:ceiling(length(unique(spec.supp.df$species)) / 5)),
            each=ceiling(length(unique(spec.supp.df$species)) / ceiling(length(unique(spec.supp.df$species)) / 5)),
            length.out=length(unique(spec.supp.df$species)))

for (i in 1:length(unique(spec.supp.df$species))) {
  spec.supp.df$row[spec.supp.df$species == unique(spec.supp.df$species)[i]] <- rows[i]
}

#Create dataframe for plot labels
labels.supp.df <- unique(spec.supp.df[c(1, 3, 4, 5)])

for (i in 1:length(labels.supp.df$species)) {
  #Add field with max y position for label
  labels.supp.df$max[i] <- max(spec.supp.df$size[spec.supp.df$species == labels.supp.df$species[i] & spec.supp.df$method == labels.supp.df$method[i]])
  #Add field with sample size
  labels.supp.df$count[i] <- paste0("n=",length(spec.supp.df$size[spec.supp.df$species == labels.supp.df$species[i] & spec.supp.df$method == labels.supp.df$method[i]]))
}

#Plot boxplots for each row of 5 species
for (i in 1:ceiling(length(unique(spec.supp.df$species)) / 5)) {
  
  gg.spec.supp <- ggplot(spec.supp.df[spec.supp.df$row == i,], aes(method, size, colour=type)) +
    geom_violin(position="dodge") +
    geom_boxplot(width=0.1) +
    facet_grid(. ~ species,
               scales="free",
               space="free",
               labeller=label_wrap_gen(width=7)) +
    geom_text(data=labels.supp.df[labels.supp.df$row == i,],
              vjust=-1,
              size=2,
              aes(x=method, y=max, label=count),
              inherit.aes=FALSE) +
    labs(y="Genome size (Mbp)", col="", size=2) +
    scale_y_continuous(limits=c(0, max(spec.supp.df$size) + 5),
                       breaks=seq(20, max(spec.supp.df$size), 20),
                       expand=expansion(mult=c(0.05,0.05))) +
    scale_color_manual(labels=c("Genome assembly", "Cytometric"), values=c("black", "red")) +
    theme(strip.text=element_text(face="italic", size=8),
          axis.title.x=element_blank(),
          axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
          axis.text.x=element_text(size=rel(0.7), angle=45, hjust=1),
          legend.title=element_blank(),
          plot.title.position="plot",
          plot.margin=unit(c(0,0,0,0), "mm"))
  
  assign(paste0("gg.spec.supp.", i), gg.spec.supp)
}

rm(gg.spec.supp)
#Create list of boxplots
supps <- mget(objects(patter="gg.spec.supp.*"))
#Combine boxplots into one plot
gg.spec.supp <- do.call(ggarrange, c(supps, ncol=1, common.legend=TRUE, legend="top"))
#Add x axis label
gg.spec.supp <- annotate_figure(gg.spec.supp,
                                bottom=text_grob("Method of genome size inference"))

#Write to file
tiff(file=paste0("speccomparisonsuppfig_", Sys.Date(), ".tiff"), height=12, width=8, units="in", res=300)
plot(gg.spec.supp)
dev.off()