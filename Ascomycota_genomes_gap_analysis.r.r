######Ascomycota genome gap analysis###########
##This script uses some webpage scraping, so results will differ as new data becomes available##

library(ape)
library(rvest)
library(magrittr)
library(ggplot2)
library(gridExtra)
library(ggstance)
library(ggtree)
library(colorspace)
library(gtable)
library(scales)
library(grid)
library(multcompView)

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
#Create vector of method to exclude (genome assembly, unreliable or unknown methods)
exclude <- c("CS", "genomic reconstruction", "", "CHEF gel electrophoresis", "DAPI-PC", "DAPI-IC", "PFGE", "CS and PFGE", "quantitative real-time PCR", "Re-association kinetics")
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
  myc$i <- tax.df$i[match(myc$Genus, tax.df$Genus)]
  ncbi$i <- tax.df$i[match(ncbi$Genus, tax.df$Genus)]
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
branches <- data.frame(node=asc.tree$edge[,2], edge_num=1:nrow(asc.tree$edge), branch_length=rep(NA,nrow(asc.tree$edge)))
#Read in tree node labels for classes
classnodes <- read.csv("classnodes.csv", header=TRUE)
#Add class
for (i in 1:length(classnodes$class)) {
  branches$branch_length[branches$edge_num == classnodes$node[i]] <- as.character(classnodes$class[i])
}
#Order by node
branches <- branches[order(branches$node),]
#Add colour
branches$colour <- classnodes$colour[match(branches$branch_length, classnodes$class)]

#Make dataframe matching taxonomy data with the Ascomycota tree tip labels for colours
tax.df3 <- data.frame(order=asc.tree$tip.label, class=tax.df2$Class[match(asc.tree$tip.label, unique(tax.df2$Order))], data=NA)
#Add column for whether the order has CS genome size data
tax.df3$data[!is.na(match(tax.df3$order, CS.df$order))] <- "CS"
#Add column for whether the order has no CS genome size data
tax.df3$data[!is.na(match(tax.df3$order, no.CS.df2$order))] <- "noCS"
tax.df3$data[is.na(tax.df3$data)] <- "None"

#Make dataframe of class data to colour class branches
class.df <- data.frame(class=unique(tax.df3$class), data="None")
#Add genome data to class dataframe
class.df$data[match(unique(tax.df3$class[which(tax.df3$data == "CS")]), class.df$class)] <- "CS"
class.df$data[match(unique(tax.df3$class[which(tax.df3$data == "noCS")]), class.df$class)] <- "noCS"
#Add genome data to branches dataframe
branches$class_branch <- class.df$data[match(branches$branch_length, class.df$class)]

#Dataframes for vertical grid lines
box.lines <- data.frame(x=seq(0,220,10), .panel="Genome size (Mbp/1C)", stringsAsFactors=TRUE)
noCSbox.lines.out <- data.frame(x=seq(0,4000,100), .panel="Genome size (Mbp/1C)", stringsAsFactors=TRUE)
bar.lines <- data.frame(x=seq(0,1100,50), .panel="Number of genome assemblies", stringsAsFactors=TRUE)
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
counts.df <- counts.df[which(counts.df$data != "None"),]
#Add column with count of sample size for no CS/CS
counts.df$noCScount <- table(unlist(no.CS.df2$order))[match(counts.df$order, names(table(unlist(no.CS.df2$order))))]
counts.df$CScount <- table(unlist(CS.df2$order))[match(counts.df$order,names(table(unlist(CS.df2$order))))]
#Put brackets around count numbers
counts.df$noCScount <- paste0("(", counts.df$noCScount,")")
counts.df$CScount <- paste0("(", counts.df$CScount,")")

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
  counts.df[counts.df$order==i,]$noCSmax <- max(no.CS.df2[no.CS.df2$order==i & no.CS.df2$size < 210,]$size)
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
              aes(color=data)) +
  scale_colour_manual(values=c("black","red","darkgrey"))
gg.dummy2 <- gg.dummy1 %<+% branches + 
  geom_label(aes(x=branch, label=branch_length, colour=class_branch),
             fill="white",
             size=2.8,
             label.padding=unit(0.15, "lines"),
             label.size=0)
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
                  data=no.CS.df2,
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
  xlim(0,220) +
  scale_x_continuous(breaks=seq(0,220,by=20), limits=c(0,220),
                     position="top",
                     sec.axis = dup_axis()) +
  scale_y_continuous(expand=expansion(mult=0.005)) +
  scale_fill_manual(values=c("snow3","white",as.vector(classnodes$colour[classnodes$colour != ""]))) +
  theme_classic() +
  theme(panel.border=element_blank(),
        legend.position="none",
        strip.placement="outside",
        axis.text.x.top = element_text(angle=45, vjust=0.3),
        axis.text.x.bottom = element_text(angle=45, hjust=1),
        #panel.spacing=unit(0, "lines"),
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
              aes(color=data)) +
  scale_colour_manual(values=c("black","red","darkgrey"))
gg.main2 <- gg.main1 %<+% branches + 
  geom_label(aes(x=branch, label=branch_length, colour=class_branch),
             fill="white",
             size=2.8,
             label.padding=unit(0.15, "lines"),
             label.size=0)
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
gg.main4 <- facet_plot(gg.main3,
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
                     position="top",
                     sec.axis = dup_axis()) +
  scale_y_continuous(expand=expansion(mult=0.005)) +
  scale_fill_manual(values=c("snow3","white",as.vector(classnodes$colour[classnodes$colour != ""]))) +
  theme_classic() +
  theme(panel.border=element_blank(),
        legend.position="none",
        strip.placement="outside",
        axis.text.x.top = element_text(angle=45, vjust=0.3),
        axis.text.x.bottom = element_text(angle=45, hjust=1),
        #strip.text.x=element_blank(),
        #panel.spacing=unit(0, "lines"),
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

gg.supp <- gg.supp + geom_tree()
gg.supp1 <- gg.supp %<+% tax.df3 +
  geom_tiplab(size=2.8,
              aes(color=data)) +
  scale_colour_manual(values=c("black","red","darkgrey"))
gg.supp2 <- gg.supp1 +
  geom_vline(data=noCSbox.lines.out,
             color="grey94",
             aes(xintercept=x))
gg.supp3 <- facet_plot(gg.supp2,
                       panel="Genome size (Mbp/1C)",
                       data=no.CS.df2,
                       geom=geom_boxploth,
                       #outlier.shape=NA,
                       outlier.size=1,
                       aes(x=size, group=label, fill=class))
gg.supp3 <- facet_plot(gg.supp3,
                       panel="Genome size (Mbp/1C)",
                       data=CS.df2,
                       geom=geom_boxploth,
                       #outlier.shape=NA,
                       outlier.size=1,
                       colour="darkgrey",
                       linetype="dotted",
                       alpha=0.3,
                       aes(x=size, group=label, fill=class))
gg.supp4 <- gg.supp3 %<+% branches + 
  geom_label(aes(x=branch, label=branch_length, colour=class_branch),
             fill="white",
             size=2.8,
             label.padding = unit(0.15, "lines"),
             label.size = 0) +
  scale_x_continuous(breaks=pretty_breaks(10),
                     position="top",
                     sec.axis = dup_axis()) +
  scale_y_continuous(expand=expansion(mult=0.005)) +
  scale_fill_manual(values=c("snow3","white",as.vector(classnodes$colour[classnodes$colour != ""]))) +
  theme_classic() +
  theme(panel.border=element_blank(),
        legend.position="none",
        strip.placement="outside",
        strip.background=element_blank())

#Remove unwanted axes on tree panel
gg.supp.tab <- ggplotGrob(gg.supp4)
elements <- gg.supp.tab$layout$name
gg.supp5 <- gtable_filter_remove(gg.supp.tab, name=elements[c(4,6,8,10)], trim=FALSE)

tiff(file=paste0("genomesizeplot_outliers_", Sys.Date(), ".tiff"), height=15, width=10, units="in", res=300)

grid.draw(gg.supp5)

dev.off()



##SPECIES-LEVEL COMPARISON FIGURE##

##Identify case study species with both cytometric and assembly-based measurements
#Add name field to genome size dataframe 
no.CS.df$name <- paste0(no.CS.df$GENUS," ",no.CS.df$SPECIES)

#Print species with >3 cytometric estimations
for (i in 1:length(no.CS.df$name)) {
  if (length(no.CS.df$name[no.CS.df$name == no.CS.df$name[i]]) > 3) {
    print(paste0(no.CS.df$name[i]," ",no.CS.df$METHOD[i]))
  }
}

#Make vector of case study species
spec <- c("Aspergillus niger", "Paracoccidioides brasiliensis", "Sporothrix schenckii", "Venturia inaequalis")
#Make vector of abbreviations
abb <- substring(spec, 1, 3)

#Make directory for assembly reports
dir.create("assembly_reports")

#Create files of ftp links for ncbi assembly reports
for (i in 1:length(spec)) {
  temp <- assembly.sum[grepl(spec[i], assembly.sum$organism_name),]
  write(paste0(temp$ftp_path,"/",temp$X..assembly_accession,"_",temp$asm_name,"_assembly_report.txt"), file=paste0("assembly_reports/",abb[i],"_links"), ncolumns=1)
}

#Add assembly methods to myc
myc$Methods <- NA
myc$Methods[grep("Asplac1", myc$Portal)] <- "Velvet"
myc$Methods[grep("Aspph1", myc$Portal)] <- "Velvet"
myc$Methods[grep("Aspni_NRRL3_1", myc$Portal)] <- "Velvet"
myc$Methods[grep("Aspni_bvT_1", myc$Portal)] <- "Velvet"
myc$Methods[grep("Aspni7", myc$Portal)] <- "JAZZ"
myc$Methods[grep("Aspni_DSM_1", myc$Portal)] <- "Phrap"
myc$Methods[grep("Parbr1", myc$Portal)] <- "Arachne"
myc$Methods[grep("Parbra1", myc$Portal)] <- "Arachne"

#Download assembly reports in command line and save as *species abbreviation*_methods

#For each species...
for (i in 1:length(spec)) {

  #Read in assembly reports
  temp.methods <- read.csv(paste0("assembly_reports/",abb[i],"_methods"), sep = '#', header = FALSE)
  #Modify ncbi assembly methods file
  temp.methods$V1 <- gsub("_assembly_report.txt:", "", temp.methods$V1)
  temp.methods$V1 <- substr(temp.methods$V1, 0, 15)
  temp.methods$V2 <- sub("Assembly method: ", "", temp.methods$V2)
  temp.methods$V2 <- sub("Spades", "SPAdes", temp.methods$V2)
  temp.methods$V2 <- sub("allpaths", "ALLPATHS", temp.methods$V2)
  temp.methods$V2 <- trimws(temp.methods$V2)
  temp.methods$V2 <- gsub(" .*", "", temp.methods$V2)
  
  #Create ncbi subset dataframe
  temp.ncbi <- ncbi[grepl(spec[i], ncbi$X.Organism.Name),]
  #Add assembly method
  temp.ncbi["methods"] <- temp.methods$V2[match(temp.ncbi$Assembly.Accession, temp.methods$V1)]
  #Add CS type coding
  temp.ncbi <- data.frame(species=temp.ncbi$X.Organism.Name, size=temp.ncbi$Size..Mb., method=temp.ncbi$methods, type="CS")
  #Create Mycocosm subset dataframe
  temp.myc <- myc[grepl(spec[i], myc$Name),]
  
  if (nrow(temp.myc) != 0) {
    #Add CS type coding
    temp.myc <- data.frame(species=temp.myc$Name, size=temp.myc$Assembly.Length, method=temp.myc$Methods, type="CS")
  }
  
  #Create fungal genome size database subset dataframe
  temp.no.CS <- data.frame(species=no.CS.df$name[grepl(spec[i], no.CS.df$name)], size=no.CS.df$X1C.in.Mbp[grepl(spec[i], no.CS.df$name)], method=no.CS.df$METHOD[grepl(spec[i], no.CS.df$name)], type="noCS")
  
  #Combine dataframes
  if (nrow(temp.myc) != 0) {
    temp.df <- rbind(temp.ncbi, temp.myc, temp.no.CS)
  } else {
    temp.df <- rbind(temp.ncbi, temp.no.CS)
  }
  #Format species name to be the same
  temp.df$species <- spec[i]
  #Remove any rows without methods
  temp.df <- temp.df[!is.na(temp.df$method),]
  temp.df <- temp.df[!temp.df$method == "",]
  #Remove hyphens for Tukey testing
  temp.df$method <- gsub("-", " ", temp.df$method)
  
  #Tukey significance testing
  temp.tukey <- TukeyHSD(aov(lm(size ~ method, data=temp.df)))
  #Make dataframe for ggplot with tukey groups
  temp.sig.df <- data.frame(multcompLetters(temp.tukey[["method"]][,4])["Letters"])
  temp.sig.df <- data.frame(Treatment=rownames(temp.sig.df), Letters=temp.sig.df$Letters)
  
  if (length(temp.sig.df$Treatment) > 0) {
    temp.sig.df$species <- spec[i]
  }
  
  #Rename with abbreviation
  assign(paste0(abb[i],"_methods"), temp.methods)
  assign(paste0(abb[i],".df"), temp.df)
  assign(paste0(abb[i],".tukey"), temp.tukey)
  assign(paste0(abb[i],".sig.df"), temp.sig.df)
}

#Combine all species dataframes
spec.df <- rbind(Asp.df, Spo.df, Par.df, Ven.df)
spec.df$size <- as.numeric(spec.df$size)

#Create dataframe for Tukey labels
tukey.df <- rbind(Asp.sig.df, Spo.sig.df, Par.sig.df, Ven.sig.df)

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

#Facet boxplot of genome sizes for each species
gg.spec <- ggplot(spec.df, aes(method, size)) +
  geom_hline(data=mean.df, 
             linetype="dashed", 
             aes(yintercept=mean)) +
  geom_violin(aes(color=type),
              show.legend=FALSE) +
  geom_boxplot(aes(color=type),
               width=0.1,
               position="dodge") +
  facet_grid(. ~ species,
             scales="free",
             space="free",
             labeller=label_wrap_gen()) +
  geom_text(data=tukey.df,
            aes(x=Treatment, y=Inf, label=Letters),
            family="mono",
            vjust=1.5,
            hjust=0.5,
            size=3) +
  geom_text(data=labels.df,
            vjust=-1,
            size=3,
            aes(x=method, y=max, label=count)) +
  scale_y_continuous(expand=expansion(mult= c(0.05,0.2))) +
  labs(x="", y="Genome size (Mbp)", col="", size=2) +
  scale_color_manual(labels=c("Genome assembly", "Cytometric"), values=c("black", "red")) +
  theme(strip.text=element_text(face="italic", size=10),
        axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)),
        axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
        axis.text.x=element_text(size=rel(0.7), angle=45, hjust=1),
        legend.position="top",
        plot.title.position="plot",
        plot.margin=unit(c(0,0,-5,0), "mm")) +
  labs(subtitle=expression(bold("a")))

gg.spec <- ggplot(spec.df, aes(method, size)) +
  geom_hline(data=mean.df, 
             linetype="dashed", 
             aes(yintercept=mean)) +
  geom_violin(aes(color=type),
              show.legend=FALSE) +
  geom_boxplot(aes(color=type),
               width=0.1,
               position="dodge") +
  facet_wrap(. ~ species,
             scales="free",
             ncol=2) +
  geom_text(data=tukey.df,
            aes(x=Treatment, y=Inf, label=Letters),
            family="mono",
            vjust=1.5,
            hjust=0.5,
            size=3) +
  geom_text(data=labels.df,
            vjust=-1,
            size=3,
            aes(x=method, y=max, label=count)) +
  scale_y_continuous(expand=expansion(mult= c(0.05,0.3))) +
  labs(x="", y="Genome size (Mbp)", col="", size=2) +
  scale_color_manual(labels=c("Genome assembly", "Cytometric"), values=c("black", "red")) +
  theme(strip.text=element_text(face="italic", size=10),
        axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)),
        axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
        axis.text.x=element_text(size=rel(0.7), angle=45, hjust=1),
        legend.position="top",
        plot.title.position="plot",
        plot.margin=unit(c(0,0,-5,0), "mm")) +
  labs(subtitle=expression(bold("a")))


#Identify case study species with a large range of sizes from assemblies
assemblies <- list()

for (i in 1:length(unique(ncbi$X.Organism.Name))) {
  if (length(ncbi$X.Organism.Name[ncbi$X.Organism.Name == ncbi$X.Organism.Name[i]]) > 2) {
    assemblies[[ncbi$X.Organism.Name[i]]] <- as.vector(ncbi$Assembly.Accession[i])
    }
}

assemblies.df <- data.frame(species=names(assemblies), assembly=unlist(assemblies))
assemblies.df <- assemblies.df[order(assemblies.df$species),]

sizes <- list()

for (i in 1:length(assemblies.df$species)) {
  sizes[[as.vector(assemblies.df$species[i])]] <- (max(as.numeric(ncbi$Size..Mb.[ncbi$X.Organism.Name == assemblies.df$species[i]])) - min(as.numeric(ncbi$Size..Mb.[ncbi$X.Organism.Name == assemblies.df$species[i]])))
}

#Create vector of species with extreme genome size disparity between different assemblies
spec.x <- c("Calonectria pseudonaviculata", "Cercospora sojina", "Fusarium proliferatum", "Hortaea werneckii", "Macrophomina phaseolina", "Ophiognomonia clavigignenti-juglandacearum")
#Create abbrevation
abb.x <- substring(spec.x, 1, 3)

#Create files of ftp links for ncbi assembly reports
for (i in 1:length(spec.x)) {
  temp <- assembly.sum[grepl(spec.x[i], assembly.sum$organism_name),]
  write(paste0(temp$ftp_path,"/",temp$X..assembly_accession,"_",temp$asm_name,"_assembly_report.txt"), file=paste0("assembly_reports/",abb.x[i],"_links"), ncolumns=1)
}

#Download assembly reports in command line and save as *species abbreviation*_methods

#For each species...
for (i in 1:length(spec.x)) {
  
  #Read in assembly reports
  temp.methods <- read.csv(paste0("assembly_reports/",abb.x[i],"_methods"), sep = '#', header = FALSE)
  #Modify methods file
  temp.methods$V1 <- gsub("_assembly_report.txt:", "", temp.methods$V1)
  temp.methods$V1 <- substr(temp.methods$V1, 0, 15)
  temp.methods$V2 <- sub("Assembly method: ", "", temp.methods$V2)
  temp.methods$V2 <- gsub("Spades", "SPAdes", temp.methods$V2)
  temp.methods$V2 <- trimws(temp.methods$V2)
  temp.methods$V2 <- gsub(" .*", "", temp.methods$V2)
  
  #Create ncbi subset dataframe
  temp.ncbi <- ncbi[grepl(spec.x[i], ncbi$X.Organism.Name),]
  #Add assembly method
  temp.ncbi["methods"] <- temp.methods$V2[match(temp.ncbi$Assembly.Accession, temp.methods$V1)]
  temp.ncbi <- data.frame(species=temp.ncbi$X.Organism.Name, size=temp.ncbi$Size..Mb., method=temp.ncbi$methods)
  #Add field for species
  temp.ncbi$species <- spec.x[i]
  #Remove any rows without methods
  temp.ncbi <- temp.ncbi[!is.na(temp.ncbi$method),]
  temp.ncbi <- temp.ncbi[!temp.ncbi$method == "",]
  
  if (length(unique(temp.ncbi$method)) > 2 & length(temp.ncbi$method) != length(unique(temp.ncbi$method))) {
    #Remove hyphens for Tukey testing
    temp.ncbi$method <- gsub("-", " ", temp.ncbi$method)
    #Tukey significance testing
    temp.tukey.x <- TukeyHSD(aov(lm(size ~ method, data=temp.ncbi)))
    #Make dataframe for ggplot with tukey groups
    temp.sig.x.df <- data.frame(multcompLetters(temp.tukey.x[["method"]][,4])["Letters"])
    temp.sig.x.df <- data.frame(Treatment=rownames(temp.sig.x.df), Letters=temp.sig.x.df$Letters)
    
    if (length(temp.sig.x.df$Treatment) > 0) {
      temp.sig.x.df$species <- spec.x[i]
    }
    
    #Rename with abbreviation
    assign(paste0(abb.x[i],".tukey.x"), temp.tukey.x)
    assign(paste0(abb.x[i],".sig.x.df"), temp.sig.x.df)
    
  }
  
  #Rename with abbreviation
  assign(paste0(abb.x[i],"_methods"), temp.methods)
  assign(paste0(abb.x[i],".ncbi"), temp.ncbi)
  
}

#Combine all species dataframes
spec.x.df <- rbind(Oph.ncbi, Mac.ncbi, Cer.ncbi, Cal.ncbi, Hor.ncbi, Fus.ncbi)
spec.x.df$size <- as.numeric(spec.x.df$size)

#Create dataframe for Tukey labels
tukey.x.df <- rbind(Mac.sig.x.df, Hor.sig.x.df, Fus.sig.x.df)

#Create dataframe for plot labels
labels.x.df <- unique(spec.x.df[c(1,3)])

for (i in 1:length(labels.x.df$species)) {
  #Add field with max y position for label
  labels.x.df$max[i] <- max(spec.x.df$size[spec.x.df$species == labels.x.df$species[i] & spec.x.df$method == labels.x.df$method[i]])
  #Add field with sample size
  labels.x.df$count[i] <- paste0("n=",length(spec.x.df$size[spec.x.df$species == labels.x.df$species[i] & spec.x.df$method == labels.x.df$method[i]]))
}
  
#Plot boxplots
gg.spec.x <- ggplot(spec.x.df, aes(method, size)) +
  geom_violin(position="dodge") +
  geom_boxplot(width=0.1) +
    facet_grid(. ~ species,
               scales="free",
               space="free",
               labeller=label_wrap_gen()) +
    geom_text(data=tukey.x.df,
              aes(x=Treatment, y=Inf, label=Letters),
              family="mono",
              vjust=1.5,
              hjust=0.5,
            size=3) +
  geom_text(data=labels.x.df,
            vjust=-1,
            size=3,
            aes(x=method, y=max, label=count)) +
  scale_y_continuous(breaks=seq(15,115,15),
                     expand=expansion(mult=c(0.05,0.15))) +
  labs(x="Method of genome size inference", y="Genome size (Mbp)", col="", size=2) +
  theme(strip.text=element_text(face="italic", size=5),
        axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)),
        axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
        axis.text.x=element_text(size=rel(0.7), angle=45, hjust=1),
        plot.title.position="plot",
        plot.margin=unit(c(0,0,0,0), "mm")) +
  labs(subtitle=expression(bold("b")))

gg.spec.x <- ggplot(spec.x.df, aes(method, size)) +
  geom_violin(position="dodge") +
  geom_boxplot(width=0.1) +
  facet_wrap(. ~ species,
             ncol=2,
             scales="free") +
  geom_text(data=tukey.x.df,
            aes(x=Treatment, y=Inf, label=Letters),
            family="mono",
            vjust=1.5,
            hjust=0.5,
            size=3) +
  geom_text(data=labels.x.df,
            vjust=-1,
            size=3,
            aes(x=method, y=max, label=count)) +
  scale_y_continuous(breaks=seq(15,115,15),
                     expand=expansion(mult=c(0.05,0.3))) +
  labs(x="Method of genome size inference", y="Genome size (Mbp)", col="", size=2) +
  theme(strip.text=element_text(face="italic", size=10),
        axis.title.x=element_text(margin=margin(t=20, r=0, b=0, l=0)),
        axis.title.y=element_text(margin=margin(t=0, r=20, b=0, l=0)),
        axis.text.x=element_text(size=rel(0.7), angle=45, hjust=1),
        plot.title.position="plot",
        plot.margin=unit(c(0,0,0,0), "mm")) +
  labs(subtitle=expression(bold("b")))

#Plot species-level comparison together
tiff(file=paste0("speccomparisonfig_", Sys.Date(), ".tiff"), height=12, width=8, units="in", res=300)
grid.arrange(gg.spec, gg.spec.x, heights=c(1.5,2))
dev.off()

tiff(file=paste0("speccomparisonfig_", Sys.Date(), ".tiff"), height=9, width=8, units="in", res=300)
grid.arrange(gg.spec, gg.spec.x)
dev.off()





#Number/proportion of orders without cytometric genome size data
paste0(length(tax.df3$order[tax.df3$data != "noCS"]), "/",length(tax.df3$order), ", ", (length(tax.df3$order[tax.df3$data != "noCS"])) / length(tax.df3$order) * 100,"%")

#Number/proportion of classes without cytometric genome size data
paste0(length(class.df[-20,]$class[class.df[-20,]$data != "noCS"]), "/",length(class.df[-20,]$class), ", ", (length(class.df[-20,]$class[class.df[-20,]$data != "noCS"])) / length(class.df[-20,]$class) * 100,"%")

