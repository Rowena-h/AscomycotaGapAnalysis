######Ascomycota genome gap analysis###########

setwd("~/Kew/Ascomycota/Gap analyses")

library(ape)
library(ggplot2)
require(gridExtra)
library(ggstance)
library(ggtree)
library(colorspace)
library(gtable)
library(scales)
library(grid)

##Genome size data (not from assemblies)##

#Read in genome size data
df <- read.csv("Fungi genome sizes.csv")
#Subset genome size dataframe for just Ascomycota
asc.df <- subset(df, PHYLUM == "Ascomycota")
#Remove genome size data based on genome assembly or unknown methods
no.CS.df <- asc.df[!(asc.df$METHOD == "CS"),]
no.CS.df <- no.CS.df[!(no.CS.df$METHOD == ""),]
#Remove rows with no order classification
#no.CS.df <- no.CS.df[!(no.CS.df$ORDER == ""),]
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
no.CS.res <- no.CS.res[lapply(no.CS.res, length) >= 5]
#Convert results list to vector
no.CS.taxa <- setNames(unlist(no.CS.res, recursive=TRUE, use.names=FALSE), rep(names(no.CS.res), lengths(no.CS.res)))
#Make dataframe of genome sizes for orders
no.CS.df2 <- cbind(read.table(text=names(no.CS.taxa)), no.CS.taxa)
colnames(no.CS.df2) <- c("order", "size")
#Identify outliers
outliers <- boxplot(no.CS.df2$size, plot=FALSE)$out
no.CS.df2.out <- no.CS.df2
#Make dataframe excluding outliers
no.CS.df2 <- no.CS.df2[-which(no.CS.df2$size %in% outliers),]


##Genome size data (from assemblies)##

ncbi <- read.csv("NCBI_genomes2.csv")
myc <- read.csv("Mycocosm_Genomes.csv")
ncbi$Genus <- as.factor(gsub(" .*", '', ncbi$X.Organism.Name))
myc$Genus <- as.factor(gsub(" .*", '', myc$Species))
myc$Genome.Size <- myc$Genome.Size / 1000000
ncbitax <- read.csv("NCBI_asc_taxonomy.csv", header=FALSE)
colnames(ncbitax) <- c("ID", "Class", "Order", "Family", "Genus", "Name")
myc$Family <- as.factor(ncbitax$Family[match(myc$Genus, ncbitax$Genus)])
ncbi$Family <- as.factor(ncbitax$Family[match(ncbi$Genus, ncbitax$Genus)])
ncbi$Order <- as.factor(ncbitax$Order[match(ncbi$Genus, ncbitax$Genus)])
ncbi$Class <- as.factor(ncbitax$Class[match(ncbi$Genus, ncbitax$Genus)])

#Make dataframe of genome sizes
CS.df <- data.frame(order=c(as.character(ncbi$Order), as.character(myc$Order)), size=c(as.character(ncbi$Size.Mb.), as.character(myc$Genome.Size)))
#Convert sizes from factor to numeric
CS.df$size <- as.numeric(levels(CS.df$size))[CS.df$size]
#Remove incertae sedis and NA orders
CS.df2 <- CS.df[!is.na(CS.df$order),]
CS.df2 <- CS.df2[!grepl("incertae sedis", CS.df2$order),]
#Make vector of orders
CS.orders <- as.vector(unique(CS.df2$order))
#Filter for orders with more than 5 genome size measurements
for (i in 1:length(CS.orders)) {
  if (length(CS.df2$order[grepl(CS.orders[i], CS.df2$order)]) < 5) {
    CS.df2 <- CS.df2[!grepl(CS.orders[i], CS.df2$order),]
  }
}
#Identify outliers
outliers <- boxplot(CS.df2$size, plot=FALSE)$out
CS.df2.out <- CS.df2
#Make dataframe excluding outliers
CS.df2 <- CS.df2[-which(CS.df2$size %in% outliers),]


##Generate Ascomycota order-level phylogeny for side by side plot##

#Read in Ascomycota taxonomy
tax.df <- read.csv("Ascomycota outline 2017.csv")
#Add phylum column
tax.df$Phylum <- as.factor("Ascomycota")
#Remove incertae sedis orders
tax.df2 <- tax.df[!grepl("incertae sedis", tax.df$Order, ignore.case=TRUE),]
#Remove duplicate orders
tax.df2 <- tax.df2[!duplicated(tax.df2$Order), ]
#Make tree with dataframe
asc.tree <- as.phylo(~Phylum/Class/Order, data=tax.df2)

#Create dataframe for branch labels
branches <- data.frame(node=asc.tree$edge[,2], edge_num=1:nrow(asc.tree$edge), branch_length=rep(NA,nrow(asc.tree$edge)))
#Read in tree node labels for classes
classnodes <- read.csv("classnodes.csv", header=TRUE)
#Add nodes to branches dataframe
for (i in 1:length(classnodes$class)) {
  branches$branch_length[classnodes$node[i]] <- as.character(classnodes$class[i])
}

#Make dataframe matching taxonomy data with the Ascomycota tree tip labels
tax.df3 <- tax.df2$Class[match(asc.tree$tip.label, tax.df2$Order)]
tax.df3 <- data.frame(order=asc.tree$tip.label, class=tax.df3, data=NA)
#Add column for whether the order has CS genome size data
tax.df3$data[!is.na(match(tax.df3$order, CS.df$order))] <- "CS"
#Add column for whether the order has no CS genome size data
tax.df3$data[!is.na(match(tax.df3$order, no.CS.df2$order))] <- "noCS"
tax.df3$data[is.na(tax.df3$data)] <- "None"

#Make vector of classes
classes <- as.vector(unique(tax.df3$class))
#Make dataframe of class data
class.df <- data.frame(class=classes, data="None", stringsAsFactors = FALSE)
#Extract classes which have genome data
truenoCS <- tax.df3[which(tax.df3$data == "noCS"),]
truenoCS <- truenoCS[!duplicated(truenoCS$class),]
trueCS <- tax.df3[which(tax.df3$data == "CS"),]
trueCS <- trueCS[!duplicated(trueCS$class),]
#Add genome data to class dataframe
class.df$data[match(trueCS$class, class.df$class)] <- "CS"
class.df$data[match(truenoCS$class, class.df$class)] <- "noCS"
#Add class data to branches dataframe
branches["class_branch"] <- class.df$data[match(branches$branch_length, class.df$class)]

#Dataframes for vertical lines
boxlines <- data.frame(x=seq(0,100,5), .panel='Genome size (Mbp)')
boxlines.out <- data.frame(x=seq(0,6000,100), .panel='Genome size (Mbp)')
barlines <- data.frame(x=seq(0,475,25), .panel='Number of genomes')
#Dataframe for Ascomycota genome size mean lines for CS/no CS with and without outliers
no.CS.mean <- data.frame(x=mean(no.CS.df2$size), .panel='Genome size (Mbp)')
no.CS.mean.out <- data.frame(x=mean(no.CS.df2.out$size), .panel='Genome size (Mbp)')
CS.mean <- data.frame(x=mean(CS.df2$size), .panel='Genome size (Mbp)')
CS.mean.out <- data.frame(x=mean(CS.df2.out$size), .panel='Genome size (Mbp)')

#Add class data to genome size dataframe (NO OUTLIERS)
no.CS.df3 <- no.CS.df2
no.CS.df3["class"] <- tax.df2$Class[match(no.CS.df3$order,tax.df2$Order)]
#Create a list of dataframes with means for each class
means <- list()
for (i in sort(truenoCS$class)) {
  means[[i]] <- data.frame(x=(mean(no.CS.df3[no.CS.df3$class==i,]$size)), .panel='Genome size (Mbp)')
}

#Add class data to genome size dataframe (INCLUDING OUTLIERS)
no.CS.df3.out <- no.CS.df2.out
no.CS.df3.out["class"] <- tax.df2$Class[match(no.CS.df3.out$order,tax.df2$Order)]
#Create a list of dataframes with means for each class
means.out <- list()
for (i in sort(truenoCS$class)) {
  means.out[[i]] <- data.frame(x=(mean(no.CS.df3.out[no.CS.df3.out$class==i,]$size)), .panel='Genome size (Mbp)')
}

#Create dataframe for adding sample size to plot
counts <- tax.df3
#Remove orders without genome size data
counts <- counts[which(counts$data != "None"),]
#Add column with count of sample size, no CS/CS with and without outliers included
counts["noCScount"] <- table(unlist(no.CS.df2$order))[match(counts$order,names(table(unlist(no.CS.df2$order))))]
counts["noCScount.out"] <- table(unlist(no.CS.df2.out$order))[match(counts$order,names(table(unlist(no.CS.df2.out$order))))]
counts["CScount"] <- table(unlist(CS.df2$order))[match(counts$order,names(table(unlist(CS.df2$order))))]
counts["CScount.out"] <- table(unlist(CS.df2.out$order))[match(counts$order,names(table(unlist(CS.df2.out$order))))]
#Put brackets around count numbers
counts[,"noCScount"] <- paste0("(", unlist(counts[,"noCScount"]),")")
counts[,"noCScount.out"] <- paste0("(", unlist(counts[,"noCScount.out"]),")")
counts[,"CScount"] <- paste0("(", unlist(counts[,"CScount"]),")")
counts[,"CScount.out"] <- paste0("(", unlist(counts[,"CScount.out"]),")")

#Add column with the upper limit per order, with and without outliers included
counts["noCSmax"] <- NA
for (i in unique(no.CS.df2$order)) {
  counts[counts$order==i,]$noCSmax <- max(no.CS.df2[no.CS.df2$order==i,]$size)
}
counts["noCSmax.out"] <- NA
for (i in unique(no.CS.df2.out$order)) {
  counts[counts$order==i,]$noCSmax.out <- max(no.CS.df2.out[no.CS.df2.out$order==i,]$size)
}
counts["CSmax"] <- NA
for (i in unique(CS.df2$order)) {
  counts[counts$order==i,]$CSmax <- max(CS.df2[CS.df2$order==i,]$size)
}
counts["CSmax.out"] <- NA
for (i in unique(CS.df2.out$order)) {
  counts[counts$order==i,]$CSmax.out <- max(CS.df2.out[CS.df2.out$order==i,]$size)
}


##Number of genomes##

#Create data frame for number of genomes per order
num.df <- data.frame(order=tax.df3$order)

#Add column with number of genomes in mycocosm and ncbi
for (i in 1:length(num.df$order)) {
  num.df$num[i] <- length(grepl(num.df$order[i], myc$Order)[grepl(num.df$order[i], myc$Order) == TRUE]) + length(grepl(num.df$order[i], ncbi$Order)[grepl(num.df$order[i], ncbi$Order) == TRUE])
}

#Delete rows with no genomes
num.df <- num.df[num.df$num > 0,]


#Plot tree, bars and boxplot side by side for assembly based sizes (CS) (NO OUTLIERS)
tiff(file="CSgenomesplot.tiff", height=15, width=15, units="in", res=300)

gg <- ggtree(asc.tree)
gg1 <- gg %<+% tax.df3 +
       #geom_tippoint(aes(color=class)) +
       #geom_text(aes(label=node), hjust=0.3) +
       geom_tiplab(size=2.8,
                   #fill="white",
                   aes(color=data))
                   #geom="label",
                   #label.padding = unit(0.15, "lines"),
                   #label.size = 0)
#for (i in seq(means.out)) {
#  gg1 <- gg1+
#         geom_vline(data=means[[i]],
#                    color=hue_pal()(9)[i],
#                    linetype="dashed",
#                    aes(xintercept=x))
#}
gg2 <- gg1 +
       geom_vline(data=boxlines,
                  color="grey",
                  aes(xintercept=x)) +
       geom_vline(data=CS.mean,
                  color="black",
                  linetype="dashed",
                  aes(xintercept=x))
gg3 <- facet_plot(gg2+xlim_tree(2.9),
                  panel="Number of genomes",
                  data=num.df,
                  geom_barh,
                  aes(x=num, fill=class),
                  stat='identity')
#scale_x_continuous(breaks=c(seq(0,100,10))) +
#theme_classic() +
#theme(panel.border=element_blank()) +
#theme(strip.background = element_rect(colour="white")) +
#theme(strip.text.x = element_blank())
#theme(legend.position="none")
gg4 <- facet_plot(gg3,
                  panel="Genome size (Mbp)",
                  data=CS.df2,
                  geom_boxploth,
                  #outlier.shape=NA,
                  outlier.colour=NULL,
                  outlier.size=1,
                  aes(x=size, group=label, col=class)) +
       #geom_text(aes(x=label, label=count, size=2.8)) +
       scale_x_continuous(scale_x_continuous(breaks=pretty_breaks(10))) +
       theme_classic() +
       theme(panel.border=element_blank()) +
       #theme(strip.background = element_rect(colour="white")) +
       theme(strip.text.x = element_blank())
       #theme(legend.position="none")
gg4 <- facet_plot(gg4,
                  panel="Genome size (Mbp)",
                  data=counts,
                  geom_text,
                  size=2,
                  nudge_x=3.5,
                  aes(x=CSmax, label=CScount))
gg5 <- gg4 %<+% branches + 
       geom_label(aes(x=branch, label=branch_length, colour=class_branch),
                  fill="white",
                  size=2.8,
                  label.padding = unit(0.15, "lines"),
                  label.size = 0) +
       scale_colour_manual(values=c("black",hue_pal()(10)[1:4],"red","grey",hue_pal()(10)[5:11])) +
       scale_fill_manual(values=c("white","white",hue_pal()(10)[1:4],"white",hue_pal()(10)[5:6],"white",hue_pal()(10)[7:10],"white"))

#Remove unwanted axis on tree panel
plot_tab <- ggplotGrob(gg5)
elements <- plot_tab$layout$name

gtable_filter_remove <- function (x, name, trim = TRUE){
  matches <- !(x$layout$name %in% name)
  x$layout <- x$layout[matches, , drop = FALSE]
  x$grobs <- x$grobs[matches]
  if (trim) 
    x <- gtable_trim(x)
  x
}

p_filtered <- gtable_filter_remove(plot_tab, name = elements[6],
                                   trim = FALSE)
p_filtered <- gtable_filter_remove(p_filtered, name = elements[8],
                                   trim = FALSE)
p_filtered <- gtable_filter_remove(p_filtered, name = elements[16],
                                   trim = FALSE)

grid.newpage()
grid.draw(p_filtered)

dev.off()


#Plot all noCS data points INCLUDING OUTLIERS
tiff(file="genomesizeplot_outliers.tiff", height=15, width=10, units="in", res=300)

gg <- ggtree(asc.tree)
gg1 <- gg %<+% tax.df3 +
       geom_tiplab(size=2.8,
                   aes(color=data)) +
       scale_colour_manual(values=c("grey",as.character(classnodes$colour[4:5]),as.character(classnodes$colour[8:9]),"red","grey",as.character(classnodes$colour[12:13]),as.character(classnodes$colour[15]),as.character(classnodes$colour[17])))
gg2 <- gg1 +
       geom_vline(data=boxlines.out,
                  color="grey",
                  aes(xintercept=x)) +
       geom_vline(data=no.CS.mean.out,
                  color="black",
                  linetype="dashed",
                  aes(xintercept=x))
gg3 <- facet_plot(gg2+xlim_tree(2.9),
                  panel="Genome size (Mbp)",
                  data=no.CS.df2.out,
                  geom_point,
                  size=1,
                  mapping=aes(x=size, colour=class)) +
       #geom_text(aes(x=label, label=count, size=2.8)) +
       scale_x_continuous(breaks=pretty_breaks(10),
                          position="top",
                          sec.axis = dup_axis()) +
       scale_y_continuous(expand=expand_scale(mult=0.005)) +
       theme_classic() +
       theme(panel.border=element_blank()) +
       theme(legend.position="none") +
       theme(strip.placement="outside") +
       #theme(strip.text.x=element_blank()) +
       theme(strip.background=element_blank())
gg3 <- facet_plot(gg3,
                  panel="Genome size (Mbp)",
                  data=counts,
                  geom_text,
                  size=2,
                  nudge_x=200,
                  aes(x=noCSmax.out, label=noCScount.out))
gg4 <- gg3 %<+% branches + 
       geom_label(aes(x=branch, label=branch_length, colour=class_branch),
                  fill="white",
                  size=2.8,
                  label.padding = unit(0.15, "lines"),
                  label.size = 0)

#Remove unwanted axis on tree panel
plot_tab <- ggplotGrob(gg4)
elements <- plot_tab$layout$name


p_filtered <- gtable_filter_remove(plot_tab, name = elements[4],
                                   trim = FALSE)
p_filtered <- gtable_filter_remove(p_filtered, name = elements[6],
                                   trim = FALSE)
p_filtered <- gtable_filter_remove(p_filtered, name = elements[8],
                                   trim = FALSE)
p_filtered <- gtable_filter_remove(p_filtered, name = elements[10],
                                   trim = FALSE)

grid.newpage()
grid.draw(p_filtered)

dev.off()


#Plot tree against number of genomes and CS/noCS genome sizes 

tiff(file="allgenomesplot.tiff", height=15, width=10, units="in", res=300)

gg <- ggtree(asc.tree)
gg1 <- gg %<+% tax.df3 +
       geom_tiplab(size=2.8,
                   aes(color=data)) +
       scale_colour_manual(values=c("black","red","darkgrey"))
gg2 <- gg1 +
       geom_vline(data=boxlines,
                  color="lightgrey",
                  aes(xintercept=x)) +
       geom_vline(data=no.CS.mean,
                  color="black",
                  linetype="dashed",
                  aes(xintercept=x)) +
       geom_vline(data=CS.mean,
                  color="darkgrey",
                  linetype="dashed",
                  aes(xintercept=x)) +
       geom_vline(data=barlines,
                  color="lightgrey",
                  aes(xintercept=x))
gg3 <- facet_plot(gg2+xlim_tree(2.9),
                  panel="Number of genomes",
                  data=num.df,
                  geom_barh,
                  aes(x=num, fill=class),
                  stat='identity')
gg4 <- facet_plot(gg3,
                  panel="Genome size (Mbp)",
                  data=no.CS.df2,
                  geom_boxploth,
                  #outlier.shape=NA,
                  outlier.size=1,
                  aes(x=size, group=label, fill=class))
gg4 <- facet_plot(gg4,
                  panel="Genome size (Mbp)",
                  data=CS.df2,
                  geom_boxploth,
                  #outlier.shape=NA,
                  outlier.size=1,
                  colour="darkgrey",
                  linetype="dotted",
                  alpha=0.3,
                  aes(x=size, group=label, fill=class))
gg4 <- facet_plot(gg4,
                  panel="Genome size (Mbp)",
                  data=counts,
                  geom_text,
                  size=2,
                  nudge_x=5,
                  aes(x=noCSmax, label=noCScount)) +
                  scale_x_continuous(breaks=pretty_breaks(10),
                                     position="top",
                                     sec.axis = dup_axis()) +
                  scale_y_continuous(expand=expand_scale(mult=0.005)) +
                  theme_classic() +
                  theme(panel.border=element_blank()) +
                  theme(legend.position="none") +
                  theme(strip.placement="outside") +
                  #theme(strip.text.x=element_blank()) +
                  theme(strip.background=element_blank())
gg5 <- gg4 %<+% branches + 
       geom_label(aes(x=branch, label=branch_length, colour=class_branch),
                  fill="white",
                  size=2.8,
                  label.padding=unit(0.15, "lines"),
                  label.size=0)

#Make plot into table
plot_tab <- ggplotGrob(gg5)
#Change widths of bar and box panels
plot_tab$widths[7] <- 0.60*plot_tab$widths[7]
plot_tab$widths[9] <- 0.70*plot_tab$widths[9]

elements <- plot_tab$layout$name

gtable_filter_remove <- function (x, name, trim = TRUE){
  matches <- !(x$layout$name %in% name)
  x$layout <- x$layout[matches, , drop = FALSE]
  x$grobs <- x$grobs[matches]
  if (trim) 
    x <- gtable_trim(x)
  x
}

#Remove unwanted axis on tree panel
p_filtered <- gtable_filter_remove(plot_tab, name = elements[5],
                                   trim = FALSE)
p_filtered <- gtable_filter_remove(p_filtered, name = elements[8],
                                   trim = FALSE)
p_filtered <- gtable_filter_remove(p_filtered, name = elements[11],
                                   trim = FALSE)
p_filtered <- gtable_filter_remove(p_filtered, name = elements[13],
                                   trim = FALSE)

grid.newpage()
grid.draw(p_filtered)

dev.off()