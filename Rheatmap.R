
# set working directory to waffle

setwd("F:/ekeim/saliva")

#need to input data from my mapping file and my OTU file


mapping <- read.delim("data_raw/040816JP27F-mapping.txt", header= TRUE, check.names =FALSE, 
                      row.names= 1, colClasses = c("#SampleID"="character"))


otu_table <- read.delim("otus/otu_table.txt",
                        header= TRUE, skip=1, check.names = FALSE, row.names=1)

##data is now input
## Now need to relate the two dataframes together note that column names in otu_table are also sample names
#that match sampleIDs in the rows of mapping


otu_table2 <- subset(otu_table, select= -c(taxonomy))

samplesums<-as.data.frame(colSums(otu_table2, na.rm=FALSE, dims=1))

mapping2 <-merge(mapping, samplesums, by=0)

#contamination

otu_contam <- subset(otu_table, 294.11==0 & 294.12==0)

#no contamination/non represented otus in samples from donor

#OTU table manipulation where taxonomical information is delegated as an object and the rest of the data is a matrix
# part above for the Moveme fxn seems redundant


taxa.names <- otu_table$taxonomy
otu_table3 <- as.matrix(otu_table[,-dim(otu_table)[2]])


#####NOT WORKING#######

# sum for each column- treating it like a "new sample"

s_abundance <-apply(otu_table3,2,sum)

s_abundance2 <- colSums(otu_table3)


#exclusion of low abundance reads threshold set at 20
bads<-otu_table3[,s_abundance2<2000]
goods<-otu_table3[,s_abundance2>20000]

ncol(otu_table3)-ncol(goods)
ncol(bads)
##Skipping doesn't see like it is working---- all the sums are >20000
#Maybe this is the point where we import the standarized OTU table
#####################

#relative abundances and transposing

otu_table3 <- scale(otu_table3, center=F, scale=colSums(otu_table3))
otu_table3 <-t(otu_table3)

#numbering taxonomic levels starting at level 1

extractNameLevel <- function(x, level){
  a <- c(unlist(strsplit(x, ';')), 'Other')
  paste(a[1:min(level, length(a))], collapse = ';')
}

#function that breaks up taxonomy names
otu2taxonomy <-  function(x,level,taxa=NULL){
  if(is.null(taxa)){
    taxa <- colnames(x)
  }
  if(length(taxa)!=dim(x)[2]){
    print("ERROR: taxonomy should have the same length 
          as the number of columns in OTU table")
    return;
  }
  level.names <-sapply(as.character(taxa), 
                       function(x)
                         extractNameLevel(x,level=level))
  t(apply(x, 1, 
          function(y) 
            tapply(y,level.names,sum)))
}


d.phylum = otu2taxonomy(otu_table3,level=2,taxa=taxa.names)
d.class = otu2taxonomy(otu_table3,level=3,taxa=taxa.names)
d.order = otu2taxonomy(otu_table3,level=4,taxa=taxa.names)
d.family = otu2taxonomy(otu_table3,level=5,taxa=taxa.names)
d.genus = otu2taxonomy(otu_table3,level=6,taxa=taxa.names)
d.species = otu2taxonomy(otu_table3,level=7,taxa=taxa.names)

#transpose, and export levels
phylum2 <-t(d.phylum)
class2 <-t(d.class)
order2 <-t(d.order)
family2 <-t(d.family)
genus2 <-t(d.genus)
species2 <-t(d.species)
write.table(phylum2, file="Routput/phyla.csv", col.names=NA,row.names=TRUE, sep=",", quote=FALSE) 
write.table(class2, file="Routput/classes.csv", col.names=NA,row.names=TRUE, sep=",", quote=FALSE)
write.table(order2, file="Routput/orders.csv", col.names=NA,row.names=TRUE, sep=",", quote=FALSE)
write.table(family2, file="Routput/families.csv", col.names=NA,row.names=TRUE, sep=",", quote=FALSE)
write.table(genus2, file="Routput/genera.csv", col.names=NA,row.names=TRUE, sep=",", quote=FALSE)
write.table(species2, file="Routput/species.csv", col.names=NA,row.names=TRUE, sep=",", quote=FALSE)


#mergin mapping file and OTU tables by family level

mergefamily <- merge(mapping, d.family, by="row.names", all.x=FALSE)

#check yourself! makes sure columns makes sense and match
dim(mapping)
dim(d.family)
dim(mergefamily)
#all have 12 columns = 12 samples!!

write.table(mergefamily, file= "Routput/mergefamily.csv", col.names=NA,row.names=TRUE,sep = ",", quote= FALSE)

#Reordering by metadata variable alphabetically-- need to think of a better way to do this

colnames(mapping)

bysource <- mergefamily[order(mergefamily$source, mergefamily$Row.names),]

#making the matrix for heatmap
## slightly worried about splitting up the mapping data check this
(OTUcol1<-ncol(mapping)+2)
(OTUcol2<-ncol(bysource))
justOTU<-bysource[,OTUcol1:OTUcol2]
justOTU[1:5,1:5]
rownames(justOTU[1:10,])
rownames(justOTU)<-bysource$Row.names
rownames(justOTU[1:10,])
justOTU2<-as.matrix(t(justOTU))
justOTU2[1:5,1:5]

#making a heatmap
install.packages("gplots")
library(gplots)
install.packages("RColorBrewer")
library(RColorBrewer)

scaleblue<-colorRampPalette(colors=c("white","blue"))(100)

#Still need to mess around with this
heatmap.2(justOTU2, Rowv=TRUE, Colv=FALSE, scale="column", 
          trace="none", col=scaleblue, xlab="Sample",
          ylab="Family", margins=c(10,15))

#
##
#####
#######
#merge by Order

mergeorder <- merge(mapping, d.order, by="row.names", all.x=FALSE)

#check yourself! makes sure columns makes sense and match
dim(mapping)
dim(d.order)
dim(mergeorder)
#all have 12 columns = 12 samples!!

write.table(mergeorder, file= "Routput/mergerorder.csv", col.names=NA,row.names=TRUE,sep = ",", quote= FALSE)

#Reordering by metadata variable alphabetically-- need to think of a better way to do this

colnames(mapping)

bysource2 <- mergeorder[order(mergeorder$source, mergeorder$Row.names),]

#making the matrix for heatmap
## slightly worried about splitting up the mapping data check this
(OTUcol1<-ncol(mapping)+2)
(OTUcol2<-ncol(bysource2))
justOTU3<-bysource2[,OTUcol1:OTUcol2]
justOTU3[1:5,1:5]
rownames(justOTU3[1:10,])
rownames(justOTU3)<-bysource$Row.names
rownames(justOTU3[1:10,])
justOTU4<-as.matrix(t(justOTU))
justOTU4[1:5,1:5]