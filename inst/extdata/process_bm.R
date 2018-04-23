library(readr)
library(tibble)
library(phyloseq)

data_rav<- read_csv("ravel_data.csv")
# Concatenating "ID_" to sample id 
for (i in 1:length(data_rav$`Sample ID`)){
    data_rav[i,1]<-paste("ID", data_rav[i,1],sep = "_")
}

# 1. Creating meta data and abundance table 
meta<- data_rav[,1:9]
col_names<- c("Sample_ID","Time_in_study","Subject_ID","Race","Age","Nugent Score","Nugent Category","Community_State_Type", "Total Read Countsd")
colnames(meta)<-col_names
meta$Race <- gsub("0", "Black", meta$Race)
meta$Race <- gsub("1", "White", meta$Race)
meta$Race<- gsub("4", "Others", meta$Race)
meta$Race <- gsub("5", "Hispanic", meta$Race)
meta$Race <- factor(meta$Race)
meta$Subject_ID <- factor(meta$Subject_ID)
row.names(meta) <- data_rav$`Sample ID`
meta$Sample_ID<-row.names(meta)

abund <- data_rav[,c(11:340)]
rownames(abund)<-meta$Sample_ID
abund<-as.data.frame(t(abund))
colnames(abund)<-data_rav$`Sample ID`
abund<-add_column(abund, organism = row.names(abund), .before = abund$ID_400_010106)

# 2. Creating and processing organisms to make OTU table
df_org<-as.data.frame(abund$organism)
colnames(df_org)<-"genus"
df_org$genus<-as.character(df_org$genus)
y<-strsplit(df_org$genus,"[.]")
df_org<-add_column(df_org,  species= "", .after = "genus")
df_org$genus<- sapply(y,"[[",1) 
for  (i in 1:330){
  if (length(y[[i]]==2)){
    df_org$species[i]<-y[[i]][2]
  }
  else{
    df_org$species[i]<-"s__"
  }
}
df_org$species<- gsub(" ", "", df_org$species)
df_org$genus[df_org$genus=="L"]<-"Lactobacillus"


# 3. Creating OTU table 
otu <- read.table("Greengenes_16S_all_2011-1.txt",sep = ";", header = FALSE)
otu[sapply(otu, is.factor)] <- lapply(otu[sapply(otu, is.factor)],as.character)
y<-strsplit(otu$V1,"\tk__")
otu<- add_column(otu, id = sapply(y,"[[",1),.before = "V1")
otu$V1<-sapply(y,"[[",2)

col_nam<-c("OTU","Domain","Phylum","Class","Order","Family","Genus","Species")
colnames(otu)<-col_nam
otu$Phylum<- gsub("p__", "", otu$Phylum)
otu$Class<- gsub("c__", "", otu$Class)
otu$Order<- gsub("o__", "", otu$Order)
otu$Family<- gsub("f__", "", otu$Family)
otu$Genus<- gsub("g__", "", otu$Genus)
otu$Species<- gsub("s__", "", otu$Species)
otu$Phylum<- gsub(" ", "", otu$Phylum)
otu$Class<- gsub(" ", "", otu$Class)
otu$Order<- gsub(" ", "", otu$Order)
otu$Family<- gsub(" ", "", otu$Family)
otu$Genus<- gsub(" ", "", otu$Genus)
otu$Species<- gsub(" ", "", otu$Species)


final_otu<-data.frame()
for (i in 1:330){
  c<- subset(otu, otu$Genus == df_org$genus[i] | otu$Family == df_org$genus[i] | otu$Order == df_org$genus[i] | otu$Class == df_org$genus[i]|otu$Phylum == df_org$genus[i]|otu$Domain == df_org$genus[i])
  final_otu<- rbind(final_otu,c[1,2:7])
}
rownames(final_otu)<-1:330

final_otu<-add_column(final_otu, Species = df_org$species, .after = "Genus")
final_otu$Species[is.na(final_otu$Species)]<-""
final_otu<-add_column(final_otu, OTU = row.names(final_otu), .before = "Domain")
final_otu$OTU<- paste("OTU",final_otu$OTU,sep = "_")
row.names(final_otu)<-final_otu$OTU
final_otu[sapply(final_otu, is.character)] <- lapply(final_otu[sapply(final_otu, is.character)],as.factor)


# 4. Creating Phyloseq object  
row.names(abund)<-final_otu$OTU
abund<-abund[,2:938]
seq.data<- as.matrix(abund+1)

seqa_otu = otu_table(seq.data, taxa_are_rows = TRUE)
metaa_sample = sample_data(meta)
taxaa_tax = tax_table(as.matrix(final_otu))
ravel_2012<- phyloseq(seqa_otu, metaa_sample, taxaa_tax)
save(ravel_2012, file = "ravel_2012.rda")
rm(data_rav,df_org,final_otu,meta,metaa_sample,otu,seq.data,y,c,abund,col_nam,col_names,i,seqa_otu,taxaa_tax)

