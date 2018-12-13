library("DOSE")
library("xlsx")
library("clusterProfiler")
library("org.Hs.eg.db")
library("STRINGdb")

########################################################################
#	Annotation script function
#		will use the data available in 
#			-org.HS.eg.db the databases with differents mappings ID for homo sapiens
#				-contains some genes with categories of interest
#				-contains id accross several database Official gene symbol/ENTREZ_ID/UNIPROT
#			-KEGG.db to get the names of the pathways
#			-Panther to get
#				-the mapping Entrez id/ pathway id
#				-the mapping Pathway name/ pathway id
#			-reactome to get
#				-the mapping Entrez id/ reactome id
#				-the mapping Pathway name/ pathway id
#		to build reference tables and will use clusterProfiler::enricher to do the statistics
#

Dir<-"./Ressources/Annotations/"
BUILD<-FALSE
if(!dir.exists(Dir)){
	dir.create(Dir,recursive = TRUE)
	BUILD<-TRUE
	library("KEGG.db")
	#library("KEGGREST")
	library("reactome.db")
	library("PANTHER.db")

}
if(BUILD==TRUE){
	ReactomeTerms<-toTable(reactome.db::reactomePATHID2NAME)
	ReactomeEntrez<-toTable(reactome.db::reactomeEXTID2PATHID)
	ReactomeEntrez$Gene_Symbol<-mapIds(org.Hs.eg.db,keys=as.character(ReactomeEntrez$gene_id),
	                                   column="SYMBOL",keytype="ENTREZID",multiVals="first")
	R2Genes<-ReactomeEntrez[,c("DB_ID","Gene_Symbol")]
	ReactomeTerms<-ReactomeTerms[ReactomeTerms$DB_ID %in% unique(ReactomeEntrez$DB_ID),]
	GoTerms<-toTable(GO.db::GOTERM)
	GoTerms.BP<-GoTerms[GoTerms$Ontology=="BP",]
	GoTerms.BP<-unique(GoTerms.BP[,c("go_id","Term")])
	GoTerms.MF<-GoTerms[GoTerms$Ontology=="MF",]
	GoTerms.MF<-unique(GoTerms.MF[,c("go_id","Term")])
	GoTerms.CC<-GoTerms[GoTerms$Ontology=="CC",]
	GoTerms.CC<-unique(GoTerms.CC[,c("go_id","Term")])
	
	Genes2Kegg<-toTable(org.Hs.egPATH)
	KeggTerms<-toTable(KEGG.db::KEGGPATHID2NAME)
	Genes2Kegg$Gene_Symbol<-mapIds(org.Hs.eg.db,keys=as.character(Genes2Kegg$gene_id),
                    column="SYMBOL",keytype="ENTREZID",multiVals="first")
	Kegg2Genes<-Genes2Kegg[,c("path_id","Gene_Symbol")]

	
	PantherInfo<-select(x=PANTHER.db,keys="HUMAN",columns=c("ENTREZ","PATHWAY_TERM","PATHWAY_ID"), keytype="SPECIES")
	PantherInfo<-PantherInfo[!is.na(PantherInfo$PATHWAY_ID),]
	PantherInfo$Gene_Symbol<-as.vector(mapIds(org.Hs.eg.db,keys=as.character(PantherInfo$ENTREZ),
	                                column="SYMBOL",keytype="ENTREZID",multiVals="first"))
	
	P2Genes<-unique(PantherInfo[,c("PATHWAY_ID","Gene_Symbol")])
	P2Genes<- apply(P2Genes,2,as.character)
	PTerms<-unique(PantherInfo[,c("PATHWAY_ID","PATHWAY_TERM")])
	PTerms<- apply(PTerms,2,as.character)
	Genes2GO<-toTable(org.Hs.egGO2ALLEGS)
	Genes2GO$Gene_Symbol<-mapIds(org.Hs.eg.db,keys=as.character(Genes2GO$gene_id),
	                             column="SYMBOL",keytype="ENTREZID",multiVals="first")
	BP2Genes<-Genes2GO[Genes2GO$Ontology=="BP",c("go_id","Gene_Symbol")]
	MF2Genes<-Genes2GO[Genes2GO$Ontology=="MF",c("go_id","Gene_Symbol")]
	CC2Genes<-Genes2GO[Genes2GO$Ontology=="CC",c("go_id","Gene_Symbol")]
	
	write.table(BP2Genes,paste(Dir,"GO.BP.GeneList.txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
	write.table(CC2Genes,paste(Dir,"GO.CC.GeneList.txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
	write.table(MF2Genes,paste(Dir,"GO.MF.GeneList.txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
	write.table(Kegg2Genes,paste(Dir,"Kegg.GeneList.txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
	write.table(P2Genes,paste(Dir,"Panther.GeneList.txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
	write.table(R2Genes,paste(Dir,"Reactome.DB.GeneList.txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
	
	write.table(GoTerms.BP,paste(Dir,"GO.BP.Names.txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
	write.table(GoTerms.CC,paste(Dir,"GO.CC.Names.txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
	write.table(GoTerms.MF,paste(Dir,"GO.MF.Names.txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
	write.table(KeggTerms,paste(Dir,"Kegg.Names.txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
	write.table(PTerms,paste(Dir,"Panther.Names.txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
	write.table(ReactomeTerms,paste(Dir,"Reactome.DB.Names.txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
}

BP2Genes<-read.table(paste(Dir,"GO.BP.GeneList.txt",sep=""),header=TRUE,sep="\t")
CC2Genes<-read.table(paste(Dir,"GO.CC.GeneList.txt",sep=""),header=TRUE,sep="\t")
MF2Genes<-read.table(paste(Dir,"GO.MF.GeneList.txt",sep=""),header=TRUE,sep="\t")
Kegg2Genes<-read.table(paste(Dir,"Kegg.GeneList.txt",sep=""),header=TRUE,sep="\t")
P2Genes<-read.table(paste(Dir,"Panther.GeneList.txt",sep=""),header=TRUE,sep="\t")
R2Genes<-read.table(paste(Dir,"Reactome.DB.GeneList.txt",sep=""),header=TRUE,sep="\t")

GoTerms.BP<-read.table(paste(Dir,"GO.BP.Names.txt",sep=""),header=TRUE,sep="\t",quote="\"")
GoTerms.CC<-read.table(paste(Dir,"GO.CC.Names.txt",sep=""),header=TRUE,sep="\t",quote="\"")
GoTerms.MF<-read.table(paste(Dir,"GO.MF.Names.txt",sep=""),header=TRUE,sep="\t",quote="\"")
KeggTerms<-read.table(paste(Dir,"Kegg.Names.txt",sep=""),header=TRUE,sep="\t",quote="\"")
PTerms<-read.table(paste(Dir,"Panther.Names.txt",sep=""),header=TRUE,sep="\t",quote="\"")
ReactomeTerms<-read.table(paste(Dir,"Reactome.DB.Names.txt",sep=""),header=TRUE,sep="\t",quote="\"")


AllTestGoKeggPanther<-function(GeneList,PopName,outputXLSX=FALSE) {
	if(outputXLSX){if(file.exists(paste(PopName,".xlsx",sep=""))){file.remove(paste(PopName,".xlsx",sep=""))}}
	BP<-enricher(GeneList,TERM2GENE=BP2Genes,TERM2NAME=GoTerms.BP,
		minGSSize=1,pAdjustMethod="fdr")
	MF<-enricher(GeneList,TERM2GENE=MF2Genes,TERM2NAME=GoTerms.MF,
		minGSSize=1,pAdjustMethod="fdr")
	CC<-enricher(GeneList,TERM2GENE=CC2Genes,TERM2NAME=GoTerms.CC,
		minGSSize=1,pAdjustMethod="fdr")
	Kegg<-enricher(GeneList,TERM2GENE=Kegg2Genes,TERM2NAME=KeggTerms,
		minGSSize=1,pAdjustMethod="fdr")
	Panther<-enricher(GeneList,TERM2GENE=P2Genes,TERM2NAME=PTerms,
		minGSSize=1,pAdjustMethod="fdr")
	Reactome<-enricher(GeneList,TERM2GENE=R2Genes,TERM2NAME=ReactomeTerms,
		minGSSize=1,pAdjustMethod="fdr")
	DESCRIPTION<-vector()
	if(dim(data.frame(BP))[1]>0){
		write.table(BP[order(BP$Count,decreasing=TRUE),],
		paste(PopName,".BP.txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
		if(outputXLSX){
			xlsx::write.xlsx(BP[order(BP$Count,decreasing=TRUE),],
				paste(PopName,".xlsx",sep=""), sheetName="Biological Process", 
				col.names=TRUE, row.names=FALSE, append=TRUE)
		}
		DESCRIPTION<-c(DESCRIPTION,as.vector(BP$Description))
	}
	if(dim(data.frame(CC))[1]>0){
		write.table(CC[order(CC$Count,decreasing=TRUE),],
			paste(PopName,".CC.txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
		if(outputXLSX){
			xlsx::write.xlsx(CC[order(CC$Count,decreasing=TRUE),],
				paste(PopName,".xlsx",sep=""), sheetName="Cellular Component", 
				col.names=TRUE, row.names=FALSE, append=TRUE)
		}
		DESCRIPTION<-c(DESCRIPTION,as.vector(CC$Description))
	}
	if(dim(data.frame(MF))[1]>0){
		write.table(MF[order(MF$Count,decreasing=TRUE),],
			paste(PopName,".MF.txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
		if(outputXLSX){
			xlsx::write.xlsx(MF[order(MF$Count,decreasing=TRUE),],
				paste(PopName,".xlsx",sep=""), sheetName="Molecular Function", 
				col.names=TRUE, row.names=FALSE, append=TRUE)
		}
		DESCRIPTION<-c(DESCRIPTION,as.vector(MF$Description))
	}
	if(dim(data.frame(Kegg))[1]>0){
		write.table(Kegg[order(Kegg$Count,decreasing=TRUE),],
			paste(PopName,".Kegg.txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
		if(outputXLSX){
			xlsx::write.xlsx(x=Kegg[order(Kegg$Count,decreasing=TRUE),],
				file=paste(PopName,".xlsx",sep=""), sheetName="KEGG", 
				col.names=TRUE, row.names=FALSE, append=TRUE)
		}
		DESCRIPTION<-c(DESCRIPTION,as.vector(Kegg$Description))
	}
	if(dim(data.frame(Panther))[1]>0){
		write.table(Panther[order(Panther$Count,decreasing=TRUE),],
			paste(PopName,".Panther.txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
		if(outputXLSX){
			xlsx::write.xlsx(Panther[order(Panther$Count,decreasing=TRUE),],
				paste(PopName,".xlsx",sep=""), sheetName="Panther", 
				col.names=TRUE, row.names=FALSE, append=TRUE)
		}
		DESCRIPTION<-c(DESCRIPTION,as.vector(Panther$Description))
	}
	if(dim(data.frame(Reactome))[1]>0){
		write.table(Reactome[order(Reactome$Count,decreasing=TRUE),],
			paste(PopName,".Reactome.txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
		if(outputXLSX){
			xlsx::write.xlsx(Reactome[order(Reactome$Count,decreasing=TRUE),],
				paste(PopName,".xlsx",sep=""), sheetName="Reactome", 
				col.names=TRUE, row.names=FALSE, append=TRUE)
		}
		DESCRIPTION<-c(DESCRIPTION,gsub("human:","",as.vector(Reactome$Description)))
	}
	DESCRIPTION<-gsub("regulation","",DESCRIPTION)
	DESCRIPTION<-gsub("process","",DESCRIPTION)
	# MArche pas 
	#system(command=COMMAND)
	#PlotWC(PopName,DESCRIPTION)0
}	

OneTestGoKeggPanther<-function(GeneList,
Ressource=c("BP","MF","CC","KEGG","reactome","Panther")){
	if(Ressource=="BP"){
		Enrich<-enricher(GeneList,TERM2GENE=BP2Genes,TERM2NAME=GoTerms.BP,
			minGSSize=1,pAdjustMethod="fdr")
	}else if(Ressource=="MF"){
		Enrich<-enricher(GeneList,TERM2GENE=MF2Genes,TERM2NAME=GoTerms.MF,
			minGSSize=1,pAdjustMethod="fdr")
	}else if(Ressource=="CC"){
		Enrich<-enricher(GeneList,TERM2GENE=CC2Genes,TERM2NAME=GoTerms.CC,
			minGSSize=1,pAdjustMethod="fdr")
	}else if(Ressource=="KEGG"){
		Enrich<-enricher(GeneList,TERM2GENE=Kegg2Genes,TERM2NAME=KeggTerms,
			minGSSize=1,pAdjustMethod="fdr")
	}else if(Ressource=="Panther"){
		Enrich<-enricher(GeneList,TERM2GENE=P2Genes,TERM2NAME=PTerms,
			minGSSize=1,pAdjustMethod="fdr")
	}else if(Ressource=="reactome"){
		Enrich<-enricher(GeneList,TERM2GENE=R2Genes,TERM2NAME=ReactomeTerms,
			minGSSize=1,pAdjustMethod="fdr")
	}
	return(Enrich)
}

PlotString<-function(GeneList){
	stringDB<-STRINGdb$new(version="10",species=9606,score_threshold=0, input_directory="" )
	Mapped <- stringDB$map( data.frame("gene"=GeneList),
	"gene", removeUnmappedRows = TRUE )
	if(dim(Mapped)[1]>400){
		Mapped<-Mapped[1:400,]
		Title<-"400 Firsts"
		stringDB$plot_network(Mapped$STRING_id)
	}
	stringDB$plot_network(Mapped$STRING_id,add_link=FALSE,add_summary=FALSE)
}

