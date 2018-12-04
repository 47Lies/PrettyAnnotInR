library("reshape2")
library("DOSE")
library("clusterProfiler")
library("PANTHER.db")
library("org.Hs.eg.db")
library("KEGG.db")
#library("KEGGREST")
library("reactome.db")
BUILD<-FALSE
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

	
	PantherKeys<-keys(PANTHER.db,keytype="PATHWAY_ID")
	PantherInfo<-select(PANTHER.db,PantherKeys, c("ENTREZ","SPECIES","PATHWAY_TERM"),"PATHWAY_ID")
	PantherInfo<-PantherInfo[PantherInfo$SPECIES=="HUMAN" & !is.na(PantherInfo$SPECIES),]
	PantherInfo$Gene_Symbol<-mapIds(org.Hs.eg.db,keys=as.character(PantherInfo$ENTREZ),
	                                column="SYMBOL",keytype="ENTREZID",multiVals="first")
	
	P2Genes<-PantherInfo[,c("PATHWAY_ID","Gene_Symbol")]
	PTerms<-PantherInfo[,c("PATHWAY_ID","PATHWAY_TERM")]  
	
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


AllTestGoKeggPanther<-function(GeneList,PopName) {
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
  
  COMMAND<-paste("ssconvert --export-type=Gnumeric_Excel:excel_biff8 --merge-to=",PopName,".Enrichement.xls ",sep="")
  DESCRIPTION<-vector()
  if(!is.null(Kegg)){
    write.table(Kegg[order(Kegg$Count,decreasing=TRUE),],
                paste(PopName,".Kegg.txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
    COMMAND<-paste(COMMAND,PopName,".Kegg.txt ",sep="")
    DESCRIPTION<-c(DESCRIPTION,as.vector(Kegg$Description))
  }
  if(!is.null(Panther)){
    write.table(Panther[order(Panther$Count,decreasing=TRUE),],
                paste(PopName,".Panther.txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
    COMMAND<-paste(COMMAND,PopName,".Panther.txt ",sep="")
    DESCRIPTION<-c(DESCRIPTION,as.vector(Panther$Description))
  }
  if(!is.null(Reactome)){
    write.table(Reactome[order(Reactome$Count,decreasing=TRUE),],
                paste(PopName,".Reactome.txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
    COMMAND<-paste(COMMAND,PopName,".Reactome.txt ",sep="")
    DESCRIPTION<-c(DESCRIPTION,gsub("human:","",as.vector(Reactome$Description)))
  }
  if(!is.null(BP)){
    write.table(BP[order(BP$Count,decreasing=TRUE),],
                paste(PopName,".BP.txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
    COMMAND<-paste(COMMAND,PopName,".BP.txt ",sep="")
    DESCRIPTION<-c(DESCRIPTION,as.vector(BP$Description))
  }
  if(!is.null(CC)){
    write.table(CC[order(CC$Count,decreasing=TRUE),],
                paste(PopName,".CC.txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
    COMMAND<-paste(COMMAND,PopName,".CC.txt ",sep="")
    DESCRIPTION<-c(DESCRIPTION,as.vector(CC$Description))
  }
  if(!is.null(MF)){
    write.table(MF[order(MF$Count,decreasing=TRUE),],
                paste(PopName,".MF.txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
    COMMAND<-paste(COMMAND,PopName,".MF.txt ",sep="")
    DESCRIPTION<-c(DESCRIPTION,as.vector(MF$Description))
  }
  DESCRIPTION<-gsub("regulation","",DESCRIPTION)
  DESCRIPTION<-gsub("process","",DESCRIPTION)
  # MArche pas 
  #system(command=COMMAND)
  #PlotWC(PopName,DESCRIPTION)0
  
}	
