############################## ARGUMENTS SECTION #############################
args <- commandArgs(TRUE)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args <- argsL;rm(argsL)

ref_genome=args$ref_genome
bam_folder=args$bam_folder
pos_file=args$pos_file #input file

if(is.null(args$sample_names)) {
  args$sample_names="FILE" #use bam file names as sample names
} else if(args$sample_names=="SAMPLE"){ #use sample names extracted from the bam files as sample names
  #Get samples names from BAM files
  if (file.exists("names.txt")==FALSE) {
    system(paste0("./create_names_txt.sh ",bam_folder))
  }
  indiv=read.table("names.txt",stringsAsFactors=F,colClasses = "character")
}
sample_names=args$sample_names


###############################################################################

############################## FUNCTIONS SECTION ##############################
#create alignments tracks
get_htTrack <- function(bams,bam_folder,plot_grtracks,UCSC,sTrack,grtrack,chr,pos,paired=FALSE){
  l=list() #store track alignment for each bam in a list
  for (i in 1:length(bams)){
    l=append(l,AlignmentsTrack(paste0(bam_folder,as.character(bams[i]),".bam"), isPaired = paired,stacking = "squish",alpha=0.95,chromosome=chr,cex.mismatch=0.5,name=as.character(bams[i]),cex.title=1.5))
  }
  if(UCSC & plot_grtracks){
    ht <- HighlightTrack(trackList = c(l,sTrack, grtrack), start = c(pos), width =0,chromosome = chr) #highlight the position of the variant
    s=c(0.05,0.1,rep(0.72/length(l),length(l)),0.05,0.08)
  }else if(UCSC==FALSE | plot_grtracks==FALSE){
    ht <- HighlightTrack(trackList = c(l, sTrack), start = c(pos), width =0,chromosome = chr)
    s=c(0.05,0.1,rep(0.8/length(l),length(l)),0.05)
  }
  
  return(list(ht,s))
}

#plot alignments
plotGviz <- function(sTrack,ref_genome,txdb,annotation,UCSC,chr,pos,bams,bam_folder,w=50,w_zoomout=1000,paired=FALSE,nb_toplot=5){
  #Define tracks common to all samples (reference sequence, chromosome representation, genome annotation)
  gtrack <- GenomeAxisTrack()
  if(UCSC){
    sTrack@chromosome <- chr
    ideoTrack <- IdeogramTrack(genome = unlist(strsplit(ref_genome,".",fixed=TRUE))[3], chromosome = chr) #chromosome representation
    #genome annotation :
    grtrack <- GeneRegionTrack(txdb,chromosome = chr,start = pos-w, end = pos-w,exonAnnotation = "exon",collapseTranscripts = "longest",shape = "arrow",showTitle=FALSE,alpha=0.95)
    displayPars(grtrack) <- list(background.title = "white")
    grtrack_zoomout <- GeneRegionTrack(txdb,chromosome = chr,start = pos-w_zoomout, end = pos+w_zoomout,transcriptAnnotation = "symbol",collapseTranscripts = "longest",alpha=0.95,showTitle=FALSE)
    if(length(gene(grtrack_zoomout))!=0){
      #check if geneIDs are found in annotation (ENTREZID)
      if( length( which( unique(gene(grtrack_zoomout)) %in% keys(annotation,keytype="ENTREZID") == TRUE))==length(unique(gene(grtrack_zoomout))) ){
          symbols <- unlist(mapIds(annotation, gene(grtrack_zoomout), "SYMBOL", "ENTREZID", multiVals = "first"))
          symbol(grtrack_zoomout) <- symbols[gene(grtrack_zoomout)]
      }
      plot_grtracks=TRUE
    }else{
      plot_grtracks=FALSE
    }
    ht_zoomout <- HighlightTrack(trackList = list(grtrack_zoomout,gtrack), start = c(pos), width =0,chromosome = chr)
  }else{ #genome annotation can not be added if non UCSC genome
    sTrack@chromosome <- chr
    ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = paste0("chr",chr)) #chromosome representation
    levels(ideoTrack@bandTable$chrom) <- sub("^chr", "", levels(ideoTrack@bandTable$chrom), ignore.case=T)
    ideoTrack@chromosome<-chr
    plot_grtracks=FALSE
  }

   
  if(UCSC & plot_grtracks){
    res=get_htTrack(bams,bam_folder,plot_grtracks,UCSC,sTrack,grtrack,chr,pos,paired) #create alignments tracks
    grid.newpage()
    pushViewport(viewport(x=0,y=1, height=0.85, width=1, just=c("left","top"))) 
    plotTracks(c(ideoTrack,gtrack,res[1]),sizes=unlist(res[2]),from = pos-w, to = pos+w,add = TRUE, add53=TRUE,min.height=4, main=paste0(chr,":",pos),title.width=0.7,littleTicks = TRUE,cex.main=1.5)
    popViewport(1)
    pushViewport(viewport(x=0,y=0, height=0.15, width=1, just=c("left","bottom")))
    plotTracks(list(ht_zoomout),chromosome = chr, add = TRUE) 
    popViewport(0)
  }else if(UCSC==FALSE | plot_grtracks==FALSE){
    res=get_htTrack(bams,bam_folder,plot_grtracks,UCSC,sTrack,grtrack,chr,pos,paired) #create alignments tracks
    plotTracks(c(ideoTrack,gtrack,res[1]),sizes=unlist(res[2]),from = pos-w, to = pos+w, add53=TRUE,min.height=4, main=paste0(chr,":",pos),title.width=0.7,littleTicks = TRUE,cex.main=1.5)
  }

  
}

#associate sample names to bam files
get_bam_file_names<- function(bams,sample_names){
  if(sample_names=="SAMPLE"){
    return( sapply(bams, function(n) indiv[,1][which(indiv[,2]==n)]) )
  }
}

###############################################################################

############################## LIBRARIES ######################################
library("Gviz")
ref_string=paste0("BSgenome.",ref_genome)
library(ref_string,character.only=TRUE)

assign("g",get(ref_string))

if(ref_genome!="Hsapiens.1000genomes.hs37d5"){
  
  #genome annotation
  library(paste0("TxDb.",ref_genome,".knownGene"),character.only=TRUE)
  assign("txdb",get(paste0("TxDb.",ref_genome,".knownGene")))
  #define SequenceTrack (reference genome)
  sTrack <- SequenceTrack(g,cex = 0.6)
  
  annotation=paste0("org.",substr(ref_genome,1,2),".eg.db")
  library(annotation,character.only=TRUE)
  assign("annotation",get(annotation))
  UCSC=TRUE #UCSC reference genome
  
}else if(ref_genome=="Hsapiens.1000genomes.hs37d5"){
  UCSC=FALSE #non UCSC reference genome
  options(ucscChromosomeNames=FALSE)
  
  #define SequenceTrack (reference genome)
  sTrack <- SequenceTrack(g,cex = 0.6)
  
}

###############################################################################
  
#read input file
graphs = read.table(pos_file,stringsAsFactors=F,colClasses = "character",header=F,fill=TRUE)

#for each row of the input file a pdf is generated 
for(i in 1:nrow(graphs)){
  chr=graphs[i,1]
  pos=as.numeric(graphs[i,2])
  #get bams to include in the alignment plot
  if(sample_names=="FILE"){ 
    bams=graphs[i,][2+which(graphs[i,3:length(graphs[i,])]!="")]
  }else if(sample_names=="SAMPLE"){
    bams=get_bam_file_names(graphs[i,][2+which(graphs[i,3:length(graphs[i,])]!="")], sample_names)
  }

  pdf(paste0(chr,":",pos,"_",paste(bams, collapse = '_'),".pdf"),11,12)
  plotGviz(sTrack,ref_genome,txdb,annotation,UCSC,chr,pos,bams,bam_folder)
  dev.off()
  
}
  
