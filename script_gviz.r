#! /usr/bin/env Rscript

# Copyright (C) 2017 IARC/WHO
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

############################## ARGUMENTS SECTION #############################
args <- commandArgs(TRUE)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args <- argsL;rm(argsL)

ref_genome=args$genome_release #name of the genome release (for the annotations)
bam_folder=args$bam_folder 
pos_file=args$pos_file #input file containing the position and samples to consider for the plot
ref=args$ref  #fasta_ref

if(is.null(args$sample_names)) {
  args$sample_names="FILE" #use bam file names as sample names
} else if(args$sample_names=="SAMPLE"){ #use sample names extracted from the bam files as sample names
  #Get samples names from BAM files
  if (file.exists("names.txt")==FALSE) {  #names.txt associates the bam file names with the sample names 
    system(paste0("./create_names_txt.sh ",bam_folder))
  }
  indiv=read.table("names.txt",stringsAsFactors=F,colClasses = "character")
}
sample_names=args$sample_names


###############################################################################

############################## FUNCTIONS SECTION ##############################
#create alignments tracks  
get_htTrack <- function(bams,bam_folder,plot_grtracks,sTrack,grtrack,chr,pos,paired=FALSE){  
  l=list() #store track alignment for each bam in a list
  for (i in 1:length(bams)){
    l=append(l,AlignmentsTrack(paste0(bam_folder,as.character(bams[i]),".bam"), isPaired = paired,stacking = "squish",alpha=0.95,chromosome=chr,cex.mismatch=0.5,name=as.character(bams[i]),cex.title=1.5))
  }
  if(plot_grtracks){
    ht <- HighlightTrack(trackList = c(l,sTrack, grtrack), start = c(pos), width =0,chromosome = chr) #highlight the position of the variant
    s=c(0.05,0.1,rep(0.72/length(l),length(l)),0.05,0.08)
  }else if(plot_grtracks==FALSE){
    ht <- HighlightTrack(trackList = c(l, sTrack), start = c(pos), width =0,chromosome = chr)
    s=c(0.05,0.1,rep(0.8/length(l),length(l)),0.05)
  }

  return(list(ht,s))
}

#plot alignments
plotGviz <- function(sTrack,ref_genome,txdb,annotation,chr,pos,bams,bam_folder,w=50,w_zoomout=1000,paired=FALSE,nb_toplot=5){
 
  chr_annotation=get_UCSC_associations(chr,ref_genome)
  if(chr_annotation==""){ #if the chromosome name is not recognized, the annotations will not be added
    only_alignment=TRUE
  }else{
    only_alignment=FALSE
  }
  
  if(only_alignment==FALSE){
    #Define tracks common to all samples (reference sequence, chromosome representation, genome annotation)
    gtrack <- GenomeAxisTrack() #genomic axis
    sTrack@sequence@ranges@NAMES <- sapply( 1:length(sTrack@sequence@ranges@NAMES), function(i) unlist(strsplit(sTrack@sequence@ranges@NAMES[i], " "))[1] )
    sTrack@chromosome <-chr
    ideoTrack <- IdeogramTrack(genome = unlist(strsplit(ref_genome,".",fixed=TRUE))[3], chromosome = chr_annotation) #chromosome representation
    levels(ideoTrack@bandTable$chrom)[which(levels(ideoTrack@bandTable$chrom)==chr_annotation)] <- chr
    ideoTrack@chromosome<-chr
    
    #annotion track
    grtrack <- GeneRegionTrack(txdb,chromosome = chr_annotation,start = pos-w_zoomout, end = pos+w_zoomout,exonAnnotation = "exon",shape = "arrow",showTitle=FALSE,alpha=0.95,cex=0.7)
    if( summary(is.na(grtrack@range@elementMetadata$gene))[2]>0 ){ #if no genes in the annotation it is not possible to use the collapseTranscripts option
      displayPars(grtrack)$collapseTranscripts <- "longest"
    }
    grtrack@range@seqinfo@seqnames<-chr
    levels(grtrack@range@seqnames)<-chr
    grtrack@chromosome<-chr
    
    #complete gene annotation
    if(length(unique(gene(grtrack)))>=6){
      sampling=sample(unique(grtrack@range@elementMetadata$gene),5)
      grtrack@range=grtrack@range[grtrack@range@elementMetadata$gene %in% sampling,]
      displayPars(grtrack)$cex <- 0.6
    }
    displayPars(grtrack) <- list(background.title = "white")
    #zoom out on the gene
    grtrack_zoomout=grtrack
    displayPars(grtrack_zoomout)$transcriptAnnotation <- "symbol"
    displayPars(grtrack_zoomout)$showExonId <- FALSE
    displayPars(grtrack_zoomout)$shape <- c("smallArrow","box")
    if(length(gene(grtrack_zoomout))!=0){
      #check if geneIDs are found in annotation (ENTREZID)
      if( length( which( unique(gene(grtrack_zoomout)) %in% keys(annotation,keytype="ENTREZID") == TRUE) )==length(unique(gene(grtrack_zoomout))) ){
        symbols <- unlist(mapIds(annotation, gene(grtrack_zoomout), "SYMBOL", "ENTREZID", multiVals = "first"))
        symbol(grtrack_zoomout) <- symbols[gene(grtrack_zoomout)]
      }
      plot_grtracks=TRUE
    }else{
      plot_grtracks=FALSE #no genes in this region
    }
    if(length(unique(gene(grtrack_zoomout)))>=6){ #to many genes to display
      grtrack_zoomout@range=grtrack_zoomout@range[grtrack_zoomout@range@elementMetadata$gene %in% sampling,]
    }
    ht_zoomout <- HighlightTrack(trackList = list(grtrack_zoomout,gtrack), start = c(pos), width =0,chromosome = chr) #highlight the position of the variant on the genomic axis and the genome annotation
    
  }else{
    plot_grtracks=FALSE
    sTrack@sequence@ranges@NAMES <- sapply( 1:length(sTrack@sequence@ranges@NAMES), function(i) unlist(strsplit(sTrack@sequence@ranges@NAMES[i], " "))[1] )
    sTrack@chromosome <-chr
  }

  if(only_alignment=="FALSE"){
    res=get_htTrack(bams,bam_folder,plot_grtracks,sTrack,grtrack,chr,pos,paired) #create alignments tracks
      if(plot_grtracks){
        grid.newpage()
        pushViewport(viewport(x=0,y=1, height=0.85, width=1, just=c("left","top")))
        plotTracks(c(ideoTrack,gtrack,res[1]),sizes=unlist(res[2]),from = pos-w, to = pos+w,add = TRUE, add53=TRUE,min.height=4, main=paste0(chr,":",pos),title.width=0.7,littleTicks = TRUE,cex.main=1.5)
        popViewport(1)
        pushViewport(viewport(x=0,y=0, height=0.15, width=1, just=c("left","bottom")))
        plotTracks(list(ht_zoomout),chromosome = chr, add = TRUE)
        popViewport(0)
      }else{
        plotTracks(c(ideoTrack,gtrack,res[1]),sizes=unlist(res[2]),from = pos-w, to = pos+w, add53=TRUE,min.height=4, main=paste0(chr,":",pos),title.width=0.7,littleTicks = TRUE,cex.main=1.5)
      }
  }else {
    res=get_htTrack(bams,bam_folder,plot_grtracks,sTrack,grtrack,chr,pos,paired) #create alignments tracks
    plotTracks(c(gtrack,res[1]),sizes=unlist(res[2])[-1],from = pos-w, to = pos+w, add53=TRUE,min.height=4, main=paste0(chr,":",pos),title.width=0.7,littleTicks = TRUE,cex.main=1.5)
  }

}

#associate sample names to bam files
get_bam_file_names<- function(bams,sample_names){
  if(sample_names=="SAMPLE"){
    return( sapply(bams, function(n) indiv[,1][which(indiv[,2]==n)]) )
  }
}

#get UCSC chromosome associations
get_UCSC_associations<- function(chr,ref_genome){
  convert_table=read.table(paste(paste0(unlist(strsplit(ref_genome,"[.]"))[3],"_chromosomeNames2UCSC.txt"),sep=""),head=FALSE,fill = TRUE)
  if(length(which(convert_table$V1==chr))==0){
    return("")
  }else{
    return(as.character(unique(convert_table[which(convert_table$V1==chr),2])))
  }
}

###############################################################################

############################## LIBRARIES ######################################

library("Gviz")
#load libraries for the annotations
library(paste0("TxDb.",ref_genome,".knownGene"),character.only=TRUE)
assign("txdb",get(paste0("TxDb.",ref_genome,".knownGene")))

annotation=paste0("org.",substr(ref_genome,1,2),".eg.db")
library(annotation,character.only=TRUE)
assign("annotation",get(annotation))

options(ucscChromosomeNames=FALSE) #in case some chromosome names are not recognized as UCSC chromosome names

#read fasta
library(Biostrings)
fasta_ref=readDNAStringSet(ref)
sTrack<-SequenceTrack(fasta_ref,cex=0.6)


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
  plotGviz(sTrack,ref_genome,txdb,annotation,chr,pos,bams,bam_folder)
  dev.off()

}
