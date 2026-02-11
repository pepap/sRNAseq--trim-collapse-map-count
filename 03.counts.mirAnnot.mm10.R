source("~/Scripts/R_func/collapsedBam2counts.R")

#load( "/storage/brno1-cerit/home/pepap/DicerX/00.stats/02.DE/mirAnnot.dt.rda",verbose=T )
load("/storage/brno12-cerit/home/pepap/brno1/DicerX/11.check_miRNAstrand/mirAnnot.pepType.dt.rda",verbose=T)
mirAnnot.dt <- mirAnnot.pepType.dt

SUFF  <- ".se.Aligned.sortedByCoord.out.bam"

#ESTOR <- "/storage/brno12-cerit/home/pepap/brno1/Valeria.Buccheri/02.mESC_mutants/BAM/"
#EBAMS <- system( command=paste0("/bin/ls ",ESTOR,"*/*",SUFF),intern=T )
#ECONS <- gsub(   SUFF,"",gsub( "^.*[/]","",EBAMS ) )
#
#FSTOR <- "/storage/brno12-cerit/home/pepap/brno1/Valeria.Buccheri/15.20231212_Buccheri/BAM/"
#FBAMS <- system( command=paste0("/bin/ls ",FSTOR,"*/*",SUFF),intern=T )
#FCONS <- gsub(   SUFF,"",gsub( "^.*[/]","",FBAMS ) )
#
#LSTOR <- "/storage/brno12-cerit/home/pepap/brno1/Valeria.Buccheri/16.mESC-DicerLoops/BAM/"
#LBAMS <- system( command=paste0("/bin/ls ",LSTOR,"*/*",SUFF),intern=T )
#LCONS <- gsub(   SUFF,"",gsub( "^.*[/]","",LBAMS ) )
#
#MSTOR <- "/storage/brno12-cerit/home/pepap/brno1/Valeria.Buccheri/17.240924_mESC_loops_mutants/BAM/"
#MBAMS <- system( command=paste0("/bin/ls ",MSTOR,"*/*",SUFF),intern=T )
#MCONS <- gsub(   SUFF,"",gsub( "^.*[/]","",MBAMS ) )

CSTOR <- "/storage/brno12-cerit/home/pepap/brno1/Valeria.Buccheri/22.mESC--Dicer-loop-mutants--TraPR--20250516/BAM/"
CBAMS <- system( command=paste0("/bin/ls ",CSTOR,"*/*",SUFF),intern=T )
CCONS <- gsub(   SUFF,"",gsub( "^.*[/]","",CBAMS ) )

#aBAMS <- c( EBAMS,FBAMS,LBAMS,MBAMS )
#aCONS <- c( ECONS,FCONS,LCONS,MCONS )
aBAMS <- c(CBAMS)
aCONS <- c(CCONS)

MINOVRL=15
#> remove spliced reads
RMV.SPL=as.logical("TRUE")
#> remove soft-clipped reads
RMV.SCL=as.logical("FALSE")
#> read-length range
RLENRAN=c(19,32)
#> allowed number of edit distances
NMRANGE=c(0,100)

#= FUN
dt2gr <- function(inDT,dt.seqnames="chr",dt.start="start",dt.end="end",dt.strand="strand",keep.metaData=T) {

 out.gr <-
  GRanges(
   seqnames = as.character( inDT[[dt.seqnames]] ),
   ranges   = IRanges(
    start   = as.integer( inDT[[dt.start]] ),
    end     = as.integer( inDT[[dt.end]] )
   ),
   strand   = as.character( inDT[[dt.strand]] )
  )

 if ( keep.metaData ) {
  mcols( out.gr ) <- inDT[ , colnames(inDT)[ !( colnames(inDT) %in% c(dt.seqnames,dt.start,dt.end,dt.strand) ) ] , with = F ]
 }

 return( out.gr )

}

#=

mirAnnot.gr <- dt2gr( inDT=mirAnnot.dt )

cat( "\n",sep="" )
frac.dt <- data.table( biotype=c(unique(mirAnnot.dt[["ID"]]),"not_annotated") )
libs.dt <- data.table()
objs.ls <- c()
for ( i in seq_along(aBAMS) ) {

 cat( " -> Loading sample : ",aCONS[i],"\n",sep="" )

 tmp.cnt     <-
  collapsedBam2counts(
   BAMFILE      = aBAMS[i],
   FRAC         = T,
   MINOVERLAP   = MINOVRL,
   IGNORE_STR   = F,
   SENSING      = NULL,
   RLENRANGE    = RLENRAN,
   NMISRANGE    = NMRANGE,
   NORMRANGE    = c(19,32),
   RMV.SPLICED  = RMV.SPL,
   RMV.SOFTCLP  = RMV.SCL,
   ANNOT        = mirAnnot.gr,
   ANNOT_BTPCOL = "ID",
   ALL.BTP      = unique(mirAnnot.dt[["ID"]])
  )
 tmp.sum     <- rowSums(as.matrix( x=tmp.cnt$EXP,rownames="biotype" ))
 tmp.cnt$EXP <- data.table( biotype=names(tmp.sum),count=as.vector(tmp.sum) )

 frac.dt                            <- merge( frac.dt,tmp.cnt[["EXP"]],by="biotype",all.x=T,sort=F )
 colnames(frac.dt)[ ncol(frac.dt) ] <- aCONS[i]
# frac.dt[[aCONS[i]]]             <- frac.dt[[aCONS[i]]]/(tmp.cnt$NORM*1e-06)
 libs.dt                            <- rbind( libs.dt,data.table( libname=aCONS[i],libsize=tmp.cnt$NORM ) )

 assign( x=paste0(aCONS[i],".frac.dt"),value=tmp.cnt )
 objs.ls                            <- c( objs.ls,paste0(aCONS[i],".frac.dt") )

 rm( list=c("tmp.sum","tmp.cnt") )
 gc( verbose=T )

 cat( "\n",sep="" )

}

save( list=objs.ls,file="counts.mirAnnot.rda" )
save( list=c("frac.dt","libs.dt","mirAnnot.dt"),file="all_objs.rda" )

