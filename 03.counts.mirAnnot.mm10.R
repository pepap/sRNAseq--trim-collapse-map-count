source("03.01.collapsedBam2counts.R")

# @pepap : load the mature miRNA annotation
load("mirAnnot.pepType.dt.rda",verbose=T)
mirAnnot.dt <- mirAnnot.pepType.dt

# @pepap : SUFFIX of the input BAM files
SUFF  <- ".se.Aligned.sortedByCoord.out.bam"

SSTOR <- "/path/to/mESC--Dicer-loop-mutants--20250411/BAM/"
SBAMS <- system( command=paste0("/bin/ls ",SSTOR,"*/*",SUFF),intern=T )
SCONS <- gsub(   SUFF,"",gsub( "^.*[/]","",SBAMS ) )

TSTOR <- "/path/to/mESC--Dicer-loop-mutants--TraPR--20250516/BAM/"
TBAMS <- system( command=paste0("/bin/ls ",TSTOR,"*/*",SUFF),intern=T )
TCONS <- gsub(   SUFF,"",gsub( "^.*[/]","",TBAMS ) )

aBAMS <- c( SBAMS,TCONS )
aCONS <- c( SCONS,TCONS )

MINOVRL=15
# @pepap : remove spliced reads
RMV.SPL=as.logical("TRUE")
# @pepap : remove soft-clipped reads
RMV.SCL=as.logical("FALSE")
# @pepap : read-length range
RLENRAN=c(19,32)
# @pepap : allowed number of edit distances
NMRANGE=c(0,100)

#  @pepap : FUN
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
 libs.dt                            <- rbind( libs.dt,data.table( libname=aCONS[i],libsize=tmp.cnt$NORM ) )

 assign( x=paste0(aCONS[i],".frac.dt"),value=tmp.cnt )
 objs.ls                            <- c( objs.ls,paste0(aCONS[i],".frac.dt") )

 rm( list=c("tmp.sum","tmp.cnt") )
 gc( verbose=T )

 cat( "\n",sep="" )

}

save( list=objs.ls,file="counts.mirAnnot.rda" )
save( list=c("frac.dt","libs.dt","mirAnnot.dt"),file="all_objs.rda" )
