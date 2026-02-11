library(data.table)
library(GenomicRanges)
library(GenomicAlignments)

cat("\n @ pepap-functions loaded : \"gr_switchStrand\", \"grUnifySeqnames\", \"collapsedBam2counts\"\n\n",sep="")

# @pepap : (1)
gr_switchStrand <- function(inpGR) {

 mcols(inpGR)["tmp_pepap_idx"] <- seq_along(inpGR)

 tmp_fwd <- inpGR[ as.character(strand(inpGR)) == "+" ]
 tmp_rev <- inpGR[ as.character(strand(inpGR)) == "-" ]

 strand(tmp_fwd) <- as.character("-")
 strand(tmp_rev) <- as.character("+")

 outGR   <-
  append(
   x      = tmp_fwd,
   values = tmp_rev
  )

 outGR   <- outGR[ order( mcols(outGR)[["tmp_pepap_idx"]] ) ]
 mcols(outGR)["tmp_pepap_idx"] <- NULL

 return( outGR )

}

# @pepap : (2)
grUnifySeqnames <- function( tmp.gr ) {

 cat(" --> erasing \"chr?_\" prefixes & \"_random\" / \"[.]?\" suffixes from seqnames\n",sep="")
 tmp.seqlevels     <- seqlevels(tmp.gr)
 tmp.seqlevels[ tmp.seqlevels=="LCMV.L_segment.NC_004291.1" ] <- "LCMV_L"
 tmp.seqlevels[ tmp.seqlevels=="LCMV.S_segment.NC_004294.1" ] <- "LCMV_S"
 seqlevels(tmp.gr) <- sub( "[.].*$","",sub( "_random$","",sub( "^chr[0-9,a-z,A-Z]*_","",tmp.seqlevels ) ) )

 return( tmp.gr )

}

# @pepap : (3)
collapsedBam2counts <-
 function(
  BAMFILE,
  XPARAM       = ScanBamParam( tag=c("NH","HI","nM"),what=c("qname") ),
  FRAC         = T,
  MINOVERLAP   = 15,
  IGNORE_STR   = T,
  SENSING      = NULL,
  RLENRANGE    = c(1,100),
  NMISRANGE    = c(0,100),
  NORMRANGE    = c(19,32),
  RMV.SPLICED  = T,
  RMV.SOFTCLP  = F,
  ANNOT,
  ANNOT_BTPCOL,
  ALL.BTP
 ) {

 if ( !is.null(SENSING) ) {
 if ( SENSING!="sense" & SENSING!="anti" ) {
  cat("\n !!! \"SENSING\"=",SENSING," !!!\n",sep="")
  cat(" !!! Acceptable values for \"SENSING\" is : \"NULL\", \"sense\", \"anti\"  !!!\n\n"); stop()
 }
 }

 cat( " ++ Reading BAM : ",BAMFILE," ++\n",sep="" )
 tmp.tga  <- readGAlignments( file=BAMFILE,param=XPARAM )
 cat( " ++ Counting reads in width range   : <",NORMRANGE[1],";",NORMRANGE[2],">\n","                   in \"nM\"  range   : <",NMISRANGE[1],";",NMISRANGE[2],"> for normalization ++\n",sep="" )
 tmp.NORM <-
  sum(as.numeric(
   gsub(
    pattern     = "^[0-9]*[-]",
    replacement = "",
    x           = mcols(tmp.tga[
                   ( mcols(tmp.tga)[["HI"]] == 1            ) &
                   ( width(tmp.tga)         >= NORMRANGE[1] ) &
                   ( width(tmp.tga)         <= NORMRANGE[2] ) &
                   ( mcols(tmp.tga)[["nM"]] >= NMISRANGE[1] ) &
                   ( mcols(tmp.tga)[["nM"]] <= NMISRANGE[2] )
                  ])[["qname"]]
   )
  ))
 if ( RMV.SPLICED ) {
 cat( " ++ Removing spliced reads ++\n",sep="" )
 print( sum( as.numeric(gsub( pattern="^[0-9]*[-]",replacement="",x=mcols(tmp.tga)[["qname"]] ))/mcols(tmp.tga)[["NH"]] ) )
 tmp.tga  <- tmp.tga[ njunc(tmp.tga)==0 ]
 print( sum( as.numeric(gsub( pattern="^[0-9]*[-]",replacement="",x=mcols(tmp.tga)[["qname"]] ))/mcols(tmp.tga)[["NH"]] ) )
 }
 if ( RMV.SOFTCLP ) {
 cat( " ++ Removing soft-clipped reads ++\n",sep="" )
 print( sum( as.numeric(gsub( pattern="^[0-9]*[-]",replacement="",x=mcols(tmp.tga)[["qname"]] ))/mcols(tmp.tga)[["NH"]] ) )
 tmp.tga  <- tmp.tga[ width(tmp.tga)==qwidth(tmp.tga) ]
 print( sum( as.numeric(gsub( pattern="^[0-9]*[-]",replacement="",x=mcols(tmp.tga)[["qname"]] ))/mcols(tmp.tga)[["NH"]] ) )
 }
 cat( " ++ Extracting reads in width range : <",RLENRANGE[1],";",RLENRANGE[2],"> ++\n",sep="" )
 tmp.gr       <- as( tmp.tga,"GRanges" )
 tmp.orig     <- tmp.gr
 tmp.check    <-
  mcols(tmp.gr[
   ( ( width(tmp.gr)         < RLENRANGE[1] ) | ( width(tmp.gr)         > RLENRANGE[2] ) ) |
   ( ( mcols(tmp.gr)[["nM"]] < NMISRANGE[1] ) | ( mcols(tmp.gr)[["nM"]] > NMISRANGE[2] ) )
  ])[["qname"]]
 tmp.gr       <-
  tmp.gr[
   ( ( width(tmp.gr)         >= RLENRANGE[1] ) & ( width(tmp.gr)         <= RLENRANGE[2] ) ) &
   ( ( mcols(tmp.gr)[["nM"]] >= NMISRANGE[1] ) & ( mcols(tmp.gr)[["nM"]] <= NMISRANGE[2] ) )
  ]
 tmp.check    <- tmp.check[ tmp.check %in% mcols(tmp.gr)[["qname"]] ]
 tmp.CHECK    <- tmp.orig[ mcols(tmp.orig)[["qname"]] %in% tmp.check ]
 tmp.CHECK    <-
  tmp.CHECK[
   ( ( width(tmp.CHECK)         < RLENRANGE[1] ) | ( width(tmp.CHECK)         > RLENRANGE[2] ) ) |
   ( ( mcols(tmp.CHECK)[["nM"]] < NMISRANGE[1] ) | ( mcols(tmp.CHECK)[["nM"]] > NMISRANGE[2] ) )
  ]
 tmp.CHECK    <- sum( as.numeric(gsub( pattern="^[0-9]*[-]",replacement="",x=mcols(tmp.CHECK)[["qname"]] ))/mcols(tmp.CHECK)[["NH"]] )
 if ( !is.null(SENSING) ) {
  if ( SENSING=="anti" ) {
   tmp.gr <- gr_switchStrand( inpGR=tmp.gr )
  }
 }
 cat( " ++ Unifying seqnames ++\n",sep="" )
 tmp.gr       <- grUnifySeqnames(tmp.gr)
 ANNOT        <- grUnifySeqnames(ANNOT)
 all.biotypes <- data.table( biotype=ALL.BTP )
 if ( FRAC ) {
 cat( " ++ Extracting fractions ++\n",sep="" )
  mcols(tmp.gr)[["OUT"]] <- as.numeric(gsub( pattern="^[0-9]*[-]",replacement="",x=mcols(tmp.gr)[["qname"]] ))
  mcols(tmp.gr)[["OUT"]] <- mcols(tmp.gr)[["OUT"]]/mcols(tmp.gr)[["NH"]]
 } else      {
 cat( " ++ Extracting counts ++\n",sep="" )
  mcols(tmp.gr)[["OUT"]] <- as.numeric(gsub( pattern="^[0-9]*[-]",replacement="",x=mcols(tmp.gr)[["qname"]] ))
 }
 cat( " ++ Extracting read-lengths ++\n",sep="" )
 mcols(tmp.gr)[["rwidth"]] <- width(tmp.gr)

 cat( " ++ Mapping statistics ++\n",sep="" )
 stats.dt  <-
  data.table(
   cnt=as.numeric(gsub( pattern="^[0-9]*[-]",replacement="",x=(mcols(tmp.gr)[["qname"]])[ !duplicated(mcols(tmp.gr)[["qname"]]) ] )),
   len=width(                                                        tmp.gr[              !duplicated(mcols(tmp.gr)[["qname"]]) ])
  )
 stats.dt  <- rbind( stats.dt,data.table( cnt=0,len=seq(RLENRANGE[1],RLENRANGE[2]) ) )
 stats.dt  <- stats.dt[ , { list( cnt=sum(cnt) ) } , by="len" ]
 stats.dt  <- stats.dt[ order(len) ]

 cat( " ++ Mapping to annotation (min.overlap=",MINOVERLAP,",ignore.strand=",IGNORE_STR,") ++\n",sep="" )
 tmp.ovl   <- findOverlaps( query=tmp.gr,subject=ANNOT,minoverlap=MINOVERLAP,ignore.strand=IGNORE_STR )
 tmp.que   <-   queryHits(tmp.ovl)
 tmp.sub   <- subjectHits(tmp.ovl)
 tmp.sub   <- tmp.sub[ !duplicated(tmp.que) ]
 tmp.que   <- tmp.que[ !duplicated(tmp.que) ]
 tmp.dt    <-
  data.table(
   biotype  = mcols(ANNOT[ tmp.sub ])[[ANNOT_BTPCOL]],
   OUT      = mcols(tmp.gr[   tmp.que ])[["OUT"]],
   rwidth   = mcols(tmp.gr[   tmp.que ])[["rwidth"]]
  )
 tmp.gr    <- tmp.gr[ !( seq_along(tmp.gr) %in% tmp.que ) ]
 not.dt    <-
  data.table(
   biotype  = "not_annotated",
   OUT      = mcols(tmp.gr)[["OUT"]],
   rwidth   = mcols(tmp.gr)[["rwidth"]]
  )
 ran.dt    <- data.table( biotype="not_annotated",OUT=0,rwidth=seq(RLENRANGE[1],RLENRANGE[2]) )
 tmp.dt    <- rbind( tmp.dt,not.dt,ran.dt )
 bcnt.dt   <- tmp.dt[ , { list( OUT=sum(OUT) ) } , by=c("biotype","rwidth") ]
 bcnt.dt   <- merge( bcnt.dt,all.biotypes,by="biotype",all=T,sort=F )
 bcnt.dt   <- bcnt.dt[ , lapply( X=.SD,FUN=function(C){ C[is.na(C)]<-0; return(C) } ) , .SDcols=colnames(bcnt.dt) ]
 bcnt.dt   <- dcast( data=bcnt.dt,formula=biotype~rwidth,fun.aggregate=sum,value.var="OUT" )
 bcnt.dt   <- bcnt.dt[ , colnames(bcnt.dt)!="0" , with=F ]

 return( list( STAT=stats.dt,EXP=bcnt.dt,RMV_BY_READLENGTH=tmp.CHECK,NORM=tmp.NORM ) )

}

