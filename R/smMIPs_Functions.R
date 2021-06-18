########################################
#### SOURCE FUNCTIONS
#######################################

########################################################################################
############## Source functions for map_smMIPs_extract_UMIs.R
########################################################################################

# Using Rsamtools to load directly from the bam file
load.bam<-function(file){
  bam <- scanBam(file)[[1]]
  df<-do.call("DataFrame",bam)
  # make sure that "chr" character will be included since different genome references may or may not include it and it must be matched with the panel file
  if(length(grep("chr",df$rname))==0){
    df$rname=paste0("chr",df$rname)
    df$mrnm=paste0("chr",df$mrnm)
  }
  return(df)
}

# Load the MIPgen output to generate a dataframe of the  panel
load.panel<-function(file){
  # load from the current MIPgen format of the panel file to generate a dataframe.
  df=as.data.frame(fread(file,header=T,sep="\t"))
  q1=apply(df[,c("ext_probe_start","ext_probe_stop","lig_probe_start","lig_probe_stop")],1,min)
  q2=apply(df[,c("ext_probe_start","ext_probe_stop","lig_probe_start","lig_probe_stop")],1,max)
  df$backbone=apply(df,1, function(x) gsub(x["ext_probe_sequence"], "", x["mip_sequence"]))
  df$backbone=apply(df,1, function(x) gsub(x["lig_probe_sequence"], "", x["backbone"]))
  df$backbone=gsub("N","",df$backbone)
  df$length.left.umi=apply(df,1, function(x) gsub(paste0(x["backbone"],".*"),"", x["mip_sequence"]))
  df$length.left.umi=nchar(df$length.left.umi)-nchar(gsub("N","", df$length.left.umi))
  df$length.right.umi=apply(df,1, function(x) gsub(paste0(".*",x["backbone"]),"", x["mip_sequence"]))
  df$length.right.umi=nchar(df$length.right.umi)-nchar(gsub("N","", df$length.right.umi))
  panel=data.frame(id=df$mip_name,chrom=df$chr,start=q1,end=q2,ext.probe=df$ext_probe_sequence,lig.probe=df$lig_probe_sequence,target_seq=df$scan_target_sequence,
                   length.left.umi=df$length.left.umi,length.right.umi=df$length.right.umi,target_start=df$mip_scan_start_position,target_stop=df$mip_scan_stop_position,probe_strand=df$probe_strand,stringsAsFactors = FALSE)
  # make sure that "chr" character will be included since different genome references may or may not include it and it must be matched with the bam file
  if(length(grep("chr",panel$chrom))==0){panel$chrom=paste0("chr","",panel$chrom)}
  panel
}

# Filter reads based on MapQ helps to remove self ligated smMIPs without target sequence amoung other bad alignments
filter.on.mappingscore<-function(d){
  cat("Filtering reads based on low mapping score\n")
  index=which(d$samtable$mapq<opt$MAPQ | is.na(d$samtable$mapq))
  d$summary$mapping.score.filtered.reads=length(index)
  if(length(index)>0){
    filtered<-d$samtable[index,]
    filtered$reason<-"mapping_score_filter"
    d$filtered<-rbind(d$filtered,filtered)
    d$samtable=d$samtable[-index,]
  }
  d
}

#Filter reads that are not properly aligned as defined by the flags below
filter.on.flag<-function(d){
  cat("Filtering reads based on bad sam flag\n")
  index=which(d$samtable$flag!=83 &
                d$samtable$flag!=163 &
                d$samtable$flag!=99 &
                d$samtable$flag!=147 &
                d$samtable$flag!=81 &
                d$samtable$flag!=161 &
                d$samtable$flag!=97 &
                d$samtable$flag!=145)
  d$summary$sam.flag.filtered.reads=length(index)
  if(length(index)>0){
    filtered<-d$samtable[index,]
    filtered$reason="sam_flag"
    d$filtered=rbind(d$filtered,filtered)
    d$samtable=d$samtable[-index,]
  }
  d
}

# Filter reads that are align too close or too far from their mate.
filter.on.mate_distance<-function(d){ #If there are smMIPs that were designed to detect known structural alternations (i.e., smMIPs arms falls on different chr), they will be filtered. If detection of such variants is required, an artificial genome reference is needed for read alignment
  cat("Filtering reads based on mate distance\n")
  a=min(nchar(as.character(d$panel$ext.probe))+nchar(as.character(d$panel$lig.probe))+nchar(as.character(d$panel$target_seq)))-5 #-+5 allows some flexibility
  b=max(nchar(as.character(d$panel$ext.probe))+nchar(as.character(d$panel$lig.probe))+nchar(as.character(d$panel$target_seq))+d$panel$length.left.umi+d$panel$length.right.umi)+5
  index=which(abs(d$samtable$isize)>b | abs(d$samtable$isize)<a)
  d$summary$unexpected.mate.distance.filtered.reads=length(index)
  if(length(index)>0){
    filtered<-d$samtable[index,]
    filtered$reason="pair_mapped_too_close_or_too_far"
    d$filtered=rbind(d$filtered,filtered)
    d$samtable=d$samtable[-index,]
  }
  d
}

# Filter reads with hard clip cigar
filter.on.hard.clip<-function(d){
  cat("Filtering hard clipped reads\n")
  index=grep("H",d$samtable$cigar)
  d$summary$hard.clip.filtered.reads=length(index)
  if(length(index)>0){
    filtered<-d$samtable[index,]
    filtered$reason="hard_clip"
    d$filtered=rbind(d$filtered,filtered)
    d$samtable=d$samtable[-index,]
  }
  d
}

#Makes sure that only paired end reads remains for the subsequent functions
filter.on.mate_missing<-function(d){
  cat("Filtering reads with missing mates\n")
  qnames=as.data.frame(table(d$samtable$qname))
  index=match(qnames$Var1[which(qnames$Freq==1)],d$samtable$qname)
  d$summary$mate_removed_filtered_reads=length(index)
  if(length(index)>0){
    filtered=d$samtable[index,]
    filtered$reason="passed_but_remove_due_to_their_pair_being_filtered"
    d$filtered=rbind(d$filtered,filtered)
    d$samtable=d$samtable[-index,]
  }
  d
}

# Mapping the smMIPs to thereads based on overalp and distance mesurements. Being called by the map.smips function
map.smip_to_site<-function(chr,l,r,panel){
  assigned<-NA
  # reduce to smmips on the same chromosome
  panel<-panel[which(panel$chrom==chr),]
  if(nrow(panel)>0){
    # calculate read pair overlap ("lpos","rpos") with every potential smmip
    panel$overlap.inserts.smmips <- apply(panel,1,function(x){ length(intersect(l:r,as.numeric(x[3]):as.numeric(x[4])))/  (as.numeric(x[4])-as.numeric(x[3]) + 1 )  } )
    # calculate overlap of each potential smmip ("lpos","rpos") with the reads
    panel$overlap.smmips.inserts <-apply(panel,1,function(x){ length(intersect(l:r,as.numeric(x[3]):as.numeric(x[4])))/  (r-l + 1 )  } )
    # consider only smmips that overlap above the user defined parameter
    panel<-panel[which(panel$overlap.inserts.smmips>=opt$OVERLAP & panel$overlap.smmips.inserts>=opt$OVERLAP),]
    if(nrow(panel)>0){
      panel$sum=panel$overlap.inserts.smmips+panel$overlap.smmips.inserts
      panel=panel[order(panel$sum,decreasing = F),]
      assigned=list(data.frame(panel$id))
    }
  }
  return(assigned)
}

#Verifing the mapping smMIPs based on local alignment scores of smMIPs extantion and ligation arms. Also determined the arms location, thus allow the extraction UMI sequence if exists.
#Being called by the map.smip function
verify_mapping<-function(key,smmip,r1,r2,panel,isize){
  sm=Y1=Y2=c()
  for(i in 1:nrow(smmip[[1]])){
    smmip.info<-panel[panel$id==smmip[[1]]$panel.id[i],]
    first.match<-match(key,r1$ukey)
    sc <- nucleotideSubstitutionMatrix(match = 1, mismatch = 0, baseOnly = TRUE)
    align.check<-0
    y2=y1=NA
    for (r in list(r1,r2)){
      sam.record<-r[first.match,]
      read<-as.character(sam.record$seq)
      ### based on information in the sam.record, determine the arm/probe (extension or ligation) and the orientation of the probe to interrogate
      if(sam.record$flag<100 & sam.record$isize>0){
        probe.to.align<-reverse(chartr("ATGC","TACG",as.character(smmip.info$lig.probe)))
        n=nchar(probe.to.align)
        align_score=matchPattern(probe.to.align, read, max.mismatch=ceiling(0.25*n)) ##### allows 25% mismatch. Might change to user-defined parameters in next version

        idx1=which(align_score@ranges@start==(smmip.info$length.left.umi+1) & align_score@ranges@width==n)
        if(length(idx1)!=0){
          y1=align_score@ranges@start[idx1]
        } else {
          idx2=which(align_score@ranges@start==(smmip.info$length.left.umi+1) & align_score@ranges@width!=n)
          if(length(idx2)!=0) {
            y1=align_score@ranges@start[idx2][which(align_score@ranges@width[idx2]==max(align_score@ranges@width[idx2]))]
          } else {
            align_score=matchPattern(probe.to.align, read, max.mismatch=ceiling(0.25*n),with.indels = T)
            idx3=which(align_score@ranges@start==(smmip.info$length.left.umi+1))
            if(length(idx3)!=0){
              y1=align_score@ranges@start[idx3]
            }
          }
        }

      } else if (sam.record$flag>100 & sam.record$isize<0){
        probe.to.align<-reverse(chartr("ATGC","TACG",as.character(smmip.info$ext.probe)))
        n=nchar(probe.to.align)
        read=gsub("N",strsplit(probe.to.align,"")[[1]][n],read) #HOW DOES IT SUPPOSED TO BE IN THE OTHERS?
        align_score=matchPattern(probe.to.align, read, max.mismatch=ceiling(0.25*n))

        idx1=which(nchar(read)-smmip.info$length.right.umi==end(align_score@ranges) & align_score@ranges@width==n)
        if(length(idx1)!=0){
          y2=end(align_score@ranges)[idx1]
        } else {
          idx2=which(nchar(read)-smmip.info$length.right.umi==end(align_score@ranges) & align_score@ranges@width!=n)
          if(length(idx2)!=0){
            y2=end(align_score@ranges)[idx2][which(align_score@ranges@width[idx2]==max(align_score@ranges@width[idx2]))]
          } else {
            align_score=matchPattern(probe.to.align, read, max.mismatch=ceiling(0.25*n),with.indels = T)
            idx3=which(nchar(read)-smmip.info$length.right.umi==end(align_score@ranges))
            if(length(idx3)!=0){
              y2=end(align_score@ranges)[idx3]
            }
          }
        }

      } else if(sam.record$flag<100 & sam.record$isize<0){
        probe.to.align<-as.character(smmip.info$lig.probe)
        n=nchar(probe.to.align)
        align_score=matchPattern(probe.to.align, read, max.mismatch=ceiling(0.25*n))

        idx1=which(nchar(read)-smmip.info$length.left.umi==end(align_score@ranges) & align_score@ranges@width==n)
        if(length(idx1)!=0){
          y1=end(align_score@ranges)[idx1]
        } else {
          idx2=which(nchar(read)-smmip.info$length.left.umi==end(align_score@ranges) & align_score@ranges@width!=n)
          if(length(idx2)!=0){
            y1=end(align_score@ranges)[idx2][which(align_score@ranges@width[idx2]==max(align_score@ranges@width[idx2]))]
          } else {
            align_score=matchPattern(probe.to.align, read, max.mismatch=ceiling(0.25*n),with.indels = T)
            idx3=which(nchar(read)-smmip.info$length.left.umi==end(align_score@ranges))
            if(length(idx3)!=0){
              y1=end(align_score@ranges)[idx3]
            }
          }
        }

      } else if(sam.record$flag>100 & sam.record$isize>0){
        probe.to.align<-as.character(smmip.info$ext.probe)
        read=gsub("N",strsplit(probe.to.align,"")[[1]][1],read)
        n=nchar(probe.to.align)
        align_score=matchPattern(probe.to.align, read, max.mismatch=ceiling(0.25*n))

        idx1=which(align_score@ranges@start==(smmip.info$length.right.umi+1) & align_score@ranges@width==n)
        if(length(idx1)!=0){
          y2=align_score@ranges@start[idx1]
        } else {
          idx2=which(align_score@ranges@start==(smmip.info$length.right.umi+1) & align_score@ranges@width!=n)
          if(length(idx2)!=0){
            y2=align_score@ranges@start[idx2][which(align_score@ranges@width[idx2]==max(align_score@ranges@width[idx2]))]
          } else {
            align_score=matchPattern(probe.to.align, read, max.mismatch=ceiling(0.25*n),with.indels = T)
            idx3=which(align_score@ranges@start==(smmip.info$length.right.umi+1))
            if(length(idx3)!=0){
              y2=align_score@ranges@start[idx3]
            }
          }
        }
      }
    }
    if (!is.na(y1) & !is.na(y2)){
      sm=append(sm,as.character(smmip[[1]]$panel.id[i]))
      Y1=append(Y1,y1)
      Y2=append(Y2,y2)
    }
  }

  if(is.null(sm)){
    sm=as.character(smmip[[1]]$panel.id[i]) #when there is an overlap but validation failed the last one tested (the one with the highest overlap) will be written
    Y1=NA
    Y2=NA
  }

  list(list(sm),list(Y1),list(Y2))
}


#Main smMIP-read mapping function.
map.smips<-function(d){
  d$samtable<-d$samtable[order(d$samtable$qname,d$samtable$flag),]   ## this order the file to allow every read to be pulled out as R1 or R2

  R1<-d$samtable[seq(1,dim(d$samtable)[1],2),]
  R2<-d$samtable[seq(2,dim(d$samtable)[1],2),]
  dnames=names(d$samtable)
  d$samtable = NULL

  # determine the unique read "family/key". This is used for mapping back the assigned smmips in the list of sites back to the file records
  R1$lpos<-apply(R1[,c("pos","mpos")],1,min)
  R1$rpos<-R1$lpos+abs(R1$isize)-1
  R1$mpos<-apply(R1[,c("pos","mpos")],1,max)
  R1$ukey<-paste0(as.character(R1$rname),":",as.character(R1$lpos),":",as.character(R1$rpos),":",as.character(R1$isize),":",as.character(R1$mpos),":",as.character(R1$cigar),":",as.character(R2$cigar))

  # for each unique "family/key" detemine what will be the smmip of which the "family/key" overap the most
  # get a list of unique sites, and map candidate smmips to each.
  # Candidates are limited the two smmips with the highest scores bases on overlap the calculated distances ("map.smip_to_site").
  s=as.data.frame(table(R1$ukey))
  sites<-as.data.frame(unique(R1[,c("rname","lpos","rpos","isize","mpos","ukey")]))
  sites$potential_smMIPs<-sites$cut2<-sites$cut1<-sites$readsR1<-NA
  sites$readsR1[match(s$Var1,sites$ukey)]=s$Freq

  # determine the best candidates for each unique "family/key"
  m = list(mclapply(1:nrow(sites) ,mc.cores = opt$threads, mc.cleanup=T, mc.silent=F ,function(i){#
    if(round(i/nrow(sites),1) %in% seq(0.1,0.9,0.1)) {
      system(paste0("printf '\\rMapping smMIPs to reads. Considering overlap and distance measurements :  ",round(100*i/nrow(sites)),"%%      '"))    }
    if(i == nrow(sites)) {
      system(paste0("printf '\\rMapping smMIPs to reads. Considering overlap and distance measurements :  100%%'"))
    }
    x<-sites[,c("rname","lpos","rpos")][i,]
    map.smip_to_site(as.character(unlist(x[1])),as.integer(x[2]),as.integer(x[3]),d$panel)[[1]]
  }))

  err.idx=which(unlist(lapply(1:length(m[[1]]), function(i) is.null(m[[1]][[i]]))=="TRUE"))
  if(length(err.idx)==0){
     system(paste0("printf '\\rMapping smMIPs to reads. Considering overlap and distance measurements :  100%%'"))
  }
  cat("\n")

  attempt=0
  while(length(err.idx)>0 & attempt<10){
    attempt=attempt+1
    if(attempt==1){
    cat("An error in one or multiply cores occurred. Trying again, please wait\n")
    }
    cat(paste("There are",length(err.idx),"missing pieces of information to fill\n"))
    m.fill=list(mclapply(err.idx ,mc.cores = opt$threads, mc.cleanup=T, mc.silent=F ,function(i){
	if(round(i/length(err.idx),1) %in% seq(0.1,1,0.1)) {
	  system(paste0("printf '\\rTrying to fill the missing data :  "," Attempt number ",attempt," ... '"))
	}
	if(round(i/length(err.idx),1) %in% seq(0,0.9,0.1)) {
          system(paste0("printf '\\rTrying to fill the missing data :  "," Attempt number ",attempt,"     '"))
        }
	x<-sites[,c("rname","lpos","rpos")][i,]
     map.smip_to_site(as.character(unlist(x[1])),as.integer(x[2]),as.integer(x[3]),d$panel)[[1]]
    }))
    m[[1]][err.idx]=m.fill[[1]]
    err.idx=which(unlist(lapply(1:length(m[[1]]), function(i) is.null(m[[1]][[i]]))=="TRUE"))
   }

  if(length(err.idx)==0 & attempt!=0){
     cat("Succeed to fill in all the missing data!\n")
  } else if (length(err.idx)!=0 & attempt==10) {
     cat("Sorry, could not fill in all the missing data. One or more cores still did not deliver a result\n")
     quit()
  }
  sites[,"potential_smMIPs"]=m

  # verification based on actual arms alignment
  idx<-which(!is.na(sites$potential_smMIPs))

  m <-mclapply(idx ,mc.cores = opt$threads, mc.cleanup=T, mc.silent=F ,function(i){
    if(round(which(idx==i)/length(idx),1) %in% seq(0.1,0.9,0.1)){
      system(paste0("printf '\\rVerifying mapping by local smMIP arms alignment :  ",round(100*which(idx==i)/length(idx)),"%%     '"))
    }
    if(i %in% idx[(length(idx)-opt$threads):length(idx)]) {
      system(paste0("printf '\\rVerifying mapping by local smMIP arms alignment : 100%%'"))
    }
    verify_mapping(sites[i,]$ukey,sites[i,]$potential_smMIPs,R1,R2,d$panel,sites[i,]$isize)
  })

  err.idx=which(unlist(lapply(1:length(m), function(i) is.null(m[[i]]))=="TRUE"))
  if(length(err.idx)==0){
    system(paste0("printf '\\rVerifying mapping by local smMIP arms alignment : 100%%'"))
  }
  cat("\n")

  attempt=0
  while(length(err.idx)>0 & attempt<10){
    attempt=attempt+1
    if(attempt==1){
      cat("An error in one or multiple cores occurred. Trying again, please wait\n")
    }
    cat(paste("There are",length(err.idx),"missing pieces of information to fill\n"))
    m.fill=list(mclapply(err.idx ,mc.cores = opt$threads, mc.cleanup=T, mc.silent=F ,function(i){
        if(round(i/length(err.idx),1) %in% seq(0.1,1,0.1)) {
          system(paste0("printf '\\rTrying to fill the missing data :  "," Attempt number ",attempt," ... '"))
        }
        if(round(i/length(err.idx),1) %in% seq(0,0.9,0.1)) {
          system(paste0("printf '\\rTrying to fill the missing data :  "," Attempt number ",attempt,"     '"))
        }
        verify_mapping(sites[idx[i],]$ukey,sites[idx[i],]$potential_smMIPs,R1,R2,d$panel,sites[idx[i],]$isize)
    }))
    m[err.idx]=m.fill[[1]]
    err.idx=which(unlist(lapply(1:length(m), function(i) is.null(m[[i]]))=="TRUE"))
  }

  if(length(err.idx)==0 & attempt!=0){
     cat("Succeed to fill in all the missing data!\n")
  } else if (length(err.idx)!=0 & attempt==10) {
     cat("Sorry, could not fill in all the missing data. One or more cores still did not deliver a result\n")
     quit()
  }

  sites$smMIP=NA
  sites$smMIP[idx]=unlist(lapply(1:length(m),function(i){paste(unlist(m[[i]][[1]]),collapse = ",")}))
  sites$cut1[idx]=unlist(lapply(1:length(m),function(i){paste(unlist(m[[i]][[2]]),collapse = ",")}))
  sites$cut2[idx]=unlist(lapply(1:length(m),function(i){paste(unlist(m[[i]][[3]]),collapse = ",")}))

  #when there are multiple smMIPs that pass the overlap thershold and validation, the smMIP with the most number of reads will be selected
  idx2=grep(",",sites$smMIP[idx])
  if(length(idx2)>0){
    tab=data.frame(V1=sites$readsR1,V2=sites$smMIP)
    tab$V2[which(tab$V2=="")]=NA
    s <- strsplit(as.character(tab$V2), split = ",")
    tab=as.data.table(data.frame(V1 = rep(tab$V1, sapply(s, length)), V2 = unlist(s)))
    tab=tab[,s:=sum(V1),by=V2]
    tab=unique(tab[,V1:=NULL])

    m2=lapply(idx2,function(i){
      x=tab[V2 %in% unlist(m[[i]][[1]])]
      as.character(x$V2[which(x$s==max(x$s))[1]])
    })
    sites$smMIP[idx[idx2]]=unlist(m2)
    sites$cut1[idx[idx2]]=gsub(",.*","",sites$cut1[idx[idx2]])
    sites$cut2[idx[idx2]]=gsub(",.*","",sites$cut2[idx[idx2]])
  }

  #assign the smmips names back to all read1 and read2, as well as the determined arm start points to decide whether to pass or flag the UMI sequences
  rownames(sites)<-sites$ukey
  R1$smMIP<-R2$smMIP<-sites[R1$ukey,]$smMIP
  R1$cut1<-R2$cut1<-sites[R1$ukey,]$cut1
  R1$cut2<-R2$cut2<-sites[R1$ukey,]$cut2

  d$samtable<-rbind(R1[,c(dnames,"smMIP","cut1","cut2")],R2[,c(dnames,"smMIP","cut1","cut2")])
  cat("")
  d
}

# Remove reads that were not assigned with any smMIP
filter.offtarget<-function(d){
  index<-which(is.na(d$samtable$smMIP)) #
  d$summary$offtarget_reads=length(index)
  if(length(index)>0){
    filtered<-d$samtable[index,]
    filtered$reason="unanticipated_target"
    d$filtered=rbind(d$filtered,filtered[,names(d$filtered)])
    d$summary$offtarget_reads=length(index)
    d$samtable=d$samtable[-index,]
  }
  d
}

# Extarct and flag UMIs
extract_umi<-function(d){
  # UMI location is dependant on read and orientation
  samtab<-as.data.frame(d$samtable)
  samtab$umi<-""

  idx1=which(samtab$flag<100 & samtab$isize>0)
  m=match(samtab$smMIP[idx1],d$panel$id)
  samtab$umi[idx1]=substr(samtab$seq[idx1],1, d$panel$length.left.umi[m])
  idx2=unique(c(grep("N",samtab$umi[idx1]),which(samtab$cut1[idx1]=="NA")))
  samtab$umi[idx1[idx2]]=paste0(samtab$umi[idx1[idx2]],"flagged")

  idx1=which(samtab$flag>100 & samtab$isize<0)
  m=match(samtab$smMIP[idx1],d$panel$id)
  samtab$umi[idx1]=substr(samtab$seq[idx1],nchar(samtab$seq[idx1])-d$panel$length.right.umi[m]+1,nchar(samtab$seq[idx1]))
  idx2=unique(c(grep("N",samtab$umi[idx1]),which(samtab$cut2[idx1]=="NA")))
  samtab$umi[idx1[idx2]]=paste0(samtab$umi[idx1[idx2]],"flagged")

  idx1=which(samtab$flag<100 & samtab$isize<0)
  m=match(samtab$smMIP[idx1],d$panel$id)
  samtab$umi[idx1]=substr(samtab$seq[idx1],nchar(samtab$seq[idx1])-d$panel$length.left.umi[m]+1,nchar(samtab$seq[idx1]))
  idx2=unique(c(grep("N",samtab$umi[idx1]),which(samtab$cut1[idx1]=="NA")))
  samtab$umi[idx1[idx2]]=paste0(samtab$umi[idx1[idx2]],"flagged")

  idx1=which(samtab$flag>100 & samtab$isize>0)
  m=match(samtab$smMIP[idx1],d$panel$id)
  samtab$umi[idx1]=substr(samtab$seq[idx1],1, as.numeric(d$panel$length.right.umi[m]))
  idx2=unique(c(grep("N",samtab$umi[idx1]),which(samtab$cut2[idx1]=="NA")))
  samtab$umi[idx1[idx2]]=paste0(samtab$umi[idx1[idx2]],"flagged")

  d$samtable=samtab
  d
}

# Incorporate the smMIP ID and the UMI sequences in the reads' names
adjust_readname<-function(d){

  qnames<-unique(d$samtable$qname)
  R1.index<-which(bamFlagAsBitMatrix(as.integer(c(d$samtable$flag)))[,"isFirstMateRead"]==1)
  umi.1<-d$samtable[R1.index,]$umi
  smmips.1<-d$samtable[R1.index,]$smMIP
  names(umi.1)<-names(smmips.1)<-d$samtable[R1.index,]$qname

  R2.index<-which(bamFlagAsBitMatrix(as.integer(c(d$samtable$flag)))[,"isSecondMateRead"]==1)
  umi.2<-d$samtable[R2.index,]$umi
  smmips.2<-d$samtable[R2.index,]$smMIP
  names(umi.2)<-names(smmips.2)<-d$samtable[R1.index,]$qname

  umi.12<-paste0(umi.1[qnames],"_",umi.2[qnames])
  if(unique(unique(d$samtable[R1.index,]$umi)=="") &  unique(unique(d$samtable[R2.index,]$umi)=="")) {
    umi.12<-""
  } else if (unique(unique(d$samtable[R2.index,]$umi)=="")) {
    #umi.12<-paste0(umi.1[qnames],flag.1[qnames])
    umi.12<-umi.1[qnames]
  } else if (unique(unique(d$samtable[R1.index,]$umi)=="")) {
    umi.12<-umi.2[qnames]
  }

  # check to make sure the smmips are identical
  if(length(which(smmips.1[qnames]!=smmips.2[qnames]))!=0){
    print("Warning: there are ", length(which(smmips.1[qnames]!=smmips.2[qnames]))," read pairs where one read was assigned one smMIP and the second read was assigned with another",quote = F)
  }

  # generate umi usage summary
  if(!((length(unique(data$panel$length.left.umi))==1 & 0 %in% unique(data$panel$length.left.umi)) &
       (length(unique(data$panel$length.right.umi))==1 & 0 %in% unique(data$panel$length.right.umi)))){
    suffix<-paste(smmips.1,umi.12,sep=":")
    if(length(grep("flag",suffix))>0){
      count=as.data.frame(matrix(unlist(strsplit(suffix[-grep("flag",suffix)],split = ":")), ncol = 2, byrow = TRUE))
      d$umi_usage=table(count$V2,count$V1)
    }
  } else {
    suffix<-smmips.1
  }

  names(suffix)<-qnames
  d$samtable$qname<-paste(d$samtable$qname,suffix[d$samtable$qname],sep="||")
  d$samtable=d$samtable[,c(1:3,5,7:13)]
  d

}

# Incorporate the reason of why the read was filtered to reads' names. Paired reason is being written to keep equal read pair names. This is nessesarry for viewing reads as pair in IGV
write_reason_to_flitered_read<-function(d){
  qnames<-unique(d$filtered$qname)
  R1.index<-which(bamFlagAsBitMatrix(as.integer(c(d$filtered$flag)))[,"isFirstMateRead"]==1)
  reason.1<-d$filtered$reason[R1.index]
  names(reason.1)<-d$filtered$qname[R1.index]

  R2.index<-which(bamFlagAsBitMatrix(as.integer(c(d$filtered$flag)))[,"isSecondMateRead"]==1)
  reason.2<-d$filtered$reason[R2.index]
  names(reason.2)<-d$filtered$qname[R2.index]

  reason.12<-d$filtered$qname
  names(reason.12)=d$filtered$qname

  d$filtered$qname=paste0(reason.12,"||pair_filter_reason:",reason.1[names(reason.12)],"/",reason.2[names(reason.12)])
  d$filtered=d$filtered[,c(1:3,5,7:13)]
  d
}


########################################################################################
############## Source functions QC_plots.R
########################################################################################
#Note: all the samples' 'raw_coverage_per_smMIP.txt' and 'filtered_read_counts.txt' files need to be copied to a single folder by the user
load_qc_coverage_per_smMIP <-function(){
  files=list.files(opt$dir,pattern="raw_coverage_per_smMIP.txt")
  if(length(files)==0){print("File names must be formated as Sample_ID_raw_coverage_per_smMIP.txt")}
  for (n in files) {
    f = fread(paste0(opt$dir,"/",n),header=T,sep="\t",showProgress = FALSE)
    if(n==files[1]){
      a=data.table(matrix(nrow=nrow(f),ncol=length(files)+1))
      names(a)=c("smMIPs",gsub("_raw_coverage_per_smMIP.txt","",files))
      a$smMIPs=f$smmips
    }
    n=gsub("_raw_coverage_per_smMIP.txt","",n)
    a[[n]]=f$coverage
  }
  a
}

load_qc_filtered.reads <- function(){
  files=list.files(opt$dir,pattern="filtered_read_counts.txt")
  if(length(files)==0){print("File names must be formated as Sample_ID_filtered_read_counts.txt")}
  for (n in files) {
    f = fread(paste0(opt$dir,"/",n),header=T,sep="\t",showProgress = FALSE)
    if(n==files[1]){
      a=data.frame(matrix(nrow=length(files),ncol=9))
      names(a)=c("sample_ID",names(f))
      a$sample_ID=gsub("_filtered_read_counts.txt","",files)
    }
    n=gsub("_filtered_read_counts.txt","",n)
    a[which(a$sample_ID==n),2:9]=f
  }
  for(j in 3:8){
    a[,j]=100*a[,j]/a[,2]
  }
  a
}

########################################################################################
############## Source functions for smMIP_level_raw_and_consensus_pileups.R
########################################################################################

pileup_foreach_smmip <- function(i){
  write.table(rbind(header,smmip_piles[[i]][,c("qname","flag","rname","pos","mapq","cigar","mrnm","mpos","isize","seq","qual")]),file=paste0(opt$tmp.output,"/",opt$sample.name,"_",names(smmip_piles)[i],"_tmp.sam"),col.names = F,row.names = F,quote = F,sep = '\t')
  suppressWarnings(asBam(paste0(opt$tmp.output,"/",opt$sample.name,"_",names(smmip_piles)[i],"_tmp.sam"),paste0(opt$tmp.output,"/",opt$sample.name,"_",names(smmip_piles)[i],"_tmp"),overwrite=T))
  invisible(file.remove(paste0(opt$tmp.output,"/",opt$sample.name,"_",names(smmip_piles)[i],"_tmp.sam")))
  bf <- open(BamFile(paste0(opt$tmp.output,"/",opt$sample.name,"_",names(smmip_piles)[i],"_tmp.bam")))
  sbp <- ScanBamParam(which=GRanges(data$panel$chr[data$panel$id==names(smmip_piles)[i]], IRanges(data$panel$target_start[data$panel$id==names(smmip_piles)[i]], data$panel$target_stop[data$panel$id==names(smmip_piles)[i]]))) #pile only in the target sequence, not including the arms
  res <- as.data.table(pileup(bf,scanBamParam = sbp,pileupParam=p_param))
  invisible(file.remove(paste0(opt$tmp.output,"/",opt$sample.name,"_",names(smmip_piles)[i],"_tmp.bam")))
  invisible(file.remove(paste0(opt$tmp.output,"/",opt$sample.name,"_",names(smmip_piles)[i],"_tmp.bam.bai")))
#  close(bf)
  if(dim(res)[1]>0){
    if(opt$rank=="F"){
      res[nucleotide=="+",count_no_insertion := 0]
      res[nucleotide!="+",count_no_insertion := count]  # Some bases may not pass the other pileup filters. This also may bias VAF)
      res[,coverage_at_position := sum(count_no_insertion),by=list(pos,strand)]
      res=res[,c("VAF","smMIP") := list(count/coverage_at_position, names(smmip_piles)[i])]
      res=res[,c("which_label","count_no_insertion"):=NULL]
    } else if (opt$rank=="T"){
      res[nucleotide=="+",count_no_insertion := 0]
      res[nucleotide!="+",count_no_insertion := count]  # Some bases may not pass the other pileup filters. This also may bias VAF)
      res[,coverage_at_position := sum(count_no_insertion),by=list(pos,strand)]
      res[, rank := frank(-count,ties.method = "min" ),by = list(pos,strand)]
      res=subset(res[,c("VAF","smMIP") := list(count/coverage_at_position, names(smmip_piles)[i])], rank<3) #if there are multiple alleles that rank 2nd, it will keep them
      res=res[,c("which_label","count_no_insertion","rank"):=NULL]
    }
  } else {
    res=data.table(seqnames=factor(), pos=integer(), strand=factor(), nucleotide=factor(), count=integer(), coverage_at_position=numeric(), VAF=numeric(),smMIP=character())
  }
  return(res)
}


pileup_foreach_smmip.umi <- function(i){
  write.table(rbind(header,smmip_umi_piles[[i]][,c("qname","flag","rname","pos","mapq","cigar","mrnm","mpos","isize","seq","qual")]),file=paste0(opt$tmp.output,"/",opt$sample.name,"_",names(smmip_umi_piles)[i],"_tmp.sam"),col.names = F,row.names = F,quote = F,sep = '\t')
  suppressWarnings(asBam(paste0(opt$tmp.output,"/",opt$sample.name,"_",names(smmip_umi_piles)[i],"_tmp.sam"),paste0(opt$tmp.output,"/",opt$sample.name,"_",names(smmip_umi_piles)[i],"_tmp"),overwrite=T))
  invisible(file.remove(paste0(opt$tmp.output,"/",opt$sample.name,"_",names(smmip_umi_piles)[i],"_tmp.sam")))
  bf <- open(BamFile(paste0(opt$tmp.output,"/",opt$sample.name,"_",names(smmip_umi_piles)[i],"_tmp.bam")))
  sbp <- ScanBamParam(which=GRanges(data$panel$chr[data$panel$id==smmip_umi_piles[[i]]$smmip[1]], IRanges(data$panel$target_start[data$panel$id==smmip_umi_piles[[i]]$smmip[1]], data$panel$target_stop[data$panel$id==smmip_umi_piles[[i]]$smmip[1]]))) #pile only in the target sequence, not including the arms
  res <- as.data.table(pileup(bf,scanBamParam = sbp,pileupParam=p_param))
  invisible(file.remove(paste0(opt$tmp.output,"/",opt$sample.name,"_",names(smmip_umi_piles)[i],"_tmp.bam")))
  invisible(file.remove(paste0(opt$tmp.output,"/",opt$sample.name,"_",names(smmip_umi_piles)[i],"_tmp.bam.bai")))
#  close(bf)

  if(dim(smmip_umi_piles[[i]])[1]/2>=opt$family.size){ #include families with user-defined size or higher. It is important to be aware that if not all the families are included it can bias the VAF toward either the mutated or non-mutated alleles
    res[nucleotide=="+",count_no_insertion := 0]
    res[nucleotide!="+",count_no_insertion := count] # Some bases may not pass the other pileup filters. This also may bias VAF)
    res[,coverage_at_position := sum(count_no_insertion),by=list(pos,strand)]
    res[, rank := frank(-count,ties.method = "min" ),by = list(pos,strand)]
    res[coverage_at_position==0, coverage_at_position:=count]
    if(opt$umi=="T"){
      res=subset(res[,c("VAF","smMIP","UMI") := list(count/coverage_at_position, smmip_umi_piles[[i]]$smmip[1], smmip_umi_piles[[i]]$umi[1])], rank==1)
    } else if (opt$umi=="F"){
      res=subset(res[,c("VAF","smMIP") := list(count/coverage_at_position, smmip_umi_piles[[i]]$smmip[1])], rank==1)
    }
    res=subset(res[,c("which_label","count_no_insertion","rank"):=NULL],res$VAF>=opt$consensus.cutoff)
    res$count=1
    return(res)
  }
}


########################################################################################
############## Source functions for Annotate_SNVs.R
########################################################################################
variant_annotation_using_cellbaseR<- function(i){
  anno=""
  options(annotated=TRUE)
  tryCatch(
    expr = {
      invisible(capture.output(anno <- getVariant(object=cb,
                ids=paste0(genomicBases$chr[i],":",genomicBases$pos[i],":",genomicBases$ref[i],":",genomicBases$alt[i]),
                resource="annotation",param=cbparam)))
    },
    error = function(e){
      print(e)
      options(annotated=FALSE)
    },
    finally = {
      if(isTRUE(getOption("annotated"))){
      #geneID annotation
      gene=anno$consequenceTypes[[1]]
      gene=paste(unique(gene$geneName[which(!is.na(gene$cdnaPosition))]),collapse = ",")
      if(gene==""){
        gene=anno$consequenceTypes[[1]]$geneName[1]
      }

      #Protein (AA change) annotation
      protein=anno$consequenceTypes[[1]]
      if(!is.null(protein$proteinVariantAnnotation)){
        aa=unique(paste0(protein$proteinVariantAnnotation$reference,protein$proteinVariantAnnotation$position,protein$proteinVariantAnnotation$alternate))
        idx=grep("NA",aa)
        if(length(idx)>0){
          l=length(unique(aa[-grep("NA",aa)]))
          if(l==1){
            protein=aa[-grep("NA",aa)]
          } else {
            protein=unique(paste0(protein$ensemblTranscriptId,":",protein$proteinVariantAnnotation$reference,protein$proteinVariantAnnotation$position,protein$proteinVariantAnnotation$alternate))
            protein=paste(protein[-grep("NA",protein)],collapse = "; ")
          }
        } else {
          l=length(unique(aa))
          if(l==1){
            protein=aa
          } else {
            protein=unique(paste0(protein$ensemblTranscriptId,":",protein$proteinVariantAnnotation$reference,protein$proteinVariantAnnotation$position,protein$proteinVariantAnnotation$alternate))
          }
        }
      } else {protein=NA}

      #Cosmic annotation
      u=unique(anno$variantTraitAssociation$cosmic[[1]]$primarySite)
      cosmic=paste(c(unique(anno$variantTraitAssociation$cosmic[[1]]$mutationId),u[!is.na(u)]),collapse = ",")

      #Minimal allele frequency annotation
      maf=anno$populationFrequencies[[1]]
      maf=suppressWarnings(mean(maf$altAlleleFreq[grep("ALL",maf$population)],na.rm = T))

      #Variant type annotation
      variant_type=anno$displayConsequenceType

      #Functional impact score annotation
      cadd_scaled=anno$functionalScore[[1]]$score[2]

      } else {
        gene=protein=cosmic=maf=variant_type=cadd_scaled=""
      }
      }
    )
  annotated=getOption("annotated")

  anno=as.data.frame(cbind(gene,protein,cosmic,maf,
                           variant_type,cadd_scaled,annotated))

  anno[] <- lapply(anno, as.character)

  idx1=which(anno[1,]=="")
  if(length(idx1)>0){
    anno[1,idx1]=NA
  }

  if(round(i/nrow(genomicBases),2) %in% seq(0.01,0.99,0.01)){
    system(paste0("printf '\\rDownloading SNV annotations :  ",round(100*i/nrow(genomicBases)),"%%     '"))
  } else if (i==nrow(genomicBases)){
    system(paste0("printf '\\rDownloading SNV annotations :  100%%     '"))
  }
  cat("")
  anno
}


########################################################################################
############## Source functions for calling_mutations.R
########################################################################################

load.annotation.file=function(d){
  d$annotated.panel=fread(opt$alleles,header=T,sep="\t",showProgress = FALSE)
  if(length(grep("chr",d$annotated.panel))==0){
    d$annotated.panel$chr=paste0("chr",d$annotated.panel)
  }
  # adding indels as alternative alleles at the end of the file since they are not included in the annotated file
  idx1=match(unique(paste(d$annotated.panel$smMIP,d$annotated.panel$chr,d$annotated.panel$pos,d$annotated.panel$ref)),paste(d$annotated.panel$smMIP,d$annotated.panel$chr,d$annotated.panel$pos,d$annotated.panel$ref))
  a1=d$annotated.panel[idx1]
  a1[a1$variant_type %in% c("missense_variant","stop_gained","synonymous_variant","start_lost","stop_lost","stop_retained_variant"), variant_type:=NA ]
  a1[,c("alt","cosmic","maf","protein","cadd_scaled"):=list("+",NA,NA,NA,NA)]
  a2=copy(a1)
  a2[,alt:="-"]
  d$annotated.panel=rbind(d$annotated.panel,a1,a2)
  d
}

populate=function(d){
  d$samples=fread(opt$file,header=T,sep="\t", na.strings = "",showProgress = FALSE, fill=TRUE)
  if("TRUE" %in% duplicated(d$samples$id)){
    cat("ERROR >>>>> The input file must contain unique IDs\n")
    quit()
  }
  for (n in d$samples$id) {
    skip_to_next <- FALSE
    s = tryCatch({
      fread(paste0(opt$summary,"/",n,"_raw_pileup.txt"),header=T,sep="\t",showProgress = FALSE)
    }, error = function(e) {
      if(is.na(d$samples$replicate[grep(n,d$samples$id)])){
        skip_to_next <<- TRUE
      } else {
        print(paste("File not found. Creating an empty pileup for",n))
        data.frame(chr=NA,pos=NA,strand=c("+","-"),nucleotide=NA,count=NA,coverage_at_position=NA,VAF=NA)
      }
    })
    if(skip_to_next) {
      print(paste("File not found. Skipping",n))
      next
    }

    a=data.frame(as.numeric(rep(NA,nrow(d$annotated.panel))))
    names(a)=n

    s.plus=s[s$strand=="+"]
    s.minus=s[s$strand=="-"]
    d$allele.frequency.plus[[n]]=d$allele.frequency.minus[[n]]=
      d$total.depth.plus[[n]]=d$total.depth.minus[[n]]=
      d$non.ref.counts.plus[[n]]=d$non.ref.counts.minus[[n]]=a

    #plus strand
    idx1=match(paste(d$annotated.panel$smMIP,d$annotated.panel$chr,d$annotated.panel$pos,d$annotated.panel$alt),paste(s.plus$smMIP,s.plus$chr,s.plus$pos,s.plus$nucleotide))
    d$non.ref.counts.plus[[n]][which(!is.na(idx1)),]=s.plus$count[idx1[which(!is.na(idx1))]]
    d$allele.frequency.plus[[n]][!is.na(idx1),]=s.plus$VAF[idx1[which(!is.na(idx1))]]
    idx1=match(paste(d$annotated.panel$smMIP,d$annotated.panel$chr,d$annotated.panel$pos),paste(s.plus$smMIP,s.plus$chr,s.plus$pos))
    d$total.depth.plus[[n]][!is.na(idx1),]=s.plus$coverage_at_position[idx1[which(!is.na(idx1))]]
    d$non.ref.counts.plus[[n]][is.na(d$non.ref.counts.plus[[n]])]=0
    d$allele.frequency.plus[[n]][is.na(d$allele.frequency.plus[[n]])]=0
    d$total.depth.plus[[n]][is.na(d$total.depth.plus[[n]])]=0

    #minus strand
    idx1=match(paste(d$annotated.panel$smMIP,d$annotated.panel$chr,d$annotated.panel$pos,d$annotated.panel$alt),paste(s.minus$smMIP,s.minus$chr,s.minus$pos,s.minus$nucleotide))
    d$non.ref.counts.minus[[n]][!is.na(idx1),]=s.minus$count[idx1[which(!is.na(idx1))]]
    d$allele.frequency.minus[[n]][!is.na(idx1),]=s.minus$VAF[idx1[which(!is.na(idx1))]]
    idx1=match(paste(d$annotated.panel$smMIP,d$annotated.panel$chr,d$annotated.panel$pos),paste(s.minus$smMIP,s.minus$chr,s.minus$pos))
    d$total.depth.minus[[n]][!is.na(idx1),]=s.minus$coverage_at_position[idx1[which(!is.na(idx1))]]

    d$non.ref.counts.minus[[n]][is.na(d$non.ref.counts.minus[[n]])]=0
    d$allele.frequency.minus[[n]][is.na(d$allele.frequency.minus[[n]])]=0
    d$total.depth.minus[[n]][is.na(d$total.depth.minus[[n]])]=0

    idx=which(d$samples$id==n)
    if(idx %in% seq(1,length(d$samples$id),ceiling(length(d$samples$id)/100))) {
      system(paste0("printf '\\rLoading pileup data :  ",round(100*idx/length(d$samples$id)),"%%     '"))
    } else if (idx==length(d$samples$id)){
      system(paste0("printf '\\rLoading pileup data :  100%%     '"))
    }
  }
  cat("\n")
  d$allele.frequency.plus=as.data.table(do.call("cbind",d$allele.frequency.plus))
  d$allele.frequency.minus=as.data.table(do.call("cbind",d$allele.frequency.minus))
  d$total.depth.plus=as.data.table(do.call("cbind",d$total.depth.plus))
  d$total.depth.minus=as.data.table(do.call("cbind",d$total.depth.minus))
  d$non.ref.counts.plus=as.data.table(do.call("cbind",d$non.ref.counts.plus))
  d$non.ref.counts.minus=as.data.table(do.call("cbind",d$non.ref.counts.minus))
  d
}

prior.knowledge = function(d){
  system("printf '\\rLeveraging data from Cosmic and information concerning known polymorphisms to increase sensitivity in those genomic positions...'")

  if(length(unique(d$samples$type))==1 & "case" %in% unique(d$samples$type)){
    d$control.names=d$samples$id
  } else if (length(unique(d$samples$type))==2 & "case" %in% unique(d$samples$type) & "control" %in% unique(d$samples$type)){
    d$control.names=d$samples$id[d$samples$type=="control"]
  } else {
      stop("Only 'case' or 'control' should be included in the configuration file under the column 'type'\n.Cases can be run alone or with controls", call.=FALSE)
  }

  idx_cosmic=which(!is.na(d$annotated.panel$cosmic)) #There are only two statuses. Reported or not-available
  idx_maf=which(d$annotated.panel$maf>opt$maf) #The cut-off for inclusion of polymorphism is defined by the user.
  idx=unique(idx_cosmic,idx_maf)

  # to reduce the chance that cosmic mutation will be missed from being called in cases, it looks at their error rates in the control samples and use the median allele frequency.
  # if there is only a single control sample it will use the median across all the cosmic alleles (and the SNPs) in the control sample. VAF cut-off will be applied based on the user input
  control.names=d$control.names
  if(length(idx)>0){
    d$control.allele.frequency.plus = copy(d$allele.frequency.plus[,..control.names])
    #plus
    if(length(control.names)>1){
      l=mclapply(1:length(idx), mc.cores = opt$threads, mc.cleanup=T, mc.silent=F ,function(i) {
        a=unlist(d$control.allele.frequency.plus[idx[i]])
        idx1=as.numeric(which(a!=0 & a<opt$vaf))
        idx2=as.numeric(which(a!=0 & a<=1))
        if(length(idx1)>0){
          m=median(a[idx1]) #median
        } else if (length(idx1)==0) {
          m=length(d$total.depth.plus[idx[i],..control.names])/sum(d$total.depth.plus[idx[i],..control.names]) #equal to 1 read per sample with non-reference allele
        }
        if(length(idx2)==0){
          set(copy(d$control.allele.frequency.plus[idx[i]]), i = NULL, j = as.integer(1:ncol(d$control.allele.frequency.plus)), m)
        } else {
          set(copy(d$control.allele.frequency.plus[idx[i]]), i = NULL, j = as.integer(idx2), m)
        }
      })
      set(d$control.allele.frequency.plus,idx,names(d$control.allele.frequency.plus),rbindlist(l))
    } else { #only one control sample
      a=unlist(d$control.allele.frequency.plus[idx])
      idx1=as.numeric(which(a!=0 & a<opt$vaf))
      if(length(idx1)>0){
        m=median(a[idx1]) #median
        d$control.allele.frequency.plus[idx]= m
      }
    }

    #minus
    d$control.allele.frequency.minus = copy(d$allele.frequency.minus[,..control.names])

    if(length(control.names)>1){
      l=mclapply(1:length(idx), mc.cores = opt$threads, mc.cleanup=T, mc.silent=F ,function(i) {
        a=unlist(d$control.allele.frequency.minus[idx[i]])
        idx1=as.numeric(which(a!=0 & a<opt$vaf))
        idx2=as.numeric(which(a!=0 & a<=1))
        if(length(idx1)>0){
          m=median(a[idx1]) #median
        } else if (length(idx1)==0 & length(control.names)>1) {
          m=length(d$total.depth.minus[idx[i],..control.names])/sum(d$total.depth.minus[idx[i],..control.names])
        } else if (length(idx1)==0 & length(control.names)==1) {
          m=length(d$total.depth.minus[idx[i],..control.names])/sum(d$total.depth.minus[idx[i],..control.names])
        }
        if(length(idx2)==0){
          set(copy(d$control.allele.frequency.minus[idx[i]]), i = NULL, j = as.integer(1:ncol(d$control.allele.frequency.minus)), m)
        } else {
          set(copy(d$control.allele.frequency.minus[idx[i]]), i = NULL, j = as.integer(idx2), m)
        }
      })
      set(d$control.allele.frequency.minus,idx,names(d$control.allele.frequency.minus),rbindlist(l))
    } else { #only one control sample
      a=unlist(d$control.allele.frequency.minus[idx])
      idx1=as.numeric(which(a!=0 & a<opt$vaf))
      if(length(idx1)>0){
        m=median(a[idx1]) #median
        d$control.allele.frequency.minus[idx]=m
      }
    }
  }
  cat(" DONE!\n")
  d
}

pval.calculation = function(d){
  control.names=d$control.names
  d$pval.plus=d$pval.minus=list()
  d$control.total.depth.plus=copy(d$total.depth.plus[,..control.names])
  d$control.total.depth.minus=copy(d$total.depth.minus[,..control.names])
  d$control.non.ref.counts.plus=copy(d$control.allele.frequency.plus*d$total.depth.plus[,..control.names])
  d$control.non.ref.counts.minus=copy(d$control.allele.frequency.minus*d$total.depth.minus[,..control.names])
  case.names=d$samples$id[d$samples$type=="case"]
  x=0
  for(n in case.names){
    if(length(unique(d$samples$type)=="case")==1 & "case" %in% unique(d$samples$type)){
      control.names=case.names
      if(!is.na(d$samples$replicate[d$samples$id==n])){ #If there are technical replicates in the experiment for that sample it removes the other replicate from the controls
        control.names=d$samples$id[-which(d$samples$replicate==d$samples$replicate[d$samples$id==n])]
      } else { #no replicates
        control.names=d$samples$id[d$samples$id!=n]
      }
    } else if(length(unique(d$samples$type))==2 & ("case" %in% d$samples$type) & ("control" %in% d$samples$type)=="TRUE"){
      control.names=d$samples$id[d$samples$type=="control"]
    }

    d$pval.plus[[n]]=NA
    d$pval.minus[[n]]=NA

    if(opt$binomial=="sum"){
      #plus
      d$control.non.ref.counts.plus[d$control.non.ref.counts.plus==0]=1
      ap=d$control.non.ref.counts.plus[,rowSums(.SD, na.rm=T),.SDcols=control.names]
      bp=d$total.depth.plus[,rowSums(.SD, na.rm=T),.SDcols=control.names]
      v=ap/bp
      d$pval.plus[[n]]=suppressWarnings(pbinom(d$non.ref.counts.plus[[n]], size=d$total.depth.plus[[n]], prob=v,lower.tail = F))
      idx1=which(d$non.ref.counts.plus[[n]]==0) # the p-value will be changed to 1 when there are no supporting reads
      d$pval.plus[[n]][idx1]=1

      #minus
      d$control.non.ref.counts.minus[d$control.non.ref.counts.minus==0]=1
      am=d$control.non.ref.counts.minus[,rowSums(.SD, na.rm=T),.SDcols=control.names]
      bm=d$total.depth.minus[,rowSums(.SD, na.rm=T),.SDcols=control.names]
      v=am/bm
      d$pval.minus[[n]]=suppressWarnings(pbinom(d$non.ref.counts.minus[[n]], size=d$total.depth.minus[[n]], prob=v,lower.tail = F))
      idx1=which(d$non.ref.counts.minus[[n]]==0) # pval change to 1 when there are no supporting reads
      d$pval.minus[[n]][idx1]=1
    } else if (opt$binomial=="max"){
      #plus
      d$control.non.ref.counts.plus[d$control.non.ref.counts.plus==0]=1
      d$control.allele.frequency.plus=d$control.non.ref.counts.plus/d$control.total.depth.plus
      v=d$control.allele.frequency.plus[, do.call(pmax, c(.SD,list(na.rm = TRUE))), .SDcols = control.names]
      d$pval.plus[[n]]=suppressWarnings(pbinom(d$non.ref.counts.plus[[n]], size=d$total.depth.plus[[n]], prob=v,lower.tail = F))
      idx1=which(d$non.ref.counts.plus[[n]]==0) # the p-value will be changed to 1 when there are no supporting reads
      d$pval.plus[[n]][idx1]=1

      #minus
      d$control.non.ref.counts.minus[d$control.non.ref.counts.minus==0]=1
      d$control.allele.frequency.minus=d$control.non.ref.counts.minus/d$control.total.depth.minus
      v=d$control.allele.frequency.minus[, do.call(pmax, c(.SD,list(na.rm = TRUE))), .SDcols = control.names]
      d$pval.minus[[n]]=suppressWarnings(pbinom(d$non.ref.counts.minus[[n]], size=d$total.depth.minus[[n]], prob=v,lower.tail = F))
      idx1=which(d$non.ref.counts.minus[[n]]==0) # pval change to 1 when there are no supporting reads
      d$pval.minus[[n]][idx1]=1
    }

    x=x+1
    if(round(x/nrow(d$samples[type=="case"]),2) %in% seq(0.01,0.99,0.01)){
      system(paste0("printf '\\rEstimating background error levels and calculating P-values :  ",round(100*x/nrow(d$samples[type=="case"])),"%%     '"))
    } else if (x==nrow(d$samples[type=="case"])){
      system(paste0("printf '\\rEstimating background error levels and calculating P-values :  100%%     '"))
    }

  }
  d$pval.plus=as.data.table(do.call("cbind", d$pval.plus))
  d$pval.minus=as.data.table(do.call("cbind", d$pval.minus))
  cat("\n")
  d
}

pval.correction.cdna.strand = function(d){
  system(paste0("printf '\\rAjusting P-values to account for the plus and minus replicated DNA strands...'"))
  idx=which(d$annotated.panel$alt!="+" & d$annotated.panel$alt!="-")
  l1=mclapply(names(d$pval.plus), mc.cores = opt$threads, mc.cleanup=T, mc.silent=F ,function(n) {
    idx1=which(d$total.depth.plus[idx,..n]<opt$overlap.coverage & d$total.depth.minus[idx,..n]>=opt$overlap.coverage)
    idx2=which(d$total.depth.plus[idx,..n]>=opt$overlap.coverage & d$total.depth.minus[idx,..n]<opt$overlap.coverage)
    p.plus=d$pval.plus[idx,..n]
    p.plus[idx1]=NA
    p.minus=d$pval.minus[idx,..n]
    p.minus[idx2]=NA
    p=pmax(p.plus,p.minus, na.rm = TRUE)
    p
  })
  #allowing detection of indels without validation on the other DNA strand
  l2=mclapply(names(d$pval.plus), mc.cores = opt$threads, mc.cleanup=T, mc.silent=F ,function(n) {
    p.plus=d$pval.plus[-idx,..n]
    p.minus=d$pval.minus[-idx,..n]
    p=pmin(p.plus,p.minus, na.rm = TRUE)
    p
  })
  d$pval=as.data.table(do.call("cbind", lapply(1:length(l1),function(i) rbind(l1[[i]],l2[[i]])))) #the indices for the indels always follows those for SNVs
  cat(" DONE!\n")
  d
}

pval.correction.overlapping.smmips = function(d){
  system("printf '\\rAjusting P-values to account for overlapping smMIPs...'")
  pval.before.correction=d$pval
  idx1=which(d$annotated.panel$alt=="+" | d$annotated.panel$alt=="-")

  coverage=as.data.table(cbind(paste(d$annotated.panel$chr,d$annotated.panel$pos,d$annotated.panel$ref,d$annotated.panel$alt),d$total.depth.minus+d$total.depth.plus))
  colnames(coverage)[1]="allele"
  idx2=which(duplicated(coverage$allele) | duplicated(coverage$allele, fromLast = TRUE))

  if(length(idx2)>0){
    p=as.data.table(cbind(paste(d$annotated.panel$chr,d$annotated.panel$pos,d$annotated.panel$ref,d$annotated.panel$alt),d$pval))
    colnames(p)[1]="allele"
    for(n in names(d$pval)){
      idx3=idx2[idx2 %in% which(coverage[,..n]<opt$overlap.coverage2)] #alleles with multiple smMIPs. The smMIP with the low coverage
      if(length(idx3)>0){
        ptmp=p[-idx3,c("allele",..n)]
        ptmp=ptmp[, lapply(.SD, max), by="allele"]
        idx4=match(p$allele,ptmp$allele)
        p[,n]=ptmp[idx4,..n]
      }
    }
    p=p[,"allele":=NULL]
    d$pval=p
  }
  d$pval[idx1]=pval.before.correction[idx1]  #allowing detection of indels without validation on other smMIPs
  cat(" DONE!\n")
  d
}



pval.correction.technical.replicates = function(d){
  system(paste0("printf '\\rAjusting P-values to account for technical replicates...'"))
  d$pval.before.replicate=d$pval
  for(i in unique(d$samples$replicate[!is.na(d$samples$replicate) & d$samples$type=="case"])){
    n1=names(d$pval)[grep(d$samples[replicate==i,id][1],names(d$pval))]
    n2=names(d$pval)[grep(d$samples[replicate==i,id][2],names(d$pval))]
    if(length(n1)==0 | length(n2)==0){print("Please check your sample desciption input file. If technical replicates are included, the file should indicate 2 replicates per each case sample")}

    idx1=which(d$total.depth.minus[,..n1]+d$total.depth.plus[,..n1]<opt$overlap.coverage2)
    idx2=which(d$total.depth.minus[,..n2]+d$total.depth.plus[,..n2]<opt$overlap.coverage2)
    d$pval[idx1,n1]=NA
    d$pval[idx2,n2]=NA
    m=pmax(d$pval[,..n1],d$pval[,..n2], na.rm = TRUE)
    d$pval[,n1]=m
    d$pval[,n2]=m
  }
  cat(" DONE!\n")
  d
}


pval.passed.mutations = function(d){
  system("printf '\\rCalling mutations. Applying user-define P-value cut-off...'")
  calls=mclapply(names(d$pval), mc.cores = opt$threads, mc.cleanup=T, mc.silent=F ,function(n){
    idx1=which(d$pval[,..n]<=opt$pval)
    mut=cbind(d$annotated.panel[idx1],d$non.ref.counts.plus[idx1,..n],d$non.ref.counts.minus[idx1,..n],d$total.depth.plus[idx1,..n],d$total.depth.minus[idx1,..n],
              d$allele.frequency.plus[idx1,..n],d$allele.frequency.minus[idx1,..n],d$pval[idx1,..n],d$pval[idx1,..n]*ps)
    names(mut)[names(mut)==n]=c("non.ref.counts.plus","non.ref.counts.minus","total.depth.plus","total.depth.minus",
                                "allele.frequency.plus","allele.frequency.minus","P-value","P-value.Bonferroni")
    mut$sample_ID=n
    mut$`P-value.Bonferroni`[mut$`P-value.Bonferroni`>1]=1
    mut
  })
  d$calls=rbindlist(calls)

  #Writing NA counts (not zero so it will work with the flag code) for indels where the P-value did not pass the cut-off. We allow indels to be called without validation on the other cDNA strand and we need to preserve the correct VAF
  idx1=grep("-|[+]",d$calls$alt)
  for(i in idx1){
    idx2=which(d$annotated.panel$pos==d$calls$pos[i] & d$annotated.panel$alt==d$calls$alt[i] & d$annotated.panel$smMIP==d$calls$smMIP[i])
    n=d$calls$sample_ID[i]
    if(d$pval.minus[idx2,..n]>opt$pval){
      d$calls$non.ref.counts.minus[i]=NA
      d$calls$total.depth.minus[i]=NA
    }
    if(d$pval.plus[idx2,..n]>opt$pval){
      d$calls$non.ref.counts.plus[i]=NA
      d$calls$total.depth.plus[i]=NA
    }
  }
  cat(" DONE!\n")
  d
}

adding.sscs.info = function(d){
  d$tab=as.data.frame(matrix(nrow=length(unique(paste(d$calls$chr,d$calls$pos,d$calls$alt,d$calls$smMIP))),ncol=nrow(d$samples)))
  row.names(d$tab)=unique(paste(d$calls$chr,d$calls$pos,d$calls$alt,d$calls$smMIP))
  colnames(d$tab)=d$samples$id
  x=0
  for(n in unlist(d$samples$id)){
    skip_to_next <- FALSE
    s=tryCatch(fread(paste0(opt$summary,"/",n,"_sscs_pileup.txt"),header=T,sep="\t",showProgress = FALSE),
        error = function(e) {
        skip_to_next <<- TRUE
       }
    )
    if(skip_to_next) {
      print(paste0("Skipping ",n,". SSCS pileup file not found. Please note that this will result with flag messages indicating that there is no SSCS for the alleles reported in this sample."),quote = F)
      next}

    idx1=which(d$calls$sample_ID==n)
    s.plus=s[match(paste(d$calls$smMIP[idx1],d$calls$chr[idx1],d$calls$pos[idx1],d$calls$alt[idx1],"+"),paste(s$smMIP,s$chr,s$pos,s$nucleotide,s$strand))]
    s.minus=s[match(paste(d$calls$smMIP[idx1],d$calls$chr[idx1],d$calls$pos[idx1],d$calls$alt[idx1],"-"),paste(s$smMIP,s$chr,s$pos,s$nucleotide,s$strand))]

    d$calls[idx1, c("SSCS.non.ref.counts.plus","SSCS.non.ref.counts.minus","SSCS.total.depth.plus","SSCS.total.depth.minus") := list(s.plus$count,s.minus$count,s.plus$coverage_at_position,s.minus$coverage_at_position)]
    d$calls[idx1, c("SSCS.allele.frequency.plus","SSCS.allele.frequency.minus") := list(SSCS.non.ref.counts.plus/SSCS.total.depth.plus,SSCS.non.ref.counts.minus/SSCS.total.depth.minus)]
    if(length(which(names(s.minus)=="UMI"))>0){ #if the user choose to include this information in the sscs pilups
      d$calls[idx1, c("SSCS.family.size.plus","SSCS.in.family.non.ref.vaf.plus","SSCS.UMIs.plus","SSCS.family.size.minus","SSCS.in.family.non.ref.vaf.minus","SSCS.UMIs.minus") := list(s.plus$family_sizes,s.plus$VAF_in_families,s.plus$UMI,s.minus$family_sizes,s.minus$VAF_in_families,s.minus$UMI)]
      s=s[paste(s$chr,s$pos,s$nucleotide,s$smMIP) %in% row.names(d$tab)]
      s=s[s[, .I[which.max(count)], by=list(chr,pos,nucleotide,smMIP)]$V1]
      d$tab[paste(s$chr,s$pos,s$nucleotide,s$smMIP),n]=s$UMI
    } else {
      d$calls[idx1, c("SSCS.family.size.plus","SSCS.in.family.non.ref.vaf.plus","SSCS.family.size.minus","SSCS.in.family.non.ref.vaf.minus") := list(s.plus$family_sizes,s.plus$VAF_in_families,s.minus$family_sizes,s.minus$VAF_in_families)]
    }

    x=x+1
    if(round(x/nrow(d$samples),2) %in% seq(0.01,0.99,0.01)){
      system(paste0("printf '\\rAdding consensus read information :  ",round(100*x/nrow(d$samples)),"%%     '"))
    } else if (x==nrow(d$samples)){
      system(paste0("printf '\\rAdding consensus read information :  100%%     '"))
    }
  }
  d$calls[is.na(SSCS.non.ref.counts.minus),SSCS.non.ref.counts.minus:=0]
  d$calls[is.na(SSCS.non.ref.counts.plus),SSCS.non.ref.counts.plus:=0]
  d$calls[is.na(SSCS.total.depth.minus),SSCS.total.depth.minus:=0]
  d$calls[is.na(SSCS.total.depth.plus),SSCS.total.depth.plus:=0]
  d$calls[is.na(SSCS.allele.frequency.minus),SSCS.allele.frequency.minus:=0]
  d$calls[is.na(SSCS.allele.frequency.plus),SSCS.allele.frequency.plus:=0]
  cat("\n")
  d
}

vaf.calculation.dna.strand = function(d){
  d$calls$flags=NA

  idx1=which(d$calls$total.depth.plus<opt$overlap.coverage | d$calls$total.depth.minus<opt$overlap.coverage)
  idx2=which(d$calls$total.depth.plus<opt$overlap.coverage & d$calls$total.depth.minus<opt$overlap.coverage)
  idx2=idx1[!(idx1 %in% idx2)]
  d$calls$flags[idx2]=paste0("Cannot be supported by both Read1 and Read2 generated from smMIP: ",d$calls$smMIP[idx2],". (Low coverage in sample: ",d$calls$sam[idx2],")")

  d$calls[,c("non.ref.counts","total.depth") := list(sum(non.ref.counts.plus,non.ref.counts.minus,na.rm = T),sum(total.depth.plus,total.depth.minus,na.rm = T)),by=list(smMIP,chr,pos,ref,alt,sample_ID)]
  d$calls[,c("non.ref.counts.plus","non.ref.counts.minus","total.depth.plus","total.depth.minus","allele.frequency.plus","allele.frequency.minus","allele.frequency") := list(NULL,NULL,NULL,NULL,NULL,NULL,non.ref.counts/total.depth)]

  if(length(list.files(opt$summary,pattern="_sscs_pileup.txt"))>0){
    if(length(grep("UMI",names(d$calls))>0)){ #if the user choose to include this information in the sscs pilups
      a=d$calls[,SSCS.non.ref.counts.plus]
      a[is.na(a)]=0
      b=d$calls[,SSCS.non.ref.counts.minus]
      b[is.na(b)]=0
      idx=a>=b
      idx[idx==TRUE]=1
      idx[idx==FALSE]=2
      idx=as.matrix(data.frame(row=1:length(idx),col=idx))
      a=cbind(d$calls$SSCS.UMIs.plus,d$calls$SSCS.UMIs.minus)
      d$calls$SSCS.UMIs=a[idx]
      d$calls[,c("SSCS.UMIs.plus","SSCS.UMIs.minus") := list(NULL,NULL)]
    }

    d$calls[,c("SSCS.non.ref.counts","SSCS.total.depth") := list(sum(SSCS.non.ref.counts.plus,SSCS.non.ref.counts.minus),sum(SSCS.total.depth.plus,SSCS.total.depth.minus)),by=list(smMIP,chr,pos,ref,alt,sample_ID)]
    d$calls[,c("SSCS.family.size","SSCS.in.family.non.ref.vaf") := list(apply(cbind(d$calls$SSCS.family.size.minus,d$calls$SSCS.family.size.plus),1, function(x) paste(x[1],x[2], sep=":D:")),
                                                                        apply(cbind(d$calls$SSCS.in.family.non.ref.vaf.minus,d$calls$SSCS.in.family.non.ref.vaf.plus),1, function(x) paste(x[1],x[2], sep=":D:")))]
    d$calls[,c("SSCS.non.ref.counts.plus","SSCS.non.ref.counts.minus","SSCS.total.depth.plus","SSCS.total.depth.minus","SSCS.allele.frequency.minus","SSCS.allele.frequency.plus","SSCS.family.size.minus","SSCS.family.size.plus","SSCS.in.family.non.ref.vaf.plus","SSCS.in.family.non.ref.vaf.minus","SSCS.allele.frequency") := list(NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,SSCS.non.ref.counts/SSCS.total.depth)]
  }
  d
}

adding.batch.info1 = function(d){
  system(paste0("printf '\\rWriting batch related information (part1)...'"))
  d$calls$num.pval.pass=NA
  idx1=match(paste(d$calls$smMIP,d$calls$chr,d$calls$pos,d$calls$ref,d$calls$alt),
             paste(d$annotated.panel$smMIP,d$annotated.panel$chr,d$annotated.panel$pos,d$annotated.panel$ref,d$annotated.panel$alt))
  if(length(which(!is.na(d$samples$replicate)))>0){ #experiment with replicates
    d$calls$num.pval.pass=rowSums(d$pval.before.replicate[idx1]<=opt$pval,na.rm = T)
    d$calls[,pass.pval.in.pairs:=.N/2, by=c("smMIP","chr","pos","ref","alt")]
  } else {
    d$calls$num.pval.pass=rowSums(d$pval[idx1]<=opt$pval)
  }
  cat(" DONE!\n")
  d
}

read.sample.misplacement = function(d){
  d$calls$SSCS.overlap=0
  m=mclapply(1:nrow(d$calls),mc.cores = opt$threads, mc.cleanup=T, mc.silent=F ,function(i){
    l=unlist(d$tab[paste(d$calls$chr[i],d$calls$pos[i],d$calls$alt[i],d$calls$smMIP[i]),-which(names(d$tab)==d$calls$sample_ID[i])])
    l=l[!is.na(l)]
    l=unlist(strsplit(paste(l),","))
    a=table(unlist(strsplit(d$calls$SSCS.UMIs[i],",")) %in% l)

    r=as.numeric(a[names(a)==TRUE])
    if(length(r)>0){
      r=round(100*(as.numeric(a[names(a)==TRUE])/sum(a)),  2)
    } else {
      r=0
    }
    if(round(i/nrow(d$calls),1) %in% seq(0.1,0.9,0.1)) {
      system(paste0("printf '\\rAdding UMI information concerning potential index-hopping :  ",round(100*i/nrow(d$calls)),"%%     '"))
    } else if (i==nrow(d$calls)){
      system(paste0("printf '\\rAdding UMI information concerning potential index-hopping :  100%%     '"))
    }
    r
  }
  )
  system(paste0("printf '\\rAdding UMI information concerning potential index-hopping :  100%%     '"))
  d$calls$SSCS.overlap=unlist(m)
  d$calls[,"SSCS.UMIs":=NULL]
  cat("\n")
  d
}

vaf.calculation.overlapping.smmips = function(d){
  idx1=which(duplicated(d$calls[,list(chr,pos,ref,alt,sample_ID)]) | duplicated(d$calls[,list(chr,pos,ref,alt,sample_ID)], fromLast = TRUE))
  idx2=which(d$calls$total.depth[idx1] < opt$overlap.coverage2)
  idx3=idx1[idx2]
  if(length(idx3)>0){
    d$calls$flags[idx3]=paste(d$calls$flags[idx3], paste0("Cannot be supported by all overlapping smMIPs. (Low coverage in sample:",d$calls$sample_ID[idx3],", smMIP:",d$calls$smMIP[idx3],")"),sep=", ")
  }

  d$calls$num.pval.pass=as.character(d$calls$num.pval.pass)
  if(length(which(!is.na(d$samples$replicate)))>0){ #experiment with replicates
    d$calls$pass.pval.in.pairs=as.character(d$calls$pass.pval.in.pairs)
    d$calls[,c("smMIP","non.ref.counts","total.depth","num.pval.pass","pass.pval.in.pairs","P-value","P-value.Bonferroni","flags") := list(paste(smMIP, collapse=","),sum(non.ref.counts),sum(total.depth), paste(num.pval.pass, collapse=","),paste(pass.pval.in.pairs,collapse=","),min(`P-value`),min(`P-value.Bonferroni`),paste(flags, collapse = ", ")) , by=list(chr,pos,ref,alt,sample_ID)]
  } else {
    d$calls[,c("smMIP","non.ref.counts","total.depth","num.pval.pass","P-value","P-value.Bonferroni","flags") := list(paste(smMIP, collapse=","),sum(non.ref.counts),sum(total.depth),paste(num.pval.pass, collapse=","),min(`P-value`),min(`P-value.Bonferroni`),paste(flags, collapse = ", ")), by=list(chr,pos,ref,alt,sample_ID)]
  }
  d$calls[,"allele.frequency" := list(non.ref.counts/total.depth)]

  if(length(list.files(opt$summary,pattern="_sscs_pileup.txt"))>0){
    d$calls$SSCS.overlap=as.character(d$calls$SSCS.overlap)
    d$calls[,c("SSCS.non.ref.counts","SSCS.total.depth","SSCS.overlap") := list(sum(SSCS.non.ref.counts),sum(SSCS.total.depth),paste(SSCS.overlap,collapse=",")),by=list(chr,pos,ref,alt,sample_ID)]
    d$calls[,"SSCS.allele.frequency" := list(SSCS.non.ref.counts/SSCS.total.depth)]
    d$calls[,c("SSCS.family.size","SSCS.in.family.non.ref.vaf") :=
              list(paste(SSCS.family.size, collapse=":S:"),
                   paste(SSCS.in.family.non.ref.vaf, collapse=":S:")),
            by=list(chr,pos,ref,alt,sample_ID)]
  }
  d$calls=unique(d$calls)
  d
}

vaf.calculation.technical.replicates = function(d){
  idx1=which(d$calls$total.depth < opt$overlap.coverage2)
  if(length(idx1)>0){
    d$calls$flags[idx1]=paste(d$calls$flags[idx1], paste0("Cannot be supported by both technical replicates. (Low coverage in sample:",d$calls$sample_ID[idx1],")"),sep=", ")
  }

  for(i in unique(d$samples$replicate[!is.na(d$samples$replicate) & d$samples$type=="case"])){
    n1=d$samples[replicate==i,id][1]
    n2=d$samples[replicate==i,id][2]
    d$calls[sample_ID==n1,sample_ID := paste0(n1,",",n2)]
    d$calls[sample_ID==n2,sample_ID := paste0(n1,",",n2)]
  }

  d$calls[,c("non.ref.counts","total.depth","flags") := list(sum(non.ref.counts),sum(total.depth),paste(flags, collapse = ", ")),by=list(chr,pos,ref,alt,sample_ID)]
  d$calls[,c("allele.frequency") := list(non.ref.counts/total.depth)]

  if(length(list.files(opt$summary,pattern="_sscs_pileup.txt"))>0){
    d$calls$SSCS.overlap=as.character(d$calls$SSCS.overlap)
    d$calls[,c("SSCS.non.ref.counts","SSCS.total.depth","SSCS.overlap","SSCS.family.size","SSCS.in.family.non.ref.vaf") :=
              list(sum(SSCS.non.ref.counts),sum(SSCS.total.depth),paste(SSCS.overlap,collapse = ","),
                   paste(SSCS.family.size, collapse=":R:"),
                   paste(SSCS.in.family.non.ref.vaf, collapse=":R:")),
            by=list(chr,pos,ref,alt,sample_ID)]
    d$calls[,"SSCS.allele.frequency" := SSCS.non.ref.counts/SSCS.total.depth]

  }
  d$calls=unique(d$calls)
  d
}

adding.batch.info2 = function(d){
  #Writing zero counts for indels where the P-values did not pass as we allowed indels to be called without validation on the other cDNA strand and we need to preserve the correct VAF
  idx1=grep("-|[+]",d$annotated.panel$alt)

  idx2=which(d$pval.minus[idx1]>opt$pval,arr.ind = T)
  idx2[,1]=idx1[idx2[,1]]
  non.ref.counts.minus=d$non.ref.counts.minus
  non.ref.counts.minus[idx2]=0
  total.depth.minus=d$total.depth.minus
  total.depth.minus[idx2]=0

  idx2=which(d$pval.plus[idx1]>opt$pval,arr.ind = T)
  idx2[,1]=idx1[idx2[,1]]
  non.ref.counts.plus=d$non.ref.counts.plus
  non.ref.counts.plus[idx2]=0
  total.depth.plus=d$total.depth.plus
  total.depth.plus[idx2]=0

  d$non.ref.counts=non.ref.counts.minus+non.ref.counts.plus
  d$total.depth=total.depth.minus+total.depth.plus

  #minimizing the data.tables to those alleles that were called
  rn=paste(d$annotated.panel$chr,d$annotated.panel$pos,d$annotated.panel$ref,d$annotated.panel$alt)
  rn=gsub("[+]","[+]",rn)
  mn=paste(d$calls$chr,d$calls$pos,d$calls$ref,d$calls$alt)
  mn=gsub("[+]","[+]",mn)
  idx1=which(rn %in% mn)
  d$non.ref.counts=d$non.ref.counts[idx1]
  d$total.depth=d$total.depth[idx1]
  #sum the read counts based on overlapping smmips
  n=d$samples$id[d$samples$type=="case"]
  d$non.ref.counts$by=rn[idx1]
  d$non.ref.counts=d$non.ref.counts[,lapply(.SD,sum),.SDcols=n, by=by]
  d$total.depth$by=rn[idx1]
  d$total.depth=d$total.depth[,lapply(.SD,sum),.SDcols=n, by=by]
  #sum the read counts based on technical replicates
  if(length(which(!is.na(d$samples$replicate) & d$samples$type=="case"))>0){ #experiment with replicates
    n=d$samples[!is.na(replicate) & d$samples$type=="case"]
    n1=unlist(lapply(unique(n$replicate),function(x) which(n$replicate==x)))
    n2=n1
    n1=n1[seq(1,length(n1),2)]+1
    n2=n2[seq(2,length(n2),2)]+1

    x1=d$non.ref.counts[,..n1]
    x2=d$non.ref.counts[,..n2]
    d$non.ref.counts=x1+x2
    row.names(d$non.ref.counts)=d$total.depth$by
    names(d$non.ref.counts)=paste(names(x1),names(x2),sep=",")

    x1=d$total.depth[,..n1]
    x2=d$total.depth[,..n2]
    d$total.depth=x1+x2
  } else {
	d$non.ref.counts[,by:=NULL]
        row.names(d$non.ref.counts)=d$total.depth$by
        d$total.depth[,by:=NULL]
  }

  d$allele.frequency=d$non.ref.counts/d$total.depth
  row.names(d$allele.frequency)=row.names(d$non.ref.counts)
  names(d$allele.frequency)=names(d$non.ref.counts)

  idx1=match(mn,row.names(d$allele.frequency))
  m=unlist(mclapply(1:nrow(d$calls) ,mc.cores = opt$threads, mc.cleanup=T, mc.silent=F ,function(i){
    idx2=which(d$allele.frequency[idx1[i]] > d$calls$allele.frequency[i])
    idx3=which(d$allele.frequency[idx1[i]] < d$calls$allele.frequency[i])
    x1=length(idx2)
    if(x1>0){
      x2=d$allele.frequency[idx1[i],..idx2]
      x2=round(as.numeric(unlist(x2)),digits=6)
      x2=x2[order(x2,decreasing = F)]
      x2=x2[1:3]
      x2=paste(x2[!is.na(x2)],collapse=",")
    } else {x2=NA}
    if(length(idx3)>0){
      x3=d$allele.frequency[idx1[i],..idx3]
      x3=round(as.numeric(unlist(x3)),digits=6)
      x3=x3[order(x3,decreasing = T)]
      x3=x3[1:3]
      x3=paste(x3[!is.na(x3)],collapse=",")
    } else {x3=NA}

    if(round(i/nrow(d$calls),1) %in% seq(0.1,0.9,0.1)) {
      system(paste0("printf '\\rWriting batch related information (part2) :  ",round(100*i/nrow(d$calls)),"%%     '"))
    } else if (i==nrow(d$calls)){
      system(paste0("printf '\\rWriting batch related information (part2) :  100%%     '"))
    }
    paste(x1,x2,x3)
  }))
  system(paste0("printf '\\rWriting batch related information (part2) :  100%%     '"))
  m=as.data.table(matrix(unlist(strsplit(m,split = " ")), ncol = 3, byrow = TRUE))
  d$calls$samples.with.higher.vaf=m$V1
  d$calls$higher.vaf=m$V2
  d$calls$lower.vaf=m$V3
  cat("\n")
  d
}

additional.flags = function(d){
  idx1=grep("synonymous_variant|intron_variant|splice_region_variant|non_coding_transcript_exon_variant|5_prime_UTR_variant|3_prime_UTR_variant",d$calls$variant_type)
  idx2=grep("A|G|C|T",d$calls$alt)
  if(length(idx1)>0){
    idx1=intersect(idx1,idx2)
    if(length(idx1)>0){d$calls[idx1, flags:=paste(flags,"Likely benign (based on genomic loci)",sep = ", ")]}
  }
  idx1=which(d$calls$maf>0.01)
  if(length(idx1)>0){d$calls[idx1, flags:=paste(flags,"Common SNP",sep = ", ")]}
  idx2=which(d$calls$maf>opt$maf & is.na(d$calls$cosmic))
  idx2=idx2[!(idx2 %in% idx1)]
  if(length(idx2)>0){d$calls[idx2, flags:=paste(flags,"Potentially germline based on user's MAF cut-off and lack of COSMIC ID",sep = ", ")]}

  m=mclapply(1:nrow(d$calls),mc.cores = opt$threads, mc.cleanup=T, mc.silent=F ,function(i){
    x=unlist(strsplit(d$calls$SSCS.family.size[i],":R:"))
    idx1=grep("NA:D|D:NA",x)
    if(length(idx1)>0){
      idx2=grep("NA:D:NA",x)
    }
    idx3=intersect(idx1,idx2)
    idx1=idx1[!(idx1 %in% idx2)]

    if(length(idx3)>0){
      if(length(which(!is.na(data$samples$replicate)))>0){ #experiment with replicates
        a="No SSCS support in one of the replicates, "
      } else (a="No SSCS support at all, ")
    } else if (length(idx1)>0) {
      if(length(which(!is.na(data$samples$replicate)))>0){ #experiment with replicates
        a="No SSCS support in either Read1 or Read2 in at least one of the replicates, "
      } else {
        a="No SSCS support in either Read1 or Read2, "
      }

    } else {a=""}
    a
  })
  d$calls$flags=paste0(d$calls$flags,", ",unlist(m))
  d$calls[SSCS.allele.frequency=="NaN",SSCS.allele.frequency:=0]
  idx=intersect(which(d$calls$SSCS.non.ref.counts==0),grep("No SSCS support in one of the replicates",d$calls$flags))
  d$calls[idx, flags:=gsub("No SSCS support in one of the replicates","No SSCS support at all",flags)]

  a=strsplit(d$calls$lower.vaf,split = ",")
  x=unlist(lapply(a,length))
  w=c()
  for(i in 1:length(x)){
    if(x[i]==3 | x[i]==2){
      w=append(w,unlist(a[i])[2])
    } else if (x[i]==1){
      w=append(w,unlist(a[i])[1])
    }
  }
  options(warn=-1)
  d$calls$drop=d$calls$allele.frequency/as.numeric(w)
  options(warn=0)
  d$calls$drop[is.na(d$calls$drop)]=Inf

  u=unique(paste(d$calls$chr,d$calls$pos,d$calls$ref,d$calls$alt))
  for(n in u){
    idx1=which(paste(d$calls$chr,d$calls$pos,d$calls$ref,d$calls$alt)==n)
    if(length(idx1)==1){
      if(!(d$calls$drop[idx1]>=opt$drop & d$calls$samples.with.higher.vaf[idx1]==0)){
        d$calls[idx1, flags:=paste0(flags,"VAF Warning, ")]}
    } else if (length(idx1)>1){
      idx2=which(d$calls$drop[idx1]>=opt$drop & as.numeric(d$calls$samples.with.higher.vaf[idx1]) %in% 0:length(idx1))
      if(length(idx2)==0){
        d$calls[idx1, flags:=paste0(flags,"VAF Warning, ")]
      } else {
        idx3=which(as.numeric(d$calls$samples.with.higher.vaf[idx1])<=max(as.numeric(d$calls$samples.with.higher.vaf[idx1[idx2]])))
        idx1=idx1[-idx3]
        d$calls[idx1, flags:=paste0(flags,"VAF Warning, ")]
      }
    }
  }

  idx1=which((d$calls$alt=="+" | d$calls$alt=="-") & (d$calls$variant_type=="intron_variant" | d$calls$variant_type=="splice_region_variant"))
  if(length(idx1>0)){
    d$calls[idx1, flags:=paste0(flags,"Intronic indel, ")]
  }

  if("SSCS.overlap" %in% names(d$calls)){
    idx1=grep("100",d$calls$SSCS.overlap)
    if(length(idx1)>0){
      options(warn=-1)
      m=mclapply(idx1,mc.cores = opt$threads, mc.cleanup=T, mc.silent=F ,function(i){
        length(which(as.numeric(unlist(strsplit(d$calls$SSCS.family[i],split = ":|,")))>1))
      })

      options(warn=0)
      d$calls[idx1[which(unlist(m)==0)], flags:=paste0(flags,"Potential index hopping")]
    }
  }

  d$calls$flags=gsub("NA[,][ ]","", d$calls$flags)

  d$calls$flags=gsub("^, Potential index hopping","Potential index hopping", d$calls$flags)
  d
}


########################################################################################
############## Source functions for mutation_categories.R
########################################################################################

categorize=function(){
  #Enriching the mutation list for somatic and pathogenic mutations
  enriched=calls[-grep("Common SNP", calls$flags)]
  enriched=enriched[-grep("Potentially germline", enriched$flags)]
  enriched=enriched[-grep("Likely benign", enriched$flags)]
  enriched=enriched[-grep("Intronic indel", enriched$flags)]

  enriched_indels=enriched[alt=="-" | alt=="+"]
  enriched_snvs=enriched[alt!="-" & alt!="+"]

  #filtering indels if there is a >2 fold drop in the proximity (same smMIP) of another indel with an higher vaf
  #Please Notice: if there are two clones, defined by two slightly different indels, detected with the same smMIP and one is smaller than the other by 2 or more fold it will be removed
  idx=c()
  u=unique(c(unique(paste(enriched_indels$sample_ID,gsub(".*,","",enriched_indels$smMIP),enriched_indels$alt)),
             unique(paste(enriched_indels$sample_ID,gsub(",.*","",enriched_indels$smMIP),enriched_indels$alt))))
  for(n in u){
    idx1=unique(c(which(paste(enriched_indels$sample_ID,gsub(".*,","",enriched_indels$smMIP),enriched_indels$alt)==n),
                  which(paste(enriched_indels$sample_ID,gsub(",.*","",enriched_indels$smMIP),enriched_indels$alt)==n)))
    if(length(idx1>1)){
      idx2=which(enriched_indels$allele.frequency[idx1]/max(enriched_indels$allele.frequency[idx1])<0.5)
      if(length(idx2>1)){
        idx=append(idx,idx1[idx2])
      }
    }
  }
  if(length(idx)>0){enriched_indels=enriched_indels[-unique(idx)]}

  #High Confidence calls
  high.conf.calls.snvs=enriched_snvs[-grep("Cannot be supported by both Read1 and Read2|Cannot be supported by both technical replicates|No SSCS support in one of the replicates|No SSCS support at all|No SSCS support in either Read1 or Read2 in at least one of the replicates|Potential index hopping|VAF Warning", enriched_snvs$flags)]
  high.conf.calls.indels=enriched_indels[-grep("Cannot be supported by both Read1 and Read2|Cannot be supported by both technical replicates|No SSCS support in one of the replicates|No SSCS support at all|No SSCS support in either Read1 or Read2 in at least one of the replicates|Potential index hopping|VAF Warning", enriched_indels$flags)]

  #Lower Confidence calls
  idx1=grep("Potential index hopping|VAF Warning", enriched_snvs$flags)
  idx2=grep("Low coverage|No SSCS support in either Read1 or Read2 in at least one of the replicates", enriched_snvs$flags)
  lower.conf.calls.snvs=enriched_snvs[idx2[!(idx2 %in% idx1)]]
  m=unlist(lapply(1:nrow(lower.conf.calls.snvs), function(i) {
    s=unlist(strsplit(lower.conf.calls.snvs$SSCS.family.size[i],":"))
    l1=length(grep("NA",s))
    l2=length(which(s!="D" & s!="R" & s!="S"))
    (l2-l1)/(l2)
  }))
  lower.conf.calls.snvs=lower.conf.calls.snvs[which(m>=0.5)]

  idx1=grep("Potential index hopping|VAF Warning", enriched_indels$flags)
  idx2=grep("Low coverage|No SSCS support in either Read1 or Read2 in at least one of the replicates", enriched_indels$flags)
  lower.conf.calls.indels=enriched_indels[idx2[!(idx2 %in% idx1)]]
  m=unlist(lapply(1:nrow(lower.conf.calls.indels), function(i) {
    s=unlist(strsplit(lower.conf.calls.indels$SSCS.family.size[i],":"))
    l1=length(grep("NA",s))
    l2=length(which(s!="D" & s!="R" & s!="S"))
    (l2-l1)/(l2)
  }))
  lower.conf.calls.indels=lower.conf.calls.indels[which(m>=0.5)]

  #writing seperate tables
  cat(paste("Writing seperate mutation files to",dirname(opt$input),"\n"))
  write.table(high.conf.calls.snvs,file=paste0(dirname(opt$input),"/High_Confidence_SNVs.txt"),col.names = T,row.names = F,quote = F,sep = '\t')
  write.table(high.conf.calls.indels,file=paste0(dirname(opt$input),"/High_Confidence_Indels.txt"),col.names = T,row.names = F,quote = F,sep = '\t')
  write.table(lower.conf.calls.snvs,file=paste0(dirname(opt$input),"/Lower_Confidence_SNVs.txt"),col.names = T,row.names = F,quote = F,sep = '\t')
  write.table(lower.conf.calls.indels,file=paste0(dirname(opt$input),"/Lower_Confidence_Indels.txt"),col.names = T,row.names = F,quote = F,sep = '\t')
  cat("DONE\n")
}

