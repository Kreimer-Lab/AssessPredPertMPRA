.libPaths(c(.libPaths(), "XXXX"))

setwd("XXXX")

library(dplyr)     
library(reshape2)  # melt
library(gridExtra) # arrangeGrob
library(ggpubr)    # as_ggplot
library(cowplot)   # draw_plot_label
library(ggplot2)
library(Matrix, lib.loc = "/XXXX")
#dyn.load('/home/jl2791/miniconda3/lib/libxml2.so.2'); dyn.load('/home/jl2791/miniconda3/lib/libglpk.so.40')
library(ggsci)
library(magrittr)
library(data.table)
library(parallel)
library(pbapply)
library(future)

datetag = gsub("\\-", "", Sys.Date())

# import raw data -----
## reviewer file -----
reviewer.files.names = list.files("../00.raw/", pattern = "Reviewer1")
reviewer.files_list = pblapply(reviewer.files.names, function(X)
  read.delim(
    paste0("../00.raw/", X),
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) %>% slice(-1))
# fimo files -----
fimo_encode = read.delim("../00.raw/fimo_Pert_MPRA_ENCODE.txt", check.names = FALSE, stringsAsFactors = FALSE) %>% mutate(dataset = "ENCODE")
fimo_SOX = read.delim("../00.raw/fimo_Pert_MPRA_SOX.txt", check.names = FALSE, stringsAsFactors = FALSE) %>% mutate(dataset = "SOX")
fimo_hg19 = read.delim("../00.raw/fimo_Pert_MPRA_hg19.txt", check.names = FALSE, stringsAsFactors = FALSE) %>% mutate(dataset = "hg19")
## design names  -----
design.names = read.delim("../00.raw/regions_names.txt", header = FALSE)

# data cleaning  -----
## fimo  -----
fimo = rbind(fimo_encode, fimo_SOX, fimo_hg19)
fimo$regions = stringr::str_extract(fimo$sequence_name, "chr[0-9]+\\:[0-9]+\\-[0-9]+")
fimo$regions[is.na(fimo$regions)] = stringr::str_extract(fimo$sequence_name[is.na(fimo$regions)], "chrX\\:[0-9]+\\-[0-9]+")
fimo$successful = "unsuccessful"
## fimo_WT -----
fimo_WT = fimo[grep("WT", fimo$sequence_name), ]
## design name  -----
design.names = data.frame(identifier = 1:nrow(design.names), sequence_name = design.names$V1, 
                          stringsAsFactors = FALSE, check.names = FALSE)

## clean reviewer files list  -----
reviewer.files_list = pblapply(reviewer.files_list, function(X) 
  right_join(design.names, X, by = "identifier"))
## reduce list into df
reviewer.files_df = 
  lapply(1:3, function(X){
    reviewer.files_list[[X]] %>% mutate(pertX = paste0("pert", X))
  })
## reviewer files list to df
reviewer.files_df = Reduce(rbind, reviewer.files_df)
#reviewer.files_df$pertX = stringr::str_extract(reviewer.files_df$sequence_name, "pert[0-9]+")

X = reviewer.files_df#[1:100, ]
### successful or not: merge reviewer to fimo  -----
print("successful or not")
X = merge(fimo %>% select("# motif_id", "sequence_name", "start", "stop", "strand", "successful"), 
          X, 
          by.x = c("# motif_id", "sequence_name", "start", "stop", "strand"), 
          by.y = c("TFBS", "sequence_name", "start", "end", "strand"), 
          all.y = T)
X$successful[X$identifier == 0] = "N/A"
X$successful[is.na(X$successful)] = "successful"

### "Found in fimo, but in another position"  -----
print( "Found in fimo, but in another position")
    # whether successful is found in fimo by motif and sequence name
    ## split  X.successful and X.unsuccessful
    X.successful = X[X$successful == "successful", ]
    X.unsuccessful = X[X$successful != "successful", ]

    ## found in fimo or not? 
    X.successful_motif.sequence = 
      mapply(function(A, B) paste0(A, ".", B), 
             X.successful$`# motif_id`, X.successful$sequence_name)
    fimo_motif.sequence = 
      mapply(function(A, B) paste0(A, ".", B), 
             fimo$`# motif_id`, fimo$sequence_name)

    X.successful$found_in_fimo_other.position = X.successful_motif.sequence %in% fimo_motif.sequence
    
    ## rbind back successful and unsuccessfu
    X = plyr::rbind.fill(X.successful, X.unsuccessful)

gc()

X.df = X
X = setDT(X)

### "overlapped" -- For the found in fimo, whether it's overlapped with the pertX motif -----
print("overlapped")
# for "found in fimo", if overlapped?
## split X.found and X.notfound
### get ID
found.ID = which(X$found_in_fimo_other.position)
notfound.ID = which(!X$found_in_fimo_other.position |
                      is.na(X$found_in_fimo_other.position))
### get split Xs
X.found <- X[found.ID,]
X.notfound <- X[notfound.ID,]

## overlap or not in X.found
## describe a function to calculate if overlap or not in X
overlappedOrnot = 
               function(sequence_name, MotifID, start){
                 #colnum.sequence_name = which(names(vector) == "sequence_name")
                 #colnum.Motifid = which(names(vector) == "# motif_id") 
                 #colnum.start = which(names(vector) == "start")
                 
                 SequenceName = sequence_name
                 MotifID = MotifID
                 vectorStart = start
                 
                 fimo_Motif_SequenceName = fimo[grep(MotifID, fimo$`# motif_id`),] %>% .[grep(SequenceName, .$sequence_name),]
                 fimo_Motif_SequenceName_Start = fimo_Motif_SequenceName$start %>% as.numeric
                 fimo_Motif_SequenceName_Stop = fimo_Motif_SequenceName$stop %>% as.numeric
                 
                 #overlapped = vectorStart %in% c(fimo_Motif_SequenceName_Start:fimo_Motif_SequenceName_Stop)
                 #output = ifelse(overlapped, "overlapped", "not overlapped")
                 
                 n.overlapped = 
                   mapply(function(Start, FimoStart, FimoStop){Start %in% FimoStart:FimoStop}, 
                          vectorStart, fimo_Motif_SequenceName_Start, fimo_Motif_SequenceName_Stop) %>% sum
                 
                 return(
                   c(n_Identical.FIMO.Motifs = length(fimo_Motif_SequenceName),
                     n_Identical.FIMO.Motifs_overlapped = n.overlapped,
                     n_Identical.FIMO.Motifs_not.overlapped = length(fimo_Motif_SequenceName) - n.overlapped
                   )
                 )
               }
             v_overlappedOrnot <- Vectorize(overlappedOrnot)  
             setDTthreads(15)
             system.time(
               {found_overlapped.or.not <- X.found[,v_overlappedOrnot(sequence_name, `# motif_id`, start)]
               found_overlapped.or.not <- t(found_overlapped.or.not)}
             )
    if (all.equal(rownames(found_overlapped.or.not), X.found$sequence_name)) {
      X.found = cbind(X.found, found_overlapped.or.not)
    } else {
      print("rownames of found_overlapped.or.not and X.found do not match")
    }
             ## rbind back X.found and X.notfound
             X = plyr::rbind.fill(X.found, X.notfound)
             X = setDT(X)

gc()


### "WT motifs information"  -----
print("WT motifs information")
checkWT <- 
  function(sequence_name, MotifID, regions, start, stop){
    print(sequence_name) 
    
    MotifID = MotifID
    regions = regions
    sequence_name = sequence_name
    start.i = start
    stop.i = stop
    
    # If >1 Motifs x in WT 
    fimo_WT.i = fimo_WT[fimo_WT$`# motif_id` == MotifID &
                          fimo_WT$regions == regions, ]
    ### output
    n_Identical.WT.MotifIDs = nrow(fimo_WT.i)
    more_Identical.WT.MotifIDs = n_Identical.WT.MotifIDs > 1
    
    ## get index of duplicated Motifs x
    dup.index.i = which(
      fimo_WT$`# motif_id` == MotifID &
        fimo_WT$regions == regions & fimo_WT$start != start.i
    )
    ### output
    n_Identical.duplicated.WT.MotifIDs = length(dup.index.i)
    dup.index.i_char = as.character(dup.index.i) 
    dup.index.i_char = paste0(dup.index.i, "-")
    dup.index.i_char <-  Reduce("paste", dup.index.i)
    IDs_Identical.WT.MotifIDs = ifelse(!is.null(dup.index.i), dup.index.i, NA)
    
    ## are Motifs x still in pert?
    ### find motifs after perturbation
    fimo_pert = fimo[fimo$sequence_name == sequence_name,]
    fimo_pert.motifs = fimo_pert$`# motif_id`
    ## output: does fimo_pert.motifs contain Motif x?
    n_Identical.PertX.MotifIDs.in.sequence = sum(MotifID == fimo_pert.motifs)
    contain_Identical.PertX.MotifIDs.in.sequence = n_Identical.PertX.MotifIDs.in.sequence >= 1
    ## how many of these Motif(s) x introduced or originally there?
    ### get start and stop of MotifID in pert, and in Fimo
    fimo_pert.MotifID = fimo_pert[fimo_pert$`# motif_id` == MotifID, ]
    pert.MotifID_startstop = 
      mapply(function(B, C) paste0(B, ".", C), 
             fimo_pert.MotifID$start, fimo_pert.MotifID$stop)
    fimo_WT.MotifID = fimo_WT[dup.index.i, ]
    WT.MotifID_startstop = 
      mapply(function(B, C) paste0(B, ".", C), 
             fimo_WT.MotifID$start, fimo_WT.MotifID$stop)
    #### output
    n_removed.dup.WT.MotifIDs = WT.MotifID_startstop %>% .[!. %in% pert.MotifID_startstop] %>% length
    n_kept.dup.WT.MotifIDs = WT.MotifID_startstop %>% .[. %in% pert.MotifID_startstop] %>% length
    n_introduced.MotifIDs = pert.MotifID_startstop %>% .[!. %in% WT.MotifID_startstop] %>% length
    
    # get WT Motifs in region 
    fimo_WT.regions = fimo_WT[fimo_WT$regions == regions,]
    fimo_WT.regions = fimo_WT.regions[-which(fimo_WT.regions$start == start.i & 
                                               fimo_WT.regions$stop == stop.i & 
                                               fimo_WT.regions$`# motif_id` == MotifID), ]
    ## find WT motifs that overlapped with the sequence
    nrow_fimo_WT.regions_overlapped_class1 = which(fimo_WT.regions$start >= start.i & fimo_WT.regions$start <= stop.i)
    nrow_fimo_WT.regions_overlapped_class2 = which(fimo_WT.regions$stop >= start.i & fimo_WT.regions$stop <= stop.i)
    nrow_fimo_WT.regions_overlapped_class3 = which(fimo_WT.regions$start <= start.i & fimo_WT.regions$stop >= stop.i)
    nrow_fimo_WT.regions_overlapped = c(nrow_fimo_WT.regions_overlapped_class1, 
                                        nrow_fimo_WT.regions_overlapped_class2, 
                                        nrow_fimo_WT.regions_overlapped_class3
    ) %>% unique
    ####test_nrow_fimo_WT.regions_notoverlapped = which(fimo_WT.regions$stop < start.i | fimo_WT.regions$start > stop.i)
    nrow_fimo_WT.regions_notoverlapped = setdiff(1:nrow(fimo_WT.regions), nrow_fimo_WT.regions_overlapped)
    fimo_WT.regions_overlapped = fimo_WT.regions[nrow_fimo_WT.regions_overlapped, ]
    fimo_WT.regions_notoverlapped = fimo_WT.regions[nrow_fimo_WT.regions_notoverlapped, ]
    
    # find pert motifs in the region
    fimo_pert.regions = fimo[fimo$sequence_name == sequence_name,]
    ## find pert motifs that overlapped with the sequence
    ### get nrows
    nrow_fimo_pert.regions_notoverlapped = which(fimo_pert.regions$stop < start.i | fimo_pert.regions$start > stop.i)
    nrow_fimo_pert.regions_overlapped = setdiff(1:nrow(fimo_pert.regions), nrow_fimo_pert.regions_notoverlapped)
    ### get data frames
    fimo_pert.regions_overlapped = fimo_pert.regions[nrow_fimo_pert.regions_overlapped, ]
    fimo_pert.regions_notoverlapped = fimo_pert.regions[nrow_fimo_pert.regions_notoverlapped, ]
    
    # comparing motifs
    ## wt motifs
    motif.WT = fimo_WT.regions$`# motif_id`
    motif.WT_start.end = 
      mapply(function(A, B, C) paste0(A, "_", B, ".", C), 
             fimo_WT.regions$`# motif_id`, fimo_WT.regions$start, fimo_WT.regions$stop)
    ### get motifIDs of the overlapped/not overlapped respectively
    motif.WT.overlapped_start.end = 
      mapply(function(A, B, C) paste0(A, "_", B, ".", C), 
             fimo_WT.regions_overlapped$`# motif_id`, 
             fimo_WT.regions_overlapped$start, 
             fimo_WT.regions_overlapped$stop)
    motif.WT.notoverlapped_start.end = 
      mapply(function(A, B, C) paste0(A, "_", B, ".", C), 
             fimo_WT.regions_notoverlapped$`# motif_id`, 
             fimo_WT.regions_notoverlapped$start, 
             fimo_WT.regions_notoverlapped$stop)
    ## pert motifs
    motif.pert = fimo_pert.regions$`# motif_id`
    motif.pert_start.end = 
      mapply(function(A, B, C) paste0(A, "_", B, ".", C), 
             fimo_pert.regions$`# motif_id`, fimo_pert.regions$start, fimo_pert.regions$stop)
    ### get motifIDs of the overlapped/not overlapped respectively
    motif.pert.overlapped_start.end = 
      mapply(function(A, B, C) paste0(A, "_", B, ".", C), 
             fimo_pert.regions_overlapped$`# motif_id`, 
             fimo_pert.regions_overlapped$start, 
             fimo_pert.regions_overlapped$stop)
    motif.pert.notoverlapped_start.end = 
      mapply(function(A, B, C) paste0(A, "_", B, ".", C), 
             fimo_pert.regions_notoverlapped$`# motif_id`, 
             fimo_pert.regions_notoverlapped$start, 
             fimo_pert.regions_notoverlapped$stop)
    ### output: unique_motif.pert = unique(motif.pert)
    n_WT.motifs = length(motif.WT); 
    n_WT.motifs.overlapped = nrow(fimo_WT.regions_overlapped)
    n_WT.motifs.overlapped_removed = length(motif.WT.overlapped_start.end[!motif.WT.overlapped_start.end %in% motif.pert_start.end])
    n_WT.motifs.overlapped_stillthere = length(motif.WT.overlapped_start.end[motif.WT.overlapped_start.end %in% motif.pert_start.end])
    n_WT.motifs.notoverlapped = nrow(fimo_WT.regions_notoverlapped)
    n_WT.motifs.notoverlapped_stillthere = length(motif.WT.notoverlapped_start.end[motif.WT.notoverlapped_start.end %in% motif.pert_start.end])
    #X$n_unique.WT.motifs[i] = length(unique_motif.WT)
    n_pert.motifs = length(motif.pert)
    #X$n_unique.pert.motifs[i] = length(unique_motif.pert)
    
    ## what are gained/lost
    motif.gained = motif.pert_start.end[!motif.pert_start.end %in% motif.WT_start.end]
    motif.lost = motif.WT_start.end[!motif.WT_start.end %in% motif.pert_start.end]
    ## gained/lost not MotifID
    motif.gained_no.MotifID = motif.gained[which(names(motif.gained) != MotifID)]
    motif.lost_no.MotifID = motif.lost[which(names(motif.lost) != MotifID)]
    ## get number of gained/lost
    n_gained = length(motif.gained)
    n_gained_no.MotifID = length(motif.gained_no.MotifID)
    n_gained_only.MotifID = n_gained - n_gained_no.MotifID
    
    n_lost = length(motif.lost)
    n_lost_no.MotifID = length(motif.lost_no.MotifID)
    n_lost_only.MotifID = n_lost - n_lost_no.MotifID
    
    ## gained overlapped/not
    motif.gained_overlapped = motif.pert.overlapped_start.end %>% .[! . %in% motif.WT_start.end]
    motif.gained_notoverlapped = motif.pert.notoverlapped_start.end %>% .[! . %in% motif.WT_start.end]
    n_gained.overlapped = length(motif.gained_overlapped)
    n_gained.notoverlapped = length(motif.gained_notoverlapped)
    # lost overlapped/not
    motif.lost_overlapped = motif.WT.overlapped_start.end %>% .[!. %in% motif.pert_start.end]
    motif.lost_overlapped_no.MotifID = motif.lost_overlapped[which(names(motif.lost_overlapped) != MotifID)]
    motif.lost_notoverlapped = motif.WT.notoverlapped_start.end %>% .[!. %in% motif.pert_start.end]
    n_lost.overlapped = length(motif.lost_overlapped)
    n_lost.overlapped_no.MotifID = length(motif.lost_overlapped_no.MotifID)
    n_lost.notoverlapped = length(motif.lost_notoverlapped)
    
    # calculate rate of wt removal
    ## wt overlapped
    removal_overlapped = length(motif.WT.overlapped_start.end[!motif.WT.overlapped_start.end %in% motif.pert_start.end])/length(motif.WT.overlapped_start.end)
    removal_notoverlapped = length(motif.WT.notoverlapped_start.end[!motif.WT.notoverlapped_start.end %in% motif.pert_start.end])/length(motif.WT.overlapped_start.end)
    ### output
    ratio_WT.motifs.overlapped_removal.ratio = removal_overlapped
    ratio_WT.motifs.notoverlapped_removal.ratio = removal_notoverlapped
    
    gc()
    
    return(
      c(
        n_Identical.WT.MotifIDs = n_Identical.WT.MotifIDs, 
        more_Identical.WT.MotifIDs = more_Identical.WT.MotifIDs, 
        n_Identical.duplicated.WT.MotifIDs = n_Identical.duplicated.WT.MotifIDs,
        IDs_Identical.WT.MotifIDs = IDs_Identical.WT.MotifIDs,
        n_Identical.PertX.MotifIDs.in.sequence = n_Identical.PertX.MotifIDs.in.sequence,
        contain_Identical.PertX.MotifIDs.in.sequence = contain_Identical.PertX.MotifIDs.in.sequence,
        n_removed.dup.WT.MotifIDs = n_removed.dup.WT.MotifIDs,
        n_kept.dup.WT.MotifIDs = n_kept.dup.WT.MotifIDs,
        n_introduced.MotifIDs = n_introduced.MotifIDs, 
        n_WT.motifs = n_WT.motifs, 
        n_WT.motifs.overlapped = n_WT.motifs.overlapped, 
        n_WT.motifs.overlapped_removed = n_WT.motifs.overlapped_removed, 
        n_WT.motifs.overlapped_stillthere = n_WT.motifs.overlapped_stillthere, 
        n_WT.motifs.notoverlapped = n_WT.motifs.notoverlapped, 
        n_pert.motifs = n_pert.motifs, 
        n_gained = n_gained, 
        n_gained_no.MotifID = n_gained_no.MotifID, 
        n_gained_only.MotifID = n_gained_only.MotifID,
        n_lost = n_lost, 
        n_lost_no.MotifID = n_lost_no.MotifID, 
        n_lost_only.MotifID = n_lost_only.MotifID, 
        ratio_WT.motifs.overlapped_removal.ratio = ratio_WT.motifs.overlapped_removal.ratio, 
        ratio_WT.motifs.notoverlapped_removal.ratio = ratio_WT.motifs.notoverlapped_removal.ratio, 
        n_gained.overlapped = n_gained.overlapped, 
        n_gained.notoverlapped = n_gained.notoverlapped, 
        n_lost.overlapped = n_lost.overlapped, 
        n_lost.overlapped_no.MotifID = n_lost.overlapped_no.MotifID, 
        n_lost.notoverlapped = n_lost.notoverlapped
      )
    )
    gc()
  }
v_checkWT <- Vectorize(checkWT)  

system.time(
  {X_checkWT <- X[, v_checkWT(
    MotifID = `# motif_id`, 
    regions = regions, sequence_name,  
    start = start, 
    stop = stop)]
  X_checkWT <-  t(X_checkWT)
  }
)
if (all.equal(nrow(X_checkWT), nrow(X))) {
  X = cbind(X, X_checkWT)
} else {
  print("rownames of X_checkWT and X do not match")
  saveRDS(X_checkWT, sprintf("%s-temp-X_checkWT.Rds",datetag))
  saveRDS(X, sprintf("%s-temp-X.Rds",datetag))
}

gc()

### Filter summary  -----
print("Filter summary")
X.Fs = select(X, paste0("F", 1:4))

n_F.equals.1 = rowSums(X.Fs)
n_F.equals.0 = 4 - n_F.equals.1

X$n_F.equals.1 = n_F.equals.1
X$n_F.equals.0 = n_F.equals.0

X$F_1s.summary = sprintf("Passed %.0f filters", X$n_F.equals.1)
X$F_0s.summary = sprintf("Didn't pass %.0f filters", X$n_F.equals.0)

# transform back
reviewer.files_df = as.data.frame(X)
# write out -----
saveRDS(reviewer.files_df, 
        sprintf("%s-reviewer.files_df.category.Rds",datetag))
