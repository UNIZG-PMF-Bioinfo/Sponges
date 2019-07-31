library(ggplot2)
library(data.table)

###########################################################################################################################################################
### VISUALIZE InterProScan tsv formated RESULTS (Pfam + Phobius + TBD(...?)) passed as a Data.Table structure
###########################################################################################################################################################

###########################################################################################################################################################
### LAYER GENERATION FUNCTIONS

### getPhobius
getPhobius <- function(x, idx = 0){
  trm <- x[V5 == "TRANSMEMBRANE", .(paste0("+ geom_rect(aes(xmin= ", V7, ", xmax=", V8,", ymin=",idx,"+0.5, ymax=",idx,"+1.5), alpha= 0.7, fill = '#9BD9E8')"))]
  exc <- x[V5 == "NON_CYTOPLASMIC_DOMAIN", .(paste0("+ geom_rect(aes(xmin= ", V7, ", xmax=", V8,", ymin=",idx,"+1.4, ymax=",idx,"+1.7), alpha= 0.7, fill = '#9BD9E8')"))]
  cyt <- x[V5 == "CYTOPLASMIC_DOMAIN", .(paste0("+ geom_rect(aes(xmin= ", V7, ", xmax=", V8,", ymin=",idx,"+0.6, ymax=",idx,"+0.3), alpha= 0.7, fill = '#9BD9E8')"))]
  return(gsub(", \\+", " +",toString(rbind(trm,exc,cyt)$V1)))
}

### getPfam
getPfam <- function(x, idx = 0){
  pfm <- x[V4 == "Pfam", .(paste0("+ geom_rect(aes(xmin= ", V7, ", xmax=", V8,", ymin=",idx,"+0.85, ymax=",idx,"+1.15), fill = '", getHexColor(V5),"', color = 'black', size = 1.5)"))]
  lbl <- x[V4 == "Pfam", .(paste0("+ geom_text(aes(", (V7+V8)/2, ",",idx,"+1), label = '", V5,"', size = 4, color = 'white', fontface = 'bold')"))]
  return(gsub(", \\+", " +",toString(rbind(pfm,lbl)$V1)))
}

### getName
getName <- function(x, idx = 0, pwidth){
  nmm <- paste0("+ geom_text(aes(", pwidth/2, ",",idx,"+1.9), label = '", x[1,1],"', size = 4, color = 'black', fontface = 'italic')")
  return(nmm)
}

### getLength 
getLength <- function(x, idx = 0){
  plen<- paste0("+ geom_segment(aes(x=0, xend= ", x[1,3], ", y=",idx,"+1, yend=",idx,"+1, size = 2), show.legend = FALSE, lineend = 'round')")
  return(plen)
}

### getHexColor
getHexColor <- function(x){
  hxcol <- paste0("#" ,as.hexmode((as.numeric(substring(x,5,5))+1)*16), as.hexmode((as.numeric(substring(x,6,6))+1)*16), as.hexmode((as.numeric(substring(x,7,7))+1)*16))
  return(hxcol)
}

### getAnno
getAnno <- function(x, pheigth){
  ann <- unique(x[V4 == "Pfam",.(V5,V6,V13)])
  ann <- ann[,.(paste0("+ geom_text(aes(", 0, ",", pheigth - match(V5,ann$V5) * 0.2,"), label = '", V5, " - ", V6, " - ", V13, "\n', size = 4, color = 'black', fontface = 'bold', hjust = 0)"))]
  return(gsub(", \\+", " +",toString(ann$V1)))
}

###########################################################################################################################################################
### DATA FILTERING FUNCTIONS

# FILTER BY KEYWORDS IN FEATURE NAME OR DESCRIPTION
# subset all data for those genes where the pattern (any ammount of keywords) [WORDS] mattched in its name
# or in its description [if the flag SEARCH.DESC is turned on]

getGenes <- function(d, words, search.desc = FALSE){
  
  tbl <- d[0]
  
  for(word in words){
    
    dname <- grep(word, d$V6, ignore.case = TRUE)
    if(search.desc == TRUE){
      ddesc <- grep(word, d$V13, ignore.case = TRUE)
      dname <- unique(c(dname, ddesc))
    }
    
    tmp <- d[dname]
    tmp <- tmp[V4 != "Phobius"]
    tbl <- rbind(tbl, d[V1 %in% tmp$V1])
  }
  
  tbl <- unique(tbl)
  
  if(dim(tbl)[1] > 0){
    return(tbl)
  }
  
  else{
    return(print("No results were found"))
  }
  
}

# FILTER BY A LIST OF PFAMs OF INTEREST
# subset all data for those genes associated to a PFAM ID

filterPfam <- function(d, pfam){
  
  tbl <- d[0]
  
  for(word in pfam){

    tmp <- d[V5 == word]
    tbl <- rbind(tbl, d[V1 %in% tmp$V1])
  }
  
  tbl <- unique(tbl)
  
  if(dim(tbl)[1] > 0){
    return(tbl)
  }
  
  else{
    return(print("No results were found"))
  }
  
}

###########################################################################################################################################################
### PLOTTING FUNCTION

drawPfab <- function(x, words = NULL, search.desc = FALSE, pfam = NULL){
  
  # Filters by a list of keywords or PFAM IDs
  if(length(words) > 0){
    x <- getGenes(x, words = words, search.desc = search.desc)
  }

  if(length(pfam) > 0){
    x <- filterPfam(x, pfam = pfam)
  }
  
  # Replaces single quotes from names and descriptions since they break the getAnno function
  x$V13 <-  gsub("'","\"", x$V13)
  x$V6 <-  gsub("'","\"", x$V6)
  
  pwidth <- max(x$V3)
  pheigth <- length(unique(x$V1))*1.7 + 1
  layout <- "ggplot()+ geom_blank()+ xlim(0,pwidth)+ ylim(0,pheigth) + theme_void()"
  
  body <- ""
  idx <- 0
  for(i in unique(x$V1)){
    y <- x[V1 == i]
    protlen <- getLength(y, idx)
    pfam <- getPfam(y, idx)
    phobius <- getPhobius(y, idx)
    name <- getName(y, idx, pwidth)
    body <- paste0(body, phobius, protlen, pfam, name)
    idx <- idx + 1.7
  }
  
  annotation <- getAnno(x, pheigth)
  toplot <- paste0(layout, body, annotation)
  eval(parse(text=toplot))
}

###########################################################################################################################################################
###########################################################################################################################################################
### DATA IMPORT

setwd("D:/Downloads/")
d <- fread("./out_some_test.out", sep = "\t", fill= T)

d <- fread("./Porifera_SLC12_IPscn.tsv.tsv", sep = "\t", fill= T)
d <- d[V3 >= 1000]

drawPfab(d, words = "domain", search.desc = TRUE)
