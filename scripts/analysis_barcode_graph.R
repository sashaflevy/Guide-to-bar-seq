library(tidyverse)

slurp_distance_tsv <- function(x) {
    big_table <- c("0"=0)
    lapply(
        lapply(
            str_split(readLines(x),"\\s+"),table),
            function(x){ 
                new_names <- setdiff(names(x),intersect(names(big_table),names(x))); 
                for (i in new_names) { big_table[i] <<- 0; }; 
                for (i in names(x)) { big_table[i] <<- big_table[i] + x[i] }; 
            }
        ); 
    return(
        data.frame(Distance=as.numeric(names(big_table)),
            Count=big_table,
            stringsAsFactors=F
            )
        )
}

can_coerce_to_numeric <- function (x) {
    return(suppressWarnings(x[!is.na(as.numeric(x))]))
}

datar <- list.files(pattern="*.tsv") %>% 
    tibble(filenames=.) %>%
    mutate(TabulatedDistances=map(filenames,slurp_distance_tsv)) %>%
    unnest() %>% 
    ungroup() %>%
    mutate(filenames=unlist(map(filenames,str_remove,".*\\/\\/"))) %>%
    mutate(barcode=unlist(map(filenames,str_replace,pattern="^barcode_distances_([ATCGN]+)_.*",replacement="\\1"))) %>%
    mutate(lineages=unlist(map(filenames,str_replace,pattern=".*_([0-9]+)lineages_.*",replacement="\\1"))) %>%
    mutate(barcodes_per=unlist(map(filenames,str_replace,pattern=".*_([0-9\\.]+)barcodesper_.*",replacement="\\1"))) %>%
    mutate(fixed_barcodes=unlist(map(filenames,str_replace,pattern=".*_fixedBarcodesPer((True)|(False))_.*",replacement="\\1"))) %>%
    mutate(replicate=unlist(map(filenames,str_replace,pattern=".*replicate(.+).tsv",replacement="\\1"))) %>%
    spread(Distance,Count,sep="_",fill=0) %>%
    gather(key=Distance,value=Count,starts_with("Distance_")) %>%
    mutate(Distance=as.numeric(sub("Distance_","",Distance)))

g <- datar %>% 
    group_by(filenames) %>%
    mutate(PropCount=Count/sum(Count)) %>%
    filter(PropCount>0) %>%
    ggplot()+theme_bw()+
    facet_grid(barcode~lineages+barcodes_per,scale="free_y")+
    aes(x=Distance,y=PropCount,fill=replicate)+
    geom_bar(stat="identity",position="dodge")+
    ylab("Proportion of barcodes")+
    xlab("Lev Distance")+
    theme(axis.text.x=element_text(angle=90))
g

ggsave(str_c("distance_distribution_prop_all.png"),g,width=7,height=7)


ggsave(str_c("distance_distribution_prop_lowend.png"),
    g+coord_cartesian(xlim=c(0,3),ylim=c(0,0.01)),
    width=7,height=7)


g <- datar %>% 
    group_by(filenames) %>%
    ggplot()+theme_bw()+
    facet_grid(barcode~lineages+barcodes_per,scale="free_y")+
    aes(x=Distance,y=Count,fill=replicate)+
    geom_bar(stat="identity",position="dodge")+
    ylab("Count of barcodes")+
    xlab("Lev Distance")+
    theme(axis.text.x=element_text(angle=90))
g

ggsave(str_c("distance_distribution_count_all.png"),g,width=7,height=7)


ggsave(str_c("distance_distribution_count_lowend.png"),
    g+coord_cartesian(xlim=c(0,3),ylim=c(0,10))+
        scale_y_continuous(breaks=c(seq(0,10),seq(20,100,10))),
    width=7,height=7)

