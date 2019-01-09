library(tidyverse)

datar <- list.files(pattern="*.tsv") %>% #,full.names=T,path="../work/11/cffb39227b50bc388d1d9247258584/") %>% 
    tibble(filenames=.) %>%
    mutate(RawFiles=map(filenames,read_tsv,col_types="ccnn")) %>%
    unnest() %>%
    group_by(filenames) %>%
    mutate(input_counts=input_counts/sum(input_counts),
        output_counts=output_counts/sum(output_counts))

pdatar <- full_join(
        datar %>% mutate(empirical_var_per_code=(input_counts-output_counts)/input_counts) %>%
            rename(input_prop_per_code=input_counts),
        datar %>% group_by(filenames,clone_id) %>%
            summarize(empirical_var_per_clone=(sum(input_counts)-sum(output_counts))/sum(input_counts),
                input_prop_per_clone=sum(input_counts)
            ),
        by=c("filenames","clone_id")
        ) %>% 
    ungroup() %>%
    mutate(filenames=unlist(map(filenames,str_remove,".*\\/\\/"))) %>%
    mutate(barcode=unlist(map(filenames,str_replace,pattern="^([ATCGN]+)_.*",replacement="\\1"))) %>%
    mutate(lineages=unlist(map(filenames,str_replace,pattern=".*_([0-9]+)lineages_.*",replacement="\\1"))) %>%
    mutate(barcodes_per=unlist(map(filenames,str_replace,pattern=".*_([0-9\\.]+)barcodesper_.*",replacement="\\1"))) %>%
    mutate(fixed_barcodes=unlist(map(filenames,str_replace,pattern=".*_fixedBarcodesPer((True)|(False))_.*",replacement="\\1"))) %>%
    mutate(cell_sample=unlist(map(filenames,str_replace,pattern=".*_([0-9\\.]+)_([0-9\\.]+)reads_.*",replacement="\\1"))) %>%
    mutate(reads=unlist(map(filenames,str_replace,pattern=".*_([0-9\\.]+)_([0-9\\.]+)reads_.*",replacement="\\2"))) %>%
    mutate(error=unlist(map(filenames,str_replace,pattern=".*_(.+)error_.*",replacement="\\1"))) 

g <- pdatar %>% 
    ggplot()+
    facet_grid(error~barcode+barcodes_per)+
    aes(x=lineages,y=empirical_var_per_code)+
    geom_boxplot()+
#    geom_dotplot(binaxis="y",stackdir="center",binwidth=0.0001)+
    theme(axis.text.x=element_text(angle=90))
g

ggsave(str_c("var_per_code.png"),width=7,height=5)

g <- pdatar %>% 
    ggplot()+
    facet_grid(error~barcode+barcodes_per)+
    aes(x=lineages,y=empirical_var_per_clone)+
    geom_boxplot()+
#    geom_dotplot(binaxis="y",stackdir="center",binwidth=0.0001)+
    theme(axis.text.x=element_text(angle=90))
g

ggsave(str_c("var_per_clone.png"),width=7,height=5)

pdatar %>% 
    { str_c("> ",.$clone_id,"_",.$code,"_",.$empirical_var_per_code,"\n",.$code) } %>% 
    write_lines(str_c("var_per_code.fasta")) 

