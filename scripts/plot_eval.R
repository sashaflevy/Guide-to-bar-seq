library(tidyverse)

datar <- read_tsv("eval_starcoded.txt",col_types="ccnn") 

pdatar <- full_join(
        datar %>% mutate(empirical_var_per_code=(input_counts-output_counts)/input_counts) %>%
            rename(input_counts_code=input_counts),
        datar %>% group_by(clone_id) %>%
            summarize(empirical_var_per_clone=(sum(input_counts)-sum(output_counts))/sum(input_counts),
                input_counts_clone=sum(input_counts)
            ),
        by="clone_id")

g <- pdatar %>% select(code,input_counts_code,empirical_var_per_code) %>%
    ggplot()+
    aes(x=input_counts_code,y=empirical_var_per_code)+
    geom_point()
g
ggsave("var_per_code.png",width=7,height=5)

g <- pdatar %>% select(clone_id,input_counts_clone,empirical_var_per_clone) %>%
    ggplot()+
    aes(x=input_counts_clone,y=empirical_var_per_clone)+
    geom_point()
g
ggsave("var_per_clone.png",width=7,height=5)

