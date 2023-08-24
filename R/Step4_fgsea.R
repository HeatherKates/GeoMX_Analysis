res2 <-  %>% 
  dplyr::select(SYMBOL, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(SYMBOL) %>% 
  summarize(stat=mean(stat))
res2

pathways.hallmark <- gmtPathways("h.all.v6.2.symbols.gmt")
for (i in 1:Num_Contrasts){
  fgsea.df <- POS.Result[[i]] %>% select(Gene,pvalue)
  ranks <- deframe(fgsea.df)
  fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
  if(nrow(fgseaResTidy%>%filter(padj<0.05))>0){
    ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
      geom_col(aes(fill=padj<0.05)) +
      coord_flip() +
      labs(x="Pathway", y="Normalized Enrichment Score",
           title="Hallmark pathways NES from GSEA") + 
      theme_minimal()
  }
}

