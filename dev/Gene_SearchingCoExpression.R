#### Esto responde la pregunta de que genes se co-detectan mas con mi genes de interes ----
## Este vector podria ser de mas de un gen.
gene <- c('Sox11',"Kcnq1ot1","Gria3","Tenm1",'Nell1')

vtor <- (assay(sce,'logcounts')[gene,] > 0) %>% apply(2,all)
mtx <- (assay(sce,'logcounts') > 0) %>% apply(1,function(x){
      tapply(X = x,sce$ds1clusterP,function(x){
        y <- sum(x & vtor[names(x)])
        w <- sum(vtor[names(x)])
        y/w
      })
})
mtx[is.na(mtx)] <- 0
n <- 'intesect'
clusters <- c('GCmat2')
thr <- 1
ldf <- (mtx >= thr) %>% apply(1,which)

if(n == 'intesect'){fct <-  intersect} else{fct <- union}
vtpr <- Reduce(fct,ldf[clusters])
if(!is.null(vtpr)){
  mtx[,vtpr] %>% t %>% Heatmap()
  'Sox11' %in% colnames(mtx[,vtpr])
}
mtx[,'Sox11']


#### Esto responde para un dado set de genes cual es el cluster donde mas se codetectan

## Este vector podria ser de mas de un gen.
gene <- c('Sox11',"Kcnq1ot1","Gria3","Tenm1",'Nell1')

vtor <- (assay(sce,'logcounts')[gene,] > 0) %>% apply(2,all)
vtor  %>% tapply(.,sce$ds1clusterP,function(x){sum(x)/length(x)})
