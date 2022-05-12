##===============================================================================##
##                                0.Data pre                                     ## 
##===============================================================================##
library(readxl)
library(tidyverse)
control_index = dplyr::filter(sample_index,TRT == "a") %>% dplyr::select('Tube ID') %>% unlist(use.names = F)
treatment_index = dplyr::filter(sample_index,TRT == "b") %>%  dplyr::select('Tube ID') %>% unlist(use.names = F)

sample_index =read_excel("/Users/liulihe95/Desktop/Methionine/Samples_RNA-Seq.xlsx")
raw_data_index = list.files("/Users/liulihe95/Desktop/Methionine/counts/")[-c(1,10,22)]

# pre raw data - ungrouped
data_expr_all_raw = data.frame()
for (i in seq_along(raw_data_index)){
  print(i)
  tmp_name = substr(raw_data_index[i],1,5)
  #tmp_rownames = read.csv(raw_data_index[i],sep = "\t",header = F) %>% dplyr::slice(1:(n()-5)) %>%  dplyr::select(V1)
  tmp_read = read.csv(paste0("/Users/liulihe95/Desktop/Methionine/counts/",raw_data_index[i]),sep = "\t",header = F) %>% dplyr::slice(1:(n()-5)) %>% column_to_rownames('V1')
  if (i == 1){
    data_expr_all_raw = data.frame(tmp_name = tmp_read)
    colnames(data_expr_all_raw) = tmp_name
  } else {
    data_expr_all_raw = cbind(data_expr_all_raw,tmp_read)
    colnames(data_expr_all_raw)[i] = tmp_name
  }
}
data_expr_all_with0 = data_expr_all_raw[,c(which(substr(raw_data_index,2,5) %in% control_index),
                                           which(substr(raw_data_index,2,5) %in% treatment_index))]


##====================================================##
## 1. differential expression (Volcano plot)
## 2. gene set enrichement (ORA)
## 3. differential DNA methylation (manhattan plot)
## 4. Integration
##================================================##

library(edgeR)
## init
dgList = DGEList(counts = data_expr_all_with0,
                  genes = rownames(data_expr_all_with0))

## filter
countsPerMillion = cpm(dgList)
countCheck = countsPerMillion > 1 
head(countCheck)
keep = which(rowSums(countCheck) >= 2)
dgList = dgList[keep,]

## normalization
dgList = calcNormFactors(dgList, method="TMM")

##
sampleType = c(rep("Methionine",9),rep("Control",10))
designMat = model.matrix(~ sampleType)

dgList <- estimateGLMCommonDisp(dgList, design=designMat) 
dgList <- estimateGLMTrendedDisp(dgList, design=designMat) 
dgList <- estimateGLMTagwiseDisp(dgList, design=designMat)

fit <- glmFit(dgList, designMat) 
lrt <- glmLRT(fit)

##
library(tidyverse)
edgeR_result <- topTags(lrt,n=35000)$table

write.xlsx(edgeR_result, file = 'DEGs.xlsx')

##==========================================
edgeR_result_plot = edgeR_result_new %>%
  dplyr::select(
    genes,
    logFC,
    PValue
  ) %>% 
  mutate(gene_type = case_when(logFC >= log10(2) & PValue <= 0.05 ~ "up",
                               logFC <= log10(0.5) & PValue <= 0.05 ~ "down",
                               TRUE ~ "ns")) 

cols <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey") 
sizes <- c("up" = 2, "down" = 2, "ns" = 1) 
alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)

# plot
vol_plot <- edgeR_result_plot %>%
  ggplot(aes(x = logFC,
             y = -log10(PValue),
             fill = gene_type,    
             size = gene_type,
             alpha = gene_type)) + 
  #ggtitle("Gene expression changes") +
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters    
             colour = "black") +
  geom_hline(yintercept = -log10(0.05),linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)), linetype = "dashed")+
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-18, 6, 2)),limits = c(-8, 6)) +
  labs(title = "Gene expression changes",
       x = "log2(fold change)",
       y = "-log10(adjusted P-value)",
       colour = "Expression \nchange")+
  theme_bw() +  # Select theme with a white background  
  theme(plot.title = element_text(hjust = 0.5,face="bold",size = 24),
        axis.text = element_text(hjust = 0.5,face="bold",size = 16),
        axis.title.x = element_text(size = 18, face="bold", colour = "black"),
        axis.title.y = element_text(size = 18, face="bold", colour = "black"),
        legend.title = element_text(face="bold",size = 16),
        legend.position = 'right',
        legend.text = element_text(size = 15,face="bold"),
        legend.key.size = unit(1.5, "cm"),
        legend.key.width = unit(0.5,"cm"))


### PLOT1
tiff("/Users/liulihe95/Desktop/UW-Course/21-22-Spring/CS776_Adv_Bioinfo/pj/Final/Fig1-VocanoPlot.tiff", width = 12, height = 8, units = 'in', res = 500)
print(vol_plot)
dev.off()


#vol_plot

## Functional Enrichment
sig_g = edgeR_result %>% 
  dplyr::filter(PValue <= 0.05) %>% 
  pull(genes)
all_g = edgeR_result %>% 
  dplyr::filter(PValue <= 1) %>% 
  pull(genes)

length(sig_g)
length(all_g)

GeneInfo = convertNformatID(GeneSetNames=c("all"),
                            SigGene_list = list(sig_g),
                            TotalGene_list = list(all_g),
                            IDtype = "ens") # Need to choose from c('ens','entrez','symbol')

all_db = c('go','kegg','interpro','mesh','msig','reactome')

for (db in all_db){
  HyperGEnrich(GeneSet = GeneInfo,
               Database = db, #c("go","kegg","interpro","mesh","msig","reactome")
               minOverlap = 4, # minimum overlap of pathway genes and total genes
               pvalue_thres = 0.05, # pvalue of fisher's exact test
               adj_pvalue_thres = 1, # adjusted pvalues based on multiple testing correction
               padj_method = "fdr", # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
               NewDB = F)
}

#
output_list = c(
  "go-Enrichment-4-0.05-1.rda",
  "kegg-Enrichment-4-0.05-1.rda",
  "interpro-Enrichment-4-0.05-1.rda",
  "mesh-Enrichment-4-0.05-1.rda",
  "msig-Enrichment-4-0.05-1.rda",
  "reactome-Enrichment-4-0.05-1.rda")

output_container_all = list()

for (i in seq_along(output_list)){
  #
  tmp = data.frame()
  load(output_list[i])
  
  output_container_all[[i]] = data.frame(results[[1]]) %>% 
    mutate(hitsPerc = -hitsPerc)
}

library(openxlsx)
write.xlsx(output_container_all, file = 'enrichment.xlsx')


### PLOT2
library(openxlsx)
library(tidyverse)
enrichment_plot_df = read.xlsx("/Users/liulihe95/Desktop/UW-Course/21-22-Spring/CS776_Adv_Bioinfo/pj/Final/Enrichment_Plot.xlsx") %>% 
  mutate(Tag = paste0(ID,'[',totalG,']')) %>% 
  relocate(Tag, .after = ID)
##
plot_all = enrichment_plot_df %>% 
  mutate(Category = as.factor(Category)) %>% 
  mutate(log10pvalue = -log10(pvalue)) %>% 
  group_by(Category)



plot_all$Tag =
  factor(plot_all$Tag,
         levels = rev(c(
           "GO:0007519[35]","GO:0014894[4]","GO:0035914[38]","IPR002546[4]","IPR039704[4]",# Muscle Devevelopment
           "GO:0006955[106]","GO:0002376[124]","GO:0002250[49]","GO:0045087[165]",# Immune system
           "bta00260[29]","bta00910[7]", # Metabolism
           "GO:0006813[51]","GO:0071805[50]","GO:0005324[10]","GO:0015909[9]" # Fatty acid and Ion transport
         )))


plot_all$Category =
  factor(plot_all$Category,
         levels = (c(
           "Muscle Development",
           "Immune system",
           "Metabolism",
           "Fatty acid and Ion transport"
         )))


ggEnrich = ggplot(plot_all) + coord_flip()+
  geom_bar(mapping = aes(x = Tag, y = hitsPerc,fill = Category),#
           stat = "identity", width = 0.2) +
  geom_point(mapping = aes(x = Tag, y = 10 * log10pvalue),size = 1) +
  #scale_y_continuous() 
  scale_y_continuous(name = expression(bold("Percentage of Significant Genes")),
                     sec.axis = sec_axis(~. * 0.1, breaks = c(0,5,10), 
                                         name = expression(bold("-log10(P-value)")))
                     ,limits = c(0,70), breaks = seq(0,70,by=5)
                     ,expand = c(0,1)) +
  guides(fill=guide_legend(title="",
                           keywidth = 1, nrow = 1, 
                           keyheight = 1))+
  theme(legend.position = "bottom",
        legend.text=element_text(size=14,face="bold",color = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=10,face="bold",color = "black"),
        axis.text.x = element_text(size=10,face="bold",color = "black"))
ggEnrich

tiff("/Users/liulihe95/Desktop/UW-Course/21-22-Spring/CS776_Adv_Bioinfo/pj/Final/Fig2-EnrichBar.tiff",
     width = 12, height = 8, units = 'in', res = 500)
print(ggEnrich)
dev.off()

########################
## download DNA methy
########################
load("/Users/liulihe95/Desktop/Methionine/Network_No_Crt/Net/0_776pj/DNAmeth.rda")


mhPlot_df = list()
DNAmeth$chr = gsub("chr","",DNAmeth$chr)


for (i in c(seq(1,29))){
  print(i)
  tmpdf = DNAmeth %>% 
    dplyr::filter(chr == i) %>% 
    rename(P = pvalue,CHR = chr) %>% 
    mutate(BP = seq.int(nrow(.))) %>% 
    dplyr::select(CHR,BP,P)
  mhPlot_df[[i]] = tmpdf
}

DNAmeth = as.data.frame(do.call(rbind,mhPlot_df))
DNAmeth$CHR = as.integer(DNAmeth$CHR)

DNAmeth = DNAmeth %>% 
  mutate(SNP = as.character(CHR))


tiff("Fig3-EnrichBar.tiff",
     width = 12, height = 8, units = 'in', res = 300)
manhattan(DNAmeth)
dev.off()

tiff("Fig3.1-QQ.tiff",
     width = 12, height = 8, units = 'in', res = 300)
qq(DNAmeth$P, main = "Q-Q plot of GWAS p-values", 
   xlim = c(0, 50), ylim = c(0, 50), 
   pch = 18, col = "blue4", cex = 1.5, las = 1)
dev.off()


#### Permutation
DEGs = read_xlsx('DEGs.xlsx') %>% 
  dplyr::filter(PValue <= 0.05) %>% 
  pull(genes)

gene_down = edgeR_result_plot %>% 
  dplyr::filter(gene_type == "down") %>% 
  pull(genes)
gene_up = edgeR_result_plot %>% 
  dplyr::filter(gene_type == "up") %>% 
  pull(genes)


geneMethLevel_down = read_xlsx('./MethLevel.xlsx') %>% 
  dplyr::filter((Gene %in% gene_down)) %>% 
  dplyr::select(Gene,Prop_Prpt, Prop_Body, Prop_All)

geneMethLevel_up = read_xlsx('./MethLevel.xlsx') %>% 
  dplyr::filter((Gene %in% gene_up)) %>% 
  dplyr::select(Gene,Prop_Prpt, Prop_Body, Prop_All)

geneMethLevel_other = read_xlsx('./MethLevel.xlsx') %>% 
  dplyr::filter(!(Gene %in% DEGs)) %>% 
  dplyr::select(Gene,Prop_Prpt, Prop_Body, Prop_All)

# 324
p_allregion = c()
p_promoter = c()
p_genebody = c()
type = "less"
set.seed(123)
for (i in c(1:2000)){
  print(i)
  sp = sample(c(1:nrow(geneMethLevel_other)),nrow(geneMethLevel_up))
  #
  p_allregion[i] = t.test(
    geneMethLevel_up$Prop_All
    ,
    geneMethLevel_other$Prop_All[sp],
    type,
    )$p.value
  #
  p_promoter[i] = t.test(
    geneMethLevel_up$Prop_Prpt,
    geneMethLevel_other$Prop_Prpt[sp],
    type,
  )$p.value
  #
  p_genebody[i] = t.test(
    geneMethLevel_up$Prop_Body,
    geneMethLevel_other$Prop_Body[sp],
    type,
  )$p.value
}

down_pcom_df_gt = data.frame(PValue = c(p_allregion,p_promoter,p_genebody),
                             Region = c(
                               rep("Overall",length(p_allregion)),
                               rep("Promoter",length(p_promoter)),
                               rep("Gene body",length(p_genebody)))
                             )

pcom_df_lt = data.frame(PValue = c(p_allregion,p_promoter,p_genebody),
                        Region = c(
                          rep("Overall",length(p_allregion)),
                          rep("Promoter",length(p_promoter)),
                          rep("Gene body",length(p_genebody)))
)

#
downreg_gt_compareP = ggplot(down_pcom_df_gt, aes(x=PValue)) + 
  geom_histogram(aes(y=..density..),
                 colour="black", fill="lightblue",
                 bins=150)+
  geom_density(alpha=.2, fill="lightblue") +
  facet_grid(Region ~ .) +
  labs(title = "Empirical p value distribution")+
  theme(plot.title = element_text(hjust = 0.5,face="bold",size = 20),
        axis.title.x = element_text(size=12,face="bold",color = "black"),
        axis.title.y= element_text(size=12,face="bold",color = "black"),
        axis.text.y = element_text(size=10,face="bold",color = "black"),
        axis.text.x = element_text(size=10,face="bold",color = "black"),
        strip.text.y = element_text(size = 14,face="bold",color = "black"))

tiff("/Users/liulihe95/Desktop/UW-Course/21-22-Spring/CS776_Adv_Bioinfo/pj/Final/Fig4-MethDownRegG.tiff",
     width = 12, height = 8, units = 'in', res = 500)
print(downreg_gt_compareP)
dev.off()


#
upreg_lt_compareP = ggplot(pcom_df_lt, aes(x=PValue)) + 
  geom_histogram(aes(y=..density..),
                 colour="black", fill="lightblue",
                 bins=150)+
  geom_density(alpha=.2, fill="lightblue") +
  facet_grid(Region ~ .) +
  theme(axis.title.x = element_text(size=12,face="bold",color = "black"),
        axis.title.y= element_text(size=12,face="bold",color = "black"),
        axis.text.y = element_text(size=10,face="bold",color = "black"),
        axis.text.x = element_text(size=10,face="bold",color = "black"),
        strip.text.y = element_text(size = 14,face="bold",color = "black"))

tiff("/Users/liulihe95/Desktop/UW-Course/21-22-Spring/CS776_Adv_Bioinfo/pj/Final/Fig5-MethUpRegG.tiff",
     width = 12, height = 8, units = 'in', res = 500)
print(upreg_lt_compareP)
dev.off()


## RcisTarget (TF Enrichment)
# library(EnrichKit)
# G_Info = convertNformatID(GeneSetNames=c("t"),
#                           SigGene_list = list(sig_g),
#                           TotalGene_list = list(all_g),
#                           IDtype = "ens") # Need to choose from c('ens','entrez','symbol')
# 
# all_g_symbol = unique(G_Info$t$SYMBOL_Suggested)
# all_g_symbol = all_g_symbol[!is.na(all_g_symbol)]
# length(all_g_symbol)
# write.table(all_g_symbol, file = "all_symbol.txt",
#             append = FALSE, quote = F,
#             row.names = F,
#             col.names = F)
# 
# 
# sig_g_symbol = G_Info$t %>% 
#   dplyr::filter(Sig == 1) %>% 
#   dplyr::select(SYMBOL_Suggested) %>% 
#   drop_na() %>% 
#   pull(SYMBOL_Suggested)
# length(sig_g_symbol)
# write.table(sig_g_symbol, 
#             file = "sig_symbol.txt",
#             append = FALSE, quote = F,
#             row.names = F,
#             col.names = F)
# 
# BiocManager::install("RcisTarget")
# browseVignettes("RcisTarget")



### DNA methylation

# # dataset pre
# options(stringsAsFactors = FALSE)
# data_loci = '/blue/mateescu/lihe.liu/Methylation_WGCNA'
# setwd(data_loci)
# load("methCov08Stat.rda")
# library(tidyverse)
# library(methylKit)
# 
# # select Significant Diff Cs
# chr_index = paste(rep('chr',30),c(seq(1,29),'X'),sep = "")
# Diff_C_all = getData(methCov08Stat) %>% 
#   mutate_at(vars(chr),as.character) %>% 
#   dplyr::filter(chr %in% chr_index) # 5136556
# 
# library(openxlsx)
# write.xlsx(Diff_C_all, file = 'DNAmeth.xlsx')




