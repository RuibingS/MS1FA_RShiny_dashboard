# test cor




FT_test <- featureTable.import.fun(here::here("Data", "feature_table","PA14_featureTable_XCMS_CAMERA.csv")) %>%
  dplyr::filter(rt >= 60 &  rt <= 1200)

FT_test[which(FT_test$feature_name=="FT0203"),]

colnames(FT_test)

mat <- FT_test[,9:32] %>%   #FT_test[,6:10] %>% 
  #dplyr::select(.,input$picker) %>%   # picked columns
  dplyr::mutate(across(where(is.numeric), ~na_if(., 0))) %>%  # replace 0 with NA in numeric columns
  dplyr::mutate(across(where(is.numeric), ~ifelse(is.na(.), NA_real_, log10(.)))) 


t_mat<-t(mat)

colnames(t_mat) <- FT_test$feature_name


#t_mat[,1:10]



# point-to-point correlations, skipping NAs
res.cor <- pairwiseCor(x = t_mat,method = "pearson") 


colnames(res.cor) <- FT_test$feature_name
rownames(res.cor) <- FT_test$feature_name


# res.cor[which(rownames(res.cor)%in%"FT0243"),which(colnames(res.cor)%in%"FT0245")]

# calculateCorrelation(col1 = t_mat[,which(colnames(t_mat) %in% "FT0243")], col2 = t_mat[,which(colnames(t_mat) %in% "FT0245")] , method = "pearson")


# cor.test(t_mat[,which(colnames(t_mat) %in% "FT0243")], t_mat[,which(colnames(t_mat) %in% "FT0245")])


# res.cor[which(rownames(res.cor)%in%"FT0203"),which(colnames(res.cor)%in%"FT0209")]
# 
# calculateCorrelation(col1 = t_mat[,which(colnames(t_mat) %in% "FT0203")], col2 = t_mat[,which(colnames(t_mat) %in% "FT0209")] , method = "pearson")
# 
# 
# cor.test(t_mat[,which(colnames(t_mat) %in% "FT0203")], t_mat[,which(colnames(t_mat) %in% "FT0209")])


res.cor[res.cor < 0.8 ] <- NA

res.cor[lower.tri(res.cor, diag = FALSE)] <- NA


res.cor[1:50,1:50]


adj_mat.long<-na.omit(data.frame(as.table(res.cor)))

# left join to add Var1, Var2,Freq,mz_x,mz_y,rt_x,rt_y,rt_diff,mz_diff
adj.full<-left_join_and_mutate_fun(FT_df=FT_test,long_df = adj_mat.long,rt_thr=5)


#head(adj.full)



# order by Freq
get.adj_order<-adj.full[order(adj.full$Freq,decreasing = T),]

# keep the unique

# get.adj_order_unique<-get.adj_order[which(get.adj_order$Var1!=get.adj_order$Var2),]
# # assign group index
# groupCorFeatures_res<-groupCorFeatures(adj_long =get.adj_order_unique,threshold=0.8)


# str(groupCorFeatures_res)

# groupCorFeatures_res_df <- groupCorFeatures_res[[1]]
# 
# groupCorFeatures_res_df[1:5,]
# 
# 
# groupCorFeatures_res_df[which(groupCorFeatures_res_df$feature_name%in% c("FT0243","FT0245")), ]
# 
# groupCorFeatures_res_df[which(groupCorFeatures_res_df$cor_group %in% c("corgroup25","corgroup29")),]
# 
# res.cor[which(rownames(res.cor) %in% "FT0243"), which(colnames(res.cor) %in% c("FT0245","FT0246" , "FT0244"))]


# res.cor[which(rownames(res.cor) %in% "FT0203"), which(colnames(res.cor) %in% c("FT0209","FT0202" , "FT0203", "FT0213" , "FT0212", "FT0216"  ))]
# 
# res.cor[which(rownames(res.cor) %in% "FT0202"), which(colnames(res.cor) %in% c("FT0209"))]
# 
# res.cor[which(rownames(res.cor) %in% "FT0202"), which(colnames(res.cor) %in% c("FT0203"))]













