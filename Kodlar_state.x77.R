
if(!require(cluster)) install.packages("cluster")
if(!require(factoextra)) install.packages("factoextra")
if(!require(pheatmap)) install.packages("pheatmap")

library(cluster)
library(factoextra)
library(pheatmap)


data(state)
df <- as.data.frame(state.x77)
colnames(df) <- c("Population", "Income", "Illiteracy", "Life_Exp", 
                  "Murder", "HS_Grad", "Frost", "Area")


df_scaled <- scale(df)


dosya_adi0 <- "PAM_Elbow_Method.png"
png(dosya_adi0, width = 800, height = 600)
print(fviz_nbclust(df_scaled, pam, method = "wss") +
        geom_vline(xintercept = k_sayisi, linetype = 2) + # Se??ti??in k'y?? ??izgiyle g??ster
        labs(title = "Elbow Metodu (Optimal K??me Say??s?? ????in)",
             subtitle = "Dirsek (K??r??lma) Noktas??n?? Ara"))
dev.off()
print(paste(dosya_adi0, "kaydedildi. (Analize ba??lamadan buna bak!)"))


w_16_9 <- 1600
h_16_9 <- 900
cozunurluk <- 150


dosya_adi_manhattan <- "PAM_Elbow_Method_Manhattan_16x9.png"

png(dosya_adi_manhattan, width = w_16_9, height = h_16_9, res = cozunurluk)


print(fviz_nbclust(df_scaled, pam, method = "wss", metric = "manhattan") +
        geom_vline(xintercept = 4, linetype = 2) + # Tahmini k=4 ??izgisini yine koyal??m
        labs(title = "Elbow Metodu (Manhattan Uzakl??????)", 
             subtitle = "Optimal K??me Say??s?? (L1 Normuna G??re)"))

dev.off()

print("Manhattan uzakl??????na g??re Elbow grafi??i kaydedildi.")



k_sayisi <- 5
metrik <- "manhattan"


set.seed(123)
pam_fit <- pam(df_scaled, k = k_sayisi, metric = metrik)


dosya_adi1 <- paste0("PAM_Cluster_Plot_k", k_sayisi , metrik, ".png")
png(dosya_adi1, width = 1600, height = 900)
print(fviz_cluster(pam_fit, 
                   data = df_scaled,
                   palette = "jco", 
                   ellipse.type = "convex",
                   repel = TRUE,
                   ggtheme = theme_minimal(),
                   main = paste("PAM Kume Grafigi (k =", k_sayisi, ")")))
dev.off()
print(paste(dosya_adi1, "kaydedildi."))


dosya_adi3 <- paste0("PAM_Silhouette_Plot_k", k_sayisi, metrik, ".png")
png(dosya_adi3, width = 1600, height = 900)
print(fviz_silhouette(pam_fit, 
                      palette = "jco", 
                      ggtheme = theme_minimal(),
                      main = paste("Siluet Analizi (k =", k_sayisi, ")")))
dev.off()
print(paste(dosya_adi3, "kaydedildi."))


dosya_adi4 <- paste0("PAM_Heatmap_k", k_sayisi, metrik , ".png")


df_ordered <- df_scaled[order(pam_fit$clustering), ]
annotation_row <- data.frame(Cluster = factor(sort(pam_fit$clustering)))
rownames(annotation_row) <- rownames(df_ordered)

png(dosya_adi4, width = 1600, height = 900, res = 150) # ????z??n??rl?????? biraz art??rd??m
pheatmap(df_ordered, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         annotation_row = annotation_row,
         main = paste("Eyalet Profilleri Is?? Haritasi (k =", k_sayisi, ")"),
         display_numbers = FALSE)
dev.off()
print(paste(dosya_adi4, "kaydedildi."))

print("--- T??m g??rseller ba??ar??yla olu??turuldu! ---")




if(!require(cluster)) install.packages("cluster")
if(!require(factoextra)) install.packages("factoextra")
if(!require(pheatmap)) install.packages("pheatmap")

library(cluster)
library(factoextra)
library(pheatmap)


data(state)
df <- as.data.frame(state.x77)
colnames(df) <- c("Population", "Income", "Illiteracy", "Life_Exp", 
                  "Murder", "HS_Grad", "Frost", "Area")
df_scaled <- scale(df)


pca_res <- prcomp(df_scaled, scale = TRUE)


df_pca <- pca_res$x[, 1:4] 


k_sayisi <- 5
# =======================================================


w_16_9 <- 1600
h_16_9 <- 900
cozunurluk <- 150

dosya_adi0 <- "PCA_PAM_Elbow_Method_16x9.png"
png(dosya_adi0, width = w_16_9, height = h_16_9, res = cozunurluk)


print(fviz_nbclust(df_pca, pam, method = "wss") +
        geom_vline(xintercept = k_sayisi, linetype = 2) +
        labs(title = "Elbow Metodu (PCA Sonras??)", 
             subtitle = "Optimal K??me Say??s?? (??lk 4 Bile??en ??zerinden)"))
dev.off()



set.seed(123)

pam_pca_fit <- pam(df_pca, k = k_sayisi) 



dosya_adi1 <- paste0("PCA_PAM_Cluster_Plot_k", k_sayisi, "_16x9.png")
png(dosya_adi1, width = w_16_9, height = h_16_9, res = cozunurluk)
print(fviz_cluster(pam_pca_fit, 
                   data = df_pca, # Veri olarak PCA skorlar??n?? veriyoruz
                   palette = "jco", 
                   ellipse.type = "convex", 
                   repel = TRUE,
                   ggtheme = theme_minimal(),
                   main = paste("PAM Kume Grafigi (PCA ile, k =", k_sayisi, ")")))
dev.off()



dosya_adi3 <- paste0("PCA_PAM_Silhouette_Plot_k", k_sayisi, "_16x9.png")
png(dosya_adi3, width = w_16_9, height = h_16_9, res = cozunurluk)
print(fviz_silhouette(pam_pca_fit, 
                      palette = "jco", 
                      ggtheme = theme_minimal(),
                      main = paste("Siluet Analizi (PCA Sonras??, k =", k_sayisi, ")")))
dev.off()



dosya_adi4 <- paste0("PCA_PAM_Heatmap_k", k_sayisi, "_16x9.png")


df_ordered <- df_scaled[order(pam_pca_fit$clustering), ]
annotation_row <- data.frame(Cluster = factor(sort(pam_pca_fit$clustering)))
rownames(annotation_row) <- rownames(df_ordered)

png(dosya_adi4, width = w_16_9, height = h_16_9, res = cozunurluk)
pheatmap(df_ordered, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         annotation_row = annotation_row,
         main = paste("Is?? Haritas?? (PCA Tabanl?? K??meler, k =", k_sayisi, ")"),
         display_numbers = FALSE, fontsize = 12)
dev.off()

print("--- PCA uygulanm???? PAM analizinin t??m ????kt??lar?? haz??r! ---")



if(!require(factoextra)) install.packages("factoextra")
library(factoextra)


data(state)
df <- as.data.frame(state.x77)

pca_res <- prcomp(df, scale = TRUE) 


w_16_9 <- 1600
h_16_9 <- 900
cozunurluk <- 150


dosya_adi <- "SADECE_PCA_Scree_Plot_16x9.png"

png(dosya_adi, width = w_16_9, height = h_16_9, res = cozunurluk)

# fviz_eig fonksiyonu bu i??in piri
print(fviz_eig(pca_res, 
               addlabels = TRUE,           # Y??zdeleri s??tunlar??n ??zerine yaz
               ylim = c(0, 60),            # Y ekseni aral?????? (Daha ferah dursun diye)
               barfill = "#4E84C4",        # ????k bir mavi tonu
               barcolor = "#4E84C4",
               linecolor = "red",          # K??r??lma ??izgisi k??rm??z?? olsun
               ncp = 8) +                  # T??m bile??enleri (8 tane) g??ster
        labs(title = "PCA Scree Plot (Yamac Grafigi)", 
             subtitle = "Her Bir Temel Bilesenin Ac??klad??g?? Varyans Oran??",
             x = "Temel Bilesenler (Dimensions)",
             y = "Ac??klanan Varyans Yuzdesi (%)") +
        theme_minimal() +
        theme(plot.title = element_text(size = 16, face = "bold"),
              plot.subtitle = element_text(size = 12)))

dev.off()

print(paste(dosya_adi, "ba??ar??yla kaydedildi hayat??m!"))


#----------------------------performans indexleri------------------------#

if(!require(cluster)) install.packages("cluster")
if(!require(fpc)) install.packages("fpc")           # Dunn i??in
if(!require(clusterSim)) install.packages("clusterSim") # DBI i??in

library(cluster)
library(fpc)
library(clusterSim)


data(state)
df <- as.data.frame(state.x77)
df_scaled <- scale(df)


dist_euclidean <- dist(df_scaled, method = "euclidean")
dist_manhattan <- dist(df_scaled, method = "manhattan")


results <- data.frame()


for(k in 2:5) {
  set.seed(123)
  pam_eucl <- pam(df_scaled, k = k, metric = "euclidean")
  
  
  sil_e <- mean(silhouette(pam_eucl$clustering, dist_euclidean)[, 3])
  dunn_e <- cluster.stats(dist_euclidean, pam_eucl$clustering)$dunn
  dbi_e <- index.DB(df_scaled, pam_eucl$clustering, d = dist_euclidean, centrotypes = "medoids")$DB
  
  
  set.seed(123)
  pam_manh <- pam(df_scaled, k = k, metric = "manhattan")
  
  
  sil_m <- mean(silhouette(pam_manh$clustering, dist_manhattan)[, 3])
  dunn_m <- cluster.stats(dist_manhattan, pam_manh$clustering)$dunn
  dbi_m <- index.DB(df_scaled, pam_manh$clustering, d = dist_manhattan, centrotypes = "medoids")$DB
  
  
  results <- rbind(results, 
                   data.frame(K = k, Metrik = "Oklid", Siluet = sil_e, Dunn = dunn_e, DBI = dbi_e),
                   data.frame(K = k, Metrik = "Manhattan", Siluet = sil_m, Dunn = dunn_m, DBI = dbi_m))
}


print("--- PAM Performans Analizi Sonu??lar?? ---")
print(results)



#devam??------------------------------------------------

if(!require(cluster)) install.packages("cluster")
if(!require(fpc)) install.packages("fpc")
if(!require(clusterSim)) install.packages("clusterSim")
if(!require(clValid)) install.packages("clValid")

library(cluster)
library(fpc)
library(clusterSim)
library(clValid)


data(state)
df_scaled <- scale(as.data.frame(state.x77))


dist_euc <- dist(df_scaled, method = "euclidean")
dist_man <- dist(df_scaled, method = "manhattan")


final_table <- data.frame()


for(k in 2:5) {
  
  # --- SENARYO A: ??KL??D ---
  set.seed(123)
  pam_e <- pam(df_scaled, k = k, metric = "euclidean")
  stats_e <- cluster.stats(dist_euc, pam_e$clustering)
  
  # DBI ve Connectivity 
  dbi_e <- index.DB(df_scaled, pam_e$clustering, d = dist_euc, centrotypes = "medoids")$DB
  conn_e <- connectivity(df_scaled, pam_e$clustering, method = "euclidean")
  
  row_e <- data.frame(K = k, Mesafe = "Oklid", 
                      Siluet = stats_e$avg.silwidth, 
                      Dunn = stats_e$dunn, 
                      DBI = dbi_e, 
                      CH = stats_e$ch, 
                      Conn = conn_e)
  
  # --- SENARYO B: MANHATTAN ---
  set.seed(123)
  pam_m <- pam(df_scaled, k = k, metric = "manhattan")
  stats_m <- cluster.stats(dist_man, pam_m$clustering)
  
  dbi_m <- index.DB(df_scaled, pam_m$clustering, d = dist_man, centrotypes = "medoids")$DB
  conn_m <- connectivity(df_scaled, pam_m$clustering, method = "manhattan")
  
  row_m <- data.frame(K = k, Mesafe = "Manhattan", 
                      Siluet = stats_m$avg.silwidth, 
                      Dunn = stats_m$dunn, 
                      DBI = dbi_m, 
                      CH = stats_m$ch, 
                      Conn = conn_m)
  
  final_table <- rbind(final_table, row_e, row_m)
}


formatted_table <- final_table
for(i in 3:7) {
  formatted_table[,i] <- format(round(final_table[,i], 3), decimal.mark = ",")
}


print("--- PERFORMANS KAR??ILA??TIRMA TABLOSU ---")
print(formatted_table, row.names = FALSE)



#--------------------------pca before after----------------------------------

if(!require(cluster)) install.packages("cluster")
if(!require(fpc)) install.packages("fpc")
if(!require(clusterSim)) install.packages("clusterSim")
if(!require(clValid)) install.packages("clValid")

library(cluster)
library(fpc)
library(clusterSim)
library(clValid)


data(state)
df_raw <- as.data.frame(state.x77)
df_scaled <- scale(df_raw)

# PCA Uygulama 
pca_res <- prcomp(df_scaled, scale = TRUE)
df_pca <- pca_res$x[, 1:3] 


comparison_results <- data.frame()


for(k in 2:5) {
  
  # --- PCA OLMADAN (HAM ??KL??D) ---
  set.seed(123)
  dist_raw <- dist(df_scaled, method = "euclidean")
  pam_raw <- pam(df_scaled, k = k, metric = "euclidean")
  
  stats_raw <- cluster.stats(dist_raw, pam_raw$clustering)
  dbi_raw <- index.DB(df_scaled, pam_raw$clustering, d = dist_raw, centrotypes = "medoids")$DB
  conn_raw <- connectivity(df_scaled, pam_raw$clustering, method = "euclidean")
  
  row_raw <- data.frame(K = k, Analiz = "PCA_Yok", 
                        Siluet = stats_raw$avg.silwidth, 
                        Dunn = stats_raw$dunn, 
                        DBI = dbi_raw, 
                        CH = stats_raw$ch, 
                        Conn = conn_raw)
  
  # ---PCA SONRASI (??KL??D) ---
  set.seed(123)
  dist_pca <- dist(df_pca, method = "euclidean")
  pam_pca <- pam(df_pca, k = k, metric = "euclidean")
  
  stats_pca <- cluster.stats(dist_pca, pam_pca$clustering)
  dbi_pca <- index.DB(df_pca, pam_pca$clustering, d = dist_pca, centrotypes = "medoids")$DB
  conn_pca <- connectivity(df_pca, pam_pca$clustering, method = "euclidean")
  
  row_pca <- data.frame(K = k, Analiz = "PCA_Var", 
                        Siluet = stats_pca$avg.silwidth, 
                        Dunn = stats_pca$dunn, 
                        DBI = dbi_pca, 
                        CH = stats_pca$ch, 
                        Conn = conn_pca)
  
  comparison_results <- rbind(comparison_results, row_raw, row_pca)
}


final_table <- comparison_results
for(i in 3:7) {
  final_table[,i] <- format(round(comparison_results[,i], 3), decimal.mark = ",")
}


print("--- ??KL??D MESAFES??: PCA ??NCES?? VE SONRASI KIYASLAMA ---")
print


