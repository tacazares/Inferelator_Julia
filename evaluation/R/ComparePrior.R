library(ggplot2)

prior1 <- read.table("/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Outs/TRAC_loop/Prior/Prior_sum.tsv")
prior2 <- read.table("/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Outs/Seurat/Prior/MEMT_050723_FIMOp5_b.tsv")

index <- intersect(rownames(prior1), rownames(prior2))
prior1 <- prior1[index,]
prior2 <- prior2[index,]

prior_df <- as.matrix(rep(0, length(colnames(prior1))))
prior_df <- cbind(prior_df, rep(0, length(colnames(prior2))))
rownames(prior_df) <- colnames(prior1)
colnames(prior_df) <- c("TRAC","Body")

for(i in 1:length(colnames(prior1))){
    prior1_targets <- length(which(prior1[,i] > 0))
    prior2_targets <- length(which(prior2[,i] > 0))
    prior_df[i,1] <- prior1_targets
    prior_df[i,2] <- prior2_targets
}

prior_df <- as.data.frame(prior_df)
pdf("Prior_Scatter.pdf", width = 8, height = 8)
ggplot(prior_df, aes(x = TRAC, y = Body)) + geom_point() + geom_abline(slope=1, intercept=0)
dev.off()

prior_df_hist <- data.frame(Freq = c(prior_df[,1], prior_df[,2]))
prior_df_hist$TF <- c(rownames(prior_df), rownames(prior_df))
prior_df_hist$Prior <- c(rep("Prior1", length(colnames(prior1))), rep("Prior2", length(colnames(prior2))))

pdf("Prior_hist.pdf", width = 8, height = 8)
ggplot(prior_df_hist, aes(x = Freq, color = Prior)) + geom_histogram(fill = "white", position = "identity", alpha = 0.7)
dev.off()