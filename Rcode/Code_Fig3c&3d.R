### below codes used for plotting Fig. 3c and 3d

load("./data_Fig3c&3d.RData")

library(ggplot2)
library(beanplot)
library(coin)
library(rcompanion)

pdf("./Figures/Fig.1.angio_range.pdf", useDingbats=FALSE, width=2, height=4)

beanplot(alt_range ~ hotspot, data = angio_pe_pd_clim_future_group_selected, ll = 0.00,
         main = "Elevation range (m)", ylab = "Elevation range (m)", side = "both", beanlines = "median",
         overallline = "mean",
         border = NA, col = list("#0073c2", "#efc000"))
legend("bottomleft", fill = c("#0073c2", "#efc000"),
       legend = c("hotspot", "nonhotspot" ))
dev.off()

pdf("./Figures/Fig.1.angio_bio1.pdf", useDingbats=FALSE, width=2, height=4)

beanplot(bio1_curr ~ hotspot, data = angio_pe_pd_clim_future_group_selected_selected, ll = 0.00,
         main = "MAT (oC)", ylab = "MAT (oC)", side = "both", beanlines = "median",
         overallline = "mean",
         border = NA, col = list("#0073c2", "#efc000"))
dev.off()

pdf("./Figures/Fig.1.angio_bio12.pdf", useDingbats=FALSE, width=2, height=4)
beanplot(bio12_curr ~ hotspot, data = angio_pe_pd_clim_future_group_selected, ll = 0.00,
         main = "AP (mm)", ylab = "AP (mm)", side = "both", beanlines = "median",
         # overallline = "mean",
         border = NA, col = list("#0073c2", "#efc000"))
dev.off()

pdf("./Figures/Fig.1.angio_miocene_bio1.pdf", useDingbats=FALSE, width=2, height=4)
beanplot(mio_temp_anom ~ hotspot, data = angio_pe_pd_clim_future_group_selected, ll = 0.00,
         main = "Miocene MAT anomaly", ylab = "MAT (oC)", side = "both", beanlines = "median",
         overallline = "mean",
         border = NA, col = list("#0073c2", "#efc000"))

dev.off()

pdf("./Figures/Fig.1.angio_miocene_bio12.pdf", useDingbats=FALSE, width=2, height=4)
beanplot(mio_prec_anom ~ hotspot, data = angio_pe_pd_clim_future_group_selected, ll = 0.00,
         main = "Miocene AP anomaly", ylab = "AP (mm)", side = "both", beanlines = "median",
         overallline = "mean",
         border = NA, col = list("#0073c2", "#efc000"))
dev.off()

pdf("./Figures/Fig.1.angio_lgm_bio1.pdf", useDingbats=FALSE, width=2, height=4)
beanplot(lgm_temp_anom ~ hotspot, data = angio_pe_pd_clim_future_group_selected, ll = 0.00,
         main = "LGM MAT anomaly", ylab = "MAT (oC)", side = "both", beanlines = "median",
         overallline = "mean",
         border = NA, col = list("#0073c2", "#efc000"))
dev.off()

pdf("./Figures/Fig.1.angio_lgm_bio11.pdf", useDingbats=FALSE, width=2, height=4)
beanplot(lgm_prec_anom ~ hotspot, data = angio_pe_pd_clim_future_group_selected, ll = 0.00,
         main = "LGM AP anomaly", ylab = "AP (mm)", side = "both", beanlines = "median",
         overallline = "mean",
         border = NA, col = list("#0073c2", "#efc000"))
dev.off()

### gymno

pdf("./Figures/Fig.1.gymno_range.pdf", useDingbats=FALSE, width=2, height=4)

beanplot(alt_range ~ hotspot, data = gymno_pe_clim_future_group_selected, ll = 0.00,
         main = "Elevation range (m)", ylab = "Elevation range (m)", side = "both", beanlines = "median",
         overallline = "mean",
         border = NA, col = list("#0073c2", "#efc000"))
dev.off()

pdf("./Figures/Fig.1.gymno_bio1.pdf", useDingbats=FALSE, width=2, height=4)

beanplot(bio1_curr ~ hotspot, data = gymno_pe_clim_future_group_selected, ll = 0.00,
         main = "MAT (oC)", ylab = "MAT (oC)", side = "both", beanlines = "median",
         overallline = "mean",
         border = NA, col = list("#0073c2", "#efc000"))
dev.off()

pdf("./Figures/Fig.1.gymno_bio12.pdf", useDingbats=FALSE, width=2, height=4)
beanplot(bio12_curr ~ hotspot, data = gymno_pe_clim_future_group_selected, ll = 0.00,
         main = "AP (mm)", ylab = "AP (mm)", side = "both", beanlines = "median",
         # overallline = "mean",
         border = NA, col = list("#0073c2", "#efc000"))
dev.off()

pdf("./Figures/Fig.1.gymno_miocene_bio1_2.pdf", useDingbats=FALSE, width=2, height=4)
beanplot(mio_temp_anom ~ hotspot, data = gymno_pe_clim_future_group_selected, ll = 0.00,
         main = "Miocene MAT anomaly", ylab = "MAT (oC)", side = "both", beanlines = "median",
         overallline = "mean",
         border = NA, col = list("#0073c2", "#efc000"))

dev.off()

pdf("./Figures/Fig.1.gymno_miocene_bio12.pdf", useDingbats=FALSE, width=2, height=4)
beanplot(mio_prec_anom ~ hotspot, data = gymno_pe_clim_future_group_selected, ll = 0.00,
         main = "Miocene AP anomaly", ylab = "AP (mm)", side = "both", beanlines = "median",
         overallline = "mean",
         border = NA, col = list("#0073c2", "#efc000"))
dev.off()

pdf("./Figures/Fig.1.gymno_lgm_bio1.pdf", useDingbats=FALSE, width=2, height=4)
beanplot(lgm_temp_anom ~ hotspot, data = gymno_pe_clim_future_group_selected, ll = 0.00,
         main = "LGM MAT anomaly", ylab = "MAT (oC)", side = "both", beanlines = "median",
         overallline = "mean",
         border = NA, col = list("#0073c2", "#efc000"))
dev.off()

pdf("./Figures/Fig.1.gymno_lgm_bio11.pdf", useDingbats=FALSE, width=2, height=4)
beanplot(lgm_prec_anom ~ hotspot, data = gymno_pe_clim_future_group_selected, ll = 0.00,
         main = "LGM AP anomaly", ylab = "AP (mm)", side = "both", beanlines = "median",
         overallline = "mean",
         border = NA, col = list("#0073c2", "#efc000"))
dev.off()

######## angiosperm: hotspot vs. non-hotspot significance analysis, fig. 3c p value##################

## elevation range
angio_pe_pd_clim_future_group_selected$hotspot <- as.factor(angio_pe_pd_clim_future_group_selected$hotspot)
sig_test_angio_canape_ele_range <- independence_test(alt_range ~ hotspot, data = angio_pe_pd_clim_future_group_selected,
                                                     distribution = approximate(nresample = 10000))
pvalue(sig_test_angio_canape_ele_range, method = "global")

PT_angio_canape_ele_range = pairwisePermutationTest(alt_range ~ hotspot, data = angio_pe_pd_clim_future_group_selected,
                                                    method="fdr")

PT_angio_canape_ele_range
## display letters
cldList(p.adjust ~ Comparison,
        data = PT_angio_canape_ele_range,
        threshold  = 0.05)

## MAT
sig_test_angio_canape_mat <- independence_test(bio1_curr~ hotspot, data = angio_pe_pd_clim_future_group_selected,
                                               distribution = approximate(nresample = 10000))
pvalue(sig_test_angio_canape_mat, method = "global")

PT_angio_canape_mat = pairwisePermutationTest(bio1_curr  ~ hotspot, data = angio_pe_pd_clim_future_group_selected,
                                              method="fdr")

PT_angio_canape_mat
## display letters
cldList(p.adjust ~ Comparison,
        data = PT_angio_canape_mat,
        threshold  = 0.05)

## AP
sig_test_angio_canape_ap <- independence_test(bio12_curr~ hotspot, data = angio_pe_pd_clim_future_group_selected,
                                              distribution = approximate(nresample = 10000))
pvalue(sig_test_angio_canape_ap, method = "global")

PT_angio_canape_ap = pairwisePermutationTest(bio12_curr  ~ hotspot, data = angio_pe_pd_clim_future_group_selected,
                                             method="fdr")

PT_angio_canape_ap
## display letters
cldList(p.adjust ~ Comparison,
        data = PT_angio_canape_ap,
        threshold  = 0.05)

## mio mat
sig_test_angio_canape_mio_mat <- independence_test(mio_temp_anom~ hotspot, data = angio_pe_pd_clim_future_group_selected,
                                                   distribution = approximate(nresample = 10000))
pvalue(sig_test_angio_canape_mio_mat, method = "global")

PT_angio_canape_mio_mat = pairwisePermutationTest(mio_temp_anom  ~ hotspot, data = angio_pe_pd_clim_future_group_selected,
                                                  method="fdr")

PT_angio_canape_mio_mat
## display letters
cldList(p.adjust ~ Comparison,
        data = PT_angio_canape_mio_mat,
        threshold  = 0.05)

## mio ap
sig_test_angio_canape_mio_ap <- independence_test(mio_prec_anom~ hotspot, data = angio_pe_pd_clim_future_group_selected,
                                                  distribution = approximate(nresample = 10000))
pvalue(sig_test_angio_canape_mio_ap, method = "global")

PT_angio_canape_mio_ap = pairwisePermutationTest(mio_prec_anom  ~ hotspot, data = angio_pe_pd_clim_future_group_selected,
                                                 method="fdr")

PT_angio_canape_mio_ap
## display letters
cldList(p.adjust ~ Comparison,
        data = PT_angio_canape_mio_ap,
        threshold  = 0.05)

## lgm mat
sig_test_angio_canape_lgm_mat <- independence_test(lgm_temp_anom~ hotspot, data = angio_pe_pd_clim_future_group_selected,
                                                   distribution = approximate(nresample = 10000))
pvalue(sig_test_angio_canape_lgm_mat, method = "global")

PT_angio_canape_lgm_mat = pairwisePermutationTest(lgm_temp_anom  ~ hotspot, data = angio_pe_pd_clim_future_group_selected,
                                                  method="fdr")

PT_angio_canape_lgm_mat
## display letters
cldList(p.adjust ~ Comparison,
        data = PT_angio_canape_lgm_mat,
        threshold  = 0.05)

## lgm ap
sig_test_angio_canape_lgm_ap <- independence_test(lgm_prec_anom~ hotspot, data = angio_pe_pd_clim_future_group_selected,
                                                  distribution = approximate(nresample = 10000))
pvalue(sig_test_angio_canape_lgm_ap, method = "global")

PT_angio_canape_lgm_ap = pairwisePermutationTest(lgm_prec_anom  ~ hotspot, data = angio_pe_pd_clim_future_group_selected,
                                                 method="fdr")

PT_angio_canape_lgm_ap
## display letters
cldList(p.adjust ~ Comparison,
        data = PT_angio_canape_lgm_ap,
        threshold  = 0.05)


###########gymnosperm: hotspot vs. non-hotspot significance analysis, fig. 4 p value############
##### imputation test for Fig. hotspot driver #############
### gymnosperm
## elevation range
gymno_pe_clim_future_group_selected$hotspot <- as.factor(gymno_pe_clim_future_group_selected$hotspot)
sig_test_gymno_canape_ele_range <- independence_test(alt_range ~ hotspot, data = gymno_pe_clim_future_group_selected,
                                                     distribution = approximate(nresample = 10000))
pvalue(sig_test_gymno_canape_ele_range, method = "global")

PT_gymno_canape_ele_range = pairwisePermutationTest(alt_range ~ hotspot, data = gymno_pe_clim_future_group_selected,
                                                    method="fdr")

PT_gymno_canape_ele_range
## display letters
cldList(p.adjust ~ Comparison,
        data = PT_gymno_canape_ele_range,
        threshold  = 0.05)

## MAT
sig_test_gymno_canape_mat <- independence_test(bio1_curr~ hotspot, data = gymno_pe_clim_future_group_selected,
                                               distribution = approximate(nresample = 10000))
pvalue(sig_test_gymno_canape_mat, method = "global")

PT_gymno_canape_mat = pairwisePermutationTest(bio1_curr  ~ hotspot, data = gymno_pe_clim_future_group_selected,
                                              method="fdr")

PT_gymno_canape_mat
## display letters
cldList(p.adjust ~ Comparison,
        data = PT_gymno_canape_mat,
        threshold  = 0.05)

## AP
sig_test_gymno_canape_ap <- independence_test(bio12_curr~ hotspot, data = gymno_pe_clim_future_group_selected,
                                              distribution = approximate(nresample = 10000))
pvalue(sig_test_gymno_canape_ap, method = "global")

PT_gymno_canape_ap = pairwisePermutationTest(bio12_curr  ~ hotspot, data = gymno_pe_clim_future_group_selected,
                                             method="fdr")

PT_gymno_canape_ap
## display letters
cldList(p.adjust ~ Comparison,
        data = PT_gymno_canape_ap,
        threshold  = 0.05)

## mio mat
sig_test_gymno_canape_mio_mat <- independence_test(mio_temp_anom~ hotspot, data = gymno_pe_clim_future_group_selected,
                                                   distribution = approximate(nresample = 10000))
pvalue(sig_test_gymno_canape_mio_mat, method = "global")

PT_gymno_canape_mio_mat = pairwisePermutationTest(mio_temp_anom  ~ hotspot, data = gymno_pe_clim_future_group_selected,
                                                  method="fdr")

PT_gymno_canape_mio_mat
## display letters
cldList(p.adjust ~ Comparison,
        data = PT_gymno_canape_mio_mat,
        threshold  = 0.05)

## mio ap
sig_test_gymno_canape_mio_ap <- independence_test(mio_prec_anom~ hotspot, data = gymno_pe_clim_future_group_selected,
                                                  distribution = approximate(nresample = 10000))
pvalue(sig_test_gymno_canape_mio_ap, method = "global")

PT_gymno_canape_mio_ap = pairwisePermutationTest(mio_prec_anom  ~ hotspot, data = gymno_pe_clim_future_group_selected,
                                                 method="fdr")

PT_gymno_canape_mio_ap
## display letters
cldList(p.adjust ~ Comparison,
        data = PT_gymno_canape_mio_ap,
        threshold  = 0.05)

## lgm mat
sig_test_gymno_canape_lgm_mat <- independence_test(lgm_temp_anom~ hotspot, data = gymno_pe_clim_future_group_selected,
                                                   distribution = approximate(nresample = 10000))
pvalue(sig_test_gymno_canape_lgm_mat, method = "global")

PT_gymno_canape_lgm_mat = pairwisePermutationTest(lgm_temp_anom  ~ hotspot, data = gymno_pe_clim_future_group_selected,
                                                  method="fdr")

PT_gymno_canape_lgm_mat
## display letters
cldList(p.adjust ~ Comparison,
        data = PT_gymno_canape_lgm_mat,
        threshold  = 0.05)

## lgm ap
sig_test_gymno_canape_lgm_ap <- independence_test(lgm_prec_anom~ hotspot, data = gymno_pe_clim_future_group_selected,
                                                  distribution = approximate(nresample = 10000))
pvalue(sig_test_gymno_canape_lgm_ap, method = "global")

PT_gymno_canape_lgm_ap = pairwisePermutationTest(lgm_prec_anom  ~ hotspot, data = gymno_pe_clim_future_group_selected,
                                                 method="fdr")

PT_gymno_canape_lgm_ap
## display letters
cldList(p.adjust ~ Comparison,
        data = PT_gymno_canape_lgm_ap,
        threshold  = 0.05)
