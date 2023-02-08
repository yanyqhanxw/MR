library(TwoSampleMR)
library(RadialMR)
library(MRPRESSO)

############################################################
# load data
# all of the SNPs are significantly associated with exposure
############################################################
exposure = read_exposure_data(filename          = "smoke_initiation",
                              sep               = "\t",
                              clump             = FALSE,
                              snp_col           = "SNP",
                              beta_col          = "Beta",
                              se_col            = "SE",
                              eaf_col           = "EAF",
                              effect_allele_col = "EA",
                              other_allele_col  = "OA",
                              pval_col          = "P",
                              samplesize_col    = "N",
                              chr_col           = "Chr",
                              pos_col           = "Pos_hg38")
                                      
exposure$exposure = "smoke_initiation"

outcome = read_outcome_data(filename          = "RA",
                            sep               = "\t",
                            snp_col           = "SNP",
                            beta_col          = "Beta",
                            se_col            = "SE",
                            eaf_col           = "EAF",
                            effect_allele_col = "EA",
                            other_allele_col  = "OA",
                            pval_col          = "P",
                            samplesize_col    = "N",
                            chr_col           = "Chr",
                            pos_col           = "Pos_hg38") 
outcome$outcome = "RA"


############################################################
# select instruments
############################################################
# clump
clump_dat = clump_data(exposure, clump_r2 = 0.01, clump_kb = 10)

# exclude SNPs directly associated with outcome
outcome = outcome[outcome$pval.outcome>5e-8,]

# harmonise allele
# note:
#   We dropped all palindromic SNPs in two-step MR analysis (action = 3).
#   Because there is no eaf data in GWAS summary for cytokines.
harmonised = harmonise_data(exposure_dat = clump_dat, outcome_dat = outcome, action = 2) # 协调暴露和结果之间的等位基因，不去除回文SNP
harmonised = harmonised[harmonised$mr_keep==TRUE,]

# RadialMR: exclude outlier pleiotropic SNPs
radial_dat = format_radial(harmonised$beta.exposure, harmonised$beta.outcome,
                           harmonised$se.exposure, harmonised$se.outcome,
                           harmonised$SNP)

res1 = ivw_radial(radial_dat,0.05,1,0.0001)
res2 = egger_radial(radial_dat,0.05,1)

remove1 = as.matrix(cbind(res1$outliers[1]))
remove2 = as.matrix(cbind(res2$outliers[1]))
remove_radial = as.data.frame(rbind(remove1, remove2))
remove_radial = unique(remove_radial)

MR_dat = harmonised[!(harmonised$SNP%in%remove_radial$SNP),]


############################################################
# calculate F statistics and exclude SNPs with F < 10
############################################################
tmp = 2*(MR_dat$beta.exposure**2)*MR_dat$eaf.exposure*(1-MR_dat$eaf.exposure)
MR_dat$R2 = tmp/(tmp + (MR_dat$se.exposure**2)*2*MR_dat$samplesize.exposure*MR_dat$eaf.exposure*(1-MR_dat$eaf.exposure))
# F statistics of single SNP
MR_dat$F_stat = MR_dat$R2*(MR_dat$samplesize.exposure-2)/(1-MR_dat$R2)
# total F statistics of all SNPs
MR_dat$F_stat_total = sum(MR_dat$R2)*(MR_dat$samplesize.exposure-nrow(MR_dat)-1)/nrow(MR_dat)/(1-sum(MR_dat$R2))
MR_dat = MR_dat[MR_dat$F_stat>10,]

############################################################
# MR alalysis
# method: IVW, WM, MR Egger, MRPRESSO
# reverse causal effect test: MR Steiger
############################################################
MR_result = mr(MR_dat, method_list=c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))
MR_PRESSO_res = mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                          SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                          OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                          data = MR_dat[MR_dat$mr_keep == T,])
Steiger = directionality_test(MR_dat)


############################################################
# sensitivity analysis
############################################################
MR_ple = mr_pleiotropy_test(MR_dat)
MR_he = mr_heterogeneity(MR_dat, method_list=c("mr_egger_regression","mr_ivw"))


############################################################
# calculate odds ratio
############################################################
MR_OR = generate_odds_ratios(MR_result) # 计算OR值
MR_PRESSO_OR = do.call(rbind,MR_PRESSO[,1])
MR_PRESSO_OR = MR_PRESSO_OR[-c(2,4,6,8,10),]
MR_PRESSO_OR$or = c(0)
MR_PRESSO_OR$or_lci95 = c(0)
MR_PRESSO_OR$or_uci95 = c(0)
for(i in 1:nrow(MR_PRESSO_OR)){
  MR_PRESSO_OR$or[i] = exp(MR_PRESSO_OR$`Causal Estimate`[i])
  MR_PRESSO_OR$or_lci95[i] = exp(MR_PRESSO_OR$`Causal Estimate`[i] - 1.96 * MR_PRESSO_OR$Sd[i])
  MR_PRESSO_OR$or_uci95[i] = exp(MR_PRESSO_OR$`Causal Estimate`[i] + 1.96 * MR_PRESSO_OR$Sd[i])
}

