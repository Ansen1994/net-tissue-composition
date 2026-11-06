# ==============================================================================
# Table S5c — proportional-hazards test for the joint four-gene Cox model
# Input:  Table_S3_ML_Input_Data_ER_Positive_Ultimate.csv (n = 762, 99 events)
# Output: Table_S5c_PH_Test.csv
# Scaled Schoenfeld residuals (survival::cox.zph) on the model of Table S5.
# ==============================================================================
suppressMessages(library(survival))

d0 <- read.csv("Table_S3_ML_Input_Data_ER_Positive_Ultimate.csv",
               row.names = 1, check.names = FALSE)
genes <- c("PLG", "CD79A", "MMP7", "ANGPT1")
d <- d0[, c("OS.time", "OS", genes)]
d <- d[complete.cases(d) & d$OS.time > 0, ]

mf <- coxph(Surv(OS.time, OS) ~ PLG + CD79A + MMP7 + ANGPT1, data = d)
z  <- cox.zph(mf)
tb <- z$table

out <- data.frame(Term       = rownames(tb),
                  Chi_square = round(as.numeric(tb[, "chisq"]), 4),
                  df         = as.integer(tb[, "df"]),
                  P_value    = signif(as.numeric(tb[, "p"]), 4),
                  row.names  = NULL)

write.csv(out, "Table_S5c_PH_Test.csv", row.names = FALSE, quote = FALSE)
cat("n =", nrow(d), " events =", sum(d$OS), "\n")
print(out, row.names = FALSE)
