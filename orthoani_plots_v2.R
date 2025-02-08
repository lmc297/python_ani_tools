################ load required packages
library(ggplot2)
library(ggdist)
library(stringr)
library(khroma)
library(gridExtra)
library(forcats)

################ load/prepare data

# create read.anires function
# reads data in Supplementary Tables S1-S6
read.anires <- function(f, method){
  
  # read file
  df <- read.delim(file = f, header = F, sep = "\t", stringsAsFactors = F, check.names = F, skip = 1)
  # add column names
  colnames(df) <- c("run", "genome", "ANI")
  
  # create new columns
  df$dataset <- unlist(lapply(strsplit(x = df$run, split = "_"), "[[", 1))
  df$method <- unlist(lapply(strsplit(x = df$run, split = "_"), "[[", 2))
  df$cpus <- unlist(lapply(strsplit(x = df$run, split = "_"), "[[", 3))
  df$cpus <- gsub(pattern = "thread", replacement = "", x = df$cpus)
  df$trial <- unlist(lapply(strsplit(x = df$run, split = "_"), "[[", 4))
  df$trial <- gsub(pattern = "trial", replacement = "", x = df$trial)
  df$matchcol <- paste(gsub(pattern = paste(method, "_", sep = ""), replacement = "", x = df$run), df$genome, sep = ".")
  
  # remove genomes that didn't produce an ANI value
  df <- df[!is.na(df$ANI),]
  
  # return final dataframe
  return(df)
}

# Load OrthoANI data
orthoani <- read.anires(f = "SuppTableS1_orthoani/SupplementaryTableS1_orthoani_final.tsv", method = "orthoani")
head(orthoani)
table(orthoani$trial)

# Load FastANI data
fastani <- read.anires(f = "SuppTableS2_fastani/SupplementaryTableS2_fastani_final.tsv", method = "fastani")
head(fastani)
table(fastani$trial)

# Load skani data
skani <- read.anires(f = "SuppTableS3_skani/SupplementaryTableS3_skani_final.tsv", method = "skani")
head(skani)
table(skani$trial)

# Load PyOrthoANI data
pyorthoani <- read.anires(f = "SuppTableS4_pyorthoani/SupplementaryTableS4_pyorthoani_final.tsv", method = "pyorthoani")
head(pyorthoani)
table(pyorthoani$trial)

# Load PyFastANI data
revpyfastani <- read.anires(f = "SuppTableS5_pyfastani/SupplementaryTableS5_revpyfastani_final.tsv", method = "revpyfastani")
head(revpyfastani)
table(revpyfastani$trial)

# Load Pyskani data
revpyskani <- read.anires(f = "SuppTableS6_pyskani/SupplementaryTableS6_revpyskani_final.tsv", method = "revpyskani")
head(revpyskani)
table(revpyskani$trial)

################ Compare OrthoANI vs PyOrthoANI

# remove genomes that did not produce ANI values for >= 1 method
table(orthoani$matchcol%in%pyorthoani$matchcol)
table(pyorthoani$matchcol%in%orthoani$matchcol)

pyorthoani.reduced <- pyorthoani[which(pyorthoani$matchcol%in%orthoani$matchcol),]

table(orthoani$matchcol%in%pyorthoani.reduced$matchcol)
table(pyorthoani.reduced$matchcol%in%orthoani$matchcol)

table(orthoani$matchcol==pyorthoani.reduced$matchcol)
table(pyorthoani.reduced$matchcol%in%orthoani$matchcol)

pyorthoani.reduced <- pyorthoani.reduced[match(orthoani$matchcol, pyorthoani.reduced$matchcol),]
table(orthoani$matchcol==pyorthoani.reduced$matchcol)
table(pyorthoani.reduced$matchcol%in%orthoani$matchcol)

# create data frame
compare.orthoani <- data.frame(orthoani$dataset, orthoani$matchcol, orthoani$ANI, pyorthoani.reduced$ANI)
colnames(compare.orthoani) <- c("Dataset", "Genome", "OrthoANI", "PyOrthoANI")
head(compare.orthoani)

# fit linear model
orthoani.lm <- lm(formula = compare.orthoani$PyOrthoANI ~ compare.orthoani$OrthoANI)
summary(orthoani.lm)

# plot results
plot.orthoani <-
ggplot(data = compare.orthoani, mapping = aes(x = OrthoANI, y = PyOrthoANI)) +
  geom_point(alpha = 0.1) +
  theme_bw() + 
  geom_abline(slope = orthoani.lm$coefficients[2], 
              intercept = orthoani.lm$coefficients[1],
              linetype = "dashed", color = "#4477AA") + 
  labs(title = paste("Adj R2 = ",signif(summary(orthoani.lm)$adj.r.squared, 5),
                     " P =",signif(summary(orthoani.lm)$coef[2,4], 5)))
plot.orthoani

################ Compare FastANI vs PyFastANI

# remove genomes that did not produce ANI values for >= 1 method
table(fastani$matchcol%in%revpyfastani$matchcol)
table(revpyfastani$matchcol%in%fastani$matchcol)

table(fastani$matchcol==revpyfastani$matchcol)
table(revpyfastani$matchcol==fastani$matchcol)

revpyfastani <- revpyfastani[match(fastani$matchcol, revpyfastani$matchcol),]
table(fastani$matchcol==revpyfastani$matchcol)
table(revpyfastani$matchcol==fastani$matchcol)

# create data frame
compare.fastani <- data.frame(fastani$dataset, fastani$matchcol, fastani$ANI, revpyfastani$ANI)
colnames(compare.fastani) <- c("Dataset", "Genome", "FastANI", "PyFastANI")
head(compare.fastani)

# fit linear model
fastani.lm <- lm(formula = compare.fastani$PyFastANI ~ compare.fastani$FastANI)
summary(fastani.lm)

# plot results
plot.fastani <-
  ggplot(data = compare.fastani, mapping = aes(x = FastANI, y = PyFastANI)) +
  geom_point(alpha = 0.1) +
  theme_bw() + 
  geom_abline(slope = fastani.lm$coefficients[2], 
              intercept = fastani.lm$coefficients[1],
              linetype = "dashed", color = "#AA3377") + 
  labs(title = paste("Adj R2 = ",signif(summary(fastani.lm)$adj.r.squared, 5),
                     " P =",signif(summary(fastani.lm)$coef[2,4], 5)))
plot.fastani

################ Compare skani vs RevPyskani

# remove genomes that did not produce ANI values for >= 1 method

table(skani$matchcol%in%revpyskani$matchcol)
table(revpyskani$matchcol%in%skani$matchcol)

skani.reduced <- skani[which(skani$matchcol%in%revpyskani$matchcol),]
table(skani.reduced$matchcol%in%revpyskani$matchcol)
table(revpyskani$matchcol%in%skani.reduced$matchcol)

revpyskani.reduced <- revpyskani[which(revpyskani$matchcol%in%skani.reduced$matchcol),]
table(skani.reduced$matchcol%in%revpyskani.reduced$matchcol)
table(revpyskani.reduced$matchcol%in%skani.reduced$matchcol)

table(skani.reduced$matchcol==revpyskani.reduced$matchcol)
table(revpyskani.reduced$matchcol==skani.reduced$matchcol)

revpyskani.reduced <- revpyskani.reduced[match(skani.reduced$matchcol, revpyskani.reduced$matchcol),]
table(skani.reduced$matchcol==revpyskani.reduced$matchcol)
table(revpyskani.reduced$matchcol==skani.reduced$matchcol)

# create data frame
compare.skani <- data.frame(skani.reduced$dataset, skani.reduced$matchcol, skani.reduced$ANI, (revpyskani.reduced$ANI*100))
colnames(compare.skani) <- c("Dataset", "Genome", "skani", "Pyskani")
head(compare.skani)

# fit linear model
skani.lm <- lm(formula = compare.skani$Pyskani ~ compare.skani$skani)
summary(skani.lm)

# plot results
plot.skani <-
  ggplot(data = compare.skani, mapping = aes(x = skani, y = Pyskani)) +
  geom_point(alpha = 0.1) +
  theme_bw() + 
  geom_abline(slope = skani.lm$coefficients[2], 
              intercept = skani.lm$coefficients[1],
              linetype = "dashed", color = "#228833") + 
  labs(title = paste("Adj R2 = ",signif(summary(skani.lm)$adj.r.squared, 5),
                     " P =",signif(summary(skani.lm)$coef[2,4], 5))) +
  ylab("Pyskani*100")
plot.skani

################ speed benchmark

# load Nextflow trace file
trace <- read.delim(file = "SuppTableS7_trace/SupplementaryTableS7_trace_final.tsv",
                    header = T, sep = "\t",
                    stringsAsFactors = F,
                    check.names = F, skip = 1)
head(trace)
table(trace$status)

# subset all runs that completed
trace <- trace[which(trace$status=="COMPLETED"),]
table(trace$status)

# add columns
trace$dataset <- unlist(lapply(strsplit(x = trace$run, split = "_"), "[[", 1))
trace$method <- unlist(lapply(strsplit(x = trace$run, split = "_"), "[[", 2))
trace$cpus <- unlist(lapply(strsplit(x = trace$run, split = "_"), "[[", 3))
trace$cpus <- gsub(pattern = "thread", replacement = "", x = trace$cpus)
trace$trial <- unlist(lapply(strsplit(x = trace$run, split = "_"), "[[", 4))
trace$trial <- gsub(pattern = "trial", replacement = "", x = trace$trial)
head(trace)

# remove runs with "pyskani" or "pyfastani" method
# these were done with query and reference genomes reversed to sanity check
# revpyfastani and revpyskani have the correct query and reference order and will be used
table(trace$method)
trace <- trace[which(trace$method!="pyfastani"),]
trace <- trace[which(trace$method!="pyskani"),]

# sanity checks
table(trace$dataset)
table(trace$method)
table(trace$cpus)
table(trace$trial)

# reorder factors
table(trace$cpus)
trace$cpus <- factor(trace$cpus, levels = c("1", "8", "16"))

table(trace$method)
trace$method <- factor(trace$method, levels = c("orthoani", "pyorthoani", "fastani", "revpyfastani", "skani", "revpyskani"))

# sanity check skani: memory
trace.skani <- droplevels(trace[which(trace$method=="skani"),])
table(trace.skani$method)
mem.box.skani <-
  ggplot(data = trace.skani, mapping = aes(x = cpus, y = peak_rss_MB, color = method, fill = method)) + 
  geom_violin(alpha = 0.3) +
  theme_bw() +
  scale_color_manual(values = c(skani = "#CCDDAA")) + 
  scale_fill_manual(values = (skani = "#CCDDAA")) + 
  ylab("Peak RSS (MB)") + 
  xlab("Number of CPUs")
mem.box.skani

# sanity check skani: realtime
time.box.skani <-
  ggplot(data = trace.skani, mapping = aes(x = cpus, y = realtime_seconds, color = method, fill = method)) + 
  geom_violin(alpha = 0.3) +
  theme_bw() +
  scale_color_manual(values = c(skani = "#CCDDAA")) + 
  scale_fill_manual(values = (skani = "#CCDDAA")) + 
  ylab("Real time (seconds)") + 
  xlab("Number of CPUs")
time.box.skani

# plot skani sanity check
pdf(file = "skani_sanity_memory_time.pdf", width = 8.5, height = 4)
grid.arrange(time.box.skani, mem.box.skani, ncol=2)
dev.off()

# full data set: memory violin plot (base10 log scale)
mem.box <-
  ggplot(data = trace, mapping = aes(x = cpus, y = peak_rss_MB, color = method, fill = method)) + 
  geom_violin(alpha = 0.5) +
  theme_bw() +
  scale_color_manual(values = c(pyorthoani = "#4477AA", orthoani = "#BBCCEE",
                     revpyfastani = "#AA3377", fastani = "#FFCCCC",
                     revpyskani = "#228833", skani = "#CCDDAA")) + 
  scale_fill_manual(values = c(pyorthoani = "#4477AA", orthoani = "#BBCCEE",
                                revpyfastani = "#AA3377", fastani = "#FFCCCC",
                                revpyskani = "#228833", skani = "#CCDDAA")) + 
  ylab("Peak RSS (MB)") + 
  xlab("Number of CPUs") + 
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))
mem.box

# full data set: memory ridge plot (base10 log scale)
mem.ridges <- 
  ggplot(data = trace, mapping = aes(x = peak_rss_MB, y = fct_rev(cpus), color = method, fill = method)) + 
  coord_cartesian(clip = "off") +
  scale_y_discrete(expand = c(.07, .07)) +
  scale_color_manual(values = c(pyorthoani = "#4477AA", orthoani = "#BBCCEE",
                                revpyfastani = "#AA3377", fastani = "#FFCCCC",
                                revpyskani = "#228833", skani = "#CCDDAA")) + 
  scale_fill_manual(values = c(pyorthoani = "#4477AA", orthoani = "#BBCCEE",
                                revpyfastani = "#AA3377", fastani = "#FFCCCC",
                                revpyskani = "#228833", skani = "#CCDDAA")) +

  ggridges::geom_density_ridges(
    alpha = .7) + 
  theme_bw() + 
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10)) + 
  ylab("Number of CPUs") + 
  xlab("Peak RSS (MB)") 
  
mem.ridges

# full data set: time violin plot (base10 log scale)
time.box <-
  ggplot(data = trace, mapping = aes(x = cpus, y = realtime_seconds, color = method, fill = method)) + 
  geom_violin(alpha = 0.5) +
  theme_bw() +
  scale_color_manual(values = c(pyorthoani = "#4477AA", orthoani = "#BBCCEE",
                                revpyfastani = "#AA3377", fastani = "#FFCCCC",
                                revpyskani = "#228833", skani = "#CCDDAA")) + 
  scale_fill_manual(values = c(pyorthoani = "#4477AA", orthoani = "#BBCCEE",
                               revpyfastani = "#AA3377", fastani = "#FFCCCC",
                               revpyskani = "#228833", skani = "#CCDDAA")) +
  ylab("Real time (seconds)") + 
  xlab("Number of CPUs") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))
time.box

# full data set: time ridge plot (base10 log scale)
time.ridges <- 
  ggplot(data = trace, mapping = aes(x = realtime_seconds, y = fct_rev(cpus), color = method, fill = method)) + 
  coord_cartesian(clip = "off") +
  scale_y_discrete(expand = c(.07, .07)) +
  scale_color_manual(values = c(pyorthoani = "#4477AA", orthoani = "#BBCCEE",
                                revpyfastani = "#AA3377", fastani = "#FFCCCC",
                                revpyskani = "#228833", skani = "#CCDDAA")) + 
  scale_fill_manual(values = c(pyorthoani = "#4477AA", orthoani = "#BBCCEE",
                               revpyfastani = "#AA3377", fastani = "#FFCCCC",
                               revpyskani = "#228833", skani = "#CCDDAA")) +
  ggridges::geom_density_ridges(
    alpha = .7) + 
  theme_bw() + 
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  ylab("Number of CPUs") + 
  xlab("Real time (seconds)") 

time.ridges

# plot full data set
pdf(file = "trace_full_dataset_violin.pdf", width = 11, height = 5)
grid.arrange(mem.box, time.box, ncol=2)
dev.off()

pdf(file = "trace_full_dataset_ridges.pdf", width = 11, height = 8.5)
grid.arrange(mem.ridges, time.ridges, ncol=1)
dev.off()

# 1 CPU only

# get 1 CPU entries only
trace.1cpu <- droplevels(trace[which(trace$cpus=="1"),])
table(trace.1cpu$cpus)

# 1 CPU memory violin plot
mem.box.1cpu <-
  ggplot(data = trace.1cpu, mapping = aes(x = method, y = peak_rss_MB, color = method, fill = method)) + 
  geom_violin(alpha = 0.5) +
  theme_bw() +
  scale_color_manual(values = c(pyorthoani = "#4477AA", orthoani = "#BBCCEE",
                                revpyfastani = "#AA3377", fastani = "#FFCCCC",
                                revpyskani = "#228833", skani = "#CCDDAA")) + 
  scale_fill_manual(values = c(pyorthoani = "#4477AA", orthoani = "#BBCCEE",
                               revpyfastani = "#AA3377", fastani = "#FFCCCC",
                               revpyskani = "#228833", skani = "#CCDDAA")) +
  ylab("Peak RSS (MB)") + 
  xlab("Method") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))
mem.box.1cpu

# 1 CPU memory ridge plot
mem.ridges.1cpu <- 
  ggplot(data = trace.1cpu, mapping = aes(x = peak_rss_MB, y = fct_rev(method), color = method, fill = method)) + 
  coord_cartesian(clip = "off") +
  scale_y_discrete(expand = c(.07, .07)) +
  scale_color_manual(values = c(pyorthoani = "#4477AA", orthoani = "#BBCCEE",
                                revpyfastani = "#AA3377", fastani = "#FFCCCC",
                                revpyskani = "#228833", skani = "#CCDDAA")) + 
  scale_fill_manual(values = c(pyorthoani = "#4477AA", orthoani = "#BBCCEE",
                               revpyfastani = "#AA3377", fastani = "#FFCCCC",
                               revpyskani = "#228833", skani = "#CCDDAA")) +
  ggridges::geom_density_ridges(
    alpha = .7) + 
  theme_bw() + 
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10))
mem.ridges.1cpu

# 1 CPU time violin plot
time.box.1cpu <-
  ggplot(data = trace.1cpu, mapping = aes(x = method, y = realtime_seconds, color = method, fill = method)) + 
  geom_violin(alpha = 0.5) +
  theme_bw() +
  scale_color_manual(values = c(pyorthoani = "#4477AA", orthoani = "#BBCCEE",
                                revpyfastani = "#AA3377", fastani = "#FFCCCC",
                                revpyskani = "#228833", skani = "#CCDDAA")) + 
  scale_fill_manual(values = c(pyorthoani = "#4477AA", orthoani = "#BBCCEE",
                               revpyfastani = "#AA3377", fastani = "#FFCCCC",
                               revpyskani = "#228833", skani = "#CCDDAA")) +
  ylab("Real time (seconds)") + 
  xlab("Method") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))

time.box.1cpu

# 1 CPU time ridge plot
time.ridges.1cpu <- 
  ggplot(data = trace.1cpu, mapping = aes(x = realtime_seconds, y = fct_rev(method), color = method, fill = method)) + 
  coord_cartesian(clip = "off") +
  scale_y_discrete(expand = c(.07, .07)) +
  scale_color_manual(values = c(pyorthoani = "#4477AA", orthoani = "#BBCCEE",
                                revpyfastani = "#AA3377", fastani = "#FFCCCC",
                                revpyskani = "#228833", skani = "#CCDDAA")) + 
  scale_fill_manual(values = c(pyorthoani = "#4477AA", orthoani = "#BBCCEE",
                               revpyfastani = "#AA3377", fastani = "#FFCCCC",
                               revpyskani = "#228833", skani = "#CCDDAA")) +
  ggridges::geom_density_ridges(
    alpha = .7) + 
  theme_bw() + 
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  ylab("Method") + 
  xlab("Real time (seconds)") +
  theme(legend.position = "none") 

time.ridges.1cpu

# compare original method to our Python implementations
trace$compare <- ifelse(test = grepl(pattern = "orthoani", x = trace$method), yes = "OrthoANI", no = 
                               ifelse(test = grepl(pattern = "fastani", x = trace$method), yes = "FastANI", no = "skani"))
table(trace$compare)

# compare OrthoANI to PyOrthoANI
trace.orthoani <- droplevels(trace[which(trace$compare=="OrthoANI"),])
table(trace.orthoani$method)

# memory for OrthoANI (all CPUs)
orthoani.mem <-
  ggplot(data = trace.orthoani, mapping = aes(x = method, y = peak_rss_MB, color = method, fill = method)) + 
  geom_violin(alpha = 0.5) +
  theme_bw() +
  scale_color_manual(values = c(pyorthoani = "#4477AA", orthoani = "#BBCCEE")) + 
  scale_fill_manual(values = c(pyorthoani = "#4477AA", orthoani = "#BBCCEE")) + 
  ylab("Peak RSS (MB)") + 
  xlab("Method") + 
  theme(legend.position = "none") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))
orthoani.mem

# time for OrthoANI (all CPUs)
orthoani.time <-
  ggplot(data = trace.orthoani, mapping = aes(x = method, y = realtime_seconds, color = method, fill = method)) + 
  geom_violin(alpha = 0.5) +
  theme_bw() +
  scale_color_manual(values = c(pyorthoani = "#4477AA", orthoani = "#BBCCEE")) + 
  scale_fill_manual(values = c(pyorthoani = "#4477AA", orthoani = "#BBCCEE")) + 
  ylab("Real time (seconds)") + 
  xlab("Method") +
  theme(legend.position = "none") + 
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))
orthoani.time

table(trace.orthoani$method)
summary(trace.orthoani[which(trace.orthoani$method=="orthoani"&trace.orthoani$cpus=="1"), "realtime_seconds"])
summary(trace.orthoani[which(trace.orthoani$method=="pyorthoani"&trace.orthoani$cpus=="1"), "realtime_seconds"])
summary(trace.orthoani[which(trace.orthoani$method=="orthoani"&trace.orthoani$cpus=="8"), "realtime_seconds"])
summary(trace.orthoani[which(trace.orthoani$method=="pyorthoani"&trace.orthoani$cpus=="8"), "realtime_seconds"])
summary(trace.orthoani[which(trace.orthoani$method=="orthoani"&trace.orthoani$cpus=="16"), "realtime_seconds"])
summary(trace.orthoani[which(trace.orthoani$method=="pyorthoani"&trace.orthoani$cpus=="16"), "realtime_seconds"])


# compare FastANI to PyFastANI
trace.fastani <- droplevels(trace[which(trace$compare=="FastANI"),])
table(trace.fastani$method)

# memory for FastANI (all CPUs)
fastani.mem <-
  ggplot(data = trace.fastani, mapping = aes(x = method, y = peak_rss_MB, color = method, fill = method)) + 
  geom_violin(alpha = 0.5) +
  theme_bw() +
  scale_color_manual(values = c(revpyfastani = "#AA3377", fastani = "#FFCCCC")) + 
  scale_fill_manual(values = c(revpyfastani = "#AA3377", fastani = "#FFCCCC")) + 
  ylab("Peak RSS (MB)") + 
  xlab("Method") +
  theme(legend.position = "none") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))
fastani.mem

# time for FastANI (all CPUs)
fastani.time <-
  ggplot(data = trace.fastani, mapping = aes(x = method, y = realtime_seconds, color = method, fill = method)) + 
  geom_violin(alpha = 0.5) +
  theme_bw() +
  scale_color_manual(values = c(revpyfastani = "#AA3377", fastani = "#FFCCCC")) + 
  scale_fill_manual(values = c(revpyfastani = "#AA3377", fastani = "#FFCCCC")) + 
  ylab("Real time (seconds)") + 
  xlab("Method") +
  theme(legend.position = "none") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))
fastani.time

table(trace.fastani$method)
summary(trace.fastani[which(trace.fastani$method=="fastani"&trace.fastani$cpus=="1"), "realtime_seconds"])
summary(trace.fastani[which(trace.fastani$method=="revpyfastani"&trace.fastani$cpus=="1"), "realtime_seconds"])
summary(trace.fastani[which(trace.fastani$method=="fastani"&trace.fastani$cpus=="8"), "realtime_seconds"])
summary(trace.fastani[which(trace.fastani$method=="revpyfastani"&trace.fastani$cpus=="8"), "realtime_seconds"])
summary(trace.fastani[which(trace.fastani$method=="fastani"&trace.fastani$cpus=="16"), "realtime_seconds"])
summary(trace.fastani[which(trace.fastani$method=="revpyfastani"&trace.fastani$cpus=="16"), "realtime_seconds"])

# compare skani to Pyskani
trace.skani <- droplevels(trace[which(trace$compare=="skani"),])
table(trace.skani$method)

# memory for skani (all CPUs)
skani.mem <-
  ggplot(data = trace.skani, mapping = aes(x = method, y = peak_rss_MB, color = method, fill = method)) + 
  geom_violin(alpha = 0.5) +
  theme_bw() +
  scale_color_manual(values = c(revpyskani = "#228833", skani = "#CCDDAA")) + 
  scale_fill_manual(values = c(revpyskani = "#228833", skani = "#CCDDAA")) + 
  ylab("Peak RSS (MB)") + 
  xlab("Method") +
  theme(legend.position = "none") + 
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))
skani.mem

# time for skani (all CPUs)
skani.time <-
  ggplot(data = trace.skani, mapping = aes(x = method, y = realtime_seconds, color = method, fill = method)) + 
  geom_violin(alpha = 0.5) +
  theme_bw() +
  scale_color_manual(values = c(revpyskani = "#228833", skani = "#CCDDAA")) + 
  scale_fill_manual(values = c(revpyskani = "#228833", skani = "#CCDDAA")) + 
  ylab("Real time (seconds)") + 
  xlab("Method") +
  theme(legend.position = "none") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))
skani.time

table(trace.skani$method)
summary(trace.skani[which(trace.skani$method=="skani"&trace.skani$cpus=="1"), "realtime_seconds"])
summary(trace.skani[which(trace.skani$method=="revpyskani"&trace.skani$cpus=="1"), "realtime_seconds"])

# plot main figures
pdf(file = "main_figure_violins.pdf", width = 7, height = 5)
grid.arrange(plot.orthoani, plot.fastani, plot.skani,
             orthoani.time, fastani.time, skani.time, ncol=3)
dev.off()
pdf(file = "mem_and_time_violins.pdf", width = 7, height = 5)
grid.arrange(orthoani.mem, fastani.mem, skani.mem,
             orthoani.time, fastani.time, skani.time, ncol=3)
dev.off()

pdf(file = "main_figure_ridge.pdf", width = 7, height = 5)
grid.arrange(
  grobs = list(plot.orthoani, plot.fastani, plot.skani, time.ridges.1cpu),
  widths = c(1, 1, 1),
  layout_matrix = rbind(c(1, 2, 3),
                        c(4))
)
dev.off()


