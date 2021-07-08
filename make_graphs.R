# making graphs
library(tidyverse)
library(gRandma)
library(EFGLmh)
library(cowplot)

load("after_fn_unrel_lower.rda")

# graphs

# error rates with panel sizes
snp_err <- tibble()
mh_err <- tibble()
mh_IS <- tibble() 
snp_IS <- tibble()
for(i in 1:length(fn_snp)){
	temp <- fn_snp[[i]][[1]] %>% mutate(pan = paste0("snp_", i*100)) %>% 
		left_join(fp_unrel_snp[[i]][[1]])
	snp_err <- snp_err %>% bind_rows(temp)
	
	temp <- fn_mh[[i]][[1]] %>% mutate(pan = paste0("mh_", i*100)) %>% 
		left_join(fp_unrel_mh[[i]][[1]])
	mh_err <- mh_err %>% bind_rows(temp)
	
	temp <- fp_unrel_mh_IS[[i]][[1]] %>% mutate(pan = paste0("mh_", i*100))
	mh_IS <- mh_IS %>% bind_rows(temp)

	temp <- fp_unrel_snp_IS[[i]][[1]] %>% mutate(pan = paste0("snp_", i*100))
	snp_IS <- snp_IS %>% bind_rows(temp)
}

snp_err %>%
	ggplot(aes(x = falseNeg, y = falsePosUnrel, color = pan)) + geom_point() + geom_line() +
	scale_y_continuous(trans="log10") + coord_cartesian(xlim = c(0, .05))
mh_err %>%
	ggplot(aes(x = falseNeg, y = falsePosUnrel, color = pan)) + geom_point() + geom_line() +
	scale_y_continuous(trans="log10") + coord_cartesian(xlim = c(0, .05))
mh_IS %>%
	ggplot(aes(x = falseNeg, y = falsePosUnrel, color = pan)) + geom_point() + geom_line() +
	scale_y_continuous(trans="log10") + coord_cartesian(xlim = c(0, .05))
snp_IS %>%
	ggplot(aes(x = falseNeg, y = falsePosUnrel, color = pan)) + geom_point() + geom_line() +
	scale_y_continuous(trans="log10") + coord_cartesian(xlim = c(0, .05))

# add larger snp panels
load("extra_snp_unrel.rda")
for(i in 1:length(fn_snp)){
	temp <- fn_snp[[i]][[1]] %>% mutate(pan = paste0("snp_", (4 + i)*100)) %>% 
		left_join(fp_unrel_snp[[i]][[1]])
	snp_err <- snp_err %>% bind_rows(temp)
	
	temp <- fp_unrel_snp_IS[[i]][[1]] %>% mutate(pan = paste0("snp_", (4 + i)*100))
	snp_IS <- snp_IS %>% bind_rows(temp)
}

snp_err %>%
	ggplot(aes(x = falseNeg, y = falsePosUnrel, color = pan)) + geom_point() + geom_line() +
	scale_y_continuous(trans="log10") + coord_cartesian(xlim = c(0, .05))
snp_IS %>%
	ggplot(aes(x = falseNeg, y = falsePosUnrel, color = pan)) + geom_point() + geom_line() +
	scale_y_continuous(trans="log10") + coord_cartesian(xlim = c(0, .05))



# IS vs stratified
mh_err %>% filter(llrThreshold >= 0, falseNeg > 0, falseNeg < .05) %>%
	select(llrThreshold, falsePosUnrel, pan) %>%
	rename(strat = falsePosUnrel) %>%
	inner_join(mh_IS %>% select(llrThreshold, falsePosUnrel, pan) %>% rename(is = falsePosUnrel),
				  by = c("llrThreshold", "pan")) %>%
	ggplot(aes(x = strat, y = is, color = pan)) + geom_point() + 
	scale_y_continuous(trans="log10") + scale_x_continuous(trans="log10") + 
	geom_abline(slope = 1, intercept = 0)


# different relationships

load("diffRel.rda")
relsToTest <- c("True_GAunt", "True_Unrel", "True_HGAunt", "True_GpCous", 
					 "GAunt_Unrel", "HGAunt_Unrel", "GpCous_Unrel", "GAunt", "GAunt_HGAunt", 
					 "Gaunt_GpCous", "HGAunt", "HGAunt_GpCous", "GpCous")
relsToTest <- c(relsToTest, relsToTest)
method <- c(rep("IS", length(relsToTest)/2), rep("strat", length(relsToTest)/2))
dRelRes <- tibble()
for(r in 1:length(fp_mh)){
	dRelRes <- bind_rows(dRelRes,
		fp_mh[[r]][[1]] %>% as_tibble %>% mutate(method = method[r], rel = relsToTest[r]) %>% 
		rename(falsePos = 3) %>% select(llrThreshold, falsePos, rel, method))
}

# add false negative
dRelRes <- dRelRes %>% left_join(mh_IS %>% filter(pan == "mh_300") %>% select(llrThreshold, falseNeg) 
							 %>% distinct, by = "llrThreshold")

dRelRes %>% filter(llrThreshold >= 0, falseNeg > 0, falseNeg < .05,
						 method == "strat", rel %in% c("True_Unrel", "True_GpCous")) %>% 
	group_by(rel) %>% filter(llrThreshold == max(llrThreshold))

fp_mh[[15]][[2]] %>% filter(llrThreshold == 20)
fp_mh[[17]][[2]] %>% filter(llrThreshold == 20)

dRelRes %>% filter(method == "IS", rel %in% c("True_GAunt", "True_Unrel", 
															 "GAunt_Unrel", "HGAunt_Unrel")) %>% 
	ggplot(aes(x = falseNeg, y = falsePos, color = rel)) + geom_point() + geom_line() +
	scale_y_continuous(trans="log10") + coord_cartesian(xlim = c(0, .05), ylim = c(1e-10, NA))

dRelRes %>% filter(method == "strat", rel %in% c("True_GAunt", "True_Unrel", 
															 "GAunt_Unrel", "HGAunt_Unrel")) %>% 
	ggplot(aes(x = falseNeg, y = falsePos, color = rel)) + geom_point() + geom_line() +
	scale_y_continuous(trans="log10") + coord_cartesian(xlim = c(0, .05), ylim = c(1e-10, NA))

# diff rel IS vs stratified
dRelRes %>% filter(llrThreshold >= 0, falseNeg > 0, falseNeg < .05) %>%
	select(llrThreshold, falsePos, rel, method) %>%
	spread(method, falsePos) %>%
	ggplot(aes(x = strat, y = IS, color = rel)) + geom_point() + 
	scale_y_continuous(trans="log10") + scale_x_continuous(trans="log10") + 
	geom_abline(slope = 1, intercept = 0)



# missing vs no missing data

load("noMissingData.rda")

temp <- bind_rows(mh_IS %>% filter(pan == "mh_300") %>% mutate(pan = "miss_3"),
			 fp_unrel_mh_IS_noMiss[[1]] %>% as_tibble %>% mutate(pan = "miss_0")) %>%
	filter(falseNeg < .055) %>% select(llrThreshold, falsePosUnrel, falseNeg, pan) %>%
	group_by(pan) %>% arrange(pan, desc(llrThreshold)) %>% mutate(num = order(llrThreshold, decreasing = TRUE)) %>%
	filter(num < 3)
approxfun(temp %>% filter(pan == "miss_0") %>% pull(falseNeg), temp %>% filter(pan == "miss_0") %>% pull(falsePosUnrel))(.05)
# 2.064542e-12
approxfun(temp %>% filter(pan == "miss_3") %>% pull(falseNeg), temp %>% filter(pan == "miss_3") %>% pull(falsePosUnrel))(.05)
# 1.320013e-11

bind_rows(mh_IS %>% filter(pan == "mh_300") %>% mutate(pan = "miss_3"),
	fp_unrel_mh_IS_noMiss[[1]] %>% as_tibble %>% mutate(pan = "miss_0")) %>%
	ggplot(aes(x = falseNeg, y = falsePosUnrel, color = pan)) + geom_point() + geom_line() +
	scale_y_continuous(trans="log10") + coord_cartesian(xlim = c(0, .05), ylim = c(1e-13, NA))


# graphs for ms

# comparing panels IS
snpISgraph <- snp_IS %>% mutate(pan = gsub("snp_", "", pan)) %>% filter(pan != "1000") %>%
	filter(as.numeric(substr(pan, 1,1)) %% 2 == 1) %>%
	ggplot(aes(x = falseNeg, y = falsePosUnrel, color = pan)) + geom_point() + geom_line() +
	scale_y_continuous(trans="log10") + coord_cartesian(xlim = c(0, .05), ylim = c(1e-17, 1e-1)) +
	labs(x = "False negative", y = "False positive", color = "Number of loci")

mhISgraph <- mh_IS %>% mutate(pan = gsub("mh_", "", pan)) %>%
	ggplot(aes(x = falseNeg, y = falsePosUnrel, color = pan)) + geom_point() + geom_line() +
	scale_y_continuous(trans="log10") + coord_cartesian(xlim = c(0, .05), ylim = c(1e-17, 1e-1)) +
	labs(x = "False negative", y = "False positive", color = "Number of loci")

snp_mh_IS_graph <- plot_grid(snpISgraph, mhISgraph, align = "hv", labels = "AUTO",
									  ncol = 1)

save_plot(filename = "combined_IS_graph.pdf", plot = snp_mh_IS_graph, base_height = 7, 
			 base_width = 6.5)

# for STHD manager meeting 2021 presentation
# ggsave("mh_IS_presentation.pdf", plot = mhISgraph, height = 5, width = 8)



# IS vs stratified
is_v_strat <- mh_err %>% filter(llrThreshold >= 0, falseNeg > 0, falseNeg < .05) %>%
	select(llrThreshold, falsePosUnrel, pan) %>%
	rename(strat = falsePosUnrel) %>%
	inner_join(mh_IS %>% select(llrThreshold, falsePosUnrel, pan) %>% rename(is = falsePosUnrel),
				  by = c("llrThreshold", "pan")) %>% mutate(pan = gsub("mh_", "", pan)) %>%
	ggplot(aes(x = strat, y = is, color = pan)) + geom_point() +
	scale_y_continuous(trans="log10", breaks = 10^(seq(-16, -2, 3))) + 
	scale_x_continuous(trans="log10", breaks = 10^(seq(-16, -2, 3)), limits = c(1e-14, NA)) + 
	geom_abline(slope = 1, intercept = 0) + 
	labs(x = "Stratified sampling", y = "Importance sampling", color = "Number of loci")
	
ggsave2("comp_methods.pdf", plot = is_v_strat, width = 5.5, height = 4)


# missing vs no missing
# Note this is with mh_300 panel
miss_vs_noMiss <-  bind_rows(mh_IS %>% filter(pan == "mh_300") %>% mutate(pan = "Missing"),
			 fp_unrel_mh_IS_noMiss[[1]] %>% as_tibble %>% mutate(pan = "No missing")) %>%
	ggplot(aes(x = falseNeg, y = falsePosUnrel, color = pan)) + geom_point() + geom_line() +
	scale_y_continuous(trans="log10") + coord_cartesian(xlim = c(0, .05), ylim = c(1e-12, NA)) +
	labs(x = "False negative", y = "False positive", color = "")

ggsave2("comp_miss.pdf", plot = miss_vs_noMiss, width = 5.5, height = 4)

# diff relationships
# Note these are with mh_300 panel
dRelRes <- dRelRes %>% mutate(rel = gsub("_", ", ", rel))

# for supplementary material
diff_rel <- dRelRes %>% filter(method == "strat", rel %in% c("True, GAunt", "True, Unrel", 
																 "GAunt, Unrel", "HGAunt, Unrel")) %>% 
	mutate(rel = factor(rel, levels = c("True, GAunt", "True, Unrel", "GAunt, Unrel", "HGAunt, Unrel"))) %>%
	ggplot(aes(x = falseNeg, y = falsePos, color = rel)) + geom_point() + geom_line() +
	scale_y_continuous(trans="log10") + coord_cartesian(xlim = c(0, .05), ylim = c(1e-10, NA)) +
	labs(x = "False negative", y = "False positive", color = "Relationship")
ggsave2("diff_rel_sup.pdf", plot = diff_rel, width = 5.5, height = 4)

# diff rel IS vs strat

inset <- dRelRes %>% filter(llrThreshold >= 0, falseNeg > 0, falseNeg < .05, llrThreshold %% 2 == 0) %>%
	select(llrThreshold, falsePos, rel, method) %>%
	spread(method, falsePos) %>%
	ggplot(aes(x = strat, y = IS, color = rel, shape = rel)) + geom_point() + 
	scale_y_continuous(trans="log10", limits = c(1e-2, NA)) + 
	scale_x_continuous(trans="log10", limits = c(.01, NA)) + 
	geom_abline(slope = 1, intercept = 0) + 
	labs(x = "Stratified sampling", y = "Importance", color = "Relationship", shape = "Relationship") +
	scale_shape_manual(values=c(1:13)) + theme(legend.position="none", axis.title = element_blank())


drel_is_strat <- dRelRes %>% filter(llrThreshold >= 0, falseNeg > 0, falseNeg < .05, llrThreshold %% 2 == 0) %>%
	select(llrThreshold, falsePos, rel, method) %>%
	spread(method, falsePos) %>%
	ggplot(aes(x = strat, y = IS, color = rel, shape = rel)) + geom_point() + 
	scale_y_continuous(trans="log10", breaks = 10^(seq(-16, 0, 2))) + 
	scale_x_continuous(trans="log10", breaks = 10^(seq(-16, 0, 2)), limits = c(NA, NA)) + 
	geom_abline(slope = 1, intercept = 0) + 
	labs(x = "Stratified sampling", y = "Importance sampling", color = "Relationship", shape = "Relationship") +
	scale_shape_manual(values=c(1:13)) +
	annotation_custom(ggplotGrob(inset), xmin = -10, xmax = -5, 
							  ymin = -4, ymax = 0)

ggsave2("comp_methods_drel.pdf", plot = drel_is_strat, width = 6.5, height = 5)
