library(ggplot2)
library(stringr)
library(ggthemes)

## Figure 1
#### Plot of D_{ij} over time
plot_data <- data.frame(day = seq(from = -2, to = 56, by = 1))
alpha = 0.5
plot_data$d <- 0
plot_data$d[plot_data$day == 0] <- 10
plot_data$d[plot_data$day == 28] <- 15
plot_data$D[1] <- plot_data$d[1]
for(i in 2:nrow(plot_data)){plot_data$D[i] <- plot_data$d[i] + plot_data$D[i-1]*alpha^((plot_data$day[i]-plot_data$day[i-1])/28)}

plot_D <- ggplot(plot_data, mapping = aes(x = day))+
  geom_col(mapping = aes(y=d, fill = "Dose administered (d)"), width = 2)+
  geom_line(mapping = aes(y = D, color = "Effective dose (D)"))+
  labs(x = "Days since first dosing", y = "Agent amount (mg)", color = "", fill = "")+
  ylim(c(0,22))+
  scale_x_continuous( breaks = 14*(0:4))+
  scale_fill_discrete(type = "lightblue")+
  scale_color_discrete(type = "black")+
  theme_bw()

ggsave(plot_D, filename = "falke_1_t.png", device = "png", dpi = 600, width = 9.5, height = 4.8, scale = 0.75, units = "in", path = "output/Figures/")
ggsave(plot_D, filename = "falke_1_t.tiff", device = "tiff", dpi = 600, width = 9.5, height = 4.8, scale = 0.75, units = "in", path = "output/Figures/")

## Figure 2
#### Plot of Dose Escalation
plot_data <- data.frame(
           pid = c(1,2,3,3,4,5, rep("mtd",5)),
           d1 = c(1, 1, 1, 2, 2, 2, 3, 3, 1, 1, 2),
           d2 = c(1, 1, 1, 2, 2, 2, 3, 1, 1, 2, 2),
           t1 = c(0, 10, 20, 48, 50, 60, 28, 65, 65, 76, 76), 
           t2 = c(28, 38, 48, 76, 65, 88, 65, 65, 76, 76, 90),
           status1 = c("enter", "enter", "enter", "escalate", "enter", "enter", "", "", "", "",""),
           status2 = c("done", "done", "escalate", "done", "dlt", "done", "", "", "", "", ""))

plot_data$d1 <- factor(plot_data$d1, levels = 1:3)
plot_data$d2 <- factor(plot_data$d2, levels = 1:3)
plot_data$pid <- factor(plot_data$pid, levels = c(5:1,"mtd"), labels = c(paste("Participant",5:1), "Estimated MTD"))
plot_data$status1 <- factor(plot_data$status1, 
                            levels = c("enter", "escalate", "dlt", "done"),
                            labels = c("Enter study", "Dose escalation", "Dose limiting toxicity", "Exited study"))
plot_data$status2 <- factor(plot_data$status2, 
                            levels = c("enter", "escalate", "dlt", "done"),
                            labels = c("Enter study", "Dose escalation", "Dose limiting toxicity", "Exited study"))

plot_escalate <- ggplot(plot_data, mapping = aes(y = d1))+
  geom_segment(mapping = aes(x = t1, xend = t2, y = d1, yend = d2), linewidth = 0.75)+
  geom_point(mapping = aes(x = t1, y = d1, shape = status1), size = 3)+
  geom_point(mapping = aes(x = t2, y = d2, shape = status2), size = 3)+
  scale_shape_manual(breaks = c("Enter study", "Dose escalation", "Dose limiting toxicity", "Exited study"), 
                     values = c(15, 1, 8, 20))+
  scale_y_discrete(drop = FALSE)+
  labs(shape = "Event", y = "Dose Level", x = "Time (days)")+
  facet_wrap(facets = vars(pid), ncol = 1)+
  theme_bw()+ 
  theme(panel.spacing = unit(0, "lines")) 

ggsave(plot_escalate, filename = "falke_2_t.png", device = "png", dpi = 600,  width = 12, height = 9, scale = 0.75, units = "in", path = "output/Figures/")
ggsave(plot_escalate, filename = "falke_2_t.tiff", device = "tiff", dpi = 600,  width = 12, height = 9, scale = 0.75, units = "in", path = "output/Figures/")

## Figure 3
#### Plot of DRC scenarios
logistic <- function(x, x0, k, phi = 0.3){signif(1/(1+exp(-(log(phi/(1-phi)) + k*(x-x0)))))}

drc_data <- data.frame(k = rep(1:5, 6), 
                       p = c(logistic(1:5, x0 = 0, k = 1/2),
                             logistic(1:5, x0 = 1, k = 1/2),
                             logistic(1:5, x0 = 2, k = 1/2),
                             logistic(1:5, x0 = 3, k = 1/2),
                             logistic(1:5, x0 = 4, k = 1/2),
                             logistic(1:5, x0 = 5, k = 1/2)),
                       scenario = factor(rep(0:5, each = 5)))

plot_drc <- ggplot(drc_data, mapping = aes(x = k, y = p, color = scenario, group = scenario))+
  geom_hline(aes(yintercept = 0.3),  linetype='dashed', show.legend = TRUE)+
  geom_line()+
  geom_point()+
  labs(x = "Dose Level", y = "Probability of DLT", color = "Scenario")+
  scale_y_continuous(breaks = c(0.1,0.3,0.5,0.7), limits = c(0,0.9))+
  scale_fill_discrete(name = "Target DLT Rate")+
  theme_bw()

ggsave(plot_drc, filename = "falke_3_t.png", device = "png", dpi = 600, width = 9.5, height = 4.8, scale = 0.75, units = "in", path = "output/Figures/")
ggsave(plot_drc, filename = "falke_3_t.tiff", device = "tiff", dpi = 600, width = 9.5, height = 4.8, scale = 0.75, units = "in", path = "output/Figures/")
