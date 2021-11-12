ggplot(wplot, aes(Time), ylim(-300:200)) + 
  geom_bar(data = subset(wplot, VARIABLE == "count.up"), 
           aes(y = VALUE, fill = TREATMENT), stat = "identity", position = "dodge") +
  geom_bar(data = subset(wplot, VARIABLE == "count.down"), 
           aes(y = VALUE, fill = TREATMENT), stat = "identity", position = "dodge") + 
  geom_hline(yintercept = 0,colour = "grey90")


last_plot() + 
  geom_text(data = subset(wplot, VARIABLE == "count.up"), 
            aes(Time, VALUE, group=TREATMENT, label=VALUE),
            position = position_dodge(width=0.9), vjust = -1,size=4) +
  geom_text(data = subset(wplot, VARIABLE == "count.down"), 
            aes(Time, VALUE, group=TREATMENT, label=VALUE),
            position = position_dodge(width=0.9), vjust=1,  size=4) +
  coord_cartesian(ylim = c(-500, 500))


ggplot()+geom_bar(data = dat1, aes(x=TREATMENT, y=VALUE, fill=TREATMENT, label=VALUE),stat = "Identity") +
  geom_bar(data = dat2, aes(x=TREATMENT, y=VALUE, fill=TREATMENT, label=VALUE),stat = "Identity") +
  scale_fill_brewer(type = "seq", palette = 1) +
  facet_wrap(~Time)