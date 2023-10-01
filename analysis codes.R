####################################################################################################
# Example of code written by Jia Wei (Nuffield Department of Medicine, University of Oxford)       #
# for analyses of SARS-CoV-2 reinfection in the UK general population.                             #
# Accompanying paper/preprint: Risk of SARS-CoV-2 reinfection during multiple Omicron variant      #
# waves in the UK general population.                                                              #
####################################################################################################

library(tidyverse)
library(ggsankey)
library(MASS)
library(nnet)
library(gmodels)
library(ggeffects)
library(cowplot)
library(emmeans)
library(survival)
library(survminer)
library(rstpm2)


####################################################################################################
#########################################################################sankey plots of reinfection

data_alpha=data_episode %>% 
  filter(variant=="Alpha"&episode==1)

ID=data_alpha$participant_id

data_alpha=data_episode %>% 
  filter(participant_id%in%ID)

a=data_alpha %>% 
  select(participant_id,episode,variant) %>% 
  pivot_wider(names_from = episode,values_from = c(3),names_prefix = "episode") %>% 
  select(-participant_id) %>% 
  select(episode1,episode2,episode3,episode4)

sum=a %>% group_by(episode1,episode2,episode3,episode4) %>% summarise(n=n()) 

a2=a %>% 
  make_long(episode1,episode2,episode3,episode4)

a3=a2 %>% 
  filter(!is.na(node)) %>% 
  mutate(node=factor(node,levels=c("Omicron BQ.1/CH.1.1/XBB.1.5","Omicron BA.4/5","Omicron BA.2",
                                   "Omicron BA.1","Delta","Alpha","Pre-Alpha")),
         next_node=factor(next_node,levels=c("Omicron BQ.1/CH.1.1/XBB.1.5","Omicron BA.4/5","Omicron BA.2",
                                             "Omicron BA.1","Delta","Alpha","Pre-Alpha")))

total=a3 %>% 
  group_by(x) %>% 
  summarise(TotalCount=n())

perc=a3 %>% 
  group_by(node,x) %>% 
  tally() %>% 
  left_join(total) %>% 
  mutate(pct=n/TotalCount)

a3=a3 %>% 
  left_join(perc)

cbp=c("#FFDC91FF","#6F99ADFF","#7876B1FF","#20854EFF","#E18727FF","#0072B5FF","#BC3C29FF")
ggplot(a3,aes(x=x,
                 next_x=next_x,
                 node=node,
                 next_node=next_node,
                 fill=factor(node),
                 label=paste0(node," (",round(pct*100,1),"%)")))+
  geom_sankey(flow.alpha=0.5,node.color="black",show.legend = F)+
  geom_sankey_label(size=3,color="black",fill="white",hjust=0)+
  theme_light()+
  theme(legend.position="none",
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank())+
  scale_fill_manual(values = cbp)+
  scale_x_discrete(labels=c("First infection","Second infection","Third infection","Fourth infection"))



####################################################################################################
#############################################################################linear regression of Ct 

fit_ct <- rlm(ct ~age + sex + ethnicity + hcw + lthc+ IMD + symptom + count +
               variant + time_vac + episode + pre_2 + alpha_2 + delta_2 + ba2_2 + ba45_2 + bq_2 + delta_3 + ba2_3 + ba45_3 + bq_3 + 
               pre_no + alpha_no + delta_no + ba2_no + ba45_no + bq_no + 
               alpha_14 + delta_14 + ba2_14 + ba45_14 + bq_14 + delta_90 + ba2_90 + ba45_90 + bq_90 + time_vac:episode,
             data=data_reg)
summary(fit_ct)
sum=data.frame(summary(fit_ct)$coefficients)

CI=confint.default(fit_ct,level=0.95)
sum=cbind(sum, CI) %>% rownames_to_column()

sum$p=2*pt(abs(sum$t.value), summary(fit_ct)$df[2], lower.tail = F)


newdata=data.frame(age=0,sex="Female",ethnicity="White",hcw="No",lthc="No",IMD=0,symptom="Yes",count=1,
                   variant="Pre-Alpha",time_vac="No vaccination",episode="1",
                   pre_2=0,alpha_2=0,delta_2=0,ba2_2=0,ba45_2=0,bq_2=0,delta_3=0,ba2_3=0,ba45_3=0,bq_3=0,
                   pre_no=1,alpha_no=0,delta_no=0,ba2_no=0,ba45_no=0,bq_no=0,alpha_14=0,delta_14=0,ba2_14=0,ba45_14=0,bq_14=0,delta_90=0,ba2_90=0,ba45_90=0,bq_90=0) %>% 
  add_row(variant="Pre-Alpha",time_vac="No vaccination",episode="2",pre_2=1,pre_no=1) %>% 
  add_row(variant="Alpha",time_vac="No vaccination",episode="1",alpha_no=1) %>% 
  add_row(variant="Alpha",time_vac="No vaccination",episode="2",alpha_2=1,alpha_no=1) %>% 
  add_row(variant="Alpha",time_vac="14-90",episode="1",alpha_14=1) %>% 
  add_row(variant="Alpha",time_vac="14-90",episode="2",alpha_2=1,alpha_14=1) %>% 
  
  add_row(variant="Delta",time_vac="No vaccination",episode="1",delta_no=1) %>% 
  add_row(variant="Delta",time_vac="No vaccination",episode="2",delta_2=1,delta_no=1) %>% 
  add_row(variant="Delta",time_vac="14-90",episode="1",delta_14=1) %>% 
  add_row(variant="Delta",time_vac="14-90",episode="2",delta_2=1,delta_14=1) %>% 
  add_row(variant="Delta",time_vac="90-180",episode="1",delta_90=1) %>% 
  add_row(variant="Delta",time_vac="90-180",episode="2",delta_2=1,delta_90=1) %>% 
  add_row(variant="Delta",time_vac=">180",episode="1") %>% 
  add_row(variant="Delta",time_vac=">180",episode="2",delta_2=1) %>% 
  
  add_row(variant="Omicron BA.1",time_vac="No vaccination",episode="1") %>% 
  add_row(variant="Omicron BA.1",time_vac="No vaccination",episode="2") %>% 
  add_row(variant="Omicron BA.1",time_vac="No vaccination",episode="3") %>% 
  add_row(variant="Omicron BA.1",time_vac="14-90",episode="1") %>% 
  add_row(variant="Omicron BA.1",time_vac="14-90",episode="2") %>% 
  add_row(variant="Omicron BA.1",time_vac="14-90",episode="3") %>% 
  add_row(variant="Omicron BA.1",time_vac="90-180",episode="1") %>% 
  add_row(variant="Omicron BA.1",time_vac="90-180",episode="2") %>% 
  add_row(variant="Omicron BA.1",time_vac="90-180",episode="3") %>% 
  add_row(variant="Omicron BA.1",time_vac=">180",episode="1") %>% 
  add_row(variant="Omicron BA.1",time_vac=">180",episode="2") %>% 
  add_row(variant="Omicron BA.1",time_vac=">180",episode="3") %>% 
  
  add_row(variant="Omicron BA.2",time_vac="No vaccination",episode="1",ba2_no=1) %>% 
  add_row(variant="Omicron BA.2",time_vac="No vaccination",episode="2",ba2_2=1,ba2_no=1) %>% 
  add_row(variant="Omicron BA.2",time_vac="No vaccination",episode="3",ba2_3=1,ba2_no=1) %>% 
  add_row(variant="Omicron BA.2",time_vac="14-90",episode="1",ba2_14=1) %>% 
  add_row(variant="Omicron BA.2",time_vac="14-90",episode="2",ba2_2=1,ba2_14=1) %>% 
  add_row(variant="Omicron BA.2",time_vac="14-90",episode="3",ba2_3=1,ba2_14=1) %>% 
  add_row(variant="Omicron BA.2",time_vac="90-180",episode="1",ba2_90=1) %>% 
  add_row(variant="Omicron BA.2",time_vac="90-180",episode="2",ba2_2=1,ba2_90=1) %>% 
  add_row(variant="Omicron BA.2",time_vac="90-180",episode="3",ba2_3=1,ba2_90=1) %>% 
  add_row(variant="Omicron BA.2",time_vac=">180",episode="1") %>% 
  add_row(variant="Omicron BA.2",time_vac=">180",episode="2",ba2_2=1) %>% 
  add_row(variant="Omicron BA.2",time_vac=">180",episode="3",ba2_3=1) %>% 
  
  add_row(variant="Omicron BA.4/5",time_vac="No vaccination",episode="1",ba45_no=1) %>% 
  add_row(variant="Omicron BA.4/5",time_vac="No vaccination",episode="2",ba45_2=1,ba45_no=1) %>% 
  add_row(variant="Omicron BA.4/5",time_vac="No vaccination",episode="3",ba45_3=1,ba45_no=1) %>% 
  add_row(variant="Omicron BA.4/5",time_vac="14-90",episode="1",ba45_14=1) %>% 
  add_row(variant="Omicron BA.4/5",time_vac="14-90",episode="2",ba45_2=1,ba45_14=1) %>% 
  add_row(variant="Omicron BA.4/5",time_vac="14-90",episode="3",ba45_3=1,ba45_14=1) %>% 
  add_row(variant="Omicron BA.4/5",time_vac="90-180",episode="1",ba45_90=1) %>% 
  add_row(variant="Omicron BA.4/5",time_vac="90-180",episode="2",ba45_2=1,ba45_90=1) %>% 
  add_row(variant="Omicron BA.4/5",time_vac="90-180",episode="3",ba45_3=1,ba45_90=1) %>% 
  add_row(variant="Omicron BA.4/5",time_vac=">180",episode="1") %>% 
  add_row(variant="Omicron BA.4/5",time_vac=">180",episode="2",ba45_2=1) %>% 
  add_row(variant="Omicron BA.4/5",time_vac=">180",episode="3",ba45_3=1) %>% 
  
  add_row(variant="Omicron BQ.1/CH.1.1/XBB.1.5",time_vac="No vaccination",episode="1",bq_no=1) %>% 
  add_row(variant="Omicron BQ.1/CH.1.1/XBB.1.5",time_vac="No vaccination",episode="2",bq_2=1,bq_no=1) %>% 
  add_row(variant="Omicron BQ.1/CH.1.1/XBB.1.5",time_vac="No vaccination",episode="3",bq_3=1,bq_no=1) %>% 
  add_row(variant="Omicron BQ.1/CH.1.1/XBB.1.5",time_vac="14-90",episode="1",bq_14=1) %>% 
  add_row(variant="Omicron BQ.1/CH.1.1/XBB.1.5",time_vac="14-90",episode="2",bq_2=1,bq_14=1) %>% 
  add_row(variant="Omicron BQ.1/CH.1.1/XBB.1.5",time_vac="14-90",episode="3",bq_3=1,bq_14=1) %>% 
  add_row(variant="Omicron BQ.1/CH.1.1/XBB.1.5",time_vac="90-180",episode="1",bq_90=1) %>% 
  add_row(variant="Omicron BQ.1/CH.1.1/XBB.1.5",time_vac="90-180",episode="2",bq_2=1,bq_90=1) %>% 
  add_row(variant="Omicron BQ.1/CH.1.1/XBB.1.5",time_vac="90-180",episode="3",bq_3=1,bq_90=1) %>% 
  add_row(variant="Omicron BQ.1/CH.1.1/XBB.1.5",time_vac=">180",episode="1") %>% 
  add_row(variant="Omicron BQ.1/CH.1.1/XBB.1.5",time_vac=">180",episode="2",bq_2=1) %>% 
  add_row(variant="Omicron BQ.1/CH.1.1/XBB.1.5",time_vac=">180",episode="3",bq_3=1) 

newdata$age=0
newdata$sex="Female"
newdata$ethnicity="White"
newdata$hcw="No"
newdata$lthc="No"
newdata$IMD=0
newdata$symptom="Yes"
newdata$count=1

newdata[is.na(newdata)]<-0

pred=predict(fit_ct,newdata,se.fit=T)

newdata$pred=pred$fit
newdata$lci=pred$fit-1.96*pred$se.fit
newdata$uci=pred$fit+1.96*pred$se.fit

newdata=newdata %>% 
  mutate(variant=factor(variant,levels=c("Pre-Alpha","Alpha","Delta","Omicron BA.1","Omicron BA.2","Omicron BA.4/5","Omicron BQ.1/CH.1.1/XBB.1.5")),
         time_vac=factor(time_vac, levels=c("No vaccination","14-90","90-180",">180")))

ggplot(newdata,aes(time_vac,pred,group=episode),size=1.2)+
  geom_point(aes(color=episode),position=position_dodge(width=0.5))+
  geom_errorbar(aes(ymin=lci,ymax=uci,color=episode),position=position_dodge(width=0.5),width=0.5)+
  labs(title="", x="Days from the most recent vaccination",y="Predicted mean Ct values (95% CI)")+
  facet_wrap(~variant)+
  coord_cartesian(ylim=c(15,30))+
  scale_y_continuous(limits=c(0,50), breaks = seq(0,35,5))+
  theme_light()+
  theme(legend.position = "bottom",axis.text.x = element_text())+
  scale_fill_nejm(name="Episode")+
  scale_color_nejm(name="Episode")



### model to include previous symptom and previous ct

fit_ct <- rlm(ct~age + sex + ethnicity + hcw + lthc + IMD + symptom + count + sympt_previous + ns(ct_previous,3) + 
                variant + time_vac + episode + 
                pre_no + alpha_no + delta_no + ba2_no + ba45_no + bq_no + 
                delta_14 + ba2_14 + ba45_14 + bq_14 + delta_90 + ba2_90 + ba45_90 + bq_90,
              data=data_reg3)
summary(fit_ct)
sum=data.frame(summary(fit_ct)$coefficients)

CI=confint.default(fit_ct,level=0.95)
sum=cbind(sum, CI) %>% rownames_to_column()

sum$p=2*pt(abs(sum$t.value), summary(fit_ct)$df[2], lower.tail = F)


p <- ggpredict(fit_ct, c("ct_previous [10,15,20,25,30,35]","sympt_previous"),
               condition = c(symptom="Yes",age=0,count=1,
                             pre_no=0,alpha_no=0,delta_no=0,ba2_no=0,ba45_no=0,bq_no=0,
                             delta_14=0,ba2_14=0,ba45_14=0,bq_14=0,
                             delta_90=0,ba2_90=0,ba45_90=0,bq_90=0)) 


ggplot(p,aes(x,predicted,color=group))+
  geom_point(position=position_dodge(width=1.3))+
  geom_line(position=position_dodge(width=1.3))+
  geom_errorbar(aes(ymin=conf.low,ymax=conf.high),position=position_dodge(width=1.3),width=0.7)+
  labs(title="", x="Ct values in the most recent previous episode",y="Ct values in the current reinfection episode")+
  coord_cartesian(ylim=c(18,26))+
  scale_y_continuous(limits=c(0,50), breaks = seq(0,35,2))+
  scale_x_continuous(limits=c(9,36), breaks = seq(10,35,5))+
  theme_light()+
  theme(legend.position = "bottom")+
  scale_fill_jama(name="Symptoms in the most recent previous episode")+
  scale_color_jama(name="Symptoms in the most recent previous episode")



####################################################################################################
###################################################################################### symptom model
###compare symptom Yes vs No
#with ct
m <- glm(sympt~age + sex + ethnicity + hcw + lthc + IMD + ct + 
           variant + time_vac + episode,family = binomial(),
         data=data_reg)
summary(m)

coef=coef(m)
confint=confint(m) %>% as.data.frame()
p=summary(m)$coefficients[,4]

OR=cbind(exp(cbind(coef, confint)),p) %>% 
  as.data.frame() %>% 
  round(3) %>% 
  rownames_to_column() 

#without ct
m <- glm(sympt~age + sex + ethnicity + hcw + lthc + IMD + 
           variant + time_vac + episode,family = binomial(),
         data=data_reg)
summary(m)

coef=coef(m)
confint=confint(m) %>% as.data.frame()
p=summary(m)$coefficients[,4]

OR=cbind(exp(cbind(coef, confint)),p) %>% 
  as.data.frame() %>% 
  round(3) %>% 
  rownames_to_column() 


###only include CIS symptoms
#with ct
m <- nnet::multinom(symptom_CIS~age + sex + ethnicity + hcw + lthc + IMD + ct + 
                      variant + time_vac + episode, maxit=1000,
                    data=data_reg)

coef=t(coef(m)) %>% as.data.frame()

confint=confint(m) %>% as.data.frame() %>% 
  rename(class2l=`2.5 %.Classic`,class3l=`2.5 %.Other`,
         class2u=`97.5 %.Classic`,class3u=`97.5 %.Other`)

z<-summary(m)$coefficients/summary(m)$standard.errors
p<- t((1-pnorm(abs(z),0,1))*2) %>% as.data.frame() %>% rename(p2=`Classic`,p3=`Other`)

OR=exp(cbind(coef, confint)) %>% 
  as.data.frame() %>% 
  round(2) %>% 
  rownames_to_column() %>% 
  mutate(classicCI=paste0(class2l,"-",class2u),
         OtherCI=paste0(class3l,"-",class3u)) %>% 
  cbind(p) %>% 
  mutate(p2=ifelse(p2>0.1, round(p2,1),ifelse(p2<0.001,"<0.001",round(p2,3))),
         p3=ifelse(p3>0.1, round(p3,1),ifelse(p3<0.001,"<0.001",round(p3,3)))) %>% 
  select(rowname,Classic,classicCI,p2,Other,OtherCI,p3)


#without ct
m <- nnet::multinom(symptom_CIS~age + sex + ethnicity + hcw + lthc + IMD + 
                      variant + time_vac + episode,maxit=1000,
                    data=data_reg)

coef=t(coef(m)) %>% as.data.frame()

confint=confint(m) %>% as.data.frame() %>% 
  rename(class2l=`2.5 %.Classic`,class3l=`2.5 %.Other`,
         class2u=`97.5 %.Classic`,class3u=`97.5 %.Other`)

z<-summary(m)$coefficients/summary(m)$standard.errors
p<- t((1-pnorm(abs(z),0,1))*2) %>% as.data.frame() %>% rename(p2=`Classic`,p3=`Other`)

OR=exp(cbind(coef, confint)) %>% 
  as.data.frame() %>% 
  round(2) %>% 
  rownames_to_column() %>% 
  mutate(classicCI=paste0(class2l,"-",class2u),
         OtherCI=paste0(class3l,"-",class3u)) %>% 
  cbind(p) %>% 
  mutate(p2=ifelse(p2>0.1, round(p2,1),ifelse(p2<0.001,"<0.001",round(p2,3))),
         p3=ifelse(p3>0.1, round(p3,1),ifelse(p3<0.001,"<0.001",round(p3,3)))) %>% 
  dplyr::select(rowname,Classic,classicCI,p2,Other,OtherCI,p3)



####################################################################################################
############################################################################# kaplan_meier estimation

fit_km<- survfit(Surv(time, event)~1, km)
ggsurv=ggsurvplot(fit_km, size=1, censor.size=2,censor.shape=124,
                  conf.int = T, 
                  xlab="Time from first infection (days)",
                  legend="none",
                  cumevents = F,risk.table =F,
                  fontsize=3,
                  surv.median.line = "hv",
                  ggtheme=theme_light())

ggsurv$plot=ggsurv$plot+
  scale_x_continuous(breaks=c(0,50,200,500,1000))

fit_km<- survfit(Surv(time, event)~variant, km)

ggsurv=ggsurvplot(fit_km, size=1, censor=F,
                  fun="event",
                  conf.int = T, break.time.by=90,
                  ylim=c(0,1),
                  xlim=c(0,990),
                  xlab="Time from earlier infection (days)",
                  ylab="Reinfection percentage (%)",
                  legend.title="Variant of earlier infection",
                  legend.labs=c("Pre-Alpha","Alpha","Delta","BA.1","BA.2","BA.4/5","BQ.1/CH.1.1/XBB.1.5"),
                  cumevents = T,risk.table =T,
                  fontsize=3.5,
                  ggtheme=theme_light(),
                  palette = "nejm")

ggsurv$plot=ggsurv$plot+guides(color=guide_legend(nrow=1))+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0","25","50","75","100"))+
  geom_hline(yintercept = 0.25,linetype=2)+
  geom_hline(yintercept = 0.5,linetype=2)+
  geom_vline(xintercept = 180,linetype=2,color="red")+
  geom_vline(xintercept = 360,linetype=2,color="red") 


sum=surv_summary(fit_km) %>% filter(n.risk>=100)

surv_plot=ggsurvplot_df(sum, size=1, censor=F,
                        fun="event",
                        conf.int = T, break.time.by=90,
                        ylim=c(0,1),
                        xlim=c(0,990),
                        xlab="Time from earlier infection (days)",
                        ylab="Reinfection percentage (%)",
                        legend.title="Variant of earlier infection",
                        legend.labs=c("Pre-Alpha","Alpha","Delta","BA.1","BA.2","BA.4/5","BQ.1/CH.1.1/XBB.1.5"),
                        cumevents = T,risk.table =T,
                        fontsize=3.5,
                        ggtheme=theme_light(),
                        palette = "nejm")+
  guides(color=guide_legend(nrow=1))+
  coord_cartesian(ylim=c(0,1),xlim=c(0,990))+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0","25","50","75","100"))+
  geom_hline(yintercept = 0.25,linetype=2)+
  geom_hline(yintercept = 0.5,linetype=2)+
  geom_vline(xintercept = 180,linetype=2,color="red")+
  geom_vline(xintercept = 360,linetype=2,color="red") 

surv_table=ggsurv$table
surv_table2=ggsurv$cumevents

plot_grid(surv_plot,surv_table,surv_table2,nrow=3,align = "v",rel_heights = c(1,0.5,0.5))



####################################################################################################
############################################################################ flexible survival model

### baseline hazard
fit=stpm2(Surv(tstart,tstop,event)~1,data=data,df=8)
AIC(fit)

pred=predict(fit,newdata=data.frame(expand.grid(tstop=seq(1,42,1))),
             type="haz",se.fit=T,full=T)

baseline=ggplot(pred,aes(x=tstop,y=Estimate,ymin=lower,ymax=upper))+
  geom_ribbon(alpha=0.1)+
  geom_line()+
  theme_light()+
  scale_x_continuous(limits=c(0,42),breaks=seq(0,42,7),
                     sec.axis = sec_axis(~.,breaks = seq(0,42,7),labels = c("2021-12-27","2022-01-03",
                                                                            "2022-01-10","2022-01-17","2022-01-24",
                                                                            "2022-01-31","2022-02-07")))+
  coord_cartesian(ylim=c(0,0.01))+
  scale_y_continuous(limits = c(0,0.01),breaks=seq(0,0.01,0.002))+
  labs(y="Baseline hazard",x="Time from the start of infection wave (days)",title="Omicron BA.1")


### with covariates
fit_cont=stpm2(Surv(tstart,tstop,event)~ns(age,3) + sex + ethnicity + hcw + lthc + IMD + ns(ct,3) + symptom + region + previous + any_ct_lfd_previous + 
            pre_time + alpha_time + delta_time + variant + time_vaccination + prevalence, data=data, df=8)

summary(fit_cont)
a=summary(fit_cont)

a=a@coef %>% as.data.frame() %>% rownames_to_column()

a=a %>% mutate(hr=exp(Estimate),LCI=exp(Estimate-1.96*`Std. Error`),UCI=exp(Estimate+1.96*`Std. Error`))


fit_cat=stpm2(Surv(tstart,tstop,event)~ns(age,3) + sex + ethnicity + hcw + lthc + IMD + ns(ct,3) + symptom + region + previous + any_ct_lfd_previous + 
            comb + time_vaccination + prevalence, data=data, df=8)
summary(fit_cat)

a=summary(fit_cat)

a=a@coef %>% as.data.frame() %>% rownames_to_column()

a=a %>% mutate(hr=exp(Estimate),LCI=exp(Estimate-1.96*`Std. Error`),UCI=exp(Estimate+1.96*`Std. Error`))


### plot time from previous infection
newdata = data.frame(tstop=40,age=0,sex=c("Female"),
                     ethnicity="White",hcw="No",lthc="No",IMD=0,any_ct_lfd_previous="No",
                     ct=0,symptom="No",alpha_time=1,delta_time=1,other_time=1,
                     time_vaccination="No vaccination",region="England",prevalence=0.01,
                     previous="1")
newdata <- rbind(newdata,newdata)
newdata$variant <- c("Delta","Pre-Alpha")

dat <- data.frame(var = c(4,5,6,7))
dat$or <- NA
dat$ll <- NA
dat$ul <- NA

for (i in 1:nrow(dat)) {
  temp <- getor_vaccine_type(newdata=newdata,  contrast=c(1,dat$var[i]),
                             model = fit,var="pre_time")
  dat$or[i] <- temp[1]
  dat$ll[i] <- temp[2]
  dat$ul[i] <- temp[3]
  
}


newdata = data.frame(tstop=40,age=0,sex=c("Female"),
                     ethnicity="White",hcw="No",lthc="No",IMD=0,any_ct_lfd_previous="No",
                     ct=0,symptom="No",pre_time=1,delta_time=1,other_time=1,
                     time_vaccination="No vaccination",region="England",prevalence=0.01,
                     previous="1")
newdata <- rbind(newdata,newdata)
newdata$variant <- c("Delta","Alpha")

dat2 <- data.frame(var = c(2,3,4,5))
dat2$or <- NA
dat2$ll <- NA
dat2$ul <- NA


for (i in 1:nrow(dat2)) {
  temp <- getor_vaccine_type(newdata=newdata,  contrast=c(1,dat2$var[i]),
                             model = fit,var="alpha_time")
  dat2$or[i] <- temp[1]
  dat2$ll[i] <- temp[2]
  dat2$ul[i] <- temp[3]
  
}


newdata = data.frame(tstop=40,age=0,sex=c("Female"),
                     ethnicity="White",hcw="No",lthc="No",IMD=0,any_ct_lfd_previous="No",
                     ct=0,symptom="No",pre_time=1,alpha_time=1,other_time=1,
                     time_vaccination="No vaccination",region="England",prevalence=0.01,
                     previous="1")
newdata <- rbind(newdata,newdata)
newdata$variant <- c("Delta","Delta")

dat3 <- data.frame(var = c(1,2))
dat3$or <- NA
dat3$ll <- NA
dat3$ul <- NA

for (i in 1:nrow(dat3)) {
  temp <- getor_vaccine_type(newdata=newdata,  contrast=c(1,dat3$var[i]),
                             model = fit,var="delta_time")
  dat3$or[i] <- temp[1]
  dat3$ll[i] <- temp[2]
  dat3$ul[i] <- temp[3]
  
}


dat$variant="Pre-Alpha"
dat2$variant="Alpha"
dat3$variant="Delta"


dat_all_ba11=rbind(dat,dat2,dat3) %>% 
  mutate(var=case_when(var==1~"120-180",
                       var==2~"180-240",
                       var==3~"240-300",
                       var==4~"300-360",
                       var==5~"360-420",
                       var==6~"420-480",
                       var==7~">480") %>% factor(levels=c("120-180","180-240","240-300","300-360","360-420","420-480",">480"))) %>% 
  mutate(variant=factor(variant,levels = c("Pre-Alpha","Alpha","Delta")))

ggplot(dat_all_ba11,aes(var,log(or),color=variant,group=variant))+
  geom_point(position=position_dodge(width=0.5))+
  geom_line(position=position_dodge(width=0.5))+
  coord_cartesian(ylim=c(-1,1))+
  geom_errorbar(aes(ymin=log(ll),ymax=log(ul)),width=0.3,position=position_dodge(width=0.5))+
  theme_light()+
  scale_color_nejm()+
  scale_y_continuous(breaks = c(log(0.5),log(1),log(2)),labels = c(0.5,1,2))+
  geom_hline(yintercept = 0,color="red",linetype=2)+
  theme(legend.position = "bottom")+
  labs(x="Time from previous infection (days)",color="Variant of previous infection",
       y="Hazard ratios (95% CIs) of Omicron BA.1 reinfection")


### plot time from previous vaccination
newdata = data.frame(tstop=40,age=0,sex=c("Female"),pre_time=1,alpha_time=1,delta_time=1,variant="Delta",any_ct_lfd_previous="No",
                     ethnicity="White",hcw="No",lthc="No",IMD=0,
                     ct=0,symptom="No",previous="1",prevalence=0.01,
                     region="England")
newdata <- rbind(newdata,newdata)

dat <- data.frame(var = unique(df$time_vaccination)) 
dat$or <- NA
dat$ll <- NA
dat$ul <- NA
dat=dat %>% mutate(var=factor(var,levels=c("No vaccination","0-14","14-90","90-180",">180"))) 
for (i in 1:nrow(dat)) {
  temp <- getor_vaccine_type(newdata=newdata,  contrast=c(">180",as.character(dat$var[i])),
                             model = fit,var="time_vaccination")
  dat$or[i] <- temp[1]
  dat$ll[i] <- temp[2]
  dat$ul[i] <- temp[3]
  
}

ggplot(dat,aes(var,log(or)))+
  geom_point()+
  coord_cartesian(ylim=c(-3,1))+
  geom_errorbar(aes(ymin=log(ll),ymax=log(ul)),width=0.3)+
  theme_light()+
  scale_color_nejm()+
  scale_y_continuous(breaks = c(log(0.1),log(0.2),log(0.5),log(1),log(2)),labels = c(0.1,0.2,0.5,1,2))+
  geom_hline(yintercept = 0,color="red",linetype=2)+
  theme(legend.position = "bottom")+
  labs(x="Time from most recent vaccination (days)",
       y="Hazard ratios (95% CIs) of Omicron BA.1 reinfection")



### plot ct
newdata = data.frame(tstop=40,age=0,sex=c("Female"),
                     ethnicity="White",hcw="No",lthc="No",IMD=0,any_ct_lfd_previous="No",
                     symptom="No",previous="1",pre_time=1,alpha_time=1,delta_time=1,variant="Delta",
                     time_vaccination="No vaccination",region="England",prevalence=0.01,)
newdata <- rbind(newdata,newdata)

ct_omicron1 <- data.frame(var = seq(-12,13,1))
ct_omicron1$or <- NA
ct_omicron1$ll <- NA
ct_omicron1$ul <- NA

for (i in 1:nrow(ct_omicron1)) {
  temp <- getor_vaccine_type(newdata=newdata,  contrast=c(0,ct_omicron1$var[i]),
                             model = fit,var="ct")
  ct_omicron1$or[i] <- temp[1]
  ct_omicron1$ll[i] <- temp[2]
  ct_omicron1$ul[i] <- temp[3]
  
}

ct_omicron1=ct_omicron1 %>% mutate(var=var+22)

ggplot(ct_omicron1,aes(var,log(or)))+
  geom_point(position=position_dodge(width=1))+
  geom_line(position=position_dodge(width=1))+
  coord_cartesian(ylim=c(-1,1))+
  scale_x_continuous(breaks = seq(10,35,5))+
  scale_y_continuous(breaks = c(log(0.5),log(1),log(2)),labels = c(0.5,1,2))+
  geom_errorbar(aes(ymin=log(ll),ymax=log(ul)),width=0.3,position=position_dodge(width=1))+
  theme_light()+
  geom_hline(yintercept = 0,color="red",linetype=2)+
  labs(x="Ct value of previous infection",y="Hazard Ratios (95% CIs) of Omicron BA.1 reinfection")


### plot age
newdata = data.frame(tstop=40,sex=c("Female"),
                     ethnicity="White",hcw="No",lthc="No",IMD=0,any_ct_lfd_previous="No",
                     ct=0,symptom="No",previous="1",pre_time=1,alpha_time=1,delta_time=1,variant="Delta",
                     time_vaccination="No vaccination",region="England",prevalence=0.01,)
newdata <- rbind(newdata,newdata)

age_omicron1 <- data.frame(var = seq(-2,4.5,0.5))
age_omicron1$or <- NA
age_omicron1$ll <- NA
age_omicron1$ul <- NA

for (i in 1:nrow(age_omicron1)) {
  temp <- getor_vaccine_type(newdata=newdata,  contrast=c(0,age_omicron1$var[i]),
                             model = fit,var="age")
  age_omicron1$or[i] <- temp[1]
  age_omicron1$ll[i] <- temp[2]
  age_omicron1$ul[i] <- temp[3]
  
}


age_omicron1=age_omicron1 %>% mutate(var=var*10+40)

ggplot(age_omicron1,aes(var,log(or)))+
  geom_point(position=position_dodge(width=1))+
  geom_line(position=position_dodge(width=1))+
  coord_cartesian(ylim=c(-2,1))+
  scale_x_continuous(breaks = seq(20,85,5))+
  scale_y_continuous(breaks = c(log(0.2),log(0.5),log(1),log(2)),labels = c(0.2,0.5,1,2))+
  geom_errorbar(aes(ymin=log(ll),ymax=log(ul)),width=0.3,position=position_dodge(width=1))+
  theme_light()+
  geom_hline(yintercept = 0,color="red",linetype=2)+
  labs(x="Age (years)",y="Hazard Ratios (95% CIs) of Omicron BA.1 reinfection")



