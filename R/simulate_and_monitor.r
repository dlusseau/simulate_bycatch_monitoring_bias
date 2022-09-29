
library(extraDistr)
# we are going to use truncated Poisson distribution to estimate the number of bycatch occurring once a bycatch event happens

# the first function makes the "real" fishing year structured in 365 days during which each boat gets opportunities to 
# carry out "fishing events" which might correspond to a haul for example
# we have declared number of vessels in the fleet, an average number of fishing event per boat per day
# and then a probability of a bycatch event for each fishing event with an associated randomly generated number of bycatch per bycatch event
# here generated from a mixture of truncated Poisson distribution, the mixture is dictated by the 
# probability of a 'large' bycatch event occuring


make_fishing_year<-function(mean.bycatch.event=1,mean.bycatch.large.event=20,p.large.event=0.01,
							nboat=100,mean.fishing.event.boat.day=2) {

#we first only deal with one metier at a time
#bycatch is not affected by vessel characteristics
#later we can for example introduce vessel size for each boat 
# and eg influence probabilities by vessel size
# and subsequently influence monitoring by vessel size

# initialise the fishing year by fishing the first day
fishing.day<-1:365
fleet<-1:nboat							
fishing.event.per.boat<-rpois(nboat,mean.fishing.event.boat.day)

fishing<-data.frame(fishing.day=fishing.day[i],boat=rep(fleet,fishing.event.per.boat),bycatch=rbinom(sum(fishing.event.per.boat),1,p.bycatch),nbycatch=0)

event.type<-rbinom(sum(fishing$bycatch),1,p.large.event)
# for all bycatch event, was it a 'large' event or not

fishing$nbycatch[fishing$bycatch==1]<-apply(cbind((1-event.type)*rtpois(sum(fishing$bycatch),mean.bycatch.event,a=0),event.type*rtpois(sum(fishing$bycatch),mean.bycatch.large.event,a=0)),1,max)

#now we replicate for the whole year
for (i in 2:365) {

fishing.event.per.boat<-rpois(nboat,mean.fishing.event.boat.day)
temp<-data.frame(fishing.day=fishing.day[i],boat=rep(fleet,fishing.event.per.boat),bycatch=rbinom(sum(fishing.event.per.boat),1,p.bycatch),nbycatch=0)
event.type<-rbinom(sum(temp$bycatch),1,p.large.event)
temp$nbycatch[temp$bycatch==1]<-apply(cbind((1-event.type)*rtpois(sum(temp$bycatch),mean.bycatch.event,a=0),event.type*rtpois(sum(temp$bycatch),mean.bycatch.large.event,a=0)),1,max)

fishing<-rbind(fishing,temp)
}

return(fishing)
}


#NOW WE MONITOR

#first we simply look at the effect of monitoring rate
# we sample 1000 times in order to get a robust estimate of the estimate CV

monitor_BPUE<-function(pmonitor=0.1,nsample=1000,BPUE_real=0,fishing=NA) {
BPUE_est<-array(nsample)

for (i in 1:nsample) {
monitored<-sample(c(1:dim(fishing)[1]),floor(pmonitor*dim(fishing)[1]),replace=FALSE) # sample without replacement
BPUE_est[i]<-(sum(fishing$nbycatch[monitored])/length(monitored))
}

BPUE_est_mean<-mean(BPUE_est)
BPUE_est_CV<-sd(BPUE_est)/mean(BPUE_est)

return(list(BPUE_est=BPUE_est_mean,CV=BPUE_est_CV))
}


#TO DO: 
#extension:
#only monitor a small proportion of the fleet (ie, only some boats are sample in the fishing year 'real' data


fishing1<-make_fishing_year()

## first assessment what happens when I increase monitoring rate
p_monitor<-c(seq(.01,.2,.01),seq(.25,.5,.05))


monitor_estimate<-data.frame(p_monitor=p_monitor,BPUE_real=BPUE_real,BPUE_est=NA,BPUE_est_CV=NA)

for (i in 1:length(p_monitor)) {

monitor<-monitor_BPUE(pmonitor=p_monitor[i],BPUE_real=BPUE_real,fishing=fishing1)
monitor_estimate$BPUE_est[i]<-monitor$BPUE_est
monitor_estimate$BPUE_est_CV[i]<-monitor$CV
}

## now we are going to replicate the real fishing year replication 100 times
monitor_estimate$rel.bias<-(monitor_estimate$BPUE_real-monitor_estimate$BPUE_est)/monitor_estimate$BPUE_real
monitor_estimate$year<-1

for (y in 2:100) {
fishing1<-make_fishing_year()
BPUE_real<-sum(fishing1$nbycatch)/dim(fishing1)[1]


temp<-data.frame(p_monitor=p_monitor,BPUE_real=BPUE_real,BPUE_est=NA,BPUE_est_CV=NA)

for (i in 1:length(p_monitor)) {

monitor<-monitor_BPUE(pmonitor=p_monitor[i],BPUE_real=BPUE_real,fishing=fishing1)
temp$BPUE_est[i]<-monitor$BPUE_est
temp$BPUE_est_CV[i]<-monitor$CV
}

temp$rel.bias<-(temp$BPUE_real-temp$BPUE_est)/temp$BPUE_real
temp$year<-y


monitor_estimate<-rbind(monitor_estimate,temp)

}


###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
### VERSION 2: introducing between boat variability in both fishing and monitoring
### 29 September 2022

make_fishing_year<-function(mean.bycatch.event=1,mean.bycatch.large.event=20,p.large.event=0.01,
							nboat=100,mean.fishing.event.boat.day=2,p.bycatch=0.1,stochastic=TRUE) {
#we first only deal with one metier at a time
#bycatch is not affected by vessel characteeristics
#later we can for example introduce vessel size for each boat 
# and influence probabilities by vessel size
#and subsequently influence monitoring by vessel size

fishing.day<-1:365
fleet<-1:nboat							
if (stochastic==TRUE) {

mean.fishing.event.boat.day<-rtpois(nboat,mean.fishing.event.boat.day,a=0)  #introduce stochasticity so that the mean number of events per boats vary
fishing.event.per.boat<-rpois(nboat,mean.fishing.event.boat.day)

} else {
fishing.event.per.boat<-rpois(nboat,mean.fishing.event.boat.day) #uniform fishing behaviour

}

i=1
fishing<-data.frame(fishing.day=fishing.day[i],boat=rep(fleet,fishing.event.per.boat),bycatch=rbinom(sum(fishing.event.per.boat),1,p.bycatch),nbycatch=0)
event.type<-rbinom(sum(fishing$bycatch),1,p.large.event)
fishing$nbycatch[fishing$bycatch==1]<-apply(cbind((1-event.type)*rtpois(sum(fishing$bycatch),mean.bycatch.event,a=0),event.type*rtpois(sum(fishing$bycatch),mean.bycatch.large.event,a=0)),1,max)


for (i in 2:365) {

if (stochastic==TRUE) {
mean.fishing.event.boat.day<-rtpois(nboat,mean.fishing.event.boat.day,a=0)  #introduce stochasticity so that the mean number of events per boats vary
fishing.event.per.boat<-rpois(nboat,mean.fishing.event.boat.day)
} else {
fishing.event.per.boat<-rpois(nboat,mean.fishing.event.boat.day) #uniform fishing behaviour
}

temp<-data.frame(fishing.day=fishing.day[i],boat=rep(fleet,fishing.event.per.boat),bycatch=rbinom(sum(fishing.event.per.boat),1,p.bycatch),nbycatch=0)
event.type<-rbinom(sum(temp$bycatch),1,p.large.event)
temp$nbycatch[temp$bycatch==1]<-apply(cbind((1-event.type)*rtpois(sum(temp$bycatch),mean.bycatch.event,a=0),event.type*rtpois(sum(temp$bycatch),mean.bycatch.large.event,a=0)),1,max)

fishing<-rbind(fishing,temp)

}

return(fishing)
}



monitor_BPUE<-function(pmonitor=0.5,nsample=1000,BPUE_real=0,fishing=NA, p_monitor_boat=.1,boat_samp=TRUE) {

BPUE_est<-array(nsample)

for (i in 1:nsample) {

if (boat_samp==TRUE) {
monitored<-sample(c(1:dim(fishing)[1]),floor(pmonitor*dim(fishing)[1]),replace=FALSE) # sample without replacement
BPUE_est[i]<-(sum(fishing$nbycatch[monitored])/length(monitored))
} else {
boat_sampled<-sample(1:max(fishing$boat),n=floor(max(fishing$boat)*p_monitor_boat),replace=FALSE)
fleet_sampled<-fishing[fishing$boat%in%boat_sampled,]
monitored<-sample(c(1:dim(fleet_sampled)[1]),floor(pmonitor*dim(fleet_sampled)[1]),replace=FALSE) # sample without replacement
BPUE_est[i]<-(sum(fleet_sampled$nbycatch[monitored])/length(monitored))
}

}
BPUE_est_mean<-mean(BPUE_est)
BPUE_est_CV<-sd(BPUE_est)/mean(BPUE_est)

return(list(BPUE_est=BPUE_est_mean,CV=BPUE_est_CV))
}







###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
