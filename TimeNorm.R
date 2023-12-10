scalEst<-function(Data)
{ 
  Data<-Data[!apply(Data,1,function(x) all(x==0)),]  # delete rows with all zeros
  if(length(apply(Data,1,function(x){any(x==0)}))!=0)
  {
    Data <- Data
  }else{
    Data <- Data[!apply(Data,1,function(x){any(x==0)}),]#Find common features 
  }
  nnz<-apply(Data, 2, function(x) sum(x!=0)) # the number of nonzeros in each row
  if(any(nnz==0))
  {
    Data=Data+1
    nnz<-apply(Data, 2, function(x) sum(x!=0))
  }
  p.0<-rowSums(Data)/sum(Data)
  n.1<-colSums(Data/p.0)/nnz
  return(n.1)
}

DSS<-function(mData)
{
  scalEsts <- NULL
  scalEsts<-scalEst(mData)
  scalFactors<-scalEsts/median(scalEsts)
  normal.count<-t(apply(mData,1,function(s){s/scalFactors}))
  return(list(scalEsts=scalEsts, scalFactors=scalFactors,normal.count=normal.count))
}

between.time.normalization<-function(normal.dat.1,normal.dat.2,n.lib=c(10,10),show.ref.features=TRUE, show.all.features=FALSE, mtcm="BH", zsc=0.00)
{
  t12data<-data.frame(normal.dat.1,normal.dat.2)
  raida(t12data, n.lib=n.lib, show.ref.features=show.ref.features, show.all.features=show.all.features, mtcm=mtcm, zsc=zsc)->raidatest
  cf.data<-t12data[raidatest$reference.features,] #get stable features value from dataset
  factor.data1<-sum(cf.data[,(1:n.lib[1])]) #sum the features for baseline group
  factor.data2<-sum(cf.data[,c((n.lib[1]+1):(sum(n.lib)))])  # sum the features for trt group
  constant <- factor.data2/factor.data1 # calculate normalization factor
  #normal.dat.2 <- round(normal.dat.2/constant,0) # normalized trt group by dividing normalized factor
  normal.dat.2 <- normal.dat.2/constant
  return(list(normal.out=normal.dat.2,constant=constant,ref.genes=raidatest$reference.features))
}

#method=c("DSS"),

step1Norm<-function(wide.data, t=t,rep=rep,g=2)
   {
  modReshape<-function(mydata1,t,s.rep)
  {
    mydata <- reshape(mydata1,
                      idvar = "Feature",
                      varying =names(mydata1)[-1],
                      v.names = paste0("r",1:s.rep),            
                      times = 1:t,
                      timevar = "time",
                      direction = "long")
    return(mydata)
  }
  dnewdat<-data.frame(wide.data)
  c1.dat<-data.frame(Features=paste0("Feature",1:nrow(dnewdat)),dnewdat[,1:(t*rep)])
  c2.dat<-data.frame(Features=paste0("Feature",1:nrow(dnewdat)),dnewdat[,(t*rep+1):(g*t*rep)])
  names(c1.dat)<-c("Feature",paste(rep(paste0("t",1:t),each=rep),rep(paste0("r",1:rep),t),sep=""))->names(c2.dat)
  c1.long<-modReshape(c1.dat,t=t,s.rep=rep)
  c2.long<-modReshape(c2.dat,t=t,s.rep=rep)
  long.dat<-rbind(c1.long,c2.long)
  long.dat$group<-rep(c(1:2),c(nrow(c1.long),nrow(c2.long)))

normal.dat<-list()
for (ii in unique(long.dat$time))
{
t1c1data<-subset(long.dat,time==ii&group==1)
t1c2data<-subset(long.dat,time==ii&group==2)
row.names(t1c1data)<-as.character(t1c1data$Feature)->row.names(t1c2data)

t1cond1data<-dplyr::select(t1c1data,-Feature,-time,-group)
t1cond2data<-dplyr::select(t1c2data,-Feature,-time,-group)

t1cond1data<-dplyr::select(t1c1data,-Feature,-time,-group)
t1cond2data<-dplyr::select(t1c2data,-Feature,-time,-group)

#if (method=="DSS")
#  {
tdss1<-DSS(t1cond1data)
tdss2<-DSS(t1cond2data)
tdssdata1<-data.frame(tdss1$normal.count)
names(tdssdata1)<-paste0("g1",names(tdssdata1))
tdssdata2<-data.frame(tdss2$normal.count)
names(tdssdata2)<-paste0("g2",names(tdssdata2))
tdssdata<-data.frame(tdssdata1,tdssdata2)
#}
 
names(tdssdata)<-c(paste0("g1",names(t1cond1data)),paste0("g2",names(t1cond2data)))
normal.dat[[ii]]<-tdssdata
}
names(normal.dat)<-unique(long.dat$time)
return(normal.dat)
}


######### step 2 normalization
step2Norm<-function(step1.normal.data)
  {

extract.group<-function(s,whichgroup)
  { group1data<-dplyr::select(s,starts_with(whichgroup))
    return(group1data)
  }

c1.normal.dat<-lapply(step1.normal.data,extract.group,whichgroup="g1")## after first step C1, (output: 10 lists, each => timepoint, only group 1) 500 x 10(rep) x 10(timepoint)
c2.normal.dat<-lapply(step1.normal.data,extract.group,whichgroup="g2") ## after first step C2

s2.c1.normal.dat<-c1.normal.dat
s2.c2.normal.dat<-c2.normal.dat

### within condition normalization (under condition 1, between timepoint )
## if start with group 1, means use group1 as baseline. the g1T1 not normalized. 
## here, group1 is trt, group2 is control. use g2 as baseline
 ii=1
 c2.constant<-vector()
 ref.fea<-list()
while (ii<length(s2.c2.normal.dat))
  
   { n.lib.c2<-c(ncol(s2.c2.normal.dat[[ii]]),ncol(s2.c2.normal.dat[[ii+1]])) #[1] 10 10 #ncol of s2.c1.normal.dat is #rep, only condition1
     
     between.time.normalization(s2.c2.normal.dat[[ii]],s2.c2.normal.dat[[ii+1]],n.lib=n.lib.c2,show.ref.features=TRUE, show.all.features=FALSE, mtcm="BH", zsc=0.00)->newNextT
     s2.c2.normal.dat[[ii+1]]<-newNextT$normal.out
     c2.constant<-c(c2.constant,newNextT$constant)
     ref.fea[[ii]]<-newNextT$ref.genes
     ii<-ii+1
     #print(ii)
   }

## between condition normalization
 n.lib.c12<-c(ncol(s2.c2.normal.dat[[1]]),ncol(s2.c1.normal.dat[[1]])) #[1] 10 10 ###[[1]] means timepoint 1 for two groups
 between.time.normalization(s2.c2.normal.dat[[1]],s2.c1.normal.dat[[1]],n.lib=n.lib.c12,
                            show.ref.features=TRUE, show.all.features=FALSE, mtcm="BH", zsc=0.00)->c1newNextT1  ### normalization only at time1 for two groups
 s2.c1.normal.dat[[1]]<-c1newNextT1$normal.out## condition 1 time 1 to condition 2 time 1
     
c1.constant<-c1newNextT1$constant
c1.ref.fea<-list()
c1.ref.fea[[1]]<-c1newNextT1$ref.genes
hh=1
while (hh<length(s2.c1.normal.dat))
  
   { n.lib.c1<-c(ncol(s2.c1.normal.dat[[hh]]),ncol(s2.c1.normal.dat[[hh+1]]))
     between.time.normalization(s2.c1.normal.dat[[hh]],s2.c1.normal.dat[[hh+1]],n.lib=n.lib.c1,show.ref.features=TRUE, show.all.features=FALSE, mtcm="BH", zsc=0.00)->c1newNextT
     #s2.c2.normal.dat[[hh]]<-c2newNextT$normal.out
     s2.c1.normal.dat[[hh+1]]<-c1newNextT$normal.out
      c1.constant<-c(c1.constant,c1newNextT$constant)
     c1.ref.fea[[hh]]<-c1newNextT$ref.genes
     hh<-hh+1
     #print(hh)
   }

return(list(s2.c2.normal.dat=s2.c2.normal.dat,s2.c1.normal.dat=s2.c1.normal.dat,
            c1.constant=c1.constant,c2.constant=c2.constant,
            c1.ref.fea=c1.ref.fea,c2.ref.fea=ref.fea))
}
