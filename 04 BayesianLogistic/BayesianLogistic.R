######################################################################
#### Using R2OpenBUGS to build Bayesian logistic multilevel model ####
######################################################################

library(R2OpenBUGS)
library(coda)
library(ggmcmc)


##### Data transformation

dfold = read.csv("d.csv")
# Change dummy coding: Female = 0 ; Male = 1
dfold$genderdummy = factor(dfold$genderdummy)
levels(dfold$genderdummy) = c(0, 1)

# R2OpenBUGS data format
idtrial = numeric(); numtrial = numeric()
for (i in unique(dfold$ID)){
    idtrial = append(idtrial, nrow(subset(dfold, ID==i)))
}
for (j in idtrial){
    numtrial = append(numtrial, rep(1:j))
}

cum_idtrial = cumsum(idtrial)
start_idtrial = append(cum_idtrial+1, 1, after = 0)
start_idtrial = start_idtrial[-length(start_idtrial)]

dfold = cbind(dfold, numtrial)
df = as.matrix(dfold[,-3])
Prob = dfold[,4]
Mag = dfold[,5]
Response = dfold[,6]
ID = dfold[,7]
Stimulation = dfold[,8]
Hedonism = dfold[,9]
Security = dfold[,10]
Gender = as.integer(dfold[,11])
my.data = list("df","start_idtrial","cum_idtrial","Prob","Mag",
               "Response","ID","Stimulation","Hedonism","Security","Gender")





##### Declare model

model = function(){
    for(i in 1:43){
        for(t in start_idtrial[i]:cum_idtrial[i]){
            Response[t] ~ dbern(p[t])
            logit(p[t]) <- B0 + B1*Security[t] + B2[i]*Prob[t] + 
                B3[i]*Mag[t] + B12*Security[t]*Prob[t] + 
                B13*Security[t]*Mag[t] + B23[i]*Prob[t]*Mag[t] + 
                B123*Security[t]*Prob[t]*Mag[t] + b[i] + B4*Gender[t]
        }
        B2[i] ~ dnorm(mu.b2, pre.b2) 
        B3[i] ~ dnorm(mu.b3, pre.b3) 
        B23[i] ~ dnorm(mu.b23, pre.b23)
        b[i] ~ dnorm(0,0.001)
        
    }
    B0 ~ dnorm(0,0.001)
    B1 ~ dnorm(0,0.001)
    B12 ~ dnorm(0,0.001)
    B13 ~ dnorm(0,0.001)
    B123 ~ dnorm(0,0.001)
    B4 ~ dnorm(0,0.001)
    # Hyperprior
    mu.b2 ~ dnorm(0,0.001)
    mu.b3 ~ dnorm(0,0.001)
    mu.b23 ~ dnorm(0,0.001)
    pre.b2 ~ dgamma(0.001,0.001)
    pre.b3 ~ dgamma(0.001,0.001)
    pre.b23 ~ dgamma(0.001,0.001)
}
my.model.file = "modelGW.odc" 
write.model(model, con = my.model.file)


# Initial values
inits = function() {  
    list(B0 = 0, B1 = 0, 
         B2 = rep(0,43), 
         B3 = rep(0,43),
         B4 = 0, B12 = 0, B13 = 0,
         B23 = rep(0,43), 
         B123 = 0,
         b = rep(0,43),
         mu.b2 = 0, mu.b3 = 0, mu.b23 = 0, pre.b2 = 0.0001, pre.b3 = 0.0001, pre.b23 = 0.0001
    )
}
params = c("B0","B1","B2", "B3","B4","B12","B13","B23","B123","b","mu.b2","mu.b3","mu.b23","pre.b2","pre.b3","pre.b23")


# Run the model
out = bugs(data=my.data, inits=inits, parameters.to.save=params, 
            model.file=my.model.file, codaPkg=T, n.iter=600, n.chains=1, 
            n.burnin=100, n.thin=10, debug = T,DIC = T)
out.coda = read.bugs(out)






##### Export the results

out.summary = summary(out.coda)
out.summary$stat[,1:2]
out.summary$q
HPDinterval(out.coda, prob = 0.95)


GMC = ggs(out.coda) # 先將coda的資料格式轉換成ggmcmc可讀的格式
str(GMC)
# ggmcmc(GMC, file="bayesian GW.pdf", plot=c("density", "running", "caterpillar"))
# 自動繪出MCMC圖表的pdf檔，plot指令可設定只畫哪些圖
# 以下這些指令預設都會將所有參數的圖繪出，可用family = "參數名稱"只畫出特定參數的圖。
# greek = T可自動將參數以希臘字母顯示(預設為FALSE)
GMC_hist = ggs_histogram(GMC)
GMC_den = ggs_density(GMC)
GMC_HPD = ggs_caterpillar(GMC, family = "[^deviance]") # 粗線為95%的highest posterior density，細線為90%的HPD。本處將參數deviance去除以避免scale被影響。
GMC_auto = ggs_autocorrelation(GMC)
ggs_crosscorrelation(GMC) # Plot the Cross-correlation between-chains.
ggs_pairs(GMC) # Evaluate posterior correlations among parameters.
GMC_trace = ggs_traceplot(GMC)