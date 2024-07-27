
# inputs from TRANSFIL model - costs and effectiveness (DALYs averted for lymphoedema and hydrocele, prob of elimination)
c <- vector(mode = "list", length = 12)
e<- vector(mode = "list", length = 12)
prob <- vector(mode = "list", length = 12)
names(c) <- c(" > 5 years , 0.5%", " >5 years, 1%"," > 5 years, 2%", " >5 years, 5%" ," >20 years, 0.5%",
              ">20 years, 1%",">20 years, 2%","all, 0.5%" ,"all, 1%", "all, 2%","all, 5%", "> 20 years, 5%")

names(c) <- c(" 0.5% , > 5 years" , "1%, >5 years","2%, >5 years","5%,>5 years",
              " 0.5% , > 20 years" , "1%, >20 years","2%, >20 years","5%,>20 years",
              " 0.5% , all" , "1%, all","2%, all","5%, all")


k_daly <- 288.57 # change this parameter 

c[[1]] <- cost_0.005_5_10_5
c[[2]] <- cost_0.005_10_20_20
c[[3]] <- cost_0.005_10_20_all
c[[4]] <- cost_0.01_10_20_5
c[[5]] <- cost_0.01_10_20_20
c[[6]] <- cost_0.01_10_20_all
c[[7]] <- cost_0.005_10_20_5
c[[8]] <- cost_0.02_10_20_20
c[[9]] <- cost_0.02_10_20_all
c[[10]] <- cost_0.05_10_20_5
c[[11]] <- cost_0.05_10_20_20
c[[12]] <- cost_0.05_10_20_all



# effectiveness

e[[1]] <- (l_0.005_10_20_5+h_0.005_10_20_5)*100
e[[2]] <- (l_0.005_10_20_20+h_0.005_10_20_20)*100
e[[3]] <- (l_0.005_10_20_all+h_0.005_10_20_all)*100
e[[4]] <- (l_0.01_10_20_5+h_0.01_10_20_5)*100
e[[5]] <- (l_0.01_10_20_20+h_0.01_10_20_20)*100
e[[6]] <- (l_0.01_10_20_all+h_0.01_10_20_all)*100
e[[7]] <- (l_0.02_10_20_5+h_0.02_10_20_5)*100
e[[8]] <- (l_0.02_10_20_20+h_0.02_10_20_20)*100
e[[9]] <- (l_0.02_10_20_all+h_0.02_10_20_all)*100
e[[10]] <- (l_0.05_10_20_5+h_0.05_10_20_5)*100
e[[11]] <- (l_0.05_10_20_20+h_0.05_10_20_20)*100
e[[12]] <- (l_0.05_10_20_all+h_0.05_10_20_all)*100

# elimination
prob[[1]] <- prob_elim_0.005_10_20_5
prob[[2]] <- prob_elim_0.005_10_20_15
prob[[3]] <- prob_elim_0.005_5_10_all
prob[[4]] <- prob_elim_0.01_5_10_5
prob[[5]] <- prob_elim_0.01_10_20_15
prob[[6]] <- prob_elim_0.01_10_20_all
prob[[7]] <- prob_elim_0.02_10_20_5
prob[[8]] <- prob_elim_0.02_10_20_15
prob[[9]] <- prob_elim_0.02_10_20_all
prob[[10]] <- prob_elim_0.05_10_20_5
prob[[11]] <- prob_elim_0.05_10_20_15
prob[[12]] <- prob_elim_0.05_10_20_all

# cost and effectiveness by strategy and simulation
n_samples=100
library("data.table")
ce <- data.table(sample = rep(seq(n_samples), length(e)),
                 strategy = rep(paste0("Strategy ", seq(1, 4)), 
                                each = n_samples * 3),
                 grp = rep(rep(c(">5 years", ">20 years","all"),
                               each = n_samples), 3),
                 cost = do.call("c", c), qalys = do.call("c", e))
head(ce)

library("hesim")
ktop <- 5
cea_out <-  cea(ce, k = k_daly_tanzania, sample = "sample", strategy = "strategy",
                grp = "grp", e = "qalys", c = "cost")
cea_pw_out <-  cea_pw(ce,  k = k_daly_tanzania, comparator = "Strategy 2",
                      sample = "sample", strategy = "strategy", grp = "grp",
                      e = "qalys", c = "cost")
head(cea_pw_out$delta)
library("ggplot2")


theme_set(theme_bw())
plot_ceplane(cea_pw_out, k = 500)
plot_ceac(cea_pw_out)

# WTP for elimination
n_samples=100
library("data.table")
ce <- data.table(sample = rep(seq(n_samples), length(e)),
                 strategy = rep(paste0("Strategy ", seq(1, 4)), 
                                each = n_samples * 3),
                 grp = rep(rep(c(">5 years", ">20 years","all"),
                               each = n_samples), 3),
                 cost = do.call("c", c), qalys = do.call("c", e), prob_elim=do.call("c",prob))
head(ce)
k_prob <- 50000
ce <- ce[, nmb := k_daly_tanzania*qalys+k_prob*(prob_elim)*100 - cost]

cea_elim_out <-  cea_pw(ce,  k = seq(0,k_prob,5000) , comparator = "Strategy 2",
                      sample = "sample", strategy = "strategy", grp = "grp",
                      e = "qalys", c = "cost", prob_elim="prob")


plot_ceplane(cea_elim_out)
write.csv(cea_elim_out$ceac,"icers_haiti_5_10.csv")
