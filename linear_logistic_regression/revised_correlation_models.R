# correlation models for Jilian; revised version

library(here)
library(tidyverse)
library(patchwork)

a <- read.csv(here("data/eDNA_distance_particle.csv"), row.names = 1) %>% 
  # filter(!Site %in% c("ESP", "POSITIVE")) %>% 
  mutate(Detected = ifelse(SiteMean_eDNA.copiesL.1. > 0, 1, 0))

# ecology solution to the log problem: add a small value to make zeros slightly non-zero
b <- a %>% 
  mutate(ModifiedSiteMean = SiteMean_eDNA.copiesL.1. + 1e-1) 
names(b) <- c("Site", "SiteMean", "SiteSD", "Distance", "NoWind", "Realistic", "Detected", "ModifiedSiteMean")

#note the SD scales w the mean, which is interesting
b %>% 
  filter(SiteMean < 250) %>% 
  ggplot(aes(x = SiteMean, y = SiteSD)) +
    geom_point()
    

#### Logistic Regressions

# Distance Only

logistic_distance <- glm(Detected ~ Distance, data = b, family = "binomial")
  summary(logistic_distance)

(p_logistic_distance <- b %>% 
    ggplot(aes(x = Distance, y = Detected)) +
    geom_jitter(width = .001, height = 0) +
    # geom_smooth(method = "glm", 
    #             method.args = list(family = "binomial"), 
    #             se = FALSE) +
    ylab("Probability of Detection") + xlab("Distance from Source (km)") + ggtitle("Distance Only"))


# No wind

logistic_noWind <- glm(Detected ~ NoWind, data = b, family = "binomial")
  summary(logistic_noWind)

(p_logistic_noWind <- b %>% 
    ggplot(aes(x = NoWind, y = Detected)) +
    geom_jitter(width = .001, height = 0) +
    # geom_smooth(method = "glm", 
    #             method.args = list(family = "binomial"), 
    #             se = FALSE) +
    ylab("Probability of Detection") + xlab("Particle Density (%)") + ggtitle("Particle Tracking, No Wind"))

# Realistic model 

  logistic_realistic <- glm(Detected ~ Realistic, data = b, family = "binomial")
    summary(logistic_realistic)
  
  (p_logistic_realistic <- b %>% 
      ggplot(aes(x = Realistic, y = Detected)) +
      geom_jitter(width = .001, height = 0) +
      geom_smooth(method = "glm", 
                  method.args = list(family = "binomial"), 
                  se = FALSE) +
      ylab("Probability of Detection") + xlab("Particle Density (%)") + ggtitle("Particle Tracking, with Wind"))


    ## chi-sq test isn't super helpful here, w small sample size
    
    data.frame(
      Obs = as.numeric(b$SiteMean > 0),
      Pred = as.numeric(b$Realistic > 0)
    ) %>% 
      xtabs(~Pred + Obs, data= .) %>% 
      chisq.test(simulate.p.value = T, B = 10000)
    


    #### Linear Regressions
    
    # Distance Only
    linear_distance <- lm(log(ModifiedSiteMean) ~ Distance, data = b)
      summary(linear_distance)
    
    (p_linear_distance <- b %>% 
        ggplot(aes(x = Distance, y = log(ModifiedSiteMean))) +
        geom_point() +
        # geom_jitter(width = .001, height = 0) +
        # geom_smooth(method = "lm", 
        #             method.args = list(family = "binomial"), 
        #             se = FALSE) +
        ylab("Log (Copies + 0.1 / L)") + xlab("Distance from Source (km)") + ggtitle("Distance Only"))

    # No wind

      linear_noWind <- lm(log(ModifiedSiteMean) ~ NoWind, data = b)
      summary(linear_noWind)
      
      (p_linear_noWind <- b %>% 
          ggplot(aes(x = NoWind, y = log(ModifiedSiteMean))) +
          geom_point() +
          # geom_jitter(width = .001, height = 0) +
          # geom_smooth(method = "lm", 
          #             method.args = list(family = "binomial"), 
          #             se = FALSE) +
          ylab("Log (Copies + 0.1 / L)") + xlab("Particle Density (%)") + ggtitle("Particle Tracking, No Wind"))
    
    # Realistic model 
    

      linear_realistic <- lm(log(ModifiedSiteMean) ~ Realistic, data = b)
      summary(linear_realistic)
      
      (p_linear_realistic <- b %>% 
          ggplot(aes(x = Realistic, y = log(ModifiedSiteMean))) +
          geom_point() +
          # geom_jitter(width = .001, height = 0) +
          geom_smooth(method = "lm",
                      # method.args = list(family = "binomial"),
                      se = T) +
          ylab("Log (Copies + 0.1 / L)") + xlab("Particle Density (%)") + ggtitle("Particle Tracking, with Wind"))
      

      # ASSEMBLE PLOT
      p_logistic_distance + p_logistic_noWind + p_logistic_realistic + 
        p_linear_distance + p_linear_noWind + p_linear_realistic +
        plot_layout(nrow = 2)

      ## gather relevant statistics:
      summary(linear_realistic)
      
      AIC(logistic_noWind,logistic_realistic) #here, AIC shows that these two models are very similar, with realistic slightly preferred
      