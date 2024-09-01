# The-SIRV-Model-Project

In this project,
MATLAB will be used to implement the SIRV (Susceptible-Infectious-Recovered-Vaccinated) 
model for predicting pandemics. Here,
▪ Susceptible (S): Fraction of population who are not infected but are at risk of becoming 
infected.
▪ Infected (I): Fraction of population who are currently infected and can spread the disease 
to susceptible individuals.
▪ Recovered (R): Fraction of population who were infected and have either recovered or 
died, and are assumed to be immune to the disease.
▪ Vaccinated (V): Fraction of population who are vaccinated and as a result they don’t have 
the ability to spread the disease anymore.
▪ S+I+R+V=N (Total number of populations in the place of our research. We consider 
here total population is conserved.)


The discipline of studying the trends, causes, and consequences of health and illness problems in 
specific populations is known as epidemiology, and it makes extensive use of the SIRV model.
In order to comprehend the dynamics of infectious disease transmission and evaluate the possible 
effects of interventions like vaccination, epidemiological models- such as the SIRV model- are 
essential.

3.1. Theoretical Background:

The main simplified equations related to the SIRV model are-
𝑑𝑆(𝑡)
𝑑𝑡 = −
𝛽𝑆(𝑡)𝐼(𝑡)
𝑁
− 𝜌𝑆(𝑡) … … … … (1)
𝑑𝐼(𝑡)
𝑑𝑡 =
𝛽𝑆(𝑡)𝐼(𝑡)
𝑁
− 𝛾𝐼(𝑡) … … … … (2)
𝑑𝑅(𝑡)
𝑑𝑡 = 𝛾𝐼(𝑡) … … … … (3)
𝑑𝑉(𝑡)
𝑑𝑡 = 𝜌𝑆(𝑡) … … … … (4)
Here,
a. S is the number of susceptible individuals.
b. I is the number of infected individuals.
c. R is the number of recovered individuals.
d. β is the transmission rate (probability of transmitting disease between a susceptible and 
an infected individual).
e. γ is the recovery rate (inverse of the average duration of infection).
f. 𝜌 is the vaccination rate. (The rate at which a certain amount of people is getting 
vaccinated in our researching area.)

