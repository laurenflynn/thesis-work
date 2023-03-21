library(tidyverse)

fitMedianPolished <- rlm(fullCarrierRisk ~ carrierRiskUnaffectedInfoMasked, data = summaryTable)
print(summary(fitMedianPolished))

# finally, plot the data and the fitted curve
medPolVis = ggplot(data = summaryTable, aes(x= carrierRiskUnaffectedInfoMasked, y= fullCarrierRisk)) +
  geom_point() +
  geom_line(aes(y = predict(fitMedianPolished)))+
  labs(x = "Carrier Risk Unaffected Info Masked", y= "Carrier Risk for Full Family Information", title = "Comparing Carrier Risk for Masked Info and Full Family Information with Loess Regression")

print(medPolVis)

fitLoess <- loess(fullCarrierRisk ~ carrierRiskUnaffectedInfoMasked, data = summaryTable, span = 0.2)
summary(fitLoess) #RSE is 0.00397 (seems good but I'm not sure)

# finally, plot the data and the fitted curve
loessVis = ggplot(data = summaryTable, aes(x= carrierRiskUnaffectedInfoMasked, y= fullCarrierRisk)) +
  geom_point() +
  geom_line(aes(y = predict(fitLoess)))+
  labs(x = "Carrier Risk Unaffected Info Masked", y= "Carrier Risk for Full Family Information", title = "Comparing Carrier Risk for Masked Info and Full Family Information with Loess Regression")

print(loessVis)