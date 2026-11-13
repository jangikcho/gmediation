# print output is stable (snapshot)

    Code
      print(res)
    Output
      ----------------------------------------
      Exposure ('exposed' value): ses(1)
      Final Outcome: dmftd 
      Sample Size Used: 205 
      Number of Bootstrap Samples / Number Used: 3 / 3
      Reference Group: Whole Sample, refmult: 1 
      Scale: MeanDiff 
      Confidence Interval Percentage: 95 
      g (between-correlation factor): 1 
      ----------------------------------------

---

    Code
      print(run_gmediate(bin_models(), xcorrv = 0.5, boots = 3))
    Output
      ----------------------------------------
      Causal Mediation Sensitivity Analysis 
      Exposure ('exposed' value): ses(1)
      Final Outcome: dmftd 
      Sample Size Used: 205 
      Number of Bootstrap Samples / Number Used: 3 / 3 
      Reference Group: Whole Sample, Reference Group Multiplier (refmult)= 1 
      Scale: MeanDiff 
      Confidence Interval Level: 95 
      ----------------------------------------
      
      ====================================
      Path-Specific Effects
      ====================================
                  Path               Xcorr  Estimate  95% CI Lower  95% CI Upper  p-value 
      ses -> dmftd                     -     0.1706      0.1499        0.3160     0.0000  
      ses -> seal -> dmftd             -     0.0022     -0.0047        0.0172     0.6667  
      ses -> brush -> dmftd           0.5   -0.0061     -0.0113        0.0050     0.6667  
      ses -> brush -> seal -> dmftd   0.5    0.0015     -0.0028        0.0076     0.6667  
      
      ====================================
      Total Effect
      ====================================
      Estimate  95% CI Lower  95% CI Upper  p-value 
       0.1681      0.1634        0.3135     0.0000  
      
      ====================================
      First Stage Mediators
      ====================================
      Mediator  Estimate  95% CI Lower  95% CI Upper  p-value 
       brush    -0.0046     -0.0037        0.0022     0.6667  
      
      ====================================
      Second Stage Mediators
      ====================================
      Mediator  Xcorr  Estimate  95% CI Lower  95% CI Upper  p-value 
        seal     0.5    0.0036     -0.0075        0.0248     0.6667  
      
      ====================================
      Proportion: Path-Specific Effects
      ====================================
                  Path               Xcorr  Estimate  95% CI Lower  95% CI Upper  p-value 
      ses -> dmftd                     -     1.0148      0.9157        1.0266     0.0000  
      ses -> seal -> dmftd             -     0.0128     -0.0150        0.1061     0.6667  
      ses -> brush -> dmftd           0.5   -0.0364     -0.0694        0.0154     0.6667  
      ses -> brush -> seal -> dmftd   0.5    0.0088     -0.0088        0.0465     0.6667  
      
      ====================================
      Proportion: First Stage Mediators
      ====================================
      Mediator  Est Prop  95% CI Lower  95% CI Upper  p-value 
       brush    -0.0276     -0.0228        0.0065     0.6667  
      
      ====================================
      Proportion: Second Stage Mediators
      ====================================
      Mediator  Xcorr  Est Prop  95% CI Lower  95% CI Upper  p-value 
        seal     0.5    0.0216     -0.0239        0.1527     0.6667  

