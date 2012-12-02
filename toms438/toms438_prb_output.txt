 
TOMS438_PRB
  Test TOMS algorithm 438, product-type
  two point Gauss-Legendre Simpson
  integration.
 
TEST01
  Integral of F(X) * G(X) from -1 to 1,
  with F(X) = 1, G(X) = 1/(1+x**2 )
 
       N    VINT
 
       1     5.49020    
       2     2.47843    
       3     2.90842    
       4     2.57255    
       5     2.69529    
       6     2.63329    
       7     2.66030    
       8     2.64773    
       9     2.65342    
      10     2.65082    
 
  Exact:     2.65164    
 
TEST02
  Integral of F(X) * G(X) from 0 to PI,
  with F(X) = EXP(-X), G(X) = COS(X)
 
       N    VINT
 
       1    0.390787    
       2    0.516996    
       3    0.520913    
       4    0.521413    
       5    0.521532    
       6    0.521572    
       7    0.521589    
       8    0.521596    
       9    0.521600    
      10    0.521603    
 
  Exact:    0.521607    
 
TOMS438_PRB
  Normal end of execution.
