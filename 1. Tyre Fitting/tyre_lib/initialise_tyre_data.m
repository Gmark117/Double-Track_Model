% Initialise tyre data struct
function [tyre_data] = initialise_tyre_data()

 % precode

  

 % main code

  tyre_data.FZ0             = 0; % FZ0            
  tyre_data.pCx1            = 0; % pCx1           
  tyre_data.pDx1            = 0; % pDx1           
  tyre_data.pDx2            = 0; % pDx2           
  tyre_data.pDx3            = 0; % pDx3           
  tyre_data.pEx1            = 0; % pEx1           
  tyre_data.pEx2            = 0; % pEx2           
  tyre_data.pEx3            = 0; % pEx3           
  tyre_data.pEx4            = 0; % pEx4           
  tyre_data.pHx1            = 0; % pHx1           
  tyre_data.pHx2            = 0; % pHx2           
  tyre_data.pKx1            = 0; % pKx1           
  tyre_data.pKx2            = 0; % pKx2           
  tyre_data.pKx3            = 0; % pKx3           
  tyre_data.pVx1            = 0; % pVx1           
  tyre_data.pVx2            = 0; % pVx2           
  tyre_data.pCy1            = 0; % pCy1           
  tyre_data.pDy1            = 0; % pDy1           
  tyre_data.pDy2            = 0; % pDy2           
  tyre_data.pDy3            = 0; % pDy3           
  tyre_data.pEy1            = 0; % pEy1           
  tyre_data.pEy2            = 0; % pEy2           
  tyre_data.pEy3            = 0; % pEy3           
  tyre_data.pEy4            = 0; % pEy4           
  tyre_data.pHy1            = 0; % pHy1           
  tyre_data.pHy2            = 0; % pHy2           
  tyre_data.pHy3            = 0; % pHy3           
  tyre_data.pKy1            = 0; % pKy1           
  tyre_data.pKy2            = 0; % pKy2           
  tyre_data.pKy3            = 0; % pKy3           
  tyre_data.pVy1            = 0; % pVy1           
  tyre_data.pVy2            = 0; % pVy2           
  tyre_data.pVy3            = 0; % pVy3           
  tyre_data.pVy4            = 0; % pVy4           
  tyre_data.qBz1            = 0; % qBz1           
  tyre_data.qBz10           = 0; % qBz10          
  tyre_data.qBz2            = 0; % qBz2           
  tyre_data.qBz3            = 0; % qBz3           
  tyre_data.qBz4            = 0; % qBz4           
  tyre_data.qBz5            = 0; % qBz5           
  tyre_data.qBz9            = 0; % qBz9           
  tyre_data.qCz1            = 0; % qCz1           
  tyre_data.qDz1            = 0; % qDz1           
  tyre_data.qDz2            = 0; % qDz2           
  tyre_data.qDz3            = 0; % qDz3           
  tyre_data.qDz4            = 0; % qDz4           
  tyre_data.qDz6            = 0; % qDz6           
  tyre_data.qDz7            = 0; % qDz7           
  tyre_data.qDz8            = 0; % qDz8           
  tyre_data.qDz9            = 0; % qDz9           
  tyre_data.qEz1            = 0; % qEz1           
  tyre_data.qEz2            = 0; % qEz2           
  tyre_data.qEz3            = 0; % qEz3           
  tyre_data.qEz4            = 0; % qEz4           
  tyre_data.qEz5            = 0; % qEz5           
  tyre_data.qHz1            = 0; % qHz1           
  tyre_data.qHz2            = 0; % qHz2           
  tyre_data.qHz3            = 0; % qHz3           
  tyre_data.qHz4            = 0; % qHz4           
  tyre_data.pDy1            = 0; % pDy1           
  tyre_data.pDy2            = 0; % pDy2           
  tyre_data.pDy3            = 0; % pDy3           
  tyre_data.rBx1            = 0; % rBx1           
  tyre_data.rBx2            = 0; % rBx2           
  tyre_data.rBy1            = 0; % rBy1           
  tyre_data.rBy2            = 0; % rBy2           
  tyre_data.rBy3            = 0; % rBy3           
  tyre_data.rCx1            = 0; % rCx1           
  tyre_data.rCy1            = 0; % rCy1           
  tyre_data.rHx1            = 0; % rHx1           
  tyre_data.rHy1            = 0; % rHy1           
  tyre_data.rVy1            = 0; % rVy1           
  tyre_data.rVy2            = 0; % rVy2           
  tyre_data.rVy3            = 0; % rVy3           
  tyre_data.rVy4            = 0; % rVy4           
  tyre_data.rVy5            = 0; % rVy5           
  tyre_data.rVy6            = 0; % rVy6           
  tyre_data.LCX             = 1; % LCX            
  tyre_data.LCY             = 1; % LCY            
  tyre_data.LEX             = 1; % LEX            
  tyre_data.LEY             = 1; % LEY            
  tyre_data.LFZ0            = 1; % LFZ0           
  tyre_data.LGAMMAY         = 1; % LGAMMAY        
  tyre_data.LHX             = 1; % LHX            
  tyre_data.LHY             = 1; % LHY            
  tyre_data.LKA             = 1; % LKA            
  tyre_data.LKXK            = 1; % LKXK           
  tyre_data.LKY             = 1; % LKY            
  tyre_data.LMR             = 1; % LMR            
  tyre_data.LMUX            = 1; % LMUX           
  tyre_data.LMUY            = 1; % LMUY           
  tyre_data.LT              = 1; % LT             
  tyre_data.LVX             = 1; % LVX            
  tyre_data.LVY             = 1; % LVY            
  tyre_data.LVYK            = 1; % LVYK           
  tyre_data.LXA             = 1; % LXA            
  tyre_data.LYK             = 1; % LYK            
  
 end
