                   if secondLoop == 1
                       break
                   end
                   
                   if k == 1
                       secondLoop = 1;
                       exitFlag = 0;
                   elseif k == length(controlPointsY)
                       secondLoop = 1;
                       exitFlag = 0;
                   else
                       secondLoop = 0;
                       exitFlag = 1;
                   end