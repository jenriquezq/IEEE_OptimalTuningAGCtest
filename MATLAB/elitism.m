 function isbetter= elitism(xnew,xtable,i)
 isbetter = false;
        if (xnew.feasibility && xtable.feasibility(i))&&(xnew.objective < xtable.objective(i))
            isbetter = true;
        elseif (~xnew.feasibility && ~xtable.feasibility(i))&&(xnew.fitness <xtable.fitness(i))
             isbetter = true;
        elseif ( xnew.feasibility && ~xtable.feasibility(i))
             isbetter = true;
        end
 end
 
 
 
 