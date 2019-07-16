function par_new = PAR_POS_UPDATE(par_new,par_old,dt,diff_co_eff,total_nb_particles,xStart,xEnd)
parfor (i=1:total_nb_particles,0)
        par_new(1,i)=par_old(1,i)+dt * par_old(2,i) + sqrt(2*diff_co_eff)* sqrt(dt)*randn;
        %periodic boundary condition
        if par_new(1,i) > xEnd
            par_new(1,i) = par_new(1,i)-xEnd;
        elseif par_new(1,i) < xStart
            par_new(1,i) = par_new(1,i)+xEnd;
        end
    end
end