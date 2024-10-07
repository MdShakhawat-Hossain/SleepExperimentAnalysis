function [events]=eventevokedchanges(dF,Speed_Out)
    % mean value for all resting period
    for restidx = 1:1:length(Speed_Out.Rest(:,1))
        events.rest.F465(restidx) = mean(dF.F465(Speed_Out.Rest(restidx,1):Speed_Out.Rest(restidx,2))); 
    end
    events.rest_average.F465 = mean(events.rest.F465);

    for runidx = 1:1:length(Speed_Out.Run(:,1))
        events.run.F465(runidx) = mean(dF.F465(Speed_Out.Run(runidx,1):Speed_Out.Run(runidx,2))); 
    end
    events.run_average.F465 = mean(events.run.F465);
    %      
    for restidx = 1:1:length(Speed_Out.Rest(:,1))
        events.rest.F560(restidx) = mean(dF.F560(Speed_Out.Rest(restidx,1):Speed_Out.Rest(restidx,2))); 
    end
    events.rest_average.F560 = mean(events.rest.F560);
            
    for runidx = 1:1:length(Speed_Out.Run(:,1))
        events.run.F560(runidx) = mean(dF.F560(Speed_Out.Run(runidx,1):Speed_Out.Run(runidx,2))); 
    end
    events.run_average.F560 = mean(events.run.F560);
    % plot the event triggered changes in fluoresence
    
    figure; 
    subplot(211)
    plot(ones(size(events.rest.F465)),events.rest.F465,'o');
    hold on
    plot(2*ones(size(events.run.F465)),events.run.F465,'o');
    boxplot([events.rest.F465;events.run.F465]', {'Rest','Run'})
    hold off 
    ylabel('Average \DeltaF/F');
    title('Event triggered change in 465 fluoresence');
        
    subplot(212)
    plot(ones(size(events.rest.F560)),events.rest.F560,'o');
    hold on
    plot(2*ones(size(events.run.F560)),events.run.F560,'o');
    boxplot([events.rest.F560;events.run.F560]', {'Rest','Run'})
    hold off
    xlabel('Event');
    ylabel('Average \DeltaF/F');
    title('Event triggered change in 560');