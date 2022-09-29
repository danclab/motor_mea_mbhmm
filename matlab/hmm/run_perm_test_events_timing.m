function run_perm_test_events_timing(data, model)

min_time_steps=10;
n_perms=10000;

CIFcn = @(x,p)prctile(x(~isnan(x)),abs([0,100]-(100-p)/2));

% Times of this event in all trials
go_times= [];
mo_times = [];
oc_times = [];
pl_times = [];


trial_state_onsets={};
trial_state_offsets={};
for i=1:model.n_states
    trial_state_onsets{i}=[];
    trial_state_offsets{i}=[];
end

% For every date
for t_idx=1:data.ntrials
    % Rows of state seq for this trial
    trial_rows=find((model.state_seq.trial==t_idx));
                            
    if length(trial_rows)
        % Get the bins that we used in the HMM (time>0 and up to reward)
        bin_idx=find((data.bins>=0) & (data.bins<=data.metadata.reward(t_idx)));
        sub_bin_idx=find(data.bins(bin_idx)-data.metadata.go(t_idx)>=-500);

        % Find time of alignment event in this trial
        go_times(t_idx) = data.metadata.go(t_idx)-data.metadata.go(t_idx);
        mo_times(t_idx) = data.metadata.hand_mvmt_onset(t_idx)-data.metadata.go(t_idx);
        oc_times(t_idx) = data.metadata.obj_contact(t_idx)-data.metadata.go(t_idx);
        pl_times(t_idx) = data.metadata.place(t_idx)-data.metadata.go(t_idx);

        % Save p states within this window
        for i=1:model.n_states
            mon=trial_state_onsets{i};
            moff=trial_state_offsets{i};
            % Get mapped index
            state_idx=find(model.metadata.state_labels==i);
            
            % Get state activations from most likely state sequence
            above_thresh=model.state_seq.state(trial_rows(sub_bin_idx))==state_idx;
            onsets = strfind([0 above_thresh'], [0 ones(1,min_time_steps)]);
            offsets = strfind([above_thresh' 0], [ones(1,min_time_steps) 0]);
            durations=offsets-onsets;
            m_idx=1;
            if length(durations)
                mon(t_idx)=data.bins(bin_idx(sub_bin_idx(onsets(m_idx))))-data.metadata.go(t_idx);
                moff(t_idx)=data.bins(bin_idx(sub_bin_idx(offsets(m_idx))))-data.metadata.go(t_idx);                
            else
                mon(t_idx)=NaN;
                moff(t_idx)=NaN;
            end
            trial_state_onsets{i}=mon;  
            trial_state_offsets{i}=moff;
        end
    end
end

event_times.go=go_times;
event_times.mo=mo_times;
event_times.oc=oc_times;
event_times.pl=pl_times;

events={'go','mo','oc','pl'};

for i=1:model.n_states
    mon=trial_state_onsets{i};
    moff=trial_state_offsets{i};
    
    for e=2:length(events)
        event=events{e};
        e_times=event_times.(event);
        other_event_times=[];
        for e2=2:length(events)
            if e2~=e
                other_event_times(end+1,:)=event_times.(events{e2});
            end
        end
        for j=1:model.n_states
            if j~=i
                jmon=trial_state_onsets{j};
                jmoff=trial_state_offsets{j};
                if j>1
                    other_event_times(end+1,:)=jmon;
                end                
            end
        end
        if i>1
            b_onsets=mon;
            on_diff=partialcorr(b_onsets',e_times',other_event_times',...
                'rows','complete','type','Spearman');
            shuffled_diffs=[];
            for j=1:n_perms
                p_idx=randperm(length(e_times));
                shuffled_e_times=e_times(p_idx);
                shuffled_diffs(j)=partialcorr(b_onsets',shuffled_e_times',...
                    other_event_times(:,:)', 'rows','complete','type','Spearman');
            end        
            b=sum(abs(shuffled_diffs)>=abs(on_diff));
            p = (b+1)/(length(shuffled_diffs)+1);
            disp(sprintf('State: %d, Event: %s, Onset, r=%.3f, p=%.3f', i, event, on_diff, p));
        end
    
        other_event_times=[];
        for e2=2:length(events)
            if e2~=e
                other_event_times(end+1,:)=event_times.(events{e2});
            end
        end
        for j=1:model.n_states
            if j~=i
                jmon=trial_state_onsets{j};
                jmoff=trial_state_offsets{j};
                if j<6
                    other_event_times(end+1,:)=jmoff;
                end
            end
        end
        if i<6
            b_offsets=moff;
            off_diff=partialcorr(b_offsets',e_times',other_event_times',...
                'rows','complete','type','Spearman');
            shuffled_diffs=[];
            for j=1:n_perms
                p_idx=randperm(length(e_times));
                shuffled_e_times=e_times(p_idx);
                shuffled_diffs(j)=partialcorr(b_offsets',shuffled_e_times',...
                    other_event_times(:,:)', 'rows','complete','type','Spearman');
            end
            b=sum(abs(shuffled_diffs)>=abs(off_diff));
            p = (b+1)/(length(shuffled_diffs)+1);
            disp(sprintf('State: %d, Event: %s, Offset, r=%.3f, p=%.3f', i, event, off_diff, p));
        end
    end
    
end

