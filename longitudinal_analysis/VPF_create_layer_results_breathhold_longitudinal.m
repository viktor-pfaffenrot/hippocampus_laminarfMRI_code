function out = VPF_create_layer_results_breathhold_longitudinal(layers,idx)


subfields = size(layers,2);
N = size(layers{1,1,1},2);
echoes = size(layers{1,1,1},4);
TE = [2 6 10 14 17 20]; %17 21

con_array = zeros(subfields-1,N,echoes);
for ss = 1:subfields-1
    if ss < subfields-1%all but CA4/DG
        Y = squeeze(mean(layers{1,ss}(idx{1,ss},:,:,:),1))...
            + squeeze(mean(layers{2,ss}(idx{2,ss},:,:,:),1));
    else
        % CA4 marks 'outer'. DG marks 'inner'
        % Interpolate final T and beta between inner and outer for visualization purposes
            %CA4
            tmp = (squeeze(mean(layers{1,end-1}(idx{1,end},:,:,:),1))...
                + squeeze(mean(layers{2,end-1}(idx{2,end},:,:,:),1)));
            %DG
            tmp2 = (squeeze(mean(layers{1,end},1))...
                + squeeze(mean(layers{2,end},1)));

            Y = permute(cat(3,tmp2,tmp), [3 1 2]);

            lindex = linspace(1, size(Y,1), N-9);
            Y = cat(1,nan([9,2,6]),interp1(1:size(Y,1), Y, lindex, 'linear', 0));

    end

    %the code above adds the signals of both hemispheres. Divide by 2 to average
    Y = Y./2;

    con_array(ss,:,:) = Y(:,2,:) - Y(:,1,:);

end

w = 1./TE;
w = w./sum(w);
con_array = squeeze(mean(bsxfun(@times,con_array,reshape(w,[1 1 6])),3));
con_array([1 4 5],1:9) = nan;
out = struct('con_array',con_array);
end