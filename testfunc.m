function testfunc(records)
%TESTFUNC Summary of this function goes here
%   Detailed explanation goes here
n_records = length(records);
% dur_total = 0;
dur_TN = 0;
dur_TP = 0;
dur_FN = 0;
dur_FP = 0;
for i = 1:n_records
    rec = records{i};
    tint = rec.tint;
    % dur_total = get_duration_tints(tint);
    true_SLAMS = rec.true_SLAMS;
    true_notSLAMS = subtract_tints(tint, true_SLAMS);
    predicted_SLAMS = rec.predicted_SLAMS;
    predicted_notSLAMS = subtract_tints(tint, predicted_SLAMS);
    dur_TN = dur_TN + get_duration_tints(intersect_tints(true_notSLAMS, predicted_notSLAMS));
    dur_TP = dur_TP + get_duration_tints(intersect_tints(true_SLAMS, predicted_SLAMS));
    dur_FN = dur_FN + get_duration_tints(subtract_tints(true_SLAMS, predicted_SLAMS));
    dur_FP = dur_FP + get_duration_tints(subtract_tints(predicted_SLAMS, true_SLAMS));
end
confusion_matrix = round([dur_TP, dur_FN;
                          dur_FP, dur_TN]/(1e9));
figure;
ch = confusionchart(confusion_matrix, ["SLAMS", "not SLAMS"]);
ch.RowSummary = 'row-normalized';
ch.ColumnSummary = 'column-normalized';
% acc = correct_duration/total_duration;
end

