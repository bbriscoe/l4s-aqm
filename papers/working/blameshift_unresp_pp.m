# post-process output matrix
# shift flow a's queue from output(:,2) to output(:,4) if not at head
output(:,4) = output(:,2) .* output(:,5);
output(:,2) = output(:,2) .* !output(:,5);
# produce staggered plot
output(:,3) = output(:,2) + output(:,3);
output(:,4) = output(:,3) + output(:,4);