function retrievalTool_complete(retrievalTool)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
required=retrievalTool.requiredFields;

for i = 1:size(required,2)
    assert(isfield(retrievalTool,required{i}),['Please complete the field: ' required{i} 'in retrievalTool'])
end

end

