function retrievalTool_complete(retrievalTool)
%==========================================================================
% NAME          | retrievalTool_complete(retrievalTool)
% TYPE          |
% AUTHOR(S)     |
% CREATION      |
%               |
% ABSTRACT      | 
%               | 
%               |
%               |
% ARGUMENTS     | INPUTS:
%               |
%               | OUTPUTS:
%               |
% CALLS         |
%               |
%               |
%               |

%==========================================================================
required=retrievalTool.requiredFields;

for i = 1:size(required,2)
    assert(isfield(retrievalTool,required{i}),['Please complete the field: ' required{i} 'in retrievalTool'])
end

end

