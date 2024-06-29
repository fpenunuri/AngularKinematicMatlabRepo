%vars = args2duals(varargin)
%
%If any of the varargin arguments are dual, this function will convert the 
%others to dual numbers and place them in the cell array vars.
%

function result = args2duals(varargin)
  numArgs = nargin;
  strc = cell(1,numArgs);
  result = varargin;
  for k =1:numArgs
    strc{k} = class(varargin{k});
  end

  if any(strcmp(strc,'dual4'))
    for k =1:numArgs
      result{k} = dual4(varargin{k});
    end
  elseif any(strcmp(strc,'dual3'))
    for k =1:numArgs
      result{k} = dual3(varargin{k});
    end    
  elseif any(strcmp(strc,'dual2'))
    for k =1:numArgs
      result{k} = dual2(varargin{k});
    end
  elseif any(strcmp(strc,'dual1'))
    for k =1:numArgs
      result{k} = dual1(varargin{k});
    end
  end
end
