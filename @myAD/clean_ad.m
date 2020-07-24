function cleaned_cell = clean_ad(varargin)
% by SeHyoun Ahn, July 2020

  for iter = 1:length(varargin)
    if isa(varargin{iter}, 'myAD')
      varargin{iter} = getvalues(varargin{iter});
    end
  end
  cleaned_cell = varargin;
end
