function y = loadSingleVariableMATFile(filename)
    
foo = load(filename);
whichVariables = fieldnames(foo);

if numel(whichVariables) == 1
    y = foo.(whichVariables{1});
else
    error(['there are many variables in ', filename, ' so matlab doesn''t know which to load'])
end
