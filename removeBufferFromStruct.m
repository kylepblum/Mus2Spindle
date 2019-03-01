function data = removeBufferFromStruct(data,bs)

for i = 1:numel(fieldnames(data))
    names = fieldnames(data);
    data.(names{i}) = data.(names{i})(bs+1:end);
end
end