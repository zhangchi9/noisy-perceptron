function saveallfigs(figtype)

if nargin == 0
    figtype = '.fig';
end
figHandles = findobj('Type', 'figure');
figaxes = findobj('Type', 'Axes');

for i = 1 : length(figHandles)
    
    curfig = figHandles(i);
    curaxes = figaxes(i);
    
    name = char([curaxes.Title.String]);
    
    if isempty(name)
        name = ['emffig',num2str(i)];
    end
    saveas(curfig,[name,figtype])
end