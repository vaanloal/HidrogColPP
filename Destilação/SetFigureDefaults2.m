function SetFigureDefaults2(hfig, hw_ratio, picturewidth, fontsize, ax)
    ax.XColor = [0, 0, 0];
    ax.YColor = [0, 0, 0];
    %ax.FontSize = 12;

    set(findall(hfig, '-property', 'FontSize'), 'FontSize', fontsize) % adjust fontsize to your document
    set(findall(hfig, '-property', 'Interpreter'), 'Interpreter', 'latex')
    set(findall(hfig, '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex')
    set(hfig, 'Units', 'centimeters', 'Position', [3 3 picturewidth hw_ratio * picturewidth])
    pos = get(hfig, 'Position');
    set(hfig, 'PaperPositionMode', 'Auto', 'PaperUnits', 'centimeters', 'PaperSize', [pos(3), pos(4)])

end
