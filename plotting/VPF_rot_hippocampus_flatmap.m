function g = VPF_rot_hippocampus_flatmap(g)
pause(0.2)
set(g, 'XAxisLocation', 'origin','XAxisLocation','origin');
view(-90, 90);
set(g, 'XDir','reverse');
axis tight
end