
# Small helper to formate TGraph objects more efficient and clean

Graph title, Axis titles and range can all set in one line of code.
Additionally the axis titles are centered. More features can be added.

## Just use this line after you have created the TGraohErrors object and before you draw it:


```
DrawGraph(TGraphErrors *graph, char const *name, char const *namex, char const *namey, double startx, double endx, double starty, double endy)
```

Especially the start and end values can be set by a config file and you don't need to recompile everytime you change the range.