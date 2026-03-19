#include "drawgraph.h"

void DrawGraph(TGraphErrors *graph, char const *name, char const *namex, char const *namey, double startx, double endx, double starty, double endy){
  graph->SetTitle(name);
  graph->GetXaxis()->SetTitle(namex);
  graph->GetYaxis()->SetTitle(namey);  
  graph->GetXaxis()->SetRangeUser(startx,endx);
  graph->GetYaxis()->SetRangeUser(starty,endy);
  graph->GetXaxis()->CenterTitle(true);
    graph->GetYaxis()->CenterTitle(true);

}
