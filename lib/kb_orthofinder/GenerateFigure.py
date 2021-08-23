import os
from PyGrace.grace import Grace
from PyGrace.Styles.el import ElGraph
from PyGrace.drawing_objects import DrawText

def define_graph(graph,xlabel,ylabel,xmax,xmin,ymax,ymin,xmajor,ymajor,xminor,yminor,xprec,yprec):
    ##limits
    graph.set_world(xmin,ymin,xmax,ymax)
    ##ticks
    graph.xaxis.tick.major=xmajor
    graph.yaxis.tick.major=ymajor
    graph.xaxis.tick.onoff="on"
    graph.yaxis.tick.onoff="on"
    graph.xaxis.tick.minor_ticks=xminor
    graph.yaxis.tick.minor_ticks=yminor
    graph.xaxis.tick.inout="out"
    graph.yaxis.tick.inout="out"
    graph.xaxis.tick.place="normal"
    graph.yaxis.tick.place="normal"
    ##labels
    graph.xaxis.label.text=xlabel
    graph.yaxis.label.text=ylabel
    graph.xaxis.ticklabel.prec=xprec
    graph.yaxis.ticklabel.prec=yprec
    graph.yaxis.label.font=6
    graph.yaxis.label.char_size=1.65
    graph.yaxis.ticklabel.font=6
    graph.yaxis.ticklabel.char_size=1.4
    graph.xaxis.label.font=6
    graph.xaxis.label.char_size=1.65
    graph.xaxis.ticklabel.font=6
    graph.xaxis.ticklabel.char_size=1.4

def log(self,message, prefix_newline=False):
    time_str = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(time.time()))
    print(('\n' if prefix_newline else '') + time_str + ': ' + message)

class GenerateFigure:

    def __init__(self,data):

        grace = Grace()
        graph = grace.add_graph()
        graph.copy_format(ElGraph)
        graph.legend.loctype="view"
        graph.legend.loc=(0.2,0.4)
        define_graph(graph,"Sequence Identity Threshold","Fraction of Curated Functions",1.0,0.0,1.0,0.0,0.2,0.2,1,1,1,1)
        grace.set_fonts('Helvetica')

        Data_Colors={"Brassicaceae":2, #Red
                     "Eudicot":11, #Orange
                     "Liliopsida":12, #Purple
                     "Embryophyta":15, #Green
                     "Chlorophyta":4} #Blue

        for key in "Brassicaceae","Eudicot","Liliopsida","Embryophyta","Chlorophyta":
            Dataset=graph.add_dataset(data[key])
            Dataset.line.color=Data_Colors[key]
            Dataset.symbol.fill_color=Data_Colors[key]
            Dataset.symbol.shape=1
            Dataset.symbol.color=1
            Dataset.legend=key

        self.grace=grace
        self.graph=graph
        pass

    def generate_figure(self, figure_path, data_point):

        Dataset=self.graph.add_dataset([(data_point['threshold'],data_point['fraction'])])
        Dataset.line.type=0
        Dataset.symbol.shape=8
        Dataset.symbol.color=1
        Dataset.symbol.fill_color=1
        Dataset.symbol.size=1.0
        Dataset.symbol.linewidth=2.5
        Dataset.legend=data_point['id'].encode("ascii").decode("ascii")

        self.grace.write_file(os.path.join(figure_path,"Annotation_Threshold_Figure.eps"))
        self.grace.write_file(os.path.join(figure_path,"Annotation_Threshold_Figure.png"))
        self.grace.write_file(os.path.join(figure_path,"Annotation_Threshold_Figure.agr"))
        pass
