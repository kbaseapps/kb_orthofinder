from bokeh.plotting import figure, output_file, save
from bokeh.models import NumeralTickFormatter

import time
import json

class GenerateFigureImpl:

    def log(self,message, prefix_newline=False):
        time_str = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(time.time()))
        print(('\n' if prefix_newline else '') + time_str + ': ' + message)

    def __init__(self,root):
        data = dict()
        Reference_Results_File = "Reference_Phytozome_Threshold.txt"
        with open(root+"/"+Reference_Results_File) as reference_results_handle:
            for line in reference_results_handle.readlines():
                line=line.strip()
                (x,group,y)=line.split('\t')
                y=float(y)
                x=float(x)

                if(group not in data):
                    data[group]={'x':[],'y':[]}
                data[group]['x'].append(x)
                data[group]['y'].append(y)
        self.data=data

    def generate_figure(self,x,y):

        bokeh_fig = figure(x_range = (0.0,1.0),
                           y_range = (0.0,1.0))
        bokeh_fig.xaxis.axis_label = "Sequence Identity Threshold"
        bokeh_fig.yaxis.axis_label = "Fraction of Curated Enzymes"
        bokeh_fig.xaxis.formatter = NumeralTickFormatter(format="0.0")
        bokeh_fig.yaxis.formatter = NumeralTickFormatter(format="0.0")

        data_order = ["Brassicaceae","Eudicot","Liliopsida","Embryophyta","Chlorophyta"]
        colors={"Brassicaceae":'red',
                "Eudicot":'orange',
                "Liliopsida":'purple',
                "Embryophyta":'green',
                "Chlorophyta":'blue'}

        for key in data_order:
            line_fig = bokeh_fig.line(self.data[key]['x'], self.data[key]['y'],
                                      line_width=2, color = colors[key],
                                      legend_label = key)
    
        point_fig = bokeh_fig.cross(x, y, line_width=2, color='black', size=12)
        bokeh_fig.legend.location = "bottom_left"
        return bokeh_fig

def main():

    figure_generator = GenerateFigureImpl("../../../data/")
    figure = figure_generator.generate_figure([0.5],[0.5])

    # Save output
    file_name='test_seqid_output.html'
    output_file(file_name)
    save(figure)

if(__name__ == "__main__"):
    main()
