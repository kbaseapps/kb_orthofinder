from urllib.request import urlopen
import time
import json
import sys
import os
import re

class GenerateTableImpl:

    def log(self,message, prefix_newline=False):
        time_str = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(time.time()))
        print(('\n' if prefix_newline else '') + time_str + ': ' + message)

    def __init__(self):
        pass

    def generate_table(self, functions):

        html_lines=list()        
        html_lines.append('<table class="table table-bordered table-striped">')

        header_list = ["Enzyme","Curated Feature","Predicted Ortholog","Potential Orthologs"]

        html_lines.append('<thead>')            
        internal_header_line = "</td><td>".join(header_list)
        html_lines.append('<tr><td>'+internal_header_line+'</td></tr>')
        html_lines.append('</thead>')

        html_lines.append("<tbody>")

        # Collecting hits that come from the same family
        # regardless of function
        family_hit_dict=dict()
        for function in functions:
            for family in functions[function]:

                if(family not in family_hit_dict):
                    family_hit_dict[family]=list()

                for hit in functions[function][family]['hits']:
                    family_hit_dict[family].append(hit['feature'])


        ############################################################################
        # Each row in the table is unique to the combination of enzymatic annotation
        # curated Arabidopsis enzyme, for which we show the retrieved predicted ortholog
        # but we only show the potential orthologs from the same protein family
        for function in sorted(list(functions.keys())):
            if(function == "Uncurated"):
                continue

            for family in functions[function]:
                for ortholog in sorted(functions[function][family]['orthologs']):

                    hit_list=list()
                    for hit in functions[function][family]['hits']:
                        hit_list.append(hit['feature'])

                    cls_list=list()
                    for cls_ftr in functions[function][family]['cluster'].keys():
                        for cls_olg in functions[function][family]['cluster'][cls_ftr]:
                            if(cls_olg == ortholog):
                                seqid=functions[function][family]['cluster'][cls_ftr][cls_olg]
                                
                                # reformat seqid to two decimal points
                                seqid = "{0:.2f}".format(float(seqid))

                                ftr_seqid_str = cls_ftr+":"+seqid
                                if(cls_ftr in hit_list):
                                    # color green
                                    ftr_seqid_str = "<font color=\"red\">"+ftr_seqid_str+"</font>"
                                elif(cls_ftr in family_hit_dict[family]):
                                     # color red
                                     ftr_seqid_str = "<font color=\"light blue\">"+ftr_seqid_str+"</font>"
                                else:
                                    # default color
                                    pass

                                if(cls_ftr in family_hit_dict):
                                    ftr_seqid_str = "<b>"+ftr_seqid_str+"</b>"
                                cls_list.append([float(seqid),ftr_seqid_str])

                    cls_sorted_list = list()
                    for entry in sorted(cls_list, key = lambda x: x[0], reverse=True):
                        cls_sorted_list.append(entry[1])

                    html_lines.append("<tr>")
                    internal_row_line = "</td><td>".join([function,ortholog,
                                                          ", ".join(sorted(hit_list)),
                                                          ", ".join(cls_sorted_list)])
                    html_lines.append("<td>"+internal_row_line+"</td>")
                    html_lines.append("</tr>")

        html_lines.append("</tbody>")
        html_lines.append("</table>")

        return "\n".join(html_lines)

def main():

    table = GenerateTableImpl()

    with open("annotation_output.json") as fh:
        annotation_output = json.load(fh)

    ##########################################################
    # Generate role table
    table_html_string = table.generate_table(annotation_output)

    with open(os.path.join('../../../data','app_report_templates',
                           'report_table_template.html')) as report_template_file:
        report_template_string = report_template_file.read()

    # Insert html table
    table_report_string = report_template_string.replace('*TABLES*', table_html_string)

    table_html_file="role_table_output.html"
    with open(os.path.join('../../../data','app_report_templates',
                           table_html_file),'w') as table_file:
        table_file.write(table_report_string)

    table.log(message="Role table written to file: "+table_html_file)

if(__name__ == "__main__"):
    main()
