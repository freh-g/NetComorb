#!/usr/bin/env python
# coding: utf-8

import ast
from bs4 import BeautifulSoup as bs
import pandas as pd
import shlex
import subprocess

# # return bold text 

def bold(str):
    return('\033[1m'+str+'\033[0m')



# # Sort a list of tuples by chosen element



def tup_sorter(tup,indx,reverse=False):
    return sorted(tup, key=lambda tup: tup[indx],reverse=reverse)


# # split a string with new lines on every word (Word=true) or every character




def split_in_lines(string,n,words=False):
    if words:
        w = string.split()
        grouped_words = [' '.join(w[i: i + n]) for i in range(0, len(w), n)]
        return '\n'.join(grouped_words)
    grouped_chars=[string[i:i+n] for i in range(0, len(string), n)]
    return '\n'.join(grouped_chars)
    

    
## split a list in chuncks of a given length ###    
def chunker(seq, size):
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))
        


# # Average of a list

def average_list(lst):
    return sum(lst) / len(lst)

## Takes in input the soup of a website, finds the tables and return a Pandas Dataframe
def html_table_parser(html_soup):
    list_header = []
    data= []
    list_of_tables=[]


    header = html_soup.find_all("table")[0].find("tr")

    for items in header:
        try:
            list_header.append(items.get_text())
        except:
            continue



    # for getting the data 


    for _round in range(len(html_soup.find_all("table"))):
        HTML_data = html_soup.find_all("table")[_round].find_all("tr")[1:]

        for element in HTML_data:
            sub_data = []
            for sub_element in element:
                try:
                    sub_data.append(sub_element.get_text())
                except:
                    continue
            sub_data=[ele for ele in sub_data if ele !='\n']
            data.append(sub_data)
        list_of_tables.append(data)
    return pd.DataFrame(data = data, columns = list_header).dropna()


def run_shell_command(command):
    process = subprocess.Popen(shlex.split(command), stdout=subprocess.PIPE,shell=True)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print(output.strip().decode("utf-8"))
    rc = process.poll()
    return rc


## DATA HANDLING ##
def normalize_array(array,n=2):
    return [round((x-min(array))/(max(array)-min(array)),n) for x in array]

def normalize_df(dataframe, axis=None):
    if axis==None:
        df_max=max(dataframe.max())
        return dataframe/df_max
    elif axis == 0:
        return dataframe/dataframe.max()
    elif axis == 1:
        return (dataframe.T/dataframe.max(axis=1)).T
    
def standardize_df(dataframe, axis=0):
    if axis== 0:
        return (dataframe - dataframe.mean())/dataframe.std()
    else:
        return ((dataframe.T-dataframe.mean(axis=1))/dataframe.std(axis=1)).T

    
# Try to return an element in its original form, e.g. a string written as a list it is interpreted as a list
def original_form(element):
    try:
        return ast.literal_eval(element)
    except:
        return element
    

 









