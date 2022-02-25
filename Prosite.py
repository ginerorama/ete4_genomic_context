
import os
from Bio.ExPASy import Prosite



def convert_RE(expression):
    """
    Convierte las expresiones regulares del fichero 'prosite.dat' al formato
    reconocido por el paquete 're' de Python.
    """
    dict = {
        '-' : '',
        '{' : '[^',
        '}' : ']',
        '(' : '{',
        ')' : '}',
        'x' : '.',
        '>' : '$',
        '<' : '^'
        }

    for original, translation in dict.items():
        expression = expression.replace(original, translation)

    return str(expression)


def create_prosite_dict():
    """
    Crea un diccionario: que relaciona los patrones RE de cada dominio
    con el número de accessión del mismo, name y descripción.

    Además, si el parámetro record.cc_skip_flag == "TRUE", se excluyen del diciconario
    dominios con alta probabilidad de ocurrencia (sistios de fosforilación,
    etc.), que contienen la etiqueta '/SKIP-FLAG=TRUE>' 
    """
    dict_pattern = {}
    with open('./prosite.dat','r') as handle:
        records = Prosite.parse(handle)
        for record in records:
            if record.pattern != "": # Hay algunos records que no tienen patrón.
                #print(record.pattern)
                converted_RE = convert_RE(record.pattern)
                    # Recupera el nombre y accessión del dominio.
                if "site" not in str(record.description):
                    dict_pattern[converted_RE] = [str(record.accession),str(record.name),str(record.description)]

        #print(dict_pattern)
    return dict_pattern


def search_domains(dict_pattern):
    import re
    """
    Recorre todas las instancias de 'Protein', y busca en su secuencia de
    aminoácidos coincidencias con alguna de las expresiones regulares de
    el diccionario de Prosite.

    Almacena los matches en la lista 'self.domains' de la instancia
    correspondiente, en forma de tuplas (accession, start, end).
    """
    
    handle = open("gcontext.csv","r")
    for record in handle:
        fields = record.rstrip().split(",")
        name = fields[0]
        protein = fields[-1]
        print(name)
        for pattern in list(dict_pattern.keys()):
            ungapped_seq = protein.replace("-", "")
            matches = re.finditer(pattern, ungapped_seq)       
            for match in matches:
                'PS00538', 'CHEMOTAXIS_TRANSDUC_1', 'Bacterial chemotaxis sensory transducers signature.'
                accesion,name,description = dict_pattern[pattern]
                protein_match =accesion,name,description,match.start(),match.end()
                print(protein_match)
 
 
 

#dict_pattern = create_prosite_dict()
#search_domains(dict_pattern)